function runvp(nSteps,pauseLen)

global Param;
global State;
global Data;

if ~exist('nSteps','var') || isempty(nSteps)
    nSteps = inf;
end

if ~exist('pauseLen','var')
    pauseLen = 0; % seconds
end

Data = load_vp_si();

% Initalize Params
%===================================================
% vehicle geometry
Param.a = 3.78; % [m]
Param.b = 0.50; % [m]
Param.L = 2.83; % [m]
Param.H = 0.76; % [m]

% 2x2 process noise on control input
sigma.vc = 0.02; % [m/s]
sigma.alpha = 2*pi/180; % [rad]
Param.Qu = diag([sigma.vc, sigma.alpha].^2);

% 3x3 process noise on model error
sigma.x = 0.1; % [m]
sigma.y = 0.1; % [m]
sigma.phi = 0.5*pi/180; % [rad]
Param.Qf = diag([sigma.x, sigma.y, sigma.phi].^2);

% 2x2 observation noise
sigma.r = 0.05; % [m]
sigma.beta = 1*pi/180; % [rad]
Param.R = diag([sigma.r, sigma.beta].^2);
%===================================================

% Initialize State
%===================================================
State.Ekf.mu = [Data.Gps.x(2), Data.Gps.y(2), 36*pi/180]';
State.Ekf.Sigma = zeros(3);


global AAr;
AAr = [0:360]*pi/360;

%Take the video

makeVideo=false;
if makeVideo
    try
        votype = 'avifile';
        vo = avifile('Task3.avi', 'fps', min(5, 1/pauseLen));
    catch
        votype = 'VideoWriter';
        vo = VideoWriter('video', 'MPEG-4');
        set(vo, 'FrameRate', min(5, 1/pauseLen));
        open(vo);
    end
end

%set up parameters for time
Predition_time=[];% predition time ; real time
Update_time=[];% Update time ; real time
landmark_nm=[];% landmark number ; real time
State.Ekf.groundtruth = [];%record ground truth


figure(1); clf;
axis equal;

ci = 1; % control index
t = min(Data.Laser.time(1), Data.Control.time(1));
for k=1:min(nSteps, length(Data.Laser.time))
    
    while (Data.Control.time(ci) < Data.Laser.time(k))
       % control available
       dt = Data.Control.time(ci) - t;
       t = Data.Control.time(ci);
       u = [Data.Control.ve(ci), Data.Control.alpha(ci)]';
       %start to count time
       tic;
       ekfpredict_vp(u, dt);
       Pred_time=toc;
       Predition_time=[Predition_time,[Pred_time;t]];
       ci = ci+1;
    end
    
    % observation available
    dt = Data.Laser.time(k) - t;
    t = Data.Laser.time(k);
    z = detectTreesI16(Data.Laser.ranges(k,:));

    % start count time 
    tic;
    %update the z beta since it need to minus pi/2
    z_new=z;
    z_new(2,:)=z_new(2,:)-pi/2;
    z_new(2,:)=minimizedAngle(z_new(2,:));
    ekfupdate(z_new);
    update_t=toc;
    Update_time=[Update_time,[update_t;t]];
    landmark_nm=[landmark_nm,[State.Ekf.nL;t]];
   
    %recored ground truth
    t_g=find(Data.Gps.time<=t);
    t_g=t_g(end);
    if(~isempty(t_g))
        State.Ekf.groundtruth=[State.Ekf.groundtruth,[Data.Gps.x(t_g);Data.Gps.y(t_g)]];
    end
    
    
    
    doGraphics(z); 
    
    drawnow;
    if pauseLen > 0
        pause(pauseLen);
    end
    if makeVideo
        F = getframe(gcf);
        switch votype
          case 'avifile'
            vo = addframe(vo, F);
          case 'VideoWriter'
            writeVideo(vo, F);
          otherwise
            error('unrecognized votype');
        end
    end
    
end
%Making vedio
    
    if makeVideo
    fprintf('Writing video...');
    switch votype
      case 'avifile'
        vo = close(vo);
      case 'VideoWriter'
        close(vo);
      otherwise
        error('unrecognized votype');
    end
    fprintf('done\n');
    end
    %Plot figures
    disp('the finial map') 
    figure(2)
    hold on
    plot(Predition_time(2,:),Predition_time(1,:),'b-','DisplayName','Prediction time')
    plot(Update_time(2,:),Update_time(1,:),'r-','DisplayName','Update time')
    title('CPU time of prediction delta and update delta')
    legend('Location','NorthWest')
    xlabel('time(s)')
    ylabel('delta CPU time')
    figure(3)
    plot(landmark_nm(2,:),landmark_nm(1,:),'k-','DisplayName','number of landmarks')
    legend('Location','NorthWest')
    xlabel('time(s)')
    ylabel('number of landmarks')
    title('number of landmarks in the map for each iteration')

end

%==========================================================================
function doGraphics(z)
% Put whatever graphics you want here for visualization
%
% WARNING: this slows down your process time, so use sparingly when trying
% to crunch the whole data set!

global Param;
global State;

% plot the robot and 3-sigma covariance ellipsoid
plotbot(State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.mu(3), 'black', 1, 'blue', 1);
hold on
plotcov2d( State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma, 'blue', 0, 'blue', 0, 3);

% restrict view to a bounding box around the current pose
BB=20;
axis([[-BB,BB]+State.Ekf.mu(1), [-BB,BB]+State.Ekf.mu(2)]);

% project raw sensor detections in global frame using estimate pose
xr = State.Ekf.mu(1);
yr = State.Ekf.mu(2);
tr = State.Ekf.mu(3);
for k=1:size(z,2)
    r = z(1,k);
    b = z(2,k);
    xl = xr + r*cos(b+tr-pi/2);
    yl = yr + r*sin(b+tr-pi/2);
    plot([xr; xl], [yr; yl],'r',xl,yl,'r*');
end
    %plot robot pose
    
    plotcov2d(State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma(1:2,1:2), 'r', false, [], [], 3);
    plotmarker(State.Ekf.mu(1:2), 'r');
    
    %plot ground truth 
    plot(State.Ekf.groundtruth(1,:),State.Ekf.groundtruth(2,:),'k-')
    % Draw landmark uncertainty
    for i = 1:State.Ekf.nL
    plotcov2d(State.Ekf.mu((2*i+3-1)), State.Ekf.mu((2*i+3)), State.Ekf.Sigma((3+2*i-1):(3+2*i),(3+2*i-1):(3+2*i)), 'm',false,[],[],3);
   
    end

hold off;

end