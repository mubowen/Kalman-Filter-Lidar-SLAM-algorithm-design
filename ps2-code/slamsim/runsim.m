function varargout = runsim(stepsOrData, pauseLen)

global Param;
global Data;
global State;

if ~exist('pauseLen','var')
    pauseLen = 0.3; % seconds
end
makeVideo=false;
if makeVideo
    try
        votype = 'avifile';
        vo = avifile('video.avi', 'fps', min(5, 1/pauseLen));
    catch
        votype = 'VideoWriter';
        vo = VideoWriter('video', 'MPEG-4');
        set(vo, 'FrameRate', min(5, 1/pauseLen));
        open(vo);
    end
end

% Initalize Params
%===================================================
Param.initialStateMean = [180 50 0]';

% max number of landmark observations per timestep
Param.maxObs = 2;

% number of landmarks per sideline of field (minimum is 3)
Param.nLandmarksPerSide = 4;

% Motion noise (in odometry space, see p.134 in book).
Param.alphas = [0.05 0.001 0.05 0.01].^2; % std of noise proportional to alphas

% Standard deviation of Gaussian sensor noise (independent of distance)
Param.beta = [10, deg2rad(10)]; % [cm, rad]
Param.R = diag(Param.beta.^2);

% Step size between filter updates, can be less than 1.
Param.deltaT=0.1; % [s]

if isscalar(stepsOrData)
    % Generate a data set of motion and sensor info consistent with
    % noise models.
    numSteps = stepsOrData;
    Data = generateScript(Param.initialStateMean, numSteps, Param.maxObs, Param.alphas, Param.beta, Param.deltaT);
else
    % use a user supplied data set from a previous run
    Data = stepsOrData;
    numSteps = size(Data, 1);
    global FIELDINFO;
    FIELDINFO = getfieldinfo;
end
%===================================================

% Initialize State
%===================================================
State.Ekf.mu = Param.initialStateMean;
State.Ekf.Sigma = zeros(3);




for t = 1:numSteps
    plotsim(t);

    %=================================================
    % data available to your filter at this time step
    %=================================================
    u = getControl(t);
    z = getObservations(t);


    %=================================================
    %TODO: update your filter here based upon the
    %      motionCommand and observation
    %=================================================
    State.Ekf.t=t;
    ekfpredict_sim(u);
    ekfupdate(z);

    %=================================================
    %TODO: plot and evaluate filter results here
    %=================================================
    %plot robot uncertianty
    
    plotcov2d(State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma(1:2,1:2), 'r', false, [], [], 3);
    plotmarker(State.Ekf.mu(1:2), 'r');
    % Draw landmark uncertainty
        for i = 1:State.Ekf.nL
            plotcov2d(State.Ekf.mu((2*i+3-1)), State.Ekf.mu((2*i+3)), State.Ekf.Sigma((3+2*i-1):(3+2*i),(3+2*i-1):(3+2*i)), 'm',false,[],[],3);
        end
        
    %{
    %record the determinate of covarriance matrix
    detcov=zeros(8,1);
    for i=1:State.Ekf.nL
    matrix_a=State.Ekf.Sigma((3+2*i-1):(3+2*i),(3+2*i-1):(3+2*i));
    det_matrix=det(matrix_a);
    detcov(i,1)=log(det_matrix);
   
    end
    State.Ekf.detcov=[State.Ekf.detcov,detcov];
    
    %record the correlate coefficient between landmarks
   
    
    coor_m=corrcoef(State.Ekf.Sigma(4:end,4:end));
    corr_m_x=zeros(8,1);
    corr_m_y=zeros(8,1);
    
    for ind=1:State.Ekf.nL
        corr_m_x(ind,1)=coor_m(1,ind*2-1);
        corr_m_y(ind,1)=coor_m(1,ind*2);
    end
    
    State.Ekf.corr_x=[State.Ekf.corr_x,corr_m_x];
    State.Ekf.corr_y=[State.Ekf.corr_y,corr_m_y];
    
    corr_m_x=zeros(8,1);
    corr_m_y=zeros(8,1);
    
    for ind=1:State.Ekf.nL
        corr_m_x(ind,1)=coor_m(2,ind*2-1);
        corr_m_y(ind,1)=coor_m(2,ind*2);
    end
    State.Ekf.corr_x2=[State.Ekf.corr_x,corr_m_x];
    State.Ekf.corr_y2=[State.Ekf.corr_y,corr_m_y];
    
    corr_m_x=zeros(8,1);
    corr_m_y=zeros(8,1);
    
    for ind=1:State.Ekf.nL
        corr_m_x(ind,1)=coor_m(3,ind*2-1);
        corr_m_y(ind,1)=coor_m(3,ind*2);
    end
    State.Ekf.corr_x3=[State.Ekf.corr_x,corr_m_x];
    State.Ekf.corr_y3=[State.Ekf.corr_y,corr_m_y];
     %}
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

if nargout >= 1
    varargout{1} = Data;
end
hold off


%{ 
 fig=2;
%plot the determinate of covarriance matrix
 figure(fig)
    for i=1:8
    plot(State.Ekf.detcov(i,:))
    hold on 
    end
    hold off
    xlabel('time step')
    ylabel('determinant(log scale)')
    title('Determinant of the feature covariance matrix for all map feature')
    legend('landmark 1','landmark 2','landmark 3','landmark 4','landmark 5',...
    'landmark 6','landmark 7','landmark 8')
%plot correlation matrix
fig=fig+1;
 figure(fig)

 corr_m=corrcov(State.Ekf.Sigma(4:end,4:end));
 imagesc(corr_m)
 set(gca, 'XTick', 1:16); % center x-axis ticks on bins
 set(gca, 'YTick', 1:16); % center y-axis ticks on bins
 title('correlation matrix', 'FontSize', 14); % set title
 colormap('jet'); % set the colorscheme
 axis([1,16,1,16])
 fig=fig+1;
 figure(fig)
 for i=1:7 
     hold on 
     plot(State.Ekf.corr_x(i,:))   
     plot(State.Ekf.corr_y(i,:))
 end
    hold off
 title('correlate coefficient for first landmark')  
  

 fig=fig+1;
 figure(fig)
 for i=1:7 
     hold on 
     plot(State.Ekf.corr_x2(i,:))   
     plot(State.Ekf.corr_y2(i,:))
 end
    hold off
title('correlate coefficient for second landmark') 
    
     fig=fig+1;
 figure(fig)
 for i=1:7 
     hold on 
     plot(State.Ekf.corr_x3(i,:))   
     plot(State.Ekf.corr_y3(i,:))
 end
    hold off
    title('correlate coefficient for third landmark')
    fig=fig+1;
%}
    
end

%==========================================================================
function u = getControl(t)
global Data;
% noisefree control command
u = Data.noisefreeControl(:,t);  % 3x1 [drot1; dtrans; drot2]

end
%==========================================================================
function z = getObservations(t)
global Data;
% noisy observations
z = Data.realObservation(:,:,t); % 3xn [range; bearing; landmark id]
ii = find(~isnan(z(1,:)));
z = z(:,ii);
end
%==========================================================================
function plotsim(t)
global Data;

%--------------------------------------------------------------
% Graphics
%--------------------------------------------------------------

NOISEFREE_PATH_COL = 'green';
ACTUAL_PATH_COL = 'blue';

NOISEFREE_BEARING_COLOR = 'cyan';
OBSERVED_BEARING_COLOR = 'red';

GLOBAL_FIGURE = 1;

%=================================================
% data *not* available to your filter, i.e., known
% only by the simulator, useful for making error plots
%=================================================
% actual position (i.e., ground truth)
x = Data.Sim.realRobot(1,t);
y = Data.Sim.realRobot(2,t);
theta = Data.Sim.realRobot(3,t);

% real observation
observation = Data.realObservation(:,:,t);

% noisefree observation
noisefreeObservation = Data.Sim.noisefreeObservation(:,:,t);

%=================================================
% graphics
%=================================================
figure(GLOBAL_FIGURE); clf; hold on; plotfield(observation(3,:));

% draw actual path (i.e., ground truth)
plot(Data.Sim.realRobot(1,1:t), Data.Sim.realRobot(2,1:t), 'Color', ACTUAL_PATH_COL);
plotrobot( x, y, theta, 'black', 1, ACTUAL_PATH_COL);

% draw noise free motion command path
plot(Data.Sim.noisefreeRobot(1,1:t), Data.Sim.noisefreeRobot(2,1:t), 'Color', NOISEFREE_PATH_COL);
plot(Data.Sim.noisefreeRobot(1,t), Data.Sim.noisefreeRobot(2,t), '*', 'Color', NOISEFREE_PATH_COL);

for k=1:size(observation,2)
    rng = Data.Sim.noisefreeObservation(1,k,t);
    ang = Data.Sim.noisefreeObservation(2,k,t);
    noisy_rng = observation(1,k);
    noisy_ang = observation(2,k);

    % indicate observed range and angle relative to actual position
    plot([x x+cos(theta+noisy_ang)*noisy_rng], [y y+sin(theta+noisy_ang)*noisy_rng], 'Color', OBSERVED_BEARING_COLOR);

    % indicate ideal noise-free range and angle relative to actual position
    plot([x x+cos(theta+ang)*rng], [y y+sin(theta+ang)*rng], 'Color', NOISEFREE_BEARING_COLOR);
end
end