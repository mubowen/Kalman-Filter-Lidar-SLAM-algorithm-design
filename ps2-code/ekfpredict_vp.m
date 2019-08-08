function ekfpredict_vp(u, dt)
% EKF-SLAM prediction for Victoria Park process model

%parameter identification
% Param.a           [m]% vehicle geometry
% Param.b           [m]
% Param.L 
% Param.H
% Param.Qu          2x2 process noise on control input
% Param.Qf          3x3 process noise on model error
% Param.R           2x2 observation noise
% State.Ekf.mu      mean of state 
% State.Ekf.Sigma   Variance of state
% AAr

% ci               % control index
% u                 control
global Param;
global State;

%[m]% vehicle geometry
a=Param.a ;          
b=Param.b;        
L=Param.L; 
H=Param.H;

% 2x2 process noise on control input
Qu=Param.Qu;     
% 3x3 process noise on model erro
Qf=Param.Qf; 
% 2x2 observation noise
R=Param.R;  

% Set robot initial pose(mu,sigma) 
mu_p=State.Ekf.mu ; %mean of state 
sigma_p=State.Ekf.Sigma ;%Variance of state

% set the map of landmarks
landmark_num=State.Ekf.nL;

%Calculate F_x
F_x=[eye(3), zeros(3,2*landmark_num)];

%Using Ackerman Model update motion contorl[ve,a]
ve=u(1); alpha=u(2);
%obtain the velocity vc which is the encoder locate in the back left wheel
vc=ve/((1-tan(alpha)*H/L));

%obtain robot state
x=mu_p(1); y=mu_p(2); theta=mu_p(3);

%process the updated pose of robot

g=[dt*(vc*cos(theta)-vc/L*tan(alpha)*(a*sin(theta)+b*cos(theta)));
    dt*(vc*sin(theta)+vc/L*tan(alpha)*(a*cos(theta)-b*sin(theta)));
    dt*vc/L*tan(alpha)];
%update mu
mu=mu_p+F_x'*g;
mu(3)=minimizedAngle(mu(3));

%observation model
%motion model jacobians
Gx=[0,0,-dt*(vc*sin(theta)+vc/L*tan(alpha)*(a*cos(theta)-b*sin(theta)));
    0,0,dt*(vc*cos(theta)+vc/L*tan(alpha)*(-a*sin(theta)-b*cos(theta)));
    0,0,0];
G_t=eye(length(mu_p))+F_x'*Gx*F_x;


%update the motion noise ( jacobian of g , dvc/dt,dalpha/dt)
V=[dt*(cos(theta)-1/L*tan(alpha)*(a*sin(theta)+b*cos(theta))),dt*-vc/L*(sec(alpha))^2*(a*sin(theta)+b*cos(theta));
    dt*(sin(theta)+1/L*tan(alpha)*(a*cos(theta)-b*sin(theta))),dt*vc/L*(sec(alpha))^2*(a*cos(theta)-b*sin(theta));
    dt*1/L*tan(alpha),dt*vc/L/(cos(alpha))^2];
%obtian Rt_x
Rt_x=V*Qu*V' ;

%Update the covariance 
sigma=G_t*sigma_p*G_t'+F_x'*Rt_x*F_x+F_x'*Qf*F_x;
%Update the mean and covariance
State.Ekf.mu=mu;  %mean of state 
State.Ekf.Sigma=sigma; %Variance of state


end

