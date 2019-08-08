function ekfpredict_sim(u)
% EKF-SLAM prediction for simulator process model

global Param;
global State;
%-----------------------------------------------
%Define variables 
% -------------------------------------------
%KNOWN
% mu         mean of state 
% Sigma     Variance of state
% u,         [drot1, dtrans, drot2]' Odometry control command
    % x' = x + ˆ?trans cos(? + ˆ?rot1)
    % y' = y + ˆ?trans sin(? + ˆ?rot1)
    % ?' = ? + ˆ?rot1 + ˆ?rot2
% deltaT     Step size between filter updates, can be less than 1.
% M,         Motion noise
% z          bearing_noise observation 
% Q          Observation noise
% markerId   landmark id

%OUTPUT
% H          Jaxobian of h w.r.t location 
% K          Kalman gain
% G          Jacobian of g w,r,t location 
% V          Jacobian of g w.r.t control
% R          transfromed motion noise covariance
% zhat,     prediction measurement mean

% Set robot initial pose(mu,sigma) 
mu_p=State.Ekf.mu ; %mean of state 
sigma_p=State.Ekf.Sigma ;%Variance of state

% set the map of landmarks
landmark_num=State.Ekf.nL;

%Calculate F_x
F_x=[eye(3), zeros(3,2*landmark_num)];

%obtain drot1, dtrans, drot2 control command
drot1=u(1);dtrans=u(2);drot2=minimizedAngle(u(3));

%obtain state of robot
x=mu_p(1);y=mu_p(2);theta=minimizedAngle(mu_p(3));

%Odometry Model
odm=[dtrans*cos(theta+drot1);dtrans*sin(theta+drot1);(drot1+drot2)];

% Update EKF prediction of mean 
mu=mu_p+F_x'*odm;

% Update EKF predictioand covariance
g=[0,0,-dtrans*sin(theta+drot1);
   0,0,dtrans*cos(theta+drot1);
   0,0,0];
G_t=eye(length(mu_p))+F_x'*g*F_x;
%obtain noise a1 a2 a3 a4
a1=Param.alphas(1);
a2=Param.alphas(2);
a3=Param.alphas(3);
a4=Param.alphas(4);
%calculate M
M=[a1*abs(drot1^2)+a2*abs(dtrans^2),0,0;...
        0,a3*abs(dtrans^2)+a4*abs(drot1^2+drot2^2),0;
        0,0,a1*abs(drot2^2)+a2*abs(dtrans^2)];
%Calculate V
V=[ -dtrans*sin(drot1+theta),cos(drot1+theta),0;
      dtrans*cos(drot1+theta),sin(drot1+theta),0;
      1,0,1]; 
R=V*M*V' ;
sigma=G_t*sigma_p*G_t'+F_x'*R*F_x;
%Update the mean and covariance
State.Ekf.mu=mu;  %mean of state 
State.Ekf.Sigma=sigma; %Variance of state
end

