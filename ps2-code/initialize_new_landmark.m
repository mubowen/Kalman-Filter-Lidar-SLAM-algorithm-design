function initialize_new_landmark(z, R)

global Param;
global State;

%obtain the robot location
u_t=State.Ekf.mu(1:3);

% noisy observations [range; bearing; landmark id]
r_t=z(1); 
phi_t=z(2);

%update observed location of landmark 
m_relat=[r_t*cos(phi_t+u_t(3));r_t*sin(phi_t+u_t(3))];

%update the observed location of landmark j
u_j=[u_t(1);u_t(2)]+m_relat;

%update the state mean
State.Ekf.mu=[State.Ekf.mu;u_j];

%update the state covariance
%set up g_r g_z
g_r=[1,0,r_t*cos(u_t(3)+phi_t);
    0,1,r_t*sin(u_t(3)+phi_t)];

g_z=[cos(u_t(3)+phi_t),-r_t*sin(u_t(3)+phi_t);
     sin(u_t(3)+phi_t),r_t*cos(u_t(3)+phi_t)];

sigma=State.Ekf.Sigma;
R=Param.R;
n=State.Ekf.nL;

if (State.Ekf.nL==0)
sigma_new=[sigma,(g_r*sigma(3,3))';
          g_r*sigma(3,3),g_r*sigma(3,3)*g_r'+g_z*R*g_z'];

else

sigma_xx=sigma(1:3,1:3);
%sigma_xm=sigma(1:3,4:end);
sigma_xm=sigma(1:3,State.Ekf.iM);
sigma_mm=sigma(State.Ekf.iM,State.Ekf.iM);
%sigma_mm=sigma(4:end,4:end);

% sigma_mm=[sigma(4:end,4:end),(g_r*sigma_xm_old);...
%    (g_r*sigma_xm_old)', g_r*sigma_xx*g_r'+g_z*R*g_z'];

sigma_new=[sigma_xx,sigma_xm,sigma_xx'*g_r';
           sigma_xm',sigma_mm,sigma_xm'*g_r';
           g_r*sigma_xx,g_r*sigma_xm,g_r*sigma_xx*g_r'+g_z*R*g_z'];


end

State.Ekf.Sigma=sigma_new;
 
 
%update the scalar number of landmarks
State.Ekf.nL = State.Ekf.nL +1; 

%update the nL vector containing signatures of landmarks
State.Ekf.sL = [State.Ekf.sL,z(3)]; 


%update the 2*nL vector containing map indices
%State.Ekf.iM=[State.Ekf.iM,u_j] ;  
State.Ekf.iM=(State.Ekf.iR(end)+1):length(State.Ekf.mu);
% update nL cell array containing indices of landmark i
State.Ekf.iL{State.Ekf.nL} = u_j';   

end
