function ekfupdate(z)
% EKF-SLAM update step for both simulator and Victoria Park data set

global Param;
global State;

%KNOWN
% mu         mean of state 
% Sigma     Variance of state
% u,         [drot1, dtrans, drot2]' Odometry control command
    % x' = x + ˆ?trans cos(? + ˆ?rot1)
    % y' = y + ˆ?trans sin(? + ˆ?rot1)
    % ?' = ? + ˆ?rot1 + ˆ?rot2
% deltaT     Step size between filter updates, can be less than 1.
% M,         Motion noise
% z          3xn [range; bearing; landmark id]
% Q          Observation noise
% markerId   landmark id

%OUTPUT
% H          Jaxobian of h w.r.t location 
% K          Kalman gain
% G          Jacobian of g w,r,t location 
% V          Jacobian of g w.r.t control
% R          transfromed motion noise covariance
% zhat,     prediction measurement mean
        
      

% returns state vector indices pairing observations with landmarks
switch lower(Param.dataAssociation)
    case 'known'
        Li = da_known(z(3,:));
    case 'nn'
        Li = da_nn(z(1:2,:), Param.R);
    case 'jcbb'
        Li = da_jcbb(z(1:2,:), Param.R);
    otherwise
        error('unrecognized data association method: "%s"', Param.dataAssociation);
end

 switch lower(Param.updateMethod)
     case 'batch'
     %obtain Qt sensor observation noise
     Q_t_i=Param.R;
     %obtian the variable to be stacked
     H_t=[];
     Q_t=[];
     delta_t=[];
     for i=1:length(Li) 
         %check whether it is new landmark
         if(Li(i)==0) 
         %add a new landmark if it is not observed
         initialize_new_landmark(z(:,i), Q_t)
         Li(i)=State.Ekf.nL;
         end
       
     end
     for i=1:length(Li) 
     %obtain the pose of robot
     u_t=State.Ekf.mu(1:3);
     %obtain the pose and index of landmark
     %landmark_id=State.Ekf.sL(Li(i));
     landmark_id=Li(i);
     u_j=State.Ekf.mu(State.Ekf.iM((landmark_id*2-1):landmark_id*2));
     %State.Ekf.iM(:,Li(i));
     
     %calculate delta_x and delta_y
     delta=u_j-u_t(1:2);
     delta_x=delta(1);
     delta_y=delta(2);
     q=delta'*delta;
     %calculate the expected observation according to the current estimate
     z_t_i=[sqrt(q);atan2(delta_y,delta_x)-u_t(3)];
     z_t_i(2)=minimizedAngle(z_t_i(2));
     %calculate the jacobian of low Ht_i
     H_t_il=(1/q)*[-sqrt(q)*delta_x,-sqrt(q)*delta_y,0,sqrt(q)*delta_x,sqrt(q)*delta_y;
                  delta_y,-delta_x,-q,-delta_y,delta_x];
     %calculate F_x_j
     n=State.Ekf.nL;
     
     F_x_j=zeros(3+2,size(State.Ekf.mu,1));
     F_x_j(1:3,1:3)=eye(3);
     %F_x_j(4:5,2*i+2:2*i+3)=eye(2);
     F_x_j(4:5,State.Ekf.iM((landmark_id*2-1):landmark_id*2)) = eye(2);
     %map the low Ht_i to high dimensinal space
     H_t_i=H_t_il*F_x_j;
     %stack H_t
     H_t=[H_t;H_t_i];
     
     
     %stack Q_t
     Q_t=blkdiag(Q_t,Q_t_i);
     
     %stack delta of (z_t-u(t))
     
     delta_t_i=z(1:2,i)-z_t_i;
     delta_t_i(2)=minimizedAngle(delta_t_i(2));
     delta_t=[delta_t;delta_t_i];
     end
     %Batch update 
     S_t=H_t*State.Ekf.Sigma*H_t'+Q_t;
     K_t=State.Ekf.Sigma*H_t'*inv(S_t);
     %update mu
     mu=State.Ekf.mu+K_t*delta_t;
     State.Ekf.mu=mu;
     %Update Covariance
     sigma=(eye(size(K_t*H_t))-K_t*H_t)*State.Ekf.Sigma;
     State.Ekf.Sigma=sigma;

 case'seq'
     %obtain Qt sensor observation noise
     Q_t_i=Param.R;
     %obtian the variable
     Q_t=[];
       
     for i=1:length(Li) 
         %check whether it is new landmark
         if(Li(i)==0) 
         %add a new landmark if it is not observed
         initialize_new_landmark(z(:,i), Q_t)
         Li(i)=State.Ekf.nL;
         end
       
     %obtain the pose of robot
     u_t=State.Ekf.mu(1:3);
     %obtain the pose and ID of landmark
     %landmark_id=State.Ekf.sL(Li(i));
     landmark_id=Li(i);
     u_j=State.Ekf.mu(State.Ekf.iM((landmark_id*2-1):landmark_id*2));
     
     %calculate delta_x and delta_y
     delta=u_j-u_t(1:2);
     delta_x=delta(1);
     delta_y=delta(2);
     q=delta'*delta;
     %calculate the expected observation according to the current estimate
     z_t_i=[sqrt(q);atan2(delta_y,delta_x)-u_t(3)];
     z_t_i(2)=minimizedAngle(z_t_i(2));
     %calculate the jacobian of low Ht_i
     H_t_il=(1/q)*[-sqrt(q)*delta_x,-sqrt(q)*delta_y,0,sqrt(q)*delta_x,sqrt(q)*delta_y;
                  delta_y,-delta_x,-q,-delta_y,delta_x];
     %calculate F_x_j
     n=State.Ekf.nL;
     
     F_x_j=zeros(3+2,size(State.Ekf.mu,1));
     F_x_j(1:3,1:3)=eye(3);
     F_x_j(4:5,State.Ekf.iM((landmark_id*2-1):landmark_id*2)) = eye(2);

     %map the low Ht_i to high dimensinal space
     H_t=H_t_il*F_x_j;
     
     
     %obtain Q_t
     Q_t=Q_t_i;
     
     %obtain delta of (z_t-u(t))
      
     delta_t_i=z(1:2,i)-z_t_i;
     delta_t_i(2)=minimizedAngle(delta_t_i(2));
     delta_t=delta_t_i;
     
     %sequential Update
     
     S_t=H_t*State.Ekf.Sigma*H_t'+Q_t;
     K_t=State.Ekf.Sigma*H_t'*pinv(S_t);
     %update mu
     mu=State.Ekf.mu+K_t*delta_t;
     State.Ekf.mu=mu;
     %Update Covariance
     sigma=(eye(size(K_t*H_t))-K_t*H_t)*State.Ekf.Sigma;
     State.Ekf.Sigma=sigma;
     
     end

    
     
     


 end 
end
