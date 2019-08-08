function Li = da_nn(z, R)
% perform nearest-neighbor data association

global Param;
global State;

Li=[];
m=size(z,2);
%obtain the robot location x, y , theta
u_t=State.Ekf.mu(1:3);
pose_t=u_t(1:2);
theta=u_t(3);

%start of Neareat Neighbor
for i=1:m
   %add varible to calculating the ML correspondences
   Pi=[];
   for j=1:State.Ekf.nL
       %find the difference with the jth landmark measurement
       landmark_j=State.Ekf.mu(3+2*j-1:3+2*j);
       delta=landmark_j-pose_t;
       q=delta'*delta;
       %obtian delta_x and delta_y
       delta_x=delta(1);
       delta_y=delta(2);
       %calculate z_t_i the expected observation according to current
       %estimation
       z_t_i=[sqrt(q);atan2(delta_y,delta_x)-theta];
       z_t_i(2)=minimizedAngle(z_t_i(2));
       
       %calculate the jacobian of low Ht_i
       H_t_il=(1/q)*[-sqrt(q)*delta_x,-sqrt(q)*delta_y,0,sqrt(q)*delta_x,sqrt(q)*delta_y;
                  delta_y,-delta_x,-q,-delta_y,delta_x];
       %calculate F_x_j     
       F_x_j=zeros(3+2,size(State.Ekf.mu,1));
       F_x_j(1:3,1:3)=eye(3);
       %F_x_j(4:5,2*i+2:2*i+3)=eye(2);
       F_x_j(4:5,3+2*j-1:3+2*j) = eye(2);
       %map the low Ht_i to high dimensinal space
       H_t=H_t_il*F_x_j;
       %calculate delta 
       if(i<=size(z,2))
       delta_t_i=z(1:2,i)-z_t_i;
       delta_t_i(2)=minimizedAngle(delta_t_i(2));
       end
       %calculate Psi
       psi_k=H_t*State.Ekf.Sigma*H_t'+R;
       %calculate pi
       pi_k=delta_t_i'*pinv(psi_k)*delta_t_i;
       pi_k=sqrt(pi_k);
       Pi=[Pi,pi_k];
   end
   %calculate the chi2(di,alpha)
   di=0.82;
   %degree of freedom
   alpha=2;
   chi=chi2inv(di,alpha);
   %find the index and minimin value
   [m,index]=min(Pi);

   if (m<=chi)
       %if the measured landmark was measured before that the differece of
       %distance smaller than the chi
       Li=[Li,index];
   else
       %if it was not measured before, give 0 to show it is new landmark
       Li=[Li, 0];
   end
   
end



end
