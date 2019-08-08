function Li = da_known(z)
% EKF-SLAM data association with known correspondences

global Param;
global State;
% z          3xn [detected landmark id]
% State.Ekf.t     = 0;          % time
% State.Ekf.mu    = zeros(3,1); % robot initial pose
% State.Ekf.Sigma = zeros(3,3); % robot initial covariance
% State.Ekf.iR    = 1:3;        % 3 vector containing robot indices
% State.Ekf.iM    = [];         % 2*nL vector containing map indices
% State.Ekf.iL    = {};         % nL cell array containing indices of landmark i
% State.Ekf.sL    = [];         % nL vector containing signatures of landmarks
% State.Ekf.nL    = 0;          % scalar number of landmarks
Li=[];
for i=1:length(z)

    %current signiture of observed landmark 
    j=z(i);  
    if ~isempty(j)
    %iR=State.Ekf.iR 
    %check whether the landmark has been observed
    current_position=find(State.Ekf.sL==j);
        if (~isempty(current_position))
        Li=[Li current_position];
        
        else
        Li=[Li 0];
        end
    end
end
