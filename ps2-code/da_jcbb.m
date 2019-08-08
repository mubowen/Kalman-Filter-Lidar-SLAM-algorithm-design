function   Li = da_jcbb(z, R)
% perform joint-compatability branch and bound data association

global Param;
global State;
global Jcbb;

%set up the threshold for the Mahalanobis distance of the observation with respect to the noisy
threshold=0.01;
%size of pairing 
m=size(z,2);

%setup the best 
Jcbb.Best=zeros(1,size(z,2));

for i=1:m
    

















end
