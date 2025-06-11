function O=Generate_obs(N)
%O=Generate_obs(N);

% generate a sequence of N trials from the set {Light,Dark} 
%which will serve as observations when performing RL

% 1 <-> 'light';
% 0 <-> 'dark';

% N is the nb of trials (= nb of steps in the simulation)

oBS_white=ones(N/4,1);
oBS_black=zeros(N/4,1);
oBS=[oBS_white;oBS_black];
oBS=oBS(randperm(N/2));

oAS_white=ones(N/4,1);
oAS_black=zeros(N/4,1);
oAS=[oAS_white;oAS_black];
oAS=oAS(randperm(N/2));

O=[oBS;oAS];

       

end