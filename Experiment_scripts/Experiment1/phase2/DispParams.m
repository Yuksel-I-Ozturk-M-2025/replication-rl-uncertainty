function [rho,wstate,bstate,pwhite,n,nBS,nAS]=DispParams

% The vector pwhite contains all the pourcentages of white to be
%used each stimilus presented.  
% n represents the nb of trials 

rho=0.65; % level of confidence for the white state
wstate=0.5307; % The propotion of white for the white state corresponding to the level of confidence rho 
bstate=0.5150; % The propotion of white for the white state corresponding to the level of confidence 1-rho 
Constants=[wstate,bstate]; % the propotions of white to be used
nConstants=length(Constants); % nb of constants

nBS=100; % nb trials prior to the switching 
nAS=100; % nb trials post-switching
n=nBS+nAS; % total nb of trials 
nPerConstantBS=nBS/2; % nb of trials/constant prior to the switch
nPerConstantAS=nAS/2; % nb of trials/constant post-switch

pwhiteBS=zeros(nBS,1); % this will become our vector containing percentages of white (before switch)
pwhiteAS=zeros(nAS,1); % this will become our vector containing percentages of white (after switch)
pwhite=zeros(n,1); % this will become our vector containing percentages of white for the entire test phase

% the filling of the array of stimulus which will be presented in the test phase

%First step: filling the first vector prior to the switching and then shuffle it
for i=1:nConstants
    for j=1:(nPerConstantBS)            
        pwhiteBS(nPerConstantBS*(i-1)+j)=Constants(i);                
    end
end
pwhiteBS=pwhiteBS(randperm(nBS)); %randomize the rank of presentation of the stimulus

%Second step: filling the second vector post-switch and then shuffle it
for i=1:nConstants
    for j=1:(nPerConstantAS)            
        pwhiteAS(nPerConstantAS*(i-1)+j)=Constants(i);                 
    end
end
pwhiteAS=pwhiteAS(randperm(nAS)); %randomize the rank of presentation of the stimulus

%Third step: concatenation of these two vectors to form pwhite
pwhite=[pwhiteBS;pwhiteAS];

end

