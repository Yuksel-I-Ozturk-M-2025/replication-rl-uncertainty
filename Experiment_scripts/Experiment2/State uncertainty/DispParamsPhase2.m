function [rho,r0,wstate,bstate,pwhite,n,nBS,nAS,rate]=DispParamsPhase2(subNo)
%[rho,r0,wstate,bstate,pwhite,n,nBS,nAS,rate]=DispParamsPhase2(1);

% The vector pwhite contains the proportions of white to use for each stimilus.  
% n represents the total number of trials 
% r0 value of positive reward 
% rate is the rate to use to transform points into monetary payments (e.g., in pounds or Australian dolars)

rho=0.7; % Uncertainty level
r0=10; % value of positive reward: 10 points
rate=1/200; % Each 200 points is worth one Australian dollar  

% Read the light/dark stimulus contrasts from the results of the first phase
data=importdata(['StimulusForPhase2_',num2str(subNo),'.txt']);
vect=data.data; % extract data
wstate=vect(2); % The propotion of white for the white state corresponding to the level of confidence rho (to modify according to the results of the 1st phase)
bstate=vect(3); % The propotion of white for the white state corresponding to the level of confidence 1-rho (to modify according to the results of the 1st phase)

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

% The white dot propotions of the stimuli that will be presented in the test phase
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

