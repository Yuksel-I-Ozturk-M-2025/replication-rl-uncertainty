function [rho,r0,n,nBS,nAS,stateVect,rate]=DispParams()
%[rho,r0,n,nBS,nAS,stimVect,rate]=DispParams();

% The vector pwhite contains the proportions of white to use for each stimilus.  
% n represents the total number of trials 
% rate is the rate to use to transform points into monetary payments (e.g., in pounds or Australian dolars)

rho=0.7; % Uncertainty level
r0=10; % value of positive reward: 10 points
rate=1/200; % Each 200 points is worth one Australian dollar 

nBS=100; % nb trials prior to the switching 
nAS=100; % nb trials post-switching
n=nBS+nAS; % total nb of trials 

% Constructing the vector containing the stimui states (either 0 for black or 1 for white)
%1st step: construct the first vector prior to the switch then shuffle it
stimBS = zeros(nBS,1); % for white patchs 
stimBS(1:(nBS/2)) = 1; % for black patch
stimBS=stimBS(randperm(nBS)); %randomize the rank of presentation of the stimulus

%2nd step: construct the second vector post-switch then shuffle it
stimAS = zeros(nAS,1); % for white patchs 
stimAS(1:(nAS/2)) = 1; % for black patch
stimAS=stimAS(randperm(nAS)); %randomize the rank of presentation of the stimulus

%Third step: concatenation of these two vectors to form the vector of all
%stimuli
stateVect = [stimBS;stimAS];

end

