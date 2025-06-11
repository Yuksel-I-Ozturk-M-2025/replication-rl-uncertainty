function [rho,constants_warm,pwhite_warm,nconstants_warm,ntrialsperconstant_warm,ntrials_warm]=DispParamsPreTraining
%[rho,constants_warm,pwhite_warm,nconstants_warm,ntrialsperconstant_warm,ntrials_warm]=DispParamsPreTraining;

% The vector pwhite contains the proportion of white to use for each stimilus.  
% n represents the total number of trials 

rho=0.7; % Uncertainty level
mu=0.5; % mean of the fitted psychometric function (Gaussian)
sigma=0.04; % its standard deviation
constants_warm=[mu-2*sigma, mu-sigma, mu-0.75*sigma, mu-0.5*sigma, mu-0.25*sigma, mu+0.25*sigma, mu+0.5*sigma, mu+0.75*sigma, mu+sigma, mu+2*sigma]; % the propotions of white to be used 
nconstants_warm=length(constants_warm); % nb of constants

%%%% Warm-up trials (10*10 = 100 trials)
ntrialsperconstant_warm=10; % nb of trial per constant
ntrials_warm=nconstants_warm*ntrialsperconstant_warm; % length = (nb of constants) x (nb of trials per constant)
pwhite_warm=zeros(ntrials_warm,1); % this will become our vector containing percentages of white

% Vector defining the stimuli that will be presented in the prelinimary phase (warm-up trials)
for i=1:nconstants_warm
    for j=1:(ntrialsperconstant_warm)            
        pwhite_warm(ntrialsperconstant_warm*(i-1)+j)=constants_warm(i);                
    end
end
pwhite_warm=pwhite_warm(randperm(ntrials_warm)); %randomize the rank of presentation of the stimulus


end