function [rho,constants,pwhite,nconstants,ntrialsperconstant,ntrials,sigma_adj]=DispParamsMain(sigma)
%[rho,constants,pwhite,nconstants,ntrialsperconstant,ntrials,sigma_adj]=DispParamsMain(sigma_warm);

% The vector pwhite contains the proportion of white to use for each stimilus.  
% n represents the total number of trials 
% sigma (sd of the psychometric gaussian distribution) is estimated using the (90) pretraining trials 

rho=0.7; % Uncertainty level
mu=0.5; % mean of the fitted psychometric function
sigma_adj = max(sigma, 0.015); % adjusted sigma, which is always greater than 0.015 (to make sure to not use a very small sigma, and hence make the contrasts very difficult for the subjects)
constants=[mu-2*sigma_adj, mu-sigma_adj, mu-0.75*sigma_adj, mu-0.5*sigma_adj, mu-0.25*sigma_adj, mu+0.25*sigma_adj, mu+0.5*sigma_adj, mu+0.75*sigma_adj, mu+sigma_adj, mu+2*sigma_adj]; % the propotions of white to be used
nconstants=length(constants); % no. of constants

%%%% Main trials (10*25 = 250 trials)
ntrialsperconstant=25; % nb of trial per constant
ntrials=nconstants*ntrialsperconstant; % length = (nb of constants) x (nb of trials per constant)
pwhite=zeros(ntrials,1); % this will become our vector containing percentages of white

% Vector defining the stimuli that will be presented in the prelinimary phase (main trials)
for i=1:nconstants
    for j=1:(ntrialsperconstant)            
        pwhite(ntrialsperconstant*(i-1)+j)=constants(i);                
    end
end
pwhite=pwhite(randperm(ntrials)); %randomize the rank of presentation of the stimulus

end