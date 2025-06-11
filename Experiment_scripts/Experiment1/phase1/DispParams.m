function [rho,mu,sigma,constants,pwhite,nconstants,ntrialsperconstant,ntrials]=DispParams

% The vector pwhite contains all the pourcentages of white to be
%used each stimilus presented.  
% ntials represents the nb of trials 

rho=0.65; % Level of confidence
mu=0.5; % the mean of the fitting psychometrical function (modifiable!)
sigma=0.04; % its variance (modifiable!)
constants=[mu-2.5*sigma,mu-2*sigma, mu-sigma, mu-0.5*sigma, mu, mu+0.5*sigma, mu+sigma, mu+2*sigma,mu+2.5*sigma]; % the propotion of white to be used (modifiable!)
nconstants=length(constants); % nb of constants
ntrialsperconstant=25; % nb of trial/constant (to be chosen!)
ntrials=nconstants*ntrialsperconstant; % the lenght= (nb of constants) x (nb of trials per constant)
pwhite=zeros(ntrials,1); % this will become our vector containing percentages of white

% the filling of the array of stimulus which will be presented in the
%study phase
for i=1:nconstants
    for j=1:(ntrialsperconstant)            
        pwhite(ntrialsperconstant*(i-1)+j)=constants(i);                
    end
end
pwhite=pwhite(randperm(ntrials)); %randomize the rank of presentation of the stimulus