function [N,Ns,r0,rho,w,pr,b,T] = DispParams(paramVector)
%[N,Ns,r0,rho,w,pr,b,T] = DispParams(paramVector);

N = paramVector(1); % Total nb of trial in the experiment - N = 200
Ns = paramVector(2); % nb of simulations to use to estimate the action likelihoods - Ns = 1000
r0 = paramVector(3); % reward of correct action - r0 = 2.5
rho = paramVector(4); % level of confidence (in perception) - rho = 0.65
w = paramVector(5); % weight placed on the grand average payoff of all past experiences - w = 0.8
pr = paramVector(6); % probability of choosing the most recent trial (t-1) when we sample the 'b' recent trials - pr = 0.8
b = paramVector(7); % bound of the memory - b = 6
T = paramVector(8); % exploration rate - T = 1

end
