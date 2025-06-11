function [N,rho,alpha,T,r0] = DispParams(paramVector)
%[N,rho,alpha,T,r0] = DispParams(paramVector);

N = paramVector(1); % Total nb of trial in the experiment - N=200
rho = paramVector(2); % level of confidence (in perception) - rho=0.7
alpha = paramVector(3); % learning parameter - alpha=0.005
T = paramVector(4); % exploration parameter - T=0.05
r0 = paramVector(5); % maximal reward - r0=10

end
