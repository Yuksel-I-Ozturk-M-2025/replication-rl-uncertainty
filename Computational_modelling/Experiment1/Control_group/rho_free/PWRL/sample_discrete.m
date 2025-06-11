function sample = sample_discrete(prob, type)
%
% Generate random samples from the discrete distribution.
% Uses a property of the Gumbel distribution to avoid for loops for speed.
%
% INPUTS:
%
% prob [Q N]
%   if type=='prob', then
%     prob(:,n) is the probability vector for the n-th sample, i.e.,
%     Prob(sample(n)==q) = prob(q,n) / sum(prob(:,n))
%     prob cannot have negative entries but doesn't need to be normalized.
%   if type=='energy', then
%     prob(:,n) is specified in energy (negative log-prob) terms, i.e.,
%     Prob(sample(n)==q) = exp(-prob(q,n)) / sum(exp(-prob(:,n)))
% type (char array)
%     Can be either 'prob' [default] or 'energy'
%
% OUTPUTS:
%
% sample [1 N]
%   Random sample from the discrete distribution. Its elements take values
%   in the set {1,..., Q}.
%
% EXAMPLES:
%
% (1) Generate 1000 IID random samples from the discrete density on 4
%     variables with symbol probabilities prob = [.4 .2 .1 .3]
%
%     Q = 4; N = 1000;
%     prob = [.4 .2 .1 .3]';
%     sample = sample_discrete(repmat(prob,[1 N]));
%     f = zeros(Q,1);
%     for q=1:Q
%       f(q) = sum(sample==q)/N;
%     end
%     figure(1);
%     bar([prob f],'grouped');
%     legend('target','histogram');
% 
% (2) Generate 1000 independent random samples from discrete densities on 4
%     variables with symbol probabilities that differ for each site and are
%     specified in terms of energies.
%
%     Q = 4; N = 1000;
%     prob = [.4 .2 .1 .3]';
%     energy = bsxfun(@plus, -log(prob), randn(Q,N));
%     sample = sample_discrete(energy,'energy');
%     f = zeros(Q,1);
%     for q=1:Q
%       f(q) = sum(sample==q)/N;
%     end
%     figure(2);
%     bar([prob f],'grouped');
%     legend('target','histogram');
% 
% For the underlying method, see the appendix of the paper:
% G. Papandreou and A. Yuille, Perturb-and-MAP Random Fields, ICCV-11.
% http://www.stat.ucla.edu/~gpapan/pubs/confr/PapandreouYuille_PerturbAndMap_supmat-iccv11.pdf
%
% Author: George Papandreou, UCLA
% Date  : December 29, 2012
%

[Q N] = size(prob);

if nargin==1, type='prob'; end

switch type
  case 'prob'
    if any(prob<0), error('sample_discrete: prob cannot be negative'); end
    energy = -log(prob+eps);
  case 'energy'
    energy = prob;
  otherwise
    error('sample_discrete: invalid type');
end

u = rand([Q N]);       % u ~ Uniform((0 1])
gumbel = log(-log(u)); % gumbel ~ Gumbel distribution
energy = energy + gumbel; % Perturbed energy

% Discrete sample is the argmin of the perturbed energy
[~,sample] = min(energy,[],1);
