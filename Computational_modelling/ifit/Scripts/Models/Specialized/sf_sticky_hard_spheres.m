function y=sf_sticky_hard_spheres(varargin)
% y = sf_sticky_hard_spheres(p, x, [y]) : Sticky Hard Sphere structure factor [Baxter/Menon]
%
%   iFunc/sf_sticky_hard_spheres Sticky Hard Sphere structure factor 
%          (for e.g. spheres in dielectric liquids)
%     The 'x' wave-vector/momentum axis is usually in nm-1 or Angs-1.
%     Typical values for parameters are R=50 Angs, rho=0.04, tau=0.15
%     The Sticky Hard Sphere model converges to the Hard Sphere model [PY] with
%       increasing tau (e.g. tau > 10).
%     The model returns the S(q) structure factor.
%
%     See: S.V.G. Menon, C. Manobar, and K. Srinivasa Rao. J. Chem. Phys., 95 (1991) 9186
%          Extracted from sasfit/sasfit_sq/sasfit_sq_StickyHardSphere.c
%     I. Bressler, et al, Journal of Applied Crystallography, 2015, 48 (5), 1587-1598
%
% input:  p: sticky hard sphere model parameters (double)
%            p = [ R=hard_sphere_radius rho=Volume_Fraction tau=stickiness ]
%          or 'guess'
%         x: wave-vector/momentum axis (double, e.g. nm-1 or Angs-1)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=sf_sticky_hard_spheres([10 0.2 0.15], 0:0.01:1); or plot(sf_sticky_hard_spheres,[10 0.2 0.15], 0:0.01:1)
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/fits, iFunc/plot
% (c) E.Farhi, ILL. License: EUPL.

y.Name      = [ 'Sticky Hard Sphere S(q) (1D) [' mfilename ']' ];
y.Description='Sticky Hard Sphere scattering structure factor [Baxter/Menon]';
y.Parameters={'R hard sphere radius [1/x]', ...
              'rho hard sphere volume fraction', ...
              'tau stickiness'};
y.Expression= { ...
  'RHS=abs(p(1)); fp=max(0, min(p(2), 1)); tau=abs(p(3));q=abs(x);' ...
  'if (p(2) == 0.0) signal=ones(size(x)); else ' ...
  'kappa = 2*q*RHS; ' ...
	'epsi = tau+(fp/(1.0-fp)); ' ...
	'gama = fp*(1+fp/2.0)/(3.0*power(1.0-fp,2.0)); ' ...
	'lamb = 6.0/fp*(epsi-sqrt(epsi*epsi-gama)); ' ...
	'mu = lamb*fp*(1-fp); ' ...
	'beta = -(3*fp*power(2+fp,2.0)-2.*mu*(1+7*fp+fp*fp)+mu*mu*(2+fp)) / (2.*power(1.-fp,4)); ' ...
	'alpha = power(1+2.*fp-mu,2)/power(1-fp,4); ' ...
	'CQ =    alpha*power(kappa,3).*(sin(kappa)-kappa.*cos(kappa)); ' ...
	'CQ = CQ+beta*power(kappa,2.).*(2.*kappa.*sin(kappa)-(power(kappa,2.)-2.).*cos(kappa)-2.); ' ...
	'CQ = CQ+0.5*fp*alpha*((4*power(kappa,3)-24*kappa).*sin(kappa)-(power(kappa,4)-12*power(kappa,2.)+24).*cos(kappa)+24); ' ...
	'CQ = -24.*fp*power(kappa,-6).*CQ; ' ...
	'CQ = CQ-2*power(fp*lamb,2.)*(1.-cos(kappa)).*power(kappa,-2.) + 2.*fp*lamb./kappa.*sin(kappa); ' ...
	'signal = 1 ./ (1-CQ); end'
};
y.Guess     = @(x,signal) [ pi/sum(signal(:).*x(:))*sum(signal(:)) max(signal(:)-1) 0.15 ];
y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

