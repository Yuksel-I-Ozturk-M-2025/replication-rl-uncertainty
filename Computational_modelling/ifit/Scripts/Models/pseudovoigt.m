function y=pseudovoigt(varargin)
% y = pseudovoigt(p, x, [y]) : Pseudo Voigt
%
%   iFunc/pseudovoigt Pseudo Voigt fitting function
%     approximation of the convolution of gauss and lorz
%     	y = a * (d * (1/(1+((x-b)/c)^2)) + (1-d) * exp(-0.5 * ((x-b)/c)^2))
%
% pseudovoigt(width)          creates a model with a specified width
% pseudovoigt([ parameters ]) creates a model with specified model parameters
%
% Reference: http://en.wikipedia.org/wiki/Voigt_profile
%            P. Thompson, D.E. Cox, J.B. Hastings, J. Appl. Cryst. 1987, 20, 79.
%
% input:  p: Pseudo Voigt model parameters (double)
%            p = [ Amplitude Centre HalfWidth BackGround LorentzianRatio ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=pseudovoigt([1 0 1 1], -10:10); or plot(pseudovoigt);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/fits, iFunc/plot
% (c) E.Farhi, ILL. License: EUPL.

y.Name      = [ 'Pseudo-Voigt (1D) [' mfilename ']' ];
y.Parameters={'Amplitude','Centre','HalfWidth','Background','LorentzianRatio'};
y.Description='Single 1D Pseudo Voigt model (convolution of gauss and lorz approx.). Ref: P. Thompson, D.E. Cox, J.B. Hastings, J. Appl. Cryst. 1987, 20, 79.';
y.Expression= @(p,x) p(1) * (p(5) * (1./(1+ ((x-p(2))/p(3)).^2 )) + (1-p(5)) * exp(-0.5 * ((x-p(2))/p(3)).^2 ));

m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));

y.Guess     = @(x,s) [ NaN m1(x, s-min(s(:))) m2(x, s-min(s(:))) NaN 0.5 ];

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 0 varargin{1} 0 0.5]};
  end
  y.ParameterValues = varargin{1};
elseif length(varargin) > 1
  y = y(varargin{:});
end

