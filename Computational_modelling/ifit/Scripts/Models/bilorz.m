function y=bilorz(varargin)
% y = bilorz(p, x, [y]) : Bi-Lorentzian
%
%   iFunc/bilorz Bi-Lorentzian fitting function
%     y = p(1)*exp(-0.5*((x-p(2))/s).^2) + p(5);
%   where s = p(3) for x < p(2) and s = p(4) for x > p(2).
%
% bilorz([ w1 w2])       creates a model with specified widths
% bilorz([ parameters ]) creates a model with specified model parameters
%
% input:  p: Bi-Lorentzian model parameters (double)
%            p = [ Amplitude Centre HalfWidth1 HalfWidth2 BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=bilorz([1 0 1 1], -10:10); or plot(bilorz)
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/fits, iFunc/plot
% (c) E.Farhi, ILL. License: EUPL.

y.Name      = [ 'Bi-Lorentzian (1D) [' mfilename ']' ];
y.Description='Bi-Lorentzian/asymmetric fitting function';
y.Parameters= { 'Amplitude', 'Centre', 'HalfWidth1', 'HalfWidth2', 'BackGround' };
y.Expression = @(p,x) p(1)*exp(-0.5*((x-p(2))./ (p(3)*(x < p(2)) + p(4) * (x >= p(2))) ).^2) + p(5);
% moments of distributions
m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));

y.Guess     = @(x,s) [ NaN m1(x, s-min(s(:))) m2(x, s-min(s(:)))/1.5 m2(x, s-min(s(:)))*1.5 NaN ];

y=iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 0 varargin{1}*.75 varargin{1}*1.25 0]};
  elseif length(varargin{1}) == 2
    varargin = {[ 1 0 varargin{1} 0]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end
