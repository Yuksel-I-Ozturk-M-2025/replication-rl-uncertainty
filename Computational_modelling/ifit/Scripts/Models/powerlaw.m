function y=powerlaw(varargin)
% y = powerlaw(p, x, [y]) : power law model
%
%   iFunc/powerlaw power law fitting function
%     y = p(1)*(x - p(2))^p(3) + p(4);
%
% poisson(centre)         creates a model with a specified centre
% poisson([ parameters ]) creates a model with specified model parameters
%
% input:  p: power law model parameters (double)
%            p = [ Amplitude Centre Exponent BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=powerlaw([1 0 1 1], -10:10); or plot(powerlaw);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/fits, iFunc/plot
% (c) E.Farhi, ILL. License: EUPL.

y.Name      = [ 'Power law (1D) [' mfilename ']' ];
y.Parameters={'Amplitude','Centre','Exponent','Background'};
y.Description='power law';
y.Expression= @(p,x) real(p(1)*(x - p(2)).^p(3) + p(4));
% fill guessed information
  % ln y = ln a + c*ln (x-b)
y.Guess     = @(x,y) [ ...
  exp(subsref(strline('guess', log(x(:)),log(y(:))), struct('type','()', 'subs',{{1}}))) ...
  mean(x(:)) ...
  subsref(strline('guess', log(x(:)),log(y(:))), struct('type','()', 'subs',{{2}})) ...
  min(y(:)) ];

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 varargin{1} 2 0]};
  end
  y.ParameterValues = varargin{1};
elseif length(varargin) > 1
  y = y(varargin{:});
end

