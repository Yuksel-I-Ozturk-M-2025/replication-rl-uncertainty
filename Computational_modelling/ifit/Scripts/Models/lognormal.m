function y=lognormal(varargin)
% y = lognormal(p, x, [y]) : Log-Normal distribution function. 
%
%   iFunc/lognormal Log-Normal distribution function. 
%     y  = p(4)+ p(1)/sqrt(2)/p(3)./x .* exp( -log(x/p(2)).^2 /2/p(3)/p(3) )
%
% lognormal(centre)         creates a model with a specified centre
% lognormal([ parameters ]) creates a model with specified model parameters
%
% Reference: http://en.wikipedia.org/wiki/Log_normal
%
% input:  p: Log-Normal model parameters (double)
%            p = [ Amplitude Center Width BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=lognormal([1 0 1 1], -10:10); or plot(lognormal);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/fits, iFunc/plot
% (c) E.Farhi, ILL. License: EUPL.

y.Name      = [ 'Log-Normal distribution function (1D) [' mfilename ']' ];
y.Parameters={'Amplitude','Centre','Width','Background'};
y.Description='Log-Normal distribution function. Ref: http://en.wikipedia.org/wiki/Log_normal';
y.Expression= @(p,x) real(p(4)+ p(1)/sqrt(2)/p(3)./x .* exp( -log(x/p(2)).^2 /2/p(3)/p(3) ));
y.Guess     = @(x,y) [ (max(y(:))-min(y(:)))/2 mean(abs(x(:))) std(x(:))/2 min(y(:)) ];

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 varargin{1} varargin{1}/4 0]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

