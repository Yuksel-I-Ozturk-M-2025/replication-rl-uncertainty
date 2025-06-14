function y=tophat(varargin)
% y = tophat(p, x, [y]) : Top-Hat rectangular function
%
%   iFunc/tophat Top-Hat rectangle fitting function
%     y=0*x+p(4); y(find(p(2)-p(3) < x & x < p(2)+p(3))) = p(1);
%     and y is set to the background outside the full width.
%
% tophat(width)          creates a model with a specified width
% tophat([ parameters ]) creates a model with specified model parameters
%
% input:  p: Top-Hat rectangular model parameters (double)
%            p = [ Amplitude Centre HalfWidth BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=tophat([1 0 1 1], -10:10); or plot(tophat);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/fits, iFunc/plot, heavisde, triangl
% (c) E.Farhi, ILL. License: EUPL.

y.Name      = [ 'Top-Hat rectangular function (1D) [' mfilename ']' ];
y.Parameters={'Amplitude','Centre','HalfWidth','Background'};
y.Description='Top-Hat rectangular function';
y.Expression= @(p,x) p(1)*(p(2)-p(3) < x & x < p(2)+p(3))+p(4);

m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));

y.Guess     = @(x,s) [ NaN m1(x, s-min(s(:))) m2(x, s-min(s(:))) NaN ];

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 0 varargin{1} 0]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

