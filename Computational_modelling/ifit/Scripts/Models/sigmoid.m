function y=sigmoid(varargin)
% y = sigmoid(p, x, [y]) : Sigmoidal
%
%   iFunc/sigmoid Sigmoidal S-shaped
%     y  = A0+(A1-A0)./(1+exp(-(x-x0)/w))
%   This is a sigmoid S-shaped curve, aka logistic.
%
% sigmoid(threshold)      creates a model with specified threshold
% sigmoid([ parameters ]) creates a model with specified model parameters
%
% Ref: http://en.wikipedia.org/wiki/Sigmoid_function
%
% input:  p: Sigmoidal model parameters (double)
%            p = [ Amplitude Center Width BackGround ]
%          or action e.g. 'identify', 'guess', 'plot' (char)
%         x: axis (double)
%         y: when values are given, a guess of the parameters is performed (double)
% output: y: model value or information structure (guess, identify)
% ex:     y=sigmoid([1 0 1 1], -10:10); or plot(sigmoid);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/fits, iFunc/plot, doseresp
% (c) E.Farhi, ILL. License: EUPL.

y.Name           = [ 'Sigmoidal curve (1D) [' mfilename ']' ];
y.Parameters     = {'Amplitude','Center','Slope','Background'};
y.Description    = 'Sigmoidal curve. Ref: http://en.wikipedia.org/wiki/Sigmoid_function';
y.Expression     = @(p,x) p(4)+p(1)./(1+exp(-(x-p(2))/p(3)));
y.Dimension      = 1;
y.Guess          =  @(x,y) [ max(y(:))-min(y(:)) mean(x(:)) std(x(:))/3 min(y(:)) ];
y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 varargin{1} varargin{1}/4 0]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

