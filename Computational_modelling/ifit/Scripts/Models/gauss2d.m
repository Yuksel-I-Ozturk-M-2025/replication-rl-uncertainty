function signal=gauss2d(varargin)
% signal = gauss2d(p, x, y, {signal}) : 2D Gaussian function
%
%   iFunc/gauss2d 2D Gaussian function (fit 2D function/model)
%     x0=p(2); y0=p(3); sx=p(4); sy=p(5); theta=p(6) [given in deg]
%     a = cos(theta)^2/2/sx/sx + sin(theta)^2/2/sy/sy;
%     b =-sin(2*theta)/4/sx/sx + sin(2*theta)/4/sy/sy;
%     c = sin(theta)^2/2/sx/sx + cos(theta)^2/2/sy/sy;
%     signal = p(1)*exp(-(a*(x-x0).^2+2*b*(x-x0).*(y-y0)+c*(y-y0).^2)) + p(7);
%
%   The HalfWidth parameters are the Gaussian square root variances (Sigma). 
%   The 'true' half widths are thus 1.177*HalfWidth.
%
% gauss2d([w1 w2])        creates a model with a specified widths
% gauss2d([ parameters ]) creates a model with specified model parameters
%
% Reference: http://en.wikipedia.org/wiki/Gaussian_function
%
% input:  p: gauss2d model parameters (double)
%            p = [  'Amplitude' 'Centre_X' 'Center_Y'
%                   'HalfWidth_X' 'HalfWidth_Y' 'Angle' 'Background' ] 
%            the rotation angle is given in degrees.
%          or 'guess'
%         x: axis along rows    (double)
%         y: axis along columns (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value
% ex:     signal=gauss2d([1 2 .5 .2 .3 30 .2], -2:.1:2, -3:.1:3); or plot(gauss2d)
%
% Version: Aug. 22, 2017
% See also iData, iFunc/fits, iFunc/plot, gauss
% (c) E.Farhi, ILL. License: EUPL.

signal.Name           = [ 'Gaussian-2D function with tilt angle (2D) [' mfilename ']' ];
signal.Description    = '2D Gaussian function with tilt angle. http://en.wikipedia.org/wiki/Gaussian_function';
signal.Parameters     = {  'Amplitude' 'Centre_X' 'Center_Y' 'HalfWidth_X' 'HalfWidth_Y' 'Angle tilt [deg]' 'Background' };
signal.Dimension      = 2;         % dimensionality of input space (axes) and result
m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));
signal.Guess          = @(x,y,signal)[ max(signal(:))-min(signal(:)) m1(x,signal) m1(y, signal) m2(x,signal) m2(y,signal) 20*randn min(signal(:)) ];        % default parameters
signal.Expression     = {'x0=p(2); y0=p(3); sx=p(4); sy=p(5);', ...
  'theta = p(6)*pi/180;', ...
  'aa = cos(theta)^2/2/sx/sx + sin(theta)^2/2/sy/sy;', ...
  'bb =-sin(2*theta)/4/sx/sx + sin(2*theta)/4/sy/sy;', ...
  'cc = sin(theta)^2/2/sx/sx + cos(theta)^2/2/sy/sy;', ...
  'signal = p(1)*exp(-(aa*(x-x0).^2+2*bb*(x-x0).*(y-y0)+cc*(y-y0).^2)) + p(7);' };

signal=iFunc(signal);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 0 0 varargin{1} varargin{1} 20*randn 0]};
  elseif length(varargin{1}) == 2
    varargin = {[ 1 0 0 varargin{:} 20*randn 0]};
  end
  signal.ParameterValues = varargin{1};
elseif nargin > 1
  signal = signal(varargin{:});
end

