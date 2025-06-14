function y=twoexp(varargin)
% signal = twoexp(p, x, {signal}) : two exponential decay functions
%
%   iFunc/twoexp two exponential decay functions (fit 1D function/model)
%     y=p(1)*exp(-x/p(2))+p(3)*exp(-x/p(4)) + p(5);
%
% twoexp([tau1 tau2])    creates a model with specified decay constants
% twoexp([ parameters ]) creates a model with specified model parameters
%
% input:  p: twoexp model parameters (double)
%            p = [   'Amplitude1' 'Tau1' 'Amplitude2' 'Tau2' 'Background' ] as a numerical array
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     signal=twoexp([1 0 1 1], -10:10); or plot(twoexp)
%
% Version: Aug. 22, 2017
% See also iData, iFunc/fits, iFunc/plot, expon, sinedamp
% (c) E.Farhi, ILL. License: EUPL.

y.Name           = [ 'Bi-Exponential decay (1D) [' mfilename ']' ];
y.Parameters     = {'Amplitude1' 'Tau1 decay in inverse "x" unit' 'Amplitude2' 'Tau2 decay in inverse "x" unit' 'Background'};
y.Dimension      = 1;         % dimensionality of input space (axes) and result
y.Description    = '2 Exponential decay';
y.Expression     = @(p,x) p(1)*exp(-x/p(2))+p(3)*exp(-x/p(4)) + p(5);
y.Guess          = @(x,y)[ ...
   exp(subsref(polyfit(x(:),log(y(:)-min(y(:))+0.01*abs(min(y(:)))),1), struct('type','()', 'subs',{{2}})))*0.66 ...
   -1/(subsref(polyfit(x(:),log(y(:)-min(y(:))+0.01*abs(min(y(:)))),1), struct('type','()', 'subs',{{1}}))-(abs(subsref(polyfit(x(:),log(y(:)-min(y(:))+0.01*abs(min(y(:)))),1), struct('type','()', 'subs',{{1}}))) < 1e-2)*.1) ...
   exp(subsref(polyfit(x(:),log(y(:)-min(y(:))+0.01*abs(min(y(:)))),1), struct('type','()', 'subs',{{2}})))/3 ...
   -2/(subsref(polyfit(x(:),log(y(:)-min(y(:))+0.01*abs(min(y(:)))),1), struct('type','()', 'subs',{{1}}))-(abs(subsref(polyfit(x(:),log(y(:)-min(y(:))+0.01*abs(min(y(:)))),1), struct('type','()', 'subs',{{1}}))) < 1e-2)*.2) ...
    min(y(:)) ];
y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    tau = varargin{1};
    varargin = {[ 1 tau 1 tau 0]};
  elseif length(varargin{1}) == 2
    tau = varargin{1};
    varargin = {[ 1 tau(1) 1 tau(2) 0]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

