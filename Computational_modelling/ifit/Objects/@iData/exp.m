function a = exp(a)
% b = exp(s) : exponential value of iData object
%
%   @iData/exp function to return the exponential value of data sets
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=exp(a);
%
% Version: Aug. 22, 2017
% See also iData, iData/exp, iData/log, iData/log10, iData/sqrt

a = iData_private_unary(a, 'exp');

