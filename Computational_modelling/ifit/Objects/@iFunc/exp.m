function a = exp(a)
% b = exp(s) : exponential value of iFunc object
%
%   @iFunc/exp function to return the exponential value of data sets
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex:     b=exp(a);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/exp, iFunc/log, iFunc/log10, iFunc/sqrt

a = iFunc_private_unary(a, 'exp');

