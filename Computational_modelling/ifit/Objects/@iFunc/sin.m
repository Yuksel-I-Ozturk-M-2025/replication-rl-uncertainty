function a = sin(a)
% b = sin(s) : computes the sine of iFunc object
%
%   @iFunc/acos function to compute the sine of data sets (using radians).
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex:     b=sin(a);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/cos, iFunc/acos, iFunc/sin, iFunc/asin, iFunc/tan, iFunc/atan

a = iFunc_private_unary(a, 'sin');

