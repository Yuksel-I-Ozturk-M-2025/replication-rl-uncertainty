function a = asinh(a)
% b = asinh(s) : computes the inverse hyperbolic sine of iFunc object
%
%   @iFunc/asinh function to compute the inverse hyperbolic sine of data sets.
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex:     b=asinh(a);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/cos, iFunc/acos, iFunc/sin, iFunc/asin, iFunc/tan, iFunc/atan

a = iFunc_private_unary(a, 'asinh');

