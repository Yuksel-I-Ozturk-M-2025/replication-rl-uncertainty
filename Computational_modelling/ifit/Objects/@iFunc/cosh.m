function a = cosh(a)
% b = cosh(s) : computes the hyperbolic cosine of iFunc object
%
%   @iFunc/cosh function to compute the hyperbolic cosine of data sets.
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex:     b=cosh(a);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/cos, iFunc/acos, iFunc/sin, iFunc/asin, iFunc/tan, iFunc/atan

a = iFunc_private_unary(a, 'cosh');

