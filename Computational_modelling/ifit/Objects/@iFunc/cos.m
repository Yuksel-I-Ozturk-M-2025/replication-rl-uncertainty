function a = acos(a)
% b = cos(s) : computes the cosine of iFunc object
%
%   @iFunc/cos function to compute the cosine of data sets (using radians).
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex:     b=cos(a);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/cos, iFunc/acos, iFunc/sin, iFunc/asin, iFunc/tan, iFunc/atan

a = iFunc_private_unary(a, 'cos');

