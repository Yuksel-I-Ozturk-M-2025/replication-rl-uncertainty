function a = acos(a)
% b = acos(s) : computes the arc cosine of iFunc object
%
%   @iFunc/acos function to compute the inverse cosine of data sets (in radians).
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex:     b=acos(a);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/cos, iFunc/acos, iFunc/sin, iFunc/asin, iFunc/tan, iFunc/atan

a = iFunc_private_unary(a, 'acos');

