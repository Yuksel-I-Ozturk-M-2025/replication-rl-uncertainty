function a = atan(a)
% b = atan(s) : computes the arc tangent of iFunc object
%
%   @iFunc/atan function to compute the inverse tangent of data sets (in radians).
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex:     b=atan(a);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/cos, iFunc/acos, iFunc/sin, iFunc/asin, iFunc/tan, iFunc/atan

a = iFunc_private_unary(a, 'atan');

