function a = sin(a)
% b = sin(s) : computes the sine of iData object
%
%   @iData/acos function to compute the sine of data sets (using radians).
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=sin(a);
%
% Version: Aug. 22, 2017
% See also iData, iData/cos, iData/acos, iData/sin, iData/asin, iData/tan, iData/atan

a = iData_private_unary(a, 'sin');

