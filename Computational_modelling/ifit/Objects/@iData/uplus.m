function a = uplus(a)
% b = uplus(s) : makes a copy of iData object
%
%   @iData/uplus function to return a duplicate of data sets.
%   b=+a creates a new iData object with same content as 'a', but different Tag/ID and Date.
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=uplus(a);
%
% Version: Aug. 22, 2017
% See also iData, iData/uminus, iData/abs, iData/real, iData/imag, iData/uplus

a = iData_private_unary(a, 'uplus');

