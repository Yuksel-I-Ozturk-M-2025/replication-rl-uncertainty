function a = isinf(a)
% b = isinf(s) : True for infinite iData object elements
%
%   @iData/isinf function to return true for infinite elements
%   of 's', i.e. that are +Inf or -Inf 
%
%   To remove nan's and inf's values use: fill(s)
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=isinf(a);
%
% Version: Aug. 22, 2017
% See also iData, iData/sign, iData/isreal, iData/isfinite, iData/isnan,
%          iData/isinf, iData/isfloat, iData/isinterger,
%          iData/isnumeric, iData/islogical, iData/isscalar, 
%          iData/isvector, iData/issparse, iData/fill
a = iData_private_unary(a, 'isinf');

