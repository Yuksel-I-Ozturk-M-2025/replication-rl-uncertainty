function a = isnan(a)
% b = isnan(s) : True for NaN iData object elements
%
%   @iData/isnan function to return true for NaN elements
%   of 's', i.e. that are NaN ('not a number')
%
%   To remove nan's and inf's values use: fill(s)
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=isnan(a);
%
% Version: Aug. 22, 2017
% See also iData, iData/sign, iData/isreal, iData/isfinite, iData/isnan,
%          iData/isinf, iData/isfloat, iData/isinterger,
%          iData/isnumeric, iData/islogical, iData/isscalar, 
%          iData/isvector, iData/issparse, iData/fill

a = iData_private_unary(a, 'isnan');
