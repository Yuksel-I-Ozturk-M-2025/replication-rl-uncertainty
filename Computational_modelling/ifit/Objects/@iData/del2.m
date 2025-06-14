function a = del2(a)
% b = del2(s) : computes the Discrete Laplacian of iData object
%
%   @iData/del2 function to compute the Discrete Laplacian of data sets.
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=del2(a);
%
% Version: Aug. 22, 2017
% See also iData, iData/gradient, del2, gradient, iData/jacobian

% make sure axes are regularly binned
%a = interp(a);

a = iData_private_unary(a, 'del2');

