function a = sparse(a)
% b = sparse(s) : Convert iFunc object storage to sparse
%
%   @iFunc/sparse function to use sparse storage.
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex:     b=sparse(a);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/full, iFunc/sparse

a = iFunc_private_unary(a, 'sparse');

