function c = convn(a,b)
% c = convn(a,b) : computes the convolution of an iData object with a response function 
%
%   @iData/convn function to compute the convolution of data sets with automatic centering
%     and normalization of the filter. This is a shortcut for
%       conv(a,b, 'same pad background center normalize')
%     When used with a single scalar value, it is used as a width to build a 
%       gaussian function, with same width along all dimensions
%     When used with a vector of same length as the object dimension, a nD
%       gaussian function with width as vector elements along each diemsions
%
% input:  a: object or array, signal (iData or numeric)
%         b: object or array, filter (iData or numeric)
% output: c: object or array (iData)
% ex:     c=convn(a,b);
%
% Version: Aug. 22, 2017
% See also iData, iData/times, iData/conv, iData/fft, iData/xcorr, fconv, fconvn, fxcorr
if nargin ==1
	b=a;
end
if isscalar(b)
  b = [ 1 mean(getaxis(a,1)) double(b) 0]; % use input as a width
  b = gauss(b, getaxis(a,1));
elseif isscalar(a)
  a = [ 1 mean(getaxis(a,1)) double(a) 0]; % use input as a width
  a = gauss(a, getaxis(b,1));
elseif isa(b,'double') && numel(b) == ndims(a)
  p = [];
  g = gauss^(ndims(a));
  ax = {};
  for index=1:ndims(a)
    p = [ p 1 mean(getaxis(a,1)) double(b(index)) 0 ];
    ax{end+1} = getaxis(a, index);
  end
  g = g(p, ax{:});
  b = g;
end
c = iData_private_binary(a, b, 'conv', 'same pad background center normalize');


