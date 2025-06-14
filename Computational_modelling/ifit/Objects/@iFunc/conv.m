function c = conv(a,b, shape)
% c = conv(a,b) : computes the convolution of iFunc models
%
%   @iFunc/conv function to compute the convolution of models (FFT based).
%     A deconvolution mode is also possible.
%     When used with a single scalar value, it is used as a width to build a 
%       gaussian function.
%     when one of the argument is a character string, it is used as-is in the 
%     operator expression. 
%
%     conv(a, 'tas') convolves the 'a' model with the TAS resolution function
%       computed from ResLibCal. The dispersion is either a 2D model S(|q|,w)
%       or a 4D model S(qh,qk,ql,w). The returned model can then further be used
%       to fit its parameters to match a measurement scan.
%
% input:  a: object or array (iFunc or numeric)
%         b: object or array (iFunc or numeric)
%     shape: optional shape of the return value
%          full         Returns the full two-dimensional convolution.
%          same         Returns the central part of the convolution of the same size as a.
%          valid        Returns only those parts of the convolution that are computed
%                       without the zero-padded edges. Using this option, y has size
%                       [ma-mb+1,na-nb+1] when all(size(a) >= size(b)).
%          deconv       Performs an FFT deconvolution.
%          pad          Pads the 'a' signal by replicating its starting/ending values
%                       in order to minimize the convolution side effects
%          center       Centers the 'b' filter so that convolution does not shift
%                       the 'a' signal.
%          normalize    Normalizes the 'b' filter so that the convolution does not
%                       change the 'a' signal integral.
%          background   Remove the background from the filter 'b' (subtracts the minimal value)
%     Default shape is 'same'
%
% output: c: object or array (iFunc)
% ex:     c=conv(a,b); c=conv(a,b, 'same pad background center normalize');
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/convn, iFunc/xcorr, fconv, fconvn, fxcorr, ResLibCal
if nargin ==1
	b = a;
end
if isnumeric(b) && isscalar(b)
  g = b;
  b = gauss;
  b.Guess = [ 1 0 double(g) 0]; % use input as a width
  c = convn(a,b);
  return
elseif isnumeric(a) && isscalar(a)
  g = a;
  a = gauss;
  a.Guess = [ 1 0 double(g) 0]; % use input as a width
  c = convn(a,b);
  return
elseif ischar(b) && any(strcmpi(b, {'tas','reslibcal','rescal','reslib'}))
  c = ResLibCal(a);
  return
elseif ischar(a) && any(strcmpi(a, {'tas','reslibcal','rescal','reslib'}))
  c = ResLibCal(b);
  return
elseif isstruct(b) && isfield(b,'EXP')
  c = ResLibCal(a);
  return
elseif isstruct(a) && isfield(a,'EXP')
  c = ResLibCal(b);
  return
end
if nargin < 3, shape = 'same'; end

c = iFunc_private_binary(a, b, 'fconv', shape);

