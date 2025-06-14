function a=double(a)
% d = double(s) : convert iData into doubles
%
%   @iData/double function to convert iData object into doubles
%
% input:  s: object or array (iData)
% output: v: value of the iData Signal/Monitor (double)
% ex:     'double(iData(rand(10)))'
%
% Version: Aug. 22, 2017
% See also  iData/cell, iData/single, iData/struct, 
%           iData/char, iData/size

% EF 11/07/00 creation
% EF 23/09/07 iData implementation

if numel(a) > 1
  b = cell(size(a));
  for index=1:numel(a)
    b{index} = double(a(index));
  end
  a = b;
  return
end

a = getaxis(a, 'Signal');
a = double(a);
