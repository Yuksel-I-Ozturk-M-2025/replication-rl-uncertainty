function a = sort(a,dim,mode)
% s = sort(a,dim) : Sort iData objects axes in ascending or descending order
%
%   @iData/sort function to sort the data set on its axes
%     sort(a,dim) sorts along axis of rank dim. 
%       If dim=0, sorting is done on all axes.
%     sort(a,dim,mode) where mode='ascend' or 'descend' select sorting order
%
% input:  a: object or array (iData)
%         dim: dimension to sort (int)
%         mode: sorting order 'ascend' or 'descend'
% output: s: sorted data (iData)
% ex:     c=sort(iData(peaks));
%
% Version: Aug. 22, 2017
% See also iData, iData/plus, iData/sort, iData/unique
if ~isa(a, 'iData')
  iData_private_error(mfilename,['syntax is sort(iData, dim, mode)']);
end

if nargin < 2, dim=1; end
if nargin < 3, mode='ascend'; end

% handle input iData arrays
if numel(a) > 1
  s = zeros(iData, numel(a), 1);
  for index=1:numel(a)
    s(index) = sort(a(index), dim, mode);
    a(index)=iData;
  end
  s = reshape(s, size(a));
  return
end
cmd=a.Command;
a = copyobj(a);

sd = subsref(a,struct('type','.','subs','Signal'));
[dummy, sl] = getaxis(a, '0');  % signal definition/label
se = subsref(a,struct('type','.','subs','Error'));
sm = subsref(a,struct('type','.','subs','Monitor'));
if numel(se) > 1 && all(se(:) == se(1)), se=se(1); end
if numel(sm) > 1 && all(sm(:) == sm(1)), sm=sm(1); end

if dim > 0
  tosort=dim;
else
  tosort=1:ndims(a)
end
was_sorted=0;
myisvector = @(c)length(c) == numel(c);

for index=tosort
  x = getaxis(a, index);
  [x, sorti] = sort(x, index, mode);
  if ~isequal(sorti, 1:size(a, index)) && ~all(sorti == sorti(1))
    S.type = '()';
    S.subs = {};
    nx = ndims(x);
    if myisvector(x), nx=1; end
    if ndims(a) == nx
      S.subs={ sorti };
    else
      for j=1:ndims(a), 
        if j ~= index, S.subs{j}=':';
        else           S.subs{j}=sorti; end
      end
    end

    sd =subsref(sd, S);

    try
      se =subsref(se, S);
    catch
      se=[];
    end
    try
      sm =subsref(sm, S);
    catch
      sm=[];
    end
    setaxis(a, index, x);
    was_sorted=1;
  end
end
if was_sorted
  a = setalias(a, 'Signal', sd, [ 'sort(' sl ')' ]); clear sd
  a = setalias(a, 'Error',  se);
  a = setalias(a, 'Monitor',sm);
  a.Command=cmd;
  a = iData_private_history(a, mfilename, a, dim, mode);
end

