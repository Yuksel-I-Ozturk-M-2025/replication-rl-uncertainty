function c = char(s)
% c = char(s) : convert iData into character
%
%   @iData/char: function to convert iData objects into char
%   returns the iData title/filename
%
% input:  s: object or array (iData) 
% output: c: iData identification (char)
%
% Version: Aug. 22, 2017
% See also  iData/cell, iData/double, iData/struct, 
%           iData/char, iData/size
%

% EF 23/09/07 iData implementation

c=[];
for index=1:numel(s)
  t = s(index);
  T = t.Title;  if ~ischar(T), T=char(T); end
  if size(T, 1)~=1, T=transpose(T); T=T(:)'; end
  T   = regexprep(T,'\s+',' '); % remove duplicated spaces
  cmd = t.Command{end};
  if length(cmd) > 23, cmd = [ cmd(1:20) '...' ]; end
  % form the signal(axes) string
  labels1 = deblank(t.Alias.Labels{1});
  if strcmpi(labels1, 'Data Signal') labels1 = ''; end
  labels2 = deblank(sprintf('%s ', t.Alias.Axis{:}));
  labels=sprintf('%s(%s)', labels1, labels2);
  c = strvcat(c, [ 'iData ' cmd ' [' num2str(size(t)) '] ' labels ' "' strtrim(T) '" <' t.Source '>' ]); 
  c(~isstrprop(c,'print'))='';
end

