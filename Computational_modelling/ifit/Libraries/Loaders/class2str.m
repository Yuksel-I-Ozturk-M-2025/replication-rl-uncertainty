function str=class2str(this, data, options)
%   This function creates a string containing Matlab code describing a variable.
%
%   class2str(this,data) 
%     Create a string [ 'this = data;' ]
%   class2str(this, data, 'no comments') 
%     same as above, but removes comments from the output file
%   class2str(this, data, 'eval') 
%     creates a compact evaluable string of the initial data
%   class2str(this, data, 'flat') 
%     creates a flat text with commented data blocks, which is not an m-file, 
%     but rather a Linux-style config file.
%
% input arguments:
%   this: string containing the name of the object to describe. Use '' to name the 
%           root level from the input variable name, and ' ' to ignore root level 
%           in structures.
%   data: any data set (struct, array, cell, char, objects)
%   options: optinal argument which may contain 'flat', 'no comments','eval'
%
% output variables:
%   str: string which contains a function code to generate the data.
%
% example: 
%  data=struct('a',1,'b','a string comment','c',{'cell'});
%  str=class2str(data)          % produces a string/script that regenerates data
%  str=class2str('this', data)  % idem, but creates 'this' instead of 'data'.
%  str=class2str(data, 'flat')  % a kind of config file
%  str=class2str(data, 'eval')  % a compact evaluable string
%          
% See also: mat2str, num2str, eval, sprintf
%
% Part of: Loaders utilities (ifit.mccode.org)
% Author:  E. Farhi <farhi@ill.fr>. Aug. 22, 2017
% (c) E.Farhi, ILL. License: BSD.

if     nargin == 0
  str=[]; return
elseif nargin == 1
  data    = this;
  this    = '';  % default root level name
  options = '';
elseif nargin == 2
  if ~ischar(this) && ischar(data) % (data, options)
    options=data; data=this; this=inputname(1);
  else options = ''; end
end
if isempty(this), this = inputname(1); end
if isempty(this), this = [ class(data) '_str' ]; end

if strfind(options, 'flat')
  str = class2str_flat(this, data, options);
elseif strfind(options, 'eval')
  str = class2str_eval(data);
else
  str = class2str_m(this, data, options);
end

return

% ------------------------------------------------------------------------------
function str = class2str_m(this, data, options)
% function to create an m-file string

nocomment = strfind(options, 'no comment');
str       = '';
NL = sprintf('\n');

if ischar(data)
  str = [ this ' = ''' class2str_validstr(data) ''';' NL ];
elseif (all(isobject(data)) || all(isstruct(data))) && numel(data) > 1
  if ~nocomment, str = [ '% ' this ' (' class(data) ') array size ' mat2str(size(data)) NL ]; end
  for index=1:numel(data)
    str = [ str class2str([ this '(' num2str(index) ')' ], data(index), options) ];
  end
  str = [ str this ' = reshape(' this ', [' num2str(size(data)) ']);' NL ];
  if ~nocomment, str = [ str  '% end of ' class(data) ' array ' this NL ]; end
elseif isa(data, 'iData')
  if ~nocomment, str = [ '% ' this ' (' class(data) ') size ' num2str(size(data)) NL ]; end
  str = [ str class2str(this, struct(data), options) ];
  if ~nocomment, str = [ str NL '% handling of iData objects -------------------------------------' NL ]; end
  str = [ str 'if ~exist(''iData''), return; end' NL ];
  str = [ str this '_s=' this '; ' this ' = rmfield(' this ',''Alias''); ' this ' = iData(' this ');' NL ...
         'setalias(' this ', ' this '_s.Alias.Names, ' this '_s.Alias.Values, ' this '_s.Alias.Labels);' NL ... 
         'if ~isempty(' this '_s.Alias.Axis)' NL ...
         '  setaxis('  this ', mat2str(1:length(' this '_s.Alias.Axis)), ' this '_s.Alias.Axis);' NL ...
         'end' NL ];
  if ~nocomment, str = [ str '% end of iData ' this NL ]; end
elseif isobject(data)
  if ~nocomment, str = [ '% ' this ' (' class(data) ') size ' num2str(size(data)) NL ]; end
  str = [ str class2str(this, struct(data), options) ];
  if ~nocomment, str = [ str NL '% handling of objects -------------------------------------' NL ]; end
  str = [ str 'try; ' this ' = ' class(data) '(' this '); end' ];
  if ~nocomment, str = [ str '% end of object ' this NL ]; end
elseif isnumeric(data) | islogical(data)
  if ~nocomment, str = [ '% ' this ' numeric (' class(data) ') size ' num2str(size(data)) NL ]; end
  str = [ str this ' = ' mat2str(data(:)) ';' NL ];
  if numel(data) > 1
    str = [ str this ' = reshape(' this ', [' num2str(size(data)) ']);' NL ];
  end
elseif isstruct(data) && ~isempty(data)
  f = fieldnames(data);
  if ~nocomment, str = [ '% ' this ' (' class(data) ') length ' num2str(length(f)) NL ]; end
  for index=1:length(f)
    if isempty(deblank(this))
      str = [ str class2str([ f{index} ], getfield(data, f{index}), options) ];
    else
      str = [ str class2str([ this '.' f{index} ], getfield(data, f{index}), options) ];
    end
  end
  if ~nocomment, str = [ str '% end of struct ' this NL ]; end
elseif iscellstr(data)
  if ~nocomment, str = [ '% ' this ' (' class(data) 'str) size ' mat2str(size(data)) NL ]; end
  str = [ str this ' = { ...' NL ];
  for index=1:numel(data)
    str = [ str '  ''' class2str_validstr(data{index}) '''' ];
    if index < numel(data), str = [ str ', ' ]; end
    str = [ str ' ...' NL ];
  end
  str = [ str '}; ' NL ];
  if prod(size(data)) > 1
    str = [ str this ' = reshape(' this ', [' mat2str(size(data)) ']);' NL ];
  end
  if ~nocomment, str = [ str '% end of cellstr ' this NL ]; end
elseif iscell(data)
  if ~nocomment, str = [ '% ' this class(data) ' size ' mat2str(size(data)) NL ]; end
  str = [ str this ' = cell(' mat2str(size(data)) ');' NL ];
  for index=1:numel(data)
    str = [ str class2str([ this '{' num2str(index) '}' ], data{index}, options) ];
  end
  if prod(size(data)) > 1
    str = [ str this ' = reshape(' this ', [' mat2str(size(data)) ']);' NL ];
  end
  if ~nocomment, str = [ str '% end of cell ' this NL ]; end
elseif isa(data, 'function_handle')
  if ~nocomment, str = [ '% ' this ' function (' class(data) ')' NL ]; end
  str = [ str this ' = ' func2str(data(:)) ';' NL ];
else
  try
    % other classes
    if ~nocomment, str = [ '% ' this ' (' class(data) ') size ' num2str(size(data)) NL ]; end
    str = [ str class2str(this, struct(data), options) ];
    str = [ str 'try; ' this ' = ' class(data) '(' this '); end' ];
    if ~nocomment, str = [ str '% end of ' class(data) ' ' this NL ]; end
  catch
    warning([ mfilename ': can not save ' this ' ' class(data) '. Skipping.' ]);
  end
end

% ------------------------------------------------------------------------------
function str = class2str_flat(this, data, options)
% function to create a flat file string (no matlab code)

str       = '';
NL  = sprintf('\n');

if isempty(data), return; end
if ischar(data)
  str = [ '# ' this ': ' class2str_validstr(data) NL ];
elseif (isobject(data) || isstruct(data)) && length(data) > 1
  for index=1:numel(data)
    str = [ str class2str([ this '(' num2str(index) ')' ], data(index), options) NL ];
  end
elseif isstruct(data)
  f = fieldnames(data);
  %str = [ '# ' class(data) ' length ' num2str(length(f)) ': ' this NL ];
  str = '';
  for index=1:length(f)
    str = [ str class2str([ this ' ' f{index} ], getfield(data, f{index}), options) ];
  end
elseif iscellstr(data)
  % str = [ '# ' class(data) 'str size ' mat2str(size(data)) ': ' this NL ];
  str = '';
  for index=1:numel(data)
    str = [ str class2str(this,data{index},'flat') ];
  end
elseif iscell(data)
  % str = [ '# ' class(data) ' size ' mat2str(size(data)) ': ' this NL ];
  str = '';
  for index=1:numel(data)
    str = [ str class2str([ this '{' num2str(index) '}' ], data{index}, options) NL ];
  end
elseif isa(data, 'function_handle')
  str = [ '# ' this ' function: ' func2str(data(:)) NL ];
elseif ~isobject(data) && (isnumeric(data) || islogical(data))
  str = num2str(data);
  str(:,end+1) = sprintf('\n');
  str = str';
  str = str(:)';
  str = [ '# ' class(data) ' size ' num2str(size(data)) ': ' this NL ...
          str ];
else
  try
    % other class
    str = [ '# ' class(data) ' size ' num2str(size(data)) ': ' this NL ];
    str = [ str class2str(this, struct(data), options) ];
  catch
    warning([ mfilename ': can not save ' this ' ' class(data) '. Skipping.' ]);
  end
end

% ------------------------------------------------------------------------------
function str = class2str_eval(data)

str = '';
NL  = sprintf('\n');

if ischar(data)
  str = [ '''' class2str_validstr(data) '''' ];
elseif numel(data) == 1 && ishandle(data)
    % ignore
elseif isobject(data) && numel(data) == 1
    str = [ str class(data) '(' class2str_eval(struct(data)) ')' ];
elseif isnumeric(data)
  if ndims(data) <= 2, str = mat2str(data); end
elseif iscell(data)
  if isempty(data), str = [     '{ ' ]; end
  for index=1:numel(data)
    if index == 1, str = [     '{ ' class2str_eval(data{index}) ];
    else           str = [ str ', ' class2str_eval(data{index}) ]; end
  end
  str = [ str ' }' ];
  if size(data,1) > 1 && size(data,2) == 1, str = [ str '''' ]; end  % transpose
elseif isa(data, 'function_handle') && numel(data) == 1
  str = func2str(data);
elseif isstruct(data) && numel(data) == 1
  index=1;
  for f=fieldnames(data)'
    if index==1, str = [ 'struct(''' f{1} ''',' class2str_eval(data.(f{1})) ];
    else         str = [ str  ', ''' f{1} ''', ' class2str_eval(data.(f{1})) ]; end
    index=index+1;
  end
  str = [ str ')' ];
elseif numel(data) > 1
  for index=1:numel(data)
    if index == 1, str = [     '[ ' class2str_eval(data(index)) ];
    else           str = [ str ', ' class2str_eval(data(index)) ]; end
  end
  str = [ str ' ]' ];
  if size(data,1) > 1 && size(data,2) == 1, str = [ str '''' ]; end  % transpose
end

% ------------------------------------------------------------------------------

function str=class2str_validstr(str)
  % validate a string as a single line
  str=strrep(str(:)', sprintf('\n'), ';');
  index = find(str < 32 | str > 127);
  str(index) = ' ';
  str=strrep(str, '''', '''''');

