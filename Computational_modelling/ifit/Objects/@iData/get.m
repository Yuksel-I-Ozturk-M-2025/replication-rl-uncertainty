function [varargout] = get(this,varargin)
% [...] = get(s, 'PropertyName', ...) : get iData object properties
%
%   @iData/get function to get iData properties.
%   get(s) displays all property names and their current values for
%     the iData object 's'.
%   get(s,'PropertyName',...) returns only particular properties.
%     the PropertyName may also be an object Alias or an Axis.
%   Input 's' can be a single iData or a iData array.
%
%   To get the Signal/Monitor, rather use getaxis(s, 'Signal') or s{0}
%              Error/Monitor              getaxis(s, 'Error')
%
% input:  s: object or array (iData)
%         PropertyName: name of Property to search (char)
% output: property: property value in 's' (cell)
% ex :    get(iData) or get(iData,'Title')
%
% Version: Aug. 22, 2017
% See also iData, iData/set, iData/getalias, iData/getaxis, iData/findobj

% EF 27/07/00 creation
% EF 23/09/07 iData implementation
% ============================================================================
% calls: subsref (mainly): get is a wrapper to subsref

persistent fields
persistent method

if isempty(fields), 
  fields=fieldnames(iData);
  % remove: Warning: Could not find an exact (case-sensitive) match for 'TITLE'.
  warning('off','MATLAB:dispatcher:InexactCaseMatch');   
end
if isempty(method), method=methods('iData'); end

% handle array of objects
varargout = {};

if ~isa(this,'iData'), return; end

if numel(this) > 1
  varg = cell(1, numel(this));
  if numel(varargin) == 1 && numel(varargin{1}) == numel(this) && iscell(varargin{1})
      vargin = varargin{1};
      for index=1:numel(this)
        try
          varg{index} = get(this(index), vargin{index});
        catch
          varg{index} = [];
        end
      end
  else
      for index=1:numel(this)
        try
          varg{index} = get(this(index), varargin{:});
        catch
          varg{index} = [];
        end
      end
  end

  sz =  cellfun('prodofsize',varg);
  if all(sz == sz(1)) && all(cellfun(@isscalar,varg))
    try
      varg = cell2mat(varg);
    end
  end

  varargout{1} = varg;
  return
end

if nargin == 1
  disp(this, inputname(1));
  varargout{1} = display(this, inputname(1));
  return
end

if length(varargin) == 1 && iscell(varargin{1})
  varargin = varargin{1};
end
out = {};

% handle single object
for index=1:length(varargin)
  property = varargin{index}; % get PropertyName
  if isempty(property), continue; end
  if iscellstr(property), property = char(property{1}); end
  if ~ischar(property)
    iData_private_error(mfilename, [ 'PropertyName should be char strings in object ' inputname(1) ' ' this.Tag ' "' this.Title '" and not ' class(property) ]);
  end
  % test if this is a unique property, or a composed one
  if isvalid(property) || isvarname(property) % extract iData field/alias
    if any(strcmp(property, fields))
      b = this.(property);               % direct static field
      if isnumeric(b) && ~isempty(b) && any(strcmp(property, {'Date','ModificationDate'}))
        b = datestr(b);
      end
      out{end+1} = b;
    elseif ~any(strcmp(property, method))
      % a string containing members MAIN TIME SPENT
      s = [];
      split = textscan(property,'%s','Delimiter','.'); split=split{end};
      split = split(:)';
      for k=split
        s(end+1).type='.';
        s(end).subs=k{1};
      end
      out{end+1} = subsref(this, s);              % calls subsref directly (single subsref level)
    end
  else % this is a compound property, such as get(this,'Data.Signal')
    out{end+1} = get_eval(this, property); % calls inline (below)
  end
end

if iscell(out) && length(out) == 1
  out = out{1};
end
varargout{1} = out;

if isempty(varargout)
  varargout={[]};
end

% ------------------------------------------------------------------------------
% evaluate property in a reduced environment
function this = get_eval(this, property)
  try
    this = eval([ 'this.' property ]);  % calls subsref by eval (recursive subsref levels)
  catch
    % disp([ 'Warning: failed evaluation of this.' property ' in object ' this.Tag ' "' this.Title '".' ]);
      try
          this = eval([ 'this.Data.' property ]);
      catch
          % disp([ 'Warning: failed evaluation of this.Data.' property ' in object ' this.Tag ' "' this.Title '".' ]);
          this = eval(property);              % this is a full expression: evaluate it...
      end
  end

function TF = isvalid(property)
  TF = isstrprop(property, 'alphanum') | (property == '_') | (property == '.');
  TF = all(TF);
  
    
