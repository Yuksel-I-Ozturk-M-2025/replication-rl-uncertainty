function h = plot(a, p, varargin)
% h = plot(model, parameters, axes, ...) plot a model
%
%   @iFunc/plot applies the function 'model' using the specified parameters and axes
%     and function parameters 'pars' with optional additional parameters. A plot of
%     the model is then produced.
%     The axis rank '1' corresponds to axis usually labelled as 'y'
%     The axis rank '2' corresponds to axis usually labelled as 'x'
%
% input:  model: model function (iFunc, single or array)
%         parameters: model parameters (vector, cell or vectors) or 'guess'
%         ...: additional parameters may be passed, which are then forwarded to the model
% output: h: handle to the plotted object (handle)
%
% ex:     b=plot(gauss); plot(gauss*lorz, [1 2 3 4, 5 6 7 8]);
%
% Version: Aug. 22, 2017
% See also iFunc, iFunc/fit, iFunc/feval

% test if further arguments are iFuncs
h = [];
noassign = 0;
if nargin > 1
  if isa(p, 'iFunc'), 
      a = [ a(:) ; p ]; p = ''; noassign = 1; 
  end
  varg = varargin;
  for index=1:length(varargin)
    if isa(varargin{index}, 'iFunc'), 
        a = [ a(:) ; varargin{index} ]; varg(index) = []; noassign = 1; 
    end
  end
  varargin = varg;
end

% handle array of objects
colors = 'bgrcmk'; color0 = floor(rand*numel(colors));
if numel(a) > 1
  is = ishold;
  for index=1:numel(a)
    if iscell(p) && length(p) == numel(a)
      h(index) = feval(mfilename, a(index), p{index}, varargin{:});
    else
      h(index) = feval(mfilename, a(index), p, varargin{:});
    end
    if ndims(a) == 1 % set the color, line style
      set(h(index), 'color', colors(1+mod(index+color0, length(colors))));
    end
    hold on
  end
  if nargout == 0 && ~isempty(inputname(1)) && ~noassign % update array inplace
    assignin('caller', inputname(1), a);
  end
  if ~is, hold off; end
  return
end

% now single object
if nargin < 2, 
  p=a.ParameterValues;
end

if strcmp(p, 'guess'),   p = []; end
if strcmp(p, 'current'), p = a.ParameterValues; end
if isempty(p) && ~isempty(a.ParameterValues)
  p=a.ParameterValues;
else p=NaN; end % force evaluation of function

% evaluate the model value, and axes
[signal, a, ax, name] = feval(a, p, varargin{:});

% Parameters are stored in the updated model
if length(inputname(1))
  assignin('caller',inputname(1),a); % update in original object
end

if iscell(signal)
  ih = ishold;
  for index=1:numel(signal)
    if index > 1
      hold on
    end
    h = [ h iFunc_plot(name{index}, signal{index}, ax{index}) ];
    if ndims(a(index)) == 1 && strcmp(get(h(index),'Type'),'line')
      % change color of line
      set(h(index), 'color', colors(1+mod(index+color0, length(colors))));
    end
    h(index) = iFunc_plot_menu(h(index), a(index), name{index});
  end
  if ih == 1, hold on; else hold off; end
  return
end

% call the single plot method
try
  h = iFunc_plot(name, signal, ax);
  if ~isempty(h)
    h = iFunc_plot_menu(h, a, name);
  end
catch ME
  disp(getReport(ME))
  disp([ mfilename ': WARNING: could not plot Model ' name '. Skipping.' ])
end
  

% ------------------------------------------------------------------------------
% simple plot of the model "name" signal(ax)
function h=iFunc_plot(name, signal, ax)
% this internal function plots a single model, 1D, 2D or 3D.
h=[];
if isempty(signal)
  h = [];
  return
elseif isvector(signal)
  if isempty(ax), ax={linspace(-2,2)}; end
  if isscalar(signal), signal = signal*ones(size(ax{1})); end
  if all(~isfinite(signal)) signal = zeros(size(signal)); end
  h = plot(ax{1}, signal);
elseif ndims(signal) == 2
  h = surf(ax{2}, ax{1}, signal);
  view(3)
  set(h,'EdgeColor','None');
elseif ndims(signal) == 3
  if all(cellfun(@(x)numel(x)==length(x), ax)), [ax{:}]=ndgrid(ax{:}); end
  h =patch(isosurface(ax{2}, ax{1}, ax{3}, signal, mean(signal(:))));
  set(h,'EdgeColor','None','FaceColor','green'); alpha(0.7);
  light
  view(3)
else
  try 
    if ndims(signal) == 4
      signal=squeeze(signal(:,:,1,:));
      if all(cellfun(@(x)numel(x)==length(x), ax)), [ax{:}]=ndgrid(ax{:}); end
      x=ax{1}; y=ax{2}; z=ax{3}; t=ax{4};
      x=squeeze(x(:,:,1,:));
      y=squeeze(y(:,:,1,:));
      t=squeeze(t(:,:,1,:));
      h =patch(isosurface(y,x,t, signal, mean(signal(:))));
      set(h,'EdgeColor','None','FaceColor','green'); alpha(0.7);
      light
      view(3)
    else
      h=[];
    end
  end
end

if isempty(h)
  % we use iData plotting
  iD = iData(ax{:}, signal); sz = size(iD); sz(4:end) = 1;
  disp([ 'iFunc.plot: ' name ': Reducing ' num2str(ndims(iD)) '-th dimensional data to 3D ' mat2str(sz) ]);
  iD = resize(iD, sz);
  h=plot(iD);
end

set(h, 'DisplayName', name);

%-------------------------------------------------------------------------------
function h=iFunc_plot_menu(h, a, name)
% contextual menu for the single object being displayed
% internal functions must be avoided as it uses LOTS of memory

  % return when a Contextual Menu already exists
  % if ~isempty(get(h,   'UIContextMenu')), return; end
  
  uicm = uicontextmenu; 
  % menu About
  uimenu(uicm, 'Label', [ 'About ' a.Name ': ' num2str(a.Dimension) 'D model ...' ], ...
    'Callback', [ 'msgbox(getfield(get(get(gco,''UIContextMenu''),''UserData''),''properties''),' ...
                  '''About: Model ' name ''',' ...
                  '''custom'',getfield(getframe(gcf),''cdata''), get(gcf,''Colormap''));' ] );
  uimenu(uicm, 'Label', name) ;
  uimenu(uicm, 'Label','Show model code...', ...
    'Callback', [ 'TextEdit(getfield(get(get(gco,''UIContextMenu''),''UserData''),''Expression''), ''' a.Name ''')' ]);

  % make up title string and Properties dialog content
  properties={ [ 'Model ' a.Tag ': ' num2str(ndims(a)) 'D model' ], ...
               [ 'Name: ' name ], ...
               [ 'Description: ' a.Description ]};

  % Expression
  u=cellstr(a); ud.Expression = u;
  u = u(~strncmp('%', u, 1)); % remove comment lines 
  u=[ u{:} ];
  u(~isstrprop(u,'print'))=''; if ~isvector(u), u=u'; end
  if length(u) > 300, u = [ u(1:297) '...' ]; end
  properties{end+1} = [ 'Expression: ' u(:)' ];

  properties{end+1} = '[Parameters]';
  for p=1:length(a.Parameters)
    [pname, R] = strtok(a.Parameters{p}); % make sure we only get the first word (not following comments)
    R = strtrim(R);
    u = sprintf('  p(%3d)=%20s', p, pname);
    val  = [];
    if ~isempty(a.ParameterValues)
    try
      val = a.ParameterValues(p);
    end
    end
    if ~isempty(val), u = [ u sprintf('=%g', val) ]; end
    if ~isempty(R),   u = [ u sprintf('  %% entered as: %s', R) ]; end
    properties{end+1} = u;
    if p<10
      uimenu(uicm, 'Label', u);
    end
    if p==10 && length(a.Parameters) > 10
      uimenu(uicm, 'Label', '...');
    end
  end

  ud.properties=properties;     
  ud.handle    = h;
  ud.ParameterVales = a.ParameterValues;

  set(uicm,'UserData', ud);
  set(h,   'UIContextMenu', uicm); 
  
  % add contextual menu to the axis ============================================
  % contextual menu for the axis frame

  uicm = uicontextmenu;
  % menu Duplicate (axis frame/window)
  uimenu(uicm, 'Label', 'Duplicate View...', 'Callback', ...
     [ 'tmp_cb.g=gca;' ...
       'tmp_cb.f=figure; tmp_cb.c=copyobj(tmp_cb.g,gcf); ' ...
       'set(tmp_cb.c,''position'',[ 0.1 0.1 0.85 0.8]);' ...
       'set(gcf,''Name'',''Copy of ' a.Name '''); ' ...
       'set(gca,''XTickLabelMode'',''auto'',''XTickMode'',''auto'');' ...
       'set(gca,''YTickLabelMode'',''auto'',''YTickMode'',''auto'');' ...
       'set(gca,''ZTickLabelMode'',''auto'',''ZTickMode'',''auto'');']);
       
  if ndims(a) == 1 && ~isfield(ud,'contextual_1d')
    ud.contextual_1d = 1;
  end
  uimenu(uicm, 'Label','Toggle grid', 'Callback','grid');
  if ndims(a) >= 2 
    uimenu(uicm, 'Label','Reset Flat/3D View', 'Callback', [ ...
      '[tmp_a,tmp_e]=view; if (tmp_a==0 & tmp_e==90) view(3); else view(2); end;' ...
      'clear tmp_a tmp_e; lighting none;alpha(1);shading flat;rotate3d off;axis tight;' ]);
    uimenu(uicm, 'Label','Smooth View','Callback', 'shading interp;');
    uimenu(uicm, 'Label','Add Light','Callback', 'light;lighting phong;');
    uimenu(uicm, 'Label','Transparency','Callback', 'alpha(0.7);');
    uimenu(uicm, 'Label',[ 'Linear/Log ' strtok(a.Name) ],'Callback', 'if strcmp(get(gca,''zscale''),''linear'')  set(gca,''zscale'',''log''); else set(gca,''zscale'',''linear''); end');
    uimenu(uicm, 'Label','Linear/Log X axis', ...
    'Callback', 'if strcmp(get(gca,''xscale''),''linear'')  set(gca,''xscale'',''log''); else set(gca,''xscale'',''linear''); end');
    uimenu(uicm, 'Label','Linear/Log Y axis', ...
    'Callback', 'if strcmp(get(gca,''yscale''),''linear'')  set(gca,''yscale'',''log''); else set(gca,''yscale'',''linear''); end');
    uimenu(uicm, 'Label','Toggle Perspective','Callback', 'if strcmp(get(gca,''Projection''),''orthographic'')  set(gca,''Projection'',''perspective''); else set(gca,''Projection'',''orthographic''); end');
  else
    uimenu(uicm, 'Label','Reset View', 'Callback','view(2);lighting none;alpha(1);shading flat;axis tight;rotate3d off;');
    uimenu(uicm, 'Label',[ 'Linear/Log ' strtok(a.Name) ],'Callback', 'if strcmp(get(gca,''yscale''),''linear'')  set(gca,''yscale'',''log''); else set(gca,''yscale'',''linear''); end');
    uimenu(uicm, 'Label', 'Linear/Log axis','Callback', 'if strcmp(get(gca,''xscale''),''linear'')  set(gca,''xscale'',''log''); else set(gca,''xscale'',''linear''); end');
  end

  uimenu(uicm, 'Separator','on','Label', 'About iFit/iData', ...
    'Callback',[ 'msgbox(''' version(iData,2) sprintf('. Visit <http://ifit.mccode.org>') ''',''About iFit'',''help'')' ]);
  set(gca, 'UIContextMenu', uicm);
  
  if a.Dimension == 1
    ylabel(a.Name)
  else
    zlabel(a.Name)
  end

  title(name);

