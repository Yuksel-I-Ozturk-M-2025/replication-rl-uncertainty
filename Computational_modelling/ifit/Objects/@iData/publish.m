function a = publish(a, filename)
% b = publish(s) : export the Dataset object as an HTML document.
%
%   @iData/publish function to export an object into a readable document.
%
%   publish(a, file)  export to given file name
%   publish(a, dir)   export to given directory
%
% input:  s: object or array (iData)
%         filename: a file name or directory where to export. When not given
%           the exportation takes place in the local directory.
% output: b: generated filename
% ex:     b=publish(iData(peaks)); webbrowser(b, 'system')
%
% Version: Aug. 22, 2017
% See also iData, iData/save

if nargin < 2, filename = ''; end

a = saveas(a, filename, 'html');

if nargout == 0
  web(a)
end
