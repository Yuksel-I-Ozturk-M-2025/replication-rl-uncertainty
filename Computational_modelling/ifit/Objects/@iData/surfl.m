function h = surfl(a, option)
% h = surfl(s,option) : Plot a 2D/3D object as surface with light
%
%   @iData/surflsurf function to plot a 2D or 3D object with light
%     2D objects are shown as a surface
%     3D objects are shown as an isosurface with median value
%     The slice(a) method opens the interactive sliceomatic 3D viewer.
%
% input:  s: object or array (iData)
%         option: global option for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
% output: h: graphics object handles (cell)
% ex:     surfl(iData(peaks)); surfl(iData(flow));
%
% Version: Aug. 22, 2017
% See also iData, iData/plot, 
%          iData/slice, iData/contour3, iData/contourf, iData/mesh
%          iData/contour, iData/slice, iData/plot3, iData/surf, iData/surfc

if nargin ==1
	option='';
end
h = plot(a, [ 'surfl ' option ]);



