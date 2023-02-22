function [p_left,p_right,c] = ExampleSurfacePlotFunction(surface,vertex_id,data,cmap,data_label)

% This is just a example for how to plot the medial and lateral sides of a
% surface in the same plot

% Inputs:
%
% surface = a structure with two fields: vertices (the vertices making up 
% the surface) and faces (the faces of the surface)
%
% vertex_id = the roi id of each vertex
%
% data = either data for each individual roi or data for each vertex.
% If you don't want any data to be displayed for a roi or vertex, set that
% value to NaN. Note that this assumes a roi with a vertex_id has no data
% to plot. Additionally, if the vertex_ids are non sequential (e.g., like
% what you get from an .annot file) then data can take the form of a ROI*2
% matrix, where ROI is the number of regions, each row is a particular 
% region with the first column being the data to plot and the second being
% the region ID (should correspond to what is in vertex_id)
%
% cmap = an N*3 matrix specifying the RGB values making up the colormap to
% use
%
% data_label = (optional) the name of the data being plotted
%
% Output
%
% p_left = patch object for the surface plotted on the left
%
% p_right = patch object for the surface plotted on the right
%
% c = colorbar object

if min(size(data)) == 1
    if size(data,2) > size(data,1)
        data = data';
    end
end

figure('Position',[461   462   560   325])
ax_sub1 = axes('Position',[0.005 .33 .49 .66]);
p_left = plotSurfaceROIBoundary(surface,vertex_id,data,'midpoint',cmap,1);
camlight(80,-10);
camlight(-80,-10);
view([-90 0])

axis off
axis image

ax_sub2 = axes('Position',[.505 .33 .489 .66]);
p_right = plotSurfaceROIBoundary(surface,vertex_id,data,'midpoint',cmap,1);
camlight(80,-10);
camlight(-80,-10);
view([90 0])
axis off
axis image
colormap(cmap)
caxis([nanmin(data(:,1)) nanmax(data(:,1))])
c = colorbar('Location','southoutside');
%set(c, 'xlim', [nanmin(data) nanmax(data)],'Position',[.1 .23 .8 .05],'FontSize',20);
set(c, 'Position',[.1 .23 .8 .05],'FontSize',20);
if exist('data_label','var')
    c.Label.String = data_label;
end