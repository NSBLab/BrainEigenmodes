function [p,boundary_plot,BOUNDARY] = plotSurfaceROIBoundary(surface,vertex_id,data,boundary_method,cmap,linewidth,climits)

% This script is a wrapper for findROIboundaries and makeFaceVertexCData so
% that they smoothly work together and you don't have to spend a lot of
% your own time making them work together. The code sets up the basics of
% the patch object

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
% boundary_method = 'faces', 'midpoint', 'centroid', 'edge_vertices', or 
% 'edge_faces'. 'faces' will find the faces which exist between ROIs and
% those will be coloured black to specify the boundary. 'midpoint' finds 
% the edges that connect the vertices of two different ROIs and takes the
% midpoint of the edge and uses those coordinates to define the boundary. 
% 'centroid' finds the faces which exist between ROIs and uses the centroid
% of those to draw the coordinates that define the boundary. 
% 'edge_vertices' finds the vertices which define the boundary of the ROI 
% and uses them for the coordinates. 'edge_faces' finds the edges along
% faces which make up the boundary and uses them for the coordinates
%
% cmap = an N*3 matrix specifying the RGB values making up the colormap to
% use
%
% linewidth = the width of the boundary when using 'midpoint', or 
% 'centroid'.
%
% climits = the range to apply the colormap. This will work perfectly fine
% if the range specified is larger than the data itself or if the upper
% limit is larger. However if the lower limit is larger than the smallest
% value in data, then it may get strange. If colorUnknownGrey = 1, then
% faces/vertices with a value smaller than the limit will be coloured grey,
% or potentially black if 'faces' is used. If it is set to 0 and 'faces' is
% used, those regions will be set to black while if 'centroid' or midpoint' 
% are selected, the colormap will work appropriately). So if you really 
% need to enforce a lower limit I would suggest threshold the data in 
% advance and all should be good.
%
% Outputs:
%
% p = the patch surface object
%
% boundary_plot = a strcture containing a line object defining each
% boundary
%
% BOUNDARY = the boundary of the rois. For 'faces' this will be a logical
% where a value of 1 indicates that face is on the boundary between ROIs.
% For 'midpoint', 'centroid' or 'edges', BOUNDARY will be a cell where each 
% element contains the coordinates of the points making up the boundary, 
% which can be plotted as a line. Note that each boundary defines a 
% continuous ROI, if a ROI is made up of multiple non-continuous parts 
% (i.e., the ROI is made up of multiple unconnected sections), there will 
% be a boundary for each of those parts

% Extract the faces and vertices
vertices = surface.vertices;
faces = surface.faces;

if nargin < 6
    linewidth = 2;
end

data_orig = data;

if sum(vertex_id==0)>0
    vert0present = 1;
else
    vert0present = 0;
end

if min(size(data)) == 1

    % Because some steps require concatination in a specific dimension,
    % the input data needs to be configured such that it is an 1*N array

    if size(data,2) > size(data,1)
        data = data';
    end
    
    if length(data) ~= length(unique(vertex_id))-vert0present && length(data) ~= length(vertex_id)
        error('''data'' needs to either contain one value per roi, contain a value for each vertex, or be an N*2 array showing which data to plot to which ROI ID')
    end
    
    if length(data) ~= length(vertices)    
       ROI_ids = (1:length(data))';
    end

else
    if size(data,2) ~= 2
       error('If providing ''data'' with ROI ids, then the first column needs to be the data and the second the ROI id') 
    end
    
    ROI_ids = data(:,2);
    data = data(:,1);
    
end

if nargin < 7 
    climits = [nanmin(data) nanmax(data)];
end

% Find the boundaries of the ROIs
switch boundary_method
    case 'none'
BOUNDARY = [];        
    case {'midpoint','centroid','edge_vertices','faces','edge_faces'}
BOUNDARY = findROIboundaries(vertices,faces,vertex_id,boundary_method);
end

% Set up some options. 
switch boundary_method
    case 'faces'
    colorFaceBoundaries = 1;
% If boundaries are defined by faces, then the face
% color method for the patch object needs to be 'flat' otherwise it won't
% work
    face_color_method = 'flat';

    case {'midpoint','centroid','edge_vertices','edge_faces'}
        
    colorFaceBoundaries = 0;
    
    if length(data) == length(vertex_id)
        % If there is data for each individual vertex, use 'interp' for the 
        % face color 
        face_color_method = 'interp';
        
    else
        % If there is data for each individual ROI, use 'flat' for the 
        % face color.      
        face_color_method = 'flat';
    end
    case 'none'
       colorFaceBoundaries = 0; 
       face_color_method = 'interp';
    otherwise
        error('Unrecognised ''boundary_method'' option')
end

% Get the vertex or face data and colormap to use

% This is the old way of doing it, it generated a new colormap instead of assigning colours to faces/vertices
%[FaceVertexCData,new_cmap,new_climits,orig_data_climits] = makeFaceVertexCData_old(vertices,faces,vertex_id,data,cmap,colorFaceBoundaries,1,climits);

FaceVertexCData = makeFaceVertexCData(vertices,faces,vertex_id,data_orig,cmap,climits,colorFaceBoundaries);

% Plot the surface using some preconfigured options
p = patch(surface);
set(p,'FaceVertexCData',FaceVertexCData,'EdgeColor','none','FaceColor',face_color_method,'Clipping','off');
p.FaceLighting = 'gouraud';
material dull

hold on

% Draw the boundary if 'midpoint' or 'centroid' was used.

switch boundary_method
    case {'midpoint','centroid','edge_vertices','edge_faces'}
        for i = 1:length(BOUNDARY)
           boundary_plot.boundary(i) = plot3(BOUNDARY{i}(:,1), BOUNDARY{i}(:,2), BOUNDARY{i}(:,3), 'Color', 'k', 'LineWidth',linewidth,'Clipping','off');
        end
    otherwise
        boundary_plot = [];
end
