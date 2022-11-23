function [fig, varargout] = draw_surface_bluewhitered(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial)
% draw_surface_bluewhitered.m
%
% Draw data on surface using blue-white-red colormap for
% negative-zero-positive values
%
% Inputs: surface_to_plot : surface structure with fields
%                           vertices - vertex locations [Vx3], V = number of vertices
%                           faces - which vertices are connected [Fx3], F = number of faces 
%         data_to_plot    : data to plot [Vx1]
%         hemisphere      : which hemisphere (string)
%                           lh - left hemisphere
%                           rh - right hemisphere
%         medial_wall     : indices of the medial wall (vector)
%         with_medial     : draw medial wall view (boolean)
%
% Output: fig             : figure handle
%
% Original: James Pang, Monash University, 2022

%%
if nargin<5
    with_medial = 0;
end

if nargin<4
    medial_wall = [];
end
    
if with_medial
    data_to_plot(medial_wall) = min(data_to_plot)*1.1;
    clims = [min(data_to_plot), max(data_to_plot)];
    if clims(2)<=0
        clims(2) = 0.01;
    end

    fig = figure('Position', [200 200 600 300]);
    ax1 = axes('Position', [0.03 0.1 0.45 0.8]);
    obj1 = patch(ax1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceLighting', 'gouraud');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    caxis(clims)
    colormap(ax1,[0.5,0.5,0.5; bluewhitered])
    axis off
    axis image

    ax2 = axes('Position', [ax1.Position(1)+ax1.Position(3)*1.1 ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
    obj2 = patch(ax2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceLighting', 'gouraud');
    if strcmpi(hemisphere, 'lh')
        view([90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([-90 0]);
    end
    caxis(clims)
    colormap(ax2,[0.5,0.5,0.5; bluewhitered])
    axis off
    axis image
    
    varargout{1} = obj1;
    varargout{2} = obj2;
else
    fig = figure;
    set(fig, 'Position', get(fig, 'Position').*[0 0 0.6 0.6]+[200 200 0 0])
    ax = axes('Position', [0.01 0.01 0.98 0.98]);
    obj1 = patch(ax, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceLighting', 'gouraud');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    colormap(ax,[bluewhitered])
    axis off
    axis image
    
    varargout{1} = obj1;
end