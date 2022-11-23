function fig = draw_surface_bluewhitered_gallery(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial)
% draw_surface_bluewhitered_gallery.m
%
% Draw multiple data on surface using blue-white-red colormap for
% negative-zero-positive values
%
% Inputs: surface_to_plot : surface structure with fields
%                           vertices - vertex locations [Vx3], V = number of vertices
%                           faces - which vertices are connected [Fx3], F = number of faces 
%         data_to_plot    : data to plot [VxP]
%                           P = number of independent data
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

num_modes = size(data_to_plot,2);
       
if with_medial
    factor_x_small = 1.02;
    factor_x_big = 2.1;
    init_x = 0.005;
    init_y = 0.01;
    length_x = (1-2*init_x)/(factor_x_small + factor_x_big*(num_modes-1) + 1);
    length_y = (1-2*init_y);
    fig = figure('Position', [200 200 num_modes*400 150]);
    for mode=1:num_modes
        data_to_plot(medial_wall,mode) = min(data_to_plot(:,mode))*1.1;
        clims = [min(data_to_plot(:,mode)), max(data_to_plot(:,mode))];
        if clims(2)<=0
            clims(2) = 0.01;
        end
    
        ax1 = axes('Position', [init_x+factor_x_big*length_x*(mode-1) init_y length_x length_y]);
        obj1 = patch(ax1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,mode), ...
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
        
        ax2 = axes('Position', [init_x+factor_x_small*length_x+factor_x_big*length_x*(mode-1) init_y length_x length_y]);
        obj2 = patch(ax2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,mode), ...
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
    end
else
    factor_x = 1.05;
    init_x = 0.01;
    init_y = 0.01;
    length_x = (1-2*init_x)/(factor_x*(num_modes-1) + 1);
    length_y = (1-2*init_y);
    fig = figure('Position', [200 200 num_modes*200 150]);
    for mode=1:num_modes
        ax = axes('Position', [init_x+factor_x*length_x*(mode-1) init_y length_x length_y]);
        obj1 = patch(ax, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,mode), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceLighting', 'gouraud');
        if strcmpi(hemisphere, 'lh')
            view([-90 0]);
        elseif strcmpi(hemisphere, 'rh')
            view([90 0]);
        end
        colormap(ax,[bluewhitered])
        axis off
        axis image
    end
end