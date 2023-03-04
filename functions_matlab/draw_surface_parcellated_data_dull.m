function [fig, varargout] = draw_surface_parcellated_data_dull(surface_to_plot, data_to_plot, parc, hemisphere, medial_wall, with_medial, cmap)
% draw_surface_parcellated_data_dull.m
%
% Draw parcellated data on surface using a defined colormap with dull plot lighting
%
% Inputs: surface_to_plot : surface structure with fields
%                           vertices - vertex locations [Vx3], V = number of vertices
%                           faces - which vertices are connected [Fx3], F = number of faces
%         data_to_plot    : parcellated data [num_parcelsx1]
%         parc            : parcellation on surface [Vx1]
%         hemisphere      : which hemisphere (string)
%                           lh - left hemisphere
%                           rh - right hemisphere
%         medial_wall     : indices of the medial wall (vector)
%         with_medial     : draw medial wall view (boolean)
%         cmap            : colormap [Nx3]
%                           N = number of colors
%
% Output: fig             : figure handle
%
% Original: James Pang, Monash University, 2022

%%
parc_interest = relabel_parcellation(parc, 1);
parcels = unique(parc_interest(parc_interest>0));
num_parcels = length(parcels);

if nargin<7
    cmap = parula;
end

if nargin<6
    with_medial = 0;
end

if nargin<5
    medial_wall = [];
end

% adding gray color for zero values
cmap = [0.5 0.5 0.5; cmap];
    
% boundary style
boundary_method = 'midpoint';
BOUNDARY = findROIboundaries(surface_to_plot.vertices, surface_to_plot.faces, parc_interest, boundary_method);
          
% convert parcellated data into surface space
data_temp = data_to_plot;
data_to_plot = zeros(size(surface_to_plot.vertices,1),1);
for parcel_ind = 1:num_parcels
    parcel = parcels(parcel_ind);
    data_to_plot(parc_interest==parcel) = data_temp(parcel_ind);
end
    
if with_medial
    if min(data_to_plot)<0
        data_to_plot(medial_wall) = min(data_to_plot)*1.1;
    elseif min(data_to_plot)>0
        data_to_plot(medial_wall) = min(data_to_plot)*0.9;
    else
        data_to_plot(medial_wall) = -0.01;
    end
    clims = [min(data_to_plot), max(data_to_plot)];
    
    fig = figure('Position', [200 200 600 300]);
    ax1 = axes('Position', [0.03 0.1 0.45 0.8]);
    obj1 = patch(ax1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on;
    for i = 1:length(BOUNDARY)
        plot3(BOUNDARY{i}(:,1), BOUNDARY{i}(:,2), BOUNDARY{i}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
    end       
    hold off;       
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    caxis(clims)
    colormap(ax1, cmap)
    camlight('headlight')
%     camlight(80,-10);
%     camlight(-80,-10);
    material dull
    axis off
    axis image

    ax2 = axes('Position', [ax1.Position(1)+ax1.Position(3)*1.1 ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
    obj2 = patch(ax2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on;
    for i = 1:length(BOUNDARY)
        plot3(BOUNDARY{i}(:,1), BOUNDARY{i}(:,2), BOUNDARY{i}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
    end       
    hold off;       
    if strcmpi(hemisphere, 'lh')
        view([90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([-90 0]);
    end
    caxis(clims)
    colormap(ax2, cmap)
    camlight('headlight')
%     camlight(80,-10);
%     camlight(-80,-10);
    material dull
    axis off
    axis image
    
    varargout{1} = obj1;
    varargout{2} = obj2;
else
    fig = figure;
    set(fig, 'Position', get(fig, 'Position').*[0 0 0.6 0.6]+[200 200 0 0])
    ax = axes('Position', [0.01 0.01 0.98 0.98]);
    obj1 = patch(ax, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'flat');
    hold on;
    for i = 1:length(BOUNDARY)
        plot3(BOUNDARY{i}(:,1), BOUNDARY{i}(:,2), BOUNDARY{i}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
    end       
    hold off;
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    colormap(ax, cmap)
    camlight('headlight')
%     camlight(80,-10);
%     camlight(-80,-10);
    material dull
    axis off
    axis image
    
    varargout{1} = obj1;           
end

end

function parc_relabel = relabel_parcellation(parc, start_val, parcels)
% relabel_parcellation.m
%
% Relabel parcellation file to have consecutive values
%
% Inputs: parc         : volume (3d array) or surface (1d array) of parcellation
%         start_val    : starting value of the relabeled parcellation (integer)
%         parcels      : original parcel labels [vector]
%
% Output: parc_relabel : relabelled parcellation
%
% Original: James Pang, Monash University, 2021

%%

if nargin<3
    parcels = unique(parc(parc>0));
end
if nargin<2
    start_val = 1;
end

num_parcels = length(parcels);

parc_relabel = zeros(size(parc));

counter = start_val;

for parcel_ind = 1:num_parcels
    parcel_interest = parcels(parcel_ind);

    ind_parcel = find(parc==parcel_interest);
    
    parc_relabel(ind_parcel) = counter;
    
    counter = counter+1;
end
end