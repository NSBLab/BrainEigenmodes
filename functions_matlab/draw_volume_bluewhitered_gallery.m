function fig = draw_volume_bluewhitered_gallery(volume_to_plot, data_to_plot, camera_view, markersize)
% draw_volume_bluewhitered_gallery.m
%
% Draw multiple data on 3D volume using blue-white-red colormap for
% negative-zero-positive values
%
% Inputs: volume_to_plot : volume of mask in nifti format (3D array)
%         data_to_plot   : data to plot in nifti format (4D array)
%                          4th dimension = number of independent data
%         camera_view    : camera view angles [1x2]
%                          azimuth and elevation
%
% Output: fig            : figure handle
%
% Original: James Pang, Monash University, 2022

%%
if nargin<4
    markersize = 80;
end

if nargin<3
    camera_view = [-27.5 40];
end

num_modes = size(data_to_plot,4);

% extract mask coordinates
ind = find(volume_to_plot~=0);
[y,x,z] = ind2sub(size(volume_to_plot), ind);
coords = [x,y,z];

% plot
factor_x = 1.0;
init_x = 0.01;
init_y = 0.01;
length_x = (1-2*init_x)/(factor_x*(num_modes-1) + 1);
length_y = (1-2*init_y);
fig = figure('Position', [200 200 num_modes*200 150]);
for mode=1:num_modes
    
    % extract voxels in data_to_plot that match the mask
    data_to_plot_temp = data_to_plot(:,:,:,mode);
    data_to_plot_temp = data_to_plot_temp(ind);

    fig_temp = figure('Visible', 'off');
    imagesc(data_to_plot_temp);
    new_map = bluewhitered(size(data_to_plot_temp,1));
    close(fig_temp)

    [~, sort_ind] = sort(data_to_plot_temp, 'ascend');

    ax = axes('Position', [init_x+factor_x*length_x*(mode-1) init_y length_x length_y]);
%     obj1 = scatter3(coords(sort_ind,1), coords(sort_ind,2), coords(sort_ind,3), markersize, new_map, 'filled');
    obj1 = scatter3(coords(sort_ind,1), coords(sort_ind,2), coords(sort_ind,3), markersize, data_to_plot_temp(sort_ind), 'filled');
    set(ax, 'xlim', [min(coords(sort_ind,1)), max(coords(sort_ind,1))], ...
            'ylim', [min(coords(sort_ind,2)), max(coords(sort_ind,2))], ...
            'zlim', [min(coords(sort_ind,3)), max(coords(sort_ind,3))])
    axis square
    view(camera_view)
    axis off
end