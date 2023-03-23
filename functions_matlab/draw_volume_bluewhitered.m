function [fig, varargout] = draw_volume_bluewhitered(volume_to_plot, data_to_plot, camera_view, markersize)
% draw_volume_bluewhitered.m
%
% Draw data on 3D volume using blue-white-red colormap for
% negative-zero-positive values
%
% Inputs: volume_to_plot : volume of mask in nifti format (3D array)
%         data_to_plot   : data to plot in nifti format (3D array)
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

% extract mask coordinates
ind = find(volume_to_plot~=0);
[y,x,z] = ind2sub(size(volume_to_plot), ind);
coords = [x,y,z];

% extract voxels in data_to_plot that match the mask
data_to_plot_temp = data_to_plot(ind);
        
% plot
fig_temp = figure('Visible', 'off');
imagesc(data_to_plot_temp);
new_map = bluewhitered(size(data_to_plot_temp,1));
close(fig_temp)

[~, sort_ind] = sort(data_to_plot_temp, 'ascend');

fig = figure;
set(fig, 'Position', get(fig, 'Position').*[0 0 0.6 0.6]+[200 200 0 0])
ax = axes('Position', [0.01 0.01 0.98 0.98]);
% obj1 = scatter3(coords(sort_ind,1), coords(sort_ind,2), coords(sort_ind,3), markersize, new_map, 'filled');
obj1 = scatter3(coords(sort_ind,1), coords(sort_ind,2), coords(sort_ind,3), markersize, data_to_plot_temp(sort_ind), 'filled');
colormap(ax, bluewhitered)
set(ax, 'xlim', [min(coords(sort_ind,1)), max(coords(sort_ind,1))], ...
        'ylim', [min(coords(sort_ind,2)), max(coords(sort_ind,2))], ...
        'zlim', [min(coords(sort_ind,3)), max(coords(sort_ind,3))])
axis square
view(camera_view)
axis off
        
varargout{1} = obj1;

