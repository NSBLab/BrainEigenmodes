function fig = video_surface_activity(surface_to_plot, data_to_plot, hemisphere, t, t_interest, ...
                                      is_time_ms, medial_wall, with_medial, cmap, show_colorbar, output_filename, save_video)
% video_surface_activity.m
%
% Draw and save (optional) video of activity data on surface
%
% Inputs: surface_to_plot : surface structure with fields
%                           vertices - vertex locations [Vx3], V = number of vertices
%                           faces - which vertices are connected [Fx3], F = number of faces 
%         data_to_plot    : data to plot [VxT]
%                           T - number of time points
%         hemisphere      : which hemisphere (string)
%                           lh - left hemisphere
%                           rh - right hemisphere
%         t               : time vector [1xT]
%         t_interest      : time snapshots vector [1xTs]
%         is_time_ms      : is time vector in ms? (boolean)
%                           1 = if time is in ms
%                           0 = if time is in s
%         medial_wall     : indices of the medial wall (vector)
%         with_medial     : draw medial wall view (boolean)
%         cmap            : colormap [mx3]
%                           m - number of colors
%         show_colorbar   : show colorbar (boolean)
%         output_filename : filename of video output (string)
%         save_video      : save video (boolean)
%
% Output: fig             : figure handle
%
% Original: James Pang, Monash University, 2022

%%
if nargin<11
    save_video = 0;
end

if nargin<10
    output_filename = 'waves';
end

if nargin<9
    show_colorbar = 1;
end

if nargin<8
    cmap = bluewhitered;
end

if nargin<7
    with_medial = 0;
end

if size(t,2)~=1
    t = t';
end
if size(t_interest,2)~=1
    t_interest = t_interest';
end

if is_time_ms==1
    time_units = 'ms';
else
    time_units = 's';
end

t_interest_ind = dsearchn(t, t_interest);
data_to_plot = data_to_plot(:,t_interest_ind);
data_to_plot = data_to_plot./repmat(max(data_to_plot,[],1), size(data_to_plot,1), 1);
data_to_plot(isnan(data_to_plot)) = 0;

if save_video
    writerObj = VideoWriter(output_filename);
    writerObj.FrameRate=100;
    writerObj.Quality=100; % 1=min, 100=best
    open(writerObj);
end

if with_medial
    factor_y = 1.1;
    init_x = 0.005;
    init_y = 0.01;
    length_x = (1-2*init_x);
    length_y = (0.95-init_y)/(factor_y*(2-1) + 1);

    fig = figure('Position', [200 200 350 500], 'color', 'w');

    t_ind = 1;
    data_orig = data_to_plot;
    nan_val = min(data_orig(setdiff(1:size(data_to_plot,1), medial_wall),t_ind));
    data_to_plot(isnan(data_to_plot(:,t_ind)),t_ind) = nan_val;
    data_to_plot(medial_wall,t_ind) = nan_val-1e-2;
    clims = [min(data_to_plot(:,t_ind)), max(data_to_plot(:,t_ind))];
    if clims(2)<=0
        clims(2) = 0.01;
    end

    ax1 = axes('Position', [init_x init_y+factor_y*length_y*(2-1) length_x length_y]);
    obj1 = patch(ax1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,t_ind), ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    caxis(clims)
    camlight('headlight')
    material dull
    colormap(ax1,[0.5 0.5 0.5; cmap])
    axis off
    axis image
    title1 = title(sprintf('t = %.2f %s', t_interest(t_ind), time_units), 'fontsize', 12, 'fontweight', 'normal');

    ax2 = axes('Position', [init_x init_y+factor_y*length_y*(1-1) length_x length_y]);
    obj2 = patch(ax2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,t_ind), ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([-90 0]);
    end
    caxis(clims)
    camlight('headlight')
    material dull
    colormap(ax2,[0.5 0.5 0.5; cmap])
    axis off
    axis image
    if show_colorbar
        cbar = colorbar(ax2,'southoutside');
        set(cbar, 'fontsize', 8, 'ticklength', 0.02, 'ytick', [min(data_to_plot(:,t_ind)), max(data_to_plot(:,t_ind))], 'yticklabel', {}, ...
            'position', [ax2.Position(1)+ax2.Position(3)*0.09, ax2.Position(2)+0.05, ax2.Position(3)*0.2, 0.02])
        ylabel(cbar, 'amplitude', 'fontsize', 12)
        annotation(fig, 'textbox', [cbar.Position(1)-0.06, cbar.Position(2)*1, 0.04, 0.02], 'string', 'min', 'edgecolor', 'none', ...
               'fontsize', 10, 'horizontalalignment', 'center', 'verticalalignment', 'middle')
        annotation(fig, 'textbox', [cbar.Position(1)+cbar.Position(3)+0.03, cbar.Position(2)*1, 0.04, 0.02], 'string', 'max', 'edgecolor', 'none', ...
                   'fontsize', 10, 'horizontalalignment', 'center', 'verticalalignment', 'middle') 
    end
    
    for t_ind = 1:length(t_interest)

        nan_val = min(data_orig(setdiff(1:size(data_to_plot,1), medial_wall),t_ind));
        data_to_plot(isnan(data_to_plot(:,t_ind)),t_ind) = nan_val;
        data_to_plot(medial_wall,t_ind) = nan_val-1e-2;
        clims = [min(data_to_plot(:,t_ind)), max(data_to_plot(:,t_ind))];
        if clims(2)<=0
            clims(2) = 0.01;
        end

        set(obj1, 'FaceVertexCData', data_to_plot(:,t_ind))
        set(obj2, 'FaceVertexCData', data_to_plot(:,t_ind))
        set(ax1, 'CLim', clims, 'Colormap', [0.5 0.5 0.5; cmap])
        set(ax2, 'CLim', clims, 'Colormap', [0.5 0.5 0.5; cmap])
        set(title1, 'String', sprintf('t = %.2f %s', t_interest(t_ind), time_units))

        if save_video
            writeVideo(writerObj, getframe(fig));
        end
        drawnow
    end
else
    init_x = 0.01;
    init_y = 0.05;
    length_x = (1-2*init_x);
    length_y = (0.95-init_y);

    fig = figure('Position', [200 200 350 300], 'color', 'w');

    t_ind = 1;
    data_orig = data_to_plot;
    nan_val = min(data_orig(setdiff(1:size(data_to_plot,1), medial_wall),t_ind));
    data_to_plot(isnan(data_to_plot(:,t_ind)),t_ind) = nan_val;
    data_to_plot(medial_wall,t_ind) = nan_val-1e-2;
    clims = [min(data_to_plot(:,t_ind)), max(data_to_plot(:,t_ind))];
    if clims(2)<=0
        clims(2) = 0.01;
    end

    ax1 = axes('Position', [init_x init_y length_x length_y]);
    obj1 = patch(ax1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,t_ind), ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    caxis(clims)
    camlight('headlight')
    material dull
    colormap(ax1,[0.5 0.5 0.5; cmap])
    axis off
    axis image
    title1 = title(sprintf('t = %.2f %s', t_interest(t_ind), time_units), 'fontsize', 12, 'fontweight', 'normal');
    if  show_colorbar
        cbar = colorbar(ax1,'southoutside');
        set(cbar, 'fontsize', 8, 'ticklength', 0.02, 'ytick', [min(data_to_plot(:,t_ind)), max(data_to_plot(:,t_ind))], 'yticklabel', {}, ...
            'position', [ax1.Position(1)+ax1.Position(3)*0.09, ax1.Position(2)+0.03, ax1.Position(3)*0.2, 0.032])
        ylabel(cbar, 'amplitude', 'fontsize', 12)
        annotation(fig, 'textbox', [cbar.Position(1)-0.06, cbar.Position(2)*1, 0.04, 0.02], 'string', 'min', 'edgecolor', 'none', ...
               'fontsize', 10, 'horizontalalignment', 'center', 'verticalalignment', 'middle')
        annotation(fig, 'textbox', [cbar.Position(1)+cbar.Position(3)+0.03, cbar.Position(2)*1, 0.04, 0.02], 'string', 'max', 'edgecolor', 'none', ...
                   'fontsize', 10, 'horizontalalignment', 'center', 'verticalalignment', 'middle') 
    end
    
    for t_ind = 1:length(t_interest)

        nan_val = min(data_orig(setdiff(1:size(data_to_plot,1), medial_wall),t_ind));
        data_to_plot(isnan(data_to_plot(:,t_ind)),t_ind) = nan_val;
        data_to_plot(medial_wall,t_ind) = nan_val-1e-2;
        clims = [min(data_to_plot(:,t_ind)), max(data_to_plot(:,t_ind))];
        if clims(2)<=0
            clims(2) = 0.01;
        end

        set(obj1, 'FaceVertexCData', data_to_plot(:,t_ind))
        set(ax1, 'CLim', clims, 'Colormap', [0.5 0.5 0.5; cmap])
        set(title1, 'String', sprintf('t = %.2f %s', t_interest(t_ind), time_units))

        if save_video
            writeVideo(writerObj, getframe(fig));
        end
        drawnow
    end
end

if save_video
    writeVideo(writerObj, getframe(fig));
    clear writerObj
end
