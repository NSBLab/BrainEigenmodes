% This demonstration is to show the different way this code can plot the
% ROI boundaries

load('surface_data.mat','lh_inflated_verts','lh_verts','lh_faces','lh_HCPMMP1','lh_aparc','lh_rand200')
load('example_data.mat','lh_sulc')

surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;

figure('Position',[0 0 1920 963])

ax1 = axes('Position',[0.01 0 .3 1]);

% This just plots the ROI ID number for each ROI

plotSurfaceROIBoundary(surface,lh_rand200,1:100,'faces',jet(100),2);

% The following options set up the patch object to look pretty. This works
% well for the left hemisphere (medial and lateral). Change the inputs to 
% 'view' to see the brain from different angles ([-90 0] for left and [90 0]
% for right I find works well)

camlight(80,-10);
camlight(-80,-10);
view([-90 0])
axis off
axis tight
axis equal
axis vis3d

ax2 = axes('Position',[0.01+(1/3) 0 .3 1]);

cmap = flipud(hot(130));

% Just make up some data for the example illustration. This represents some
% value for each ROI
random_data = normpdf(1:180,100,100);

plotSurfaceROIBoundary(surface,lh_HCPMMP1,random_data,'midpoint',cmap(1:100,:),2);
camlight(80,-10);
camlight(-80,-10);

view([90 0])

axis off
axis tight
axis equal
axis vis3d

ax3 = axes('Position',[0.01+(2/3) 0 .3 1]);

% This plots sulcal depth, which is defined for each vertex

surface.vertices = lh_verts;
plotSurfaceROIBoundary(surface,lh_aparc,lh_sulc,'centroid',parula(100),4);
camlight(80,-10);
camlight(-80,-10);

view([-90 0])

axis off
axis tight
axis equal
axis vis3d

% Demonstrate different types of plots

surface.vertices = lh_inflated_verts;
boundary_type = {'faces','midpoint','centroid','edge_vertices','edge_faces',...
    'faces','midpoint','centroid','edge_vertices','edge_faces'};
linewidth = 4;


% Set up the colors for each vertex. This mirrors how the redone colour 
% assignment is performed by makeFaceVertexCData
 lh_rand200_color_map = lines(34);
 lh_rand200_ = lh_rand200;
 lh_rand200_(lh_rand200==0)=NaN;
 lh_rand200_color_ind = ceil(rescale(lh_rand200_,1,size(lh_rand200_color_map,1)));
 lh_rand200_color_ind(isnan(lh_rand200_)) = 1;
 lh_rand200_color = lh_rand200_color_map(lh_rand200_color_ind,:);
 lh_rand200_color(isnan(lh_rand200_),:) = .5;

 FaceVertexCData = makeFaceVertexCData(surface.vertices,surface.faces,lh_rand200,lh_rand200,lh_rand200_color_map);
 
for i = 1:10

figure('Position',[0 0  1680 933])
    
    % The data here is just each ROIs own ID number

        if i < 6
            data = 1:100;
            cmap = lines(34);
            savename = ['ROIS_',boundary_type{i},'_flat.png'];
        else
            data = lh_sulc;
            cmap = parula(100);
            if i == 6
            savename = ['sulc_',boundary_type{i},'_flat.png'];
            else
            savename = ['sulc_',boundary_type{i},'_interp.png'];    
            end
        end
    
    p = plotSurfaceROIBoundary(surface,lh_rand200,data,boundary_type{i},cmap,linewidth);

    camlight(80,-10);
    camlight(-80,-10);

    view([-90 0])

    axis off
    axis tight
    axis equal
    %axis vis3d
    
    % Mapping on the ROI id of each vertex to help with understanding how
    % this all works
    hold on
        
    s = scatter3(lh_inflated_verts(:,1)*1.01,lh_inflated_verts(:,2),lh_inflated_verts(:,3),40,lh_rand200_color,'filled');
    s.Clipping = 'off';
    s.MarkerEdgeColor = 'k';
    s.LineWidth = 1;
    
    % This just zooms into the area of interest    
    ylim([-25.2699   -8.7600])
    zlim([20.2174   31.3705])

    p.EdgeColor = 'k';
    p.EdgeAlpha = .5;
    
    % If you wanted to make a colorbar, this is what you would have to do:
    % colormap(cmap)
    % caxis([min(data) max(data)])
    % c = colorbar
    
    %print(['./figures/',savename],'-dpng')

end
