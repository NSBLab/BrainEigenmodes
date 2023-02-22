function [FaceVertexCData,new_cmap,new_climits,orig_data_climits] = makeFaceVertexCData_old(vertices,faces,vertex_id,data,cmap,colorFaceBoundaries,colorUnknownGrey,climits)

% This script will plot the boundaries defined by some parcellation/ROIS on
% a surface projection. It can also alter the colormap so regions that do
% not have any information are coloured grey.
%
% Inputs:
%
% vertices = the vertices making up the surface
%
% faces = the faces of the surface
%
% vertex_id = the roi id of each vertex
%
% data = either data for each individual roi or data for each vertex.
% If you don't want any data to be displayed for a roi or vertex, set that
% value to NaN.
%
% cmap = an N*3 matrix specifying the RGB values making up the colormap to
% use
%
% colorFaceBoundaries = set to 1 if you want the faces which make up the
% boundaries of each ROI to be coloured black. The code will then cofigure 
% FaceVertexCData to be a value per face instead of per vertex
%
% colorUnknownGrey = set to 1 to colour any unknown rois (i.e. vertices
% with an id of 0) grey
%
% climits = the range to apply the colormap. This will work perfectly fine
% if the range specified is larger than the data itself or if the upper
% limit is larger. However if the lower limit is larger than the smallest
% value in data, then it may get strange. If colorUnknownGrey = 1, then
% faces/vertices with a value smaller than the limit will be coloured grey.
% If it is set to 0 and 'faces' is used, those regions will be set to black
% (if 'centroid' or midpoint' are selected, the colormap will work
% appropriately). So if you really need to enforce a lower limit I would
% suggest threshold the data in advance and all should be good.
%
% Outputs:
% 
% FaceVertexCData = the new value for each vertex or face (depending on how
% colorFaceBoundaries was configured). 
%
% new_cmap = the new color map configured to work with the data in FaceVertexCData
%
% new_climits = the new caxis limits for FaceVertexCData.
%
% orig_data_climits = the limits for data. If climits is used, this is just
% that, but if not this will give you the limits for the original data. If
% making a colorbar this should be used to set the range of that.
%
% Stuart Oldham, Monash University, 2020
% Thanks to the coronavirus for giving me the time to make this script

% Set the boundary colour (black)
boundary_color = [0 0 0];

% Set the color of unknown ROIs (grey)
unknown_color = [.5 .5 .5];

if nargin < 7
    colorUnknownGrey = 1;
end

if nargin < 8
    climits = [nanmin(data) nanmax(data)];
end

if sum(vertex_id==0)>0
    vert0present = 1;
else
    vert0present = 0;
end

if length(data) ~= length(unique(vertex_id))-vert0present && length(data) ~= length(vertex_id)
    error('data needs to either contain one value per roi, or contain a value for each vertex')
end

% Because some steps require concatination in a specific dimension,
% the input data needs to be configured such that it is an 1*N array

if size(data,1) > size(data,2)
    data = data';
end

cmax = nanmax(climits);
cmin = nanmin(climits);
orig_data_climits = [cmin cmax];

if colorFaceBoundaries == 1

        % Check if the input_data is data for each ROI, or is data for each
        % vertex
        
        if length(data) ~= length(vertices)
            
            % Find the rois each face is connected to

            faces_roi_ids = vertex_id(faces);
            
            % Find the faces the exist entirely within a roi

            Faces_same_roi = ~logical(diff(faces_roi_ids,2,2));
            
            % Define a matrix specifying the value of each face. By default
            % it assume each face is on the boundary (done so with a value
            % of -1)

            FACES_ROI_DATA = ones(length(faces),1)*-inf;
            
            % For faces that exist entirely within a roi, assign them the
            % id of the roi they reside in

            FACES_ROI_DATA(Faces_same_roi) = faces_roi_ids(Faces_same_roi,1);
            
            % Make it so the colormap has value inserted that represent the
            % boundary and unknown ROIs but the presence of these values 
            % won't affect how the colormap is applied to the data

            Nrois = length(data);
            
            % The boundary will be coloured according to 'boundary_color'
            % and the faces that don't belong to any roi are coloured 
            % according to 'unknown_color'

            new_cmap = [boundary_color;unknown_color; cmap];
            cmap_length = size(cmap,1);
            
            % Get the range of values we want to use the original colormap
            % for
            crange = linspace(cmin,cmax,cmap_length);

            % Find values outside the original colormap range we can use to
            % define the boundary and unassigned roi faces
            
            boundary_val = cmin-(diff(crange(1:2))*2);
            unknown_val = cmin-diff(crange(1:2));
                      
            data(isnan(data)) = unknown_val;
                        
            % Assign faces their new values

            %FaceVertexCData = changem(FACES_ROI_DATA,[boundary_val unknown_val data(1:Nrois)],[-inf 0:Nrois]);            
            newval = [boundary_val unknown_val data(1:Nrois)];
            oldval = [-inf 0:Nrois];
            
            FaceVertexCData = FACES_ROI_DATA;
            for k = 1:numel(newval)
                FaceVertexCData(FACES_ROI_DATA == oldval(k)) = newval(k);
            end

            
            % Define the new colormap limits
            
            new_climits = [boundary_val cmax];

        else

            % Define the value of each face as the mean of the values
            % assigned to its associated vertices
            
            FACES_MEAN_ROI_DATA = nanmean(data(faces),2);
            
            % Find the rois each face is connected to

            faces_roi_ids = vertex_id(faces);
            
            % Find the faces the exist entirely within a roi

            Faces_same_roi = ~logical(diff(faces_roi_ids,2,2));
            
            % Define a matrix specifying the value of each face. By default
            % it assume each face is on the boundary (done so with a value
            % of -inf)

            FACES_ROI_DATA = ones(length(faces),1)*-inf;
            
            % For faces that exist entirely within a roi, assign them the
            % new face value we defined earlier

            FACES_ROI_DATA(Faces_same_roi) = FACES_MEAN_ROI_DATA(Faces_same_roi,1);

            if colorUnknownGrey == 0 && sum(isnan(data)) == 0

            % Make it so the colormap has a value inserted that represents 
            % the boundary and unknown ROIs but the presence of these  
            % values won't affect how the colormap is applied to the data

                 % The boundary will be coloured according to 'boundary_color'

                new_cmap = [boundary_color; cmap];

                cmap_length = size(cmap,1);
                crange = linspace(cmin,cmax,cmap_length);

                boundary_val = cmin-(diff(crange(1:2)));
                
                new_climits = [boundary_val cmax];

                %FaceVertexCData = changem(FACES_ROI_DATA,boundary_val,-inf);
                            
                FaceVertexCData = FACES_ROI_DATA;
                
                FaceVertexCData(FACES_ROI_DATA == -inf) = boundary_val;
                

            elseif colorUnknownGrey == 1 || sum(isnan(data)) ~= 0

                % Any faces which are not assigned to a roi (i.e. none of their 
                % vertices have any roi information) are given value of inf. 
                % We find such faces by taking the maximum of the roi ids of
                % the associated vertices. Those that have a max of 0 are
                % unassigned faces

                FACES_ROI_DATA(max(faces_roi_ids,[],2)==0) = inf;

                FACES_ROI_DATA(isnan(data)) = inf;

                % The boundary will be coloured aaccording to 'boundary_color' and the faces
                % they don't belong to any roi are coloured according to
                % 'unknown_color'

                new_cmap = [boundary_color;unknown_color; cmap];
                cmap_length = size(cmap,1);
                crange = linspace(cmin,cmax,cmap_length);

                boundary_val = cmin-(diff(crange(1:2))*2);
                unknown_val = cmin-diff(crange(1:2));

                new_climits = [boundary_val cmax];
                
                %FaceVertexCData = changem(FACES_ROI_DATA,[boundary_val unknown_val],[-inf inf]);

                FaceVertexCData = FACES_ROI_DATA;
                FaceVertexCData(FACES_ROI_DATA == -inf) = boundary_val;
                FaceVertexCData(FACES_ROI_DATA == inf) = unknown_val;

            end

        end

elseif colorFaceBoundaries == 0

% This will assign a value per vertex

    if length(data) == length(vertex_id)

        if colorUnknownGrey == 0 && sum(isnan(data)) == 0
            warning('input_data is already configured how you want it')
            FaceVertexCData = data;
            new_cmap = cmap;
            new_climits = climits;
        else
        
            % See above if you wanna know how this works, I ain't writing
            % it out again

        new_cmap = [unknown_color; cmap];
        cmap_length = size(cmap,1);
        crange = linspace(cmin,cmax,cmap_length);

        unknown_val = cmin-diff(crange(1:2));

        new_climits = [unknown_val cmax];
        
        FaceVertexCData = data;
        
        FaceVertexCData(isnan(FaceVertexCData)) = unknown_val;

        FaceVertexCData(vertex_id==0) = unknown_val;

        end

    else
        
        Nrois = length(data);

        new_cmap = [unknown_color; cmap];
        cmap_length = size(cmap,1);
        crange = linspace(cmin,cmax,cmap_length);

        unknown_val = cmin-diff(crange(1:2));
        
        new_climits = [unknown_val cmax];
        
        data(isnan(data)) = unknown_val;

        %FaceVertexCData = changem(vertex_id,[unknown_val data(1:Nrois)],0:1:Nrois);    
        
        newval = [unknown_val data(1:Nrois)];
        oldval = 0:Nrois;

        FaceVertexCData = vertex_id;
        for k = 1:numel(newval)
            FaceVertexCData(vertex_id == oldval(k)) = newval(k);
        end
        
    end

end

% Make FaceVertexCData N*1, as it should be

if size(FaceVertexCData,2) > size(FaceVertexCData,1)
    FaceVertexCData = FaceVertexCData';
end