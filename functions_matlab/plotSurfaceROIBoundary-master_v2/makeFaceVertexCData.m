function FaceVertexCData = makeFaceVertexCData(vertices,faces,vertex_id,data,cmap,climits,colorFaceBoundaries,unknown_color,boundary_color)

% This script will assign colours to each face/vertex by assigning each
% value in 'data' a colour from 'cmap'. This code also allows you to 
% indicate if a particular area should be displayed with any colour at
% all.
%
% Inputs:
%
% vertices = the vertices making up the surface
%
% faces = the faces of the surface
%
% vertex_id = the roi id of each vertex. Note a vertex_id = 0 means that
% region won't have any colours assigned to it
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
% colorFaceBoundaries = set to 1 if you want the faces which make up the
% boundaries of each ROI to be coloured black. The code will then cofigure 
% FaceVertexCData to be a value per face instead of per vertex
%
% unknown_color = the color to assign to all unknown regions (areas with a
% vertex_id == 0 or data == NaN
%
% boundary_color = the color for the boundary if it is being drawn
%
% Output:
% 
% FaceVertexCData = the color value for each vertex or face (depending on how
% colorFaceBoundaries was configured). 
%
% Stuart Oldham, Monash University, 2020
% Thanks to the coronavirus for giving me the time to make this script

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

if nargin < 6
    if length(data) == length(vertex_id)
        climits = [nanmin(data(vertex_id>0)) nanmax(data(vertex_id>0))];
    else
        climits = [nanmin(data) nanmax(data)];
    end
end

if nargin < 7
    colorFaceBoundaries = 0;
end

if nargin < 8
    % Set the color of unknown ROIs (grey)
    unknown_color = [.5 .5 .5];
end

if nargin < 9
    % Set the boundary colour (black)
    boundary_color = [0 0 0];    
end


cmax = nanmax(climits);
cmin = nanmin(climits);

if colorFaceBoundaries == 1

        % Check if the input_data is data for each ROI, or is data for each
        % vertex
        
        % Find the rois each face is connected to

        faces_roi_ids = vertex_id(faces);

        % Find the ROI id of each face

        face_roi_id = faces_roi_ids(:,1);

        % Find the boundary faces

        boundary = logical(diff(faces_roi_ids,2,2));
            
        if length(data) ~= length(vertices)

            % Map the data from each ROI onto each face which is part of that ROI

            Nrois = length(data);
            
            newval = [NaN; data(1:Nrois)];
            oldval = [0; ROI_ids];

            face_data = face_roi_id;
            for k = 1:numel(newval)
              face_data(face_roi_id == oldval(k)) = newval(k);
            end

        else

            % Define the value of each face as the mean of the values
            % assigned to its associated vertices
            
            face_data = nanmean(data(faces),2);
            face_data(face_roi_id==0) = NaN;
        end
         
            % Scale the data if needed

            face_data(face_data<cmin) = cmin;
            face_data(face_data>cmax) = cmax;

            Ndata = length(face_data);

            face_data(Ndata+1) = cmin;
            face_data(Ndata+2) = cmax;

            % Map to an index for the colormap

            color_ind = ceil(rescale(face_data,1,size(cmap,1)));

            % Temporarily assign NaNs (i.e., the value for unknown regions) 
            % to a value so logical indexing doesn't screw up 
            color_ind(isnan(color_ind)) = 1;
            
            % Get the color for each face

            FaceVertexCData = cmap(color_ind(1:Ndata),:);

            % Color unknown regions 

            FaceVertexCData(isnan(face_data),1) = unknown_color(1);
            FaceVertexCData(isnan(face_data),2) = unknown_color(2);
            FaceVertexCData(isnan(face_data),3) = unknown_color(3);

            % Color boundary faces

            FaceVertexCData(boundary,1) = boundary_color(1);
            FaceVertexCData(boundary,2) = boundary_color(2);
            FaceVertexCData(boundary,3) = boundary_color(3);

elseif colorFaceBoundaries == 0

% This will assign a value per vertex

    if length(data) == length(vertex_id)
        vert_data = data;
        vert_data(vertex_id==0) = NaN;
    else
    
        Nrois = length(data);

        newval = [NaN; data(1:Nrois)];

        oldval = [0; ROI_ids];

        vert_data = vertex_id;
        for k = 1:numel(newval)
            vert_data(vertex_id == oldval(k)) = newval(k);
        end
    
    end

    vert_data(vert_data<cmin) = cmin;
    vert_data(vert_data>cmax) = cmax;

    Ndata = length(vert_data);

    vert_data(Ndata+1) = cmin;
    vert_data(Ndata+2) = cmax;
    
    color_ind = ceil(rescale(vert_data,1,size(cmap,1)));

    % Temporarily assign NaNs (i.e., the value for unknown regions) 
    % to a value so logical indexing doesn't screw up 
    color_ind(isnan(color_ind)) = 1;
    
    FaceVertexCData = cmap(color_ind(1:Ndata),:);

    FaceVertexCData(isnan(vert_data),1) = unknown_color(1);
    FaceVertexCData(isnan(vert_data),2) = unknown_color(2);
    FaceVertexCData(isnan(vert_data),3) = unknown_color(3);

end

% Make FaceVertexCData N*1, as it should be

if size(FaceVertexCData,2) > size(FaceVertexCData,1)
    FaceVertexCData = FaceVertexCData';
end
