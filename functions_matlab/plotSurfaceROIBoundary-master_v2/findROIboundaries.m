function [BOUNDARY,BOUNDARY_ROI_ID,ROI_COMPONENTS] = findROIboundaries(vertices,faces,vertex_id,boundary_method)

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
% Outputs:
%
% BOUNDARY = the boundary of the rois. For 'faces' this will be a logical
% where a value of 1 indicates that face is on the boundary between ROIs.
% For 'midpoint', 'centroid' or 'edges', BOUNDARY will be a cell where each 
% element contains the coordinates of the points making up a boundary, 
% which can be plotted as a line. Note that each boundary defines a 
% continuous ROI, if a ROI is made up of multiple non-continuous parts 
% (i.e., the ROI is made up of multiple unconnected sections), there will 
% be a boundary for each of those parts
%
% BOUNDARY_ROI_ID = if 'midpoint' or 'centroid' is used, this will
% indicate to which roi each boundary belongs (as indicated by the 
% respective element in BOUNDARY e.g., BOUNDARY_ROI_ID(1) = 1 means
% BOUNDARY{1} defines a boundary for ROI 1, BOUNDARY_ROI_ID(2) = 1 means
% BOUNDARY{2} defines a boundary for ROI 1, BOUNDARY_ROI_ID(3) = 2 means
% BOUNDARY{3} defines a boundary for ROI 2 etc etc). BOUNDARY_ROI_ID will 
% be empty if 'faces' is used for boundary_method
%
% ROI_COMPONENTS = for each ROI, indicates how many componets make it up.
% In other words, this indicates if there are multiple nonspatially
% continuous parts to the ROI. This will only work if 'faces' is not 
% used
%
% Stuart Oldham, Monash University, 2020
% Thanks to the coronavirus for giving me the time to make this script

if nargin < 4
   boundary_method = 'midpoint'; 
end

switch boundary_method
    case 'faces'
           
            % Find the rois each face is connected to

            faces_roi_ids = vertex_id(faces);
            
            % Find the boundary faces 
            BOUNDARY = logical(diff(faces_roi_ids,2,2));
            
            BOUNDARY_ROI_ID = [];
            
            ROI_COMPONENTS = [];
            
    case {'midpoint','centroid','edge_vertices','edge_faces'}
    
    % Find the rois each face is connected to

    faces_roi_ids = vertex_id(faces);

    % Get the colouring (i.e., the ROI id) of the face. MATLAB
    % colours faces according to the value of the vertex in the
    % first column.

    face_color = faces_roi_ids(:,1);

    % Calculate how many rois each face belongs to
    n_rois_per_face = sum(diff(sort(faces_roi_ids,2),1,2)~=0,2)+1;
    
    boundary_faces = logical(diff(faces_roi_ids,2,2));

    % Get all possible combinations of edges

    edges = [[faces(:,1); faces(:,1); faces(:,2); faces(:,2); faces(:,3); faces(:,3)], ...
    [faces(:,3); faces(:,2); faces(:,1); faces(:,3); faces(:,2); faces(:,1)]];

    % Get the colour for each edge (an edge is assigned the colour
    % of the face it is assigned to, each edge therefore has two
    % colours)
    edge_color = repmat(face_color,6,1);

    % Get the number of faces
    Nfaces = size(faces,1);

    % Make an index to find which edge corresponds to each face)
    edge_faceID = repmat((1:Nfaces)',6,1);

    % Sort the edge vertices from low to high, allows unqiue edges
    % to be found
    [edges_sorted,edges_sorted_ind] = sortrows(sort(edges,2));

    clear edges

    % Extract all of the edges metadata into one matrix. 
    % edges_metadata(:,1:2) are the IDs of the vertices making up that
    % edge, edges_metadata(:,3:4) are the vertex ROI ids for the vertices
    % that make up the edge, edges_metadata(:,5) is the colour of the face
    % that edge is part of, edges_metadata(:,6) is the face ID of the faces
    % the edge is part of, edges_metadata(:,7) equals 1 if that edge is 
    % part of a boundary face
    edges_metadata = [edges_sorted  vertex_id(edges_sorted) edge_color(edges_sorted_ind,:) edge_faceID(edges_sorted_ind,:) boundary_faces(edge_faceID(edges_sorted_ind,:))];

    % Clean up some large matrices/arrays
    clear edges_sorted edge_faceID edge_color

    % Get the unique data for each edge. Each edge is assigned to 
    % two faces so will appear twice in this list
    edges_unqiue = unique(edges_metadata,'rows');

    % Because EDGES_unique is sorted such that each the data for
    % each edge is pair, we can split 
    edges_unique1 = edges_unqiue(1:2:end,:);
    edges_unique2 = edges_unqiue(2:2:end,:);

    % Extract all of the unique edges metadata into one matrix. 
    % unique_edges_metadata(:,1:2) are the IDs of the vertices making up 
    % that edge, unique_edges_metadata(:,3:4) are the vertex ROI ids for
    % the vertices that make up the edge, unique_edges_metadata(:,5)
    % indicates if that edge belongs to two faces of the same color, 
    % unique_edges_metadata(:,6:7) are the face IDs of the faces the edge
    % is part of, unique_edges_metadata(:,8) equals 1 if that edge is part
    % of a boundary face

    unique_edges_metadata = [edges_unique1(:,1:4) double((edges_unique1(:,5)-edges_unique2(:,5))~=0) edges_unique1(:,6) edges_unique2(:,6) max(edges_unique1(:,7),edges_unique2(:,7))];
    
    clear edges_unqiue edges_unique1 edges_unique2 
    
    switch boundary_method
    
    case {'midpoint','centroid'}
    
    % Find edges which are contained within the one roi
    
    non_boundary_edges_ind = diff(unique_edges_metadata(:,3:4),1,2)==0;
    
    % Find edges which lie between two rois (boundary edges)

    boundary_edges_vert_id = unique_edges_metadata(~non_boundary_edges_ind,1:2);

    % Get the associated roi ids for the boundary edges
    
    boundary_edges_roi_id = unique_edges_metadata(~non_boundary_edges_ind,3:4);
    
    % Get the face ids for the boundary edges

    boundary_edges_face_id = unique_edges_metadata(~non_boundary_edges_ind,6:7);
    
    % For the face ids assigned to the boundary edges, get how many rois
    % they connect

    boundary_edges_face_nrois = n_rois_per_face(boundary_edges_face_id);

    % Get all the roi ids and the number of rois
    
    roi_ids = unique(vertex_id);
    nrois = length(roi_ids);

    % Potentially a roi is not spatially contiguous, so will have multiple 
    % boundaries. The code will iteratively add in boundaries because of this
    
    nBounds = 1;
    
    ROI_COMPONENTS = zeros(nrois,1);

    for i = 1:nrois
        
        % Find the boundary edges for a given roi
        
        index = find(any(boundary_edges_roi_id == roi_ids(i), 2));
        
        % Get the face ids for the boundary edges for a given roi

        edge_roi_boundaryfaces = boundary_edges_face_id(index,:);

        % Get all the unique faces along the roi boundary and the number of 
        % such faces 
        
        roi_faces = unique(edge_roi_boundaryfaces);
        nfaces = length(roi_faces);

        % Relabel the faces from 1:nfaces so the boundary of the roi can be
        % found easily using graph theory methods
                
        newval = 1:nfaces;
        oldval = roi_faces;
        roi_boundaryfaces_ordered = edge_roi_boundaryfaces;
        
        for k = 1:numel(newval)
            roi_boundaryfaces_ordered(edge_roi_boundaryfaces == oldval(k)) = newval(k);
        end
    
        % Make an edge list using the relabeled faces as nodes. Each row
        % defines a pair of faces which are connected. The edge list needs
        % to be symmetric

        edgeList = [roi_boundaryfaces_ordered; fliplr(roi_boundaryfaces_ordered)];
        
        % Get the id of each edge
        
        edgeID = [index; index];
        
        % Make an adjacency matrix where each connection defines a pair of
        % faces which are adjacent, and the connection weight is the
        % corresponding edge id that connects those two faces

        adj = sparse(edgeList(:, 1), edgeList(:, 2), edgeID, nfaces, nfaces);

        % The boundary of a ROI should form a circle. If for some reason the 
        % roi is not spatially continuous, each of the subrois should form a 
        % cycle. If not something is very wrong with your parcellation

        % Find the components of the adjacency matrix. Each component
        % represents a spatially continuous area for the current roi
        
        [comp,ROI_COMPONENTS(i)] = graphComponents(adj);
        
        % Loop over each roi component (i.e., each boundary for a roi) 

        NodesWithDeg2 = sum(full(adj)>0)==2;
        
        for j = 1:ROI_COMPONENTS(i)
            
            % Find the "nodes" making up the component

            Nodes2Use = find((comp == j).*NodesWithDeg2);
            
            N1 = Nodes2Use(1);
            
            % Find the neighbours of the "node"

            NodeNeighbors = find(adj( N1,:));
            
            % Pick one of those neighbours to be the end "node"

            N2 = NodeNeighbors(1);
            
            % Find the connection linking the start and end "nodes". Once
            % found, record it, then remove it from the adjacency matrix 

            FinalConnection = full(adj(N1,N2));

            adj(N1,N2) = 0;

            adj(N2,N1) = 0;
            
            % The boundary should be a cycle in the graph. If you remove
            % one connection between "node" N1 and N2, if you find the 
            % shortest path between N1 and N2, you get the ordering of
            % faces which define the boundary. You can also get the order
            % of the connections traversed to make this path "E". 

            G = graph(adj);

            [~,~,E] = shortestpathtree(G,N1,N2,'OutputForm','cell');
            
            % We know the connections traversed, so we use this index to 
            % get the original edge ids which make up the boundary. To make
            % the boundary complete we add in the one edge we originally
            % removed. Each edge defines a coordinate.

            edgesINDS = [G.Edges.Weight(E{1}) ; FinalConnection];

            switch boundary_method

                case 'midpoint'
                    
    % The commented lines below define a quick and easy way to get the
    % boundary. However if done this way, where three rois meet the
    % boundaries of those rois will not intersect. This bothered me. So I
    % wrote more code to get around this you can see below. If this does
    % not bother you, you mad person, you can just set DoSimple to 1

            DoSimple = 0;
            
            if DoSimple
    
            Boundary_verts = boundary_edges_vert_id(edgesINDS,:);
                    
            midpoints = (vertices(Boundary_verts(:, 1), :) + vertices(Boundary_verts(:, 2), :))./2;
                            
            BOUNDARY{nBounds} = [midpoints; midpoints(1,:)];
        
            nBounds = nBounds + 1;
    
            else
            % For the rois boundary edges, get the number of rois the
            % corresponding faces are connected to

            Boundary_verts_nrois = boundary_edges_face_nrois(edgesINDS,:); 

            % Find the "corners" (i.e. where three rois intersect) of the
            % boundary

            boundary_corner_points = find(max(Boundary_verts_nrois,[],2)==3);

            % Check there are corners for this boundary

            if ~isempty(boundary_corner_points)

            % Prevent the boundary from starting on a corner, this messes up
            % code later. The corners should come in pairs in EdgeINDS and 
            % the code expects as such. If the first edge but not second
            % edge in a corner edge, this means the last edge is a corner
            % edge. If so we just reorder to matrix to add the starting
            % edge to the end
                
                if boundary_corner_points(1) == 1 && boundary_corner_points(2) == 1          
                    EdgeINDS_1 = edgesINDS(1,:);

                    edgesINDS(1,:) = [];

                    edgesINDS(length(edgesINDS)+1,:) = EdgeINDS_1;

                    Boundary_verts_nrois = boundary_edges_face_nrois(edgesINDS,:);

                    boundary_corner_points = find(max(Boundary_verts_nrois,[],2)==3);
                end        
            
                % Find the corresponding faces of the rois boundary edge 

                Boundary_faces = boundary_edges_face_id(edgesINDS,:);

                % Find the corresponding vertices of the rois boundary edge 

                Boundary_verts = boundary_edges_vert_id(edgesINDS,:);

                % Calculate the midpoints of the boundary edges based on the
                % vertices which make up that edge

                midpoints = (vertices(Boundary_verts(:, 1), :) + vertices(Boundary_verts(:, 2), :))./2;

                % Create a binary mask of which edges are corners

                Boundary_corner_faces_mask = Boundary_verts_nrois==3;

                % Find the faces which make up the corners

                corner_faces = Boundary_faces;

                corner_faces(~Boundary_corner_faces_mask) = 0;

                boundary_corner_faces = max(corner_faces(boundary_corner_points,:),[],2);

                % Find the vertices of the corner faces

                boundary_corner_faces_verts = faces(boundary_corner_faces,:);

                % Get the centroid of the corner faces

                corner_coords = (vertices(boundary_corner_faces_verts(:, 3), :) + vertices(boundary_corner_faces_verts(:, 2), :) + vertices(boundary_corner_faces_verts(:, 1), :))./3;

                % Multiple corner faces may be adjacent. This detects which
                % are adjacent, selects the appropriate coordinates to use
                % and the appropriate spot to place them in. This assumes
                % for a pair of adjacent corners, you take the coordinates
                % of the second of the pair but use the first to get the
                % point to insert the coordinate
                corner_coords2use_ind = find(diff(boundary_corner_points)==1);

                corner_coords2use = corner_coords(corner_coords2use_ind+1,:);

                corner_coords_insert = boundary_corner_points(corner_coords2use_ind);

                % Create the matrix defining the final number of boundary
                % coordinates

                boundary_coords = zeros(size(midpoints,1)+size(corner_coords2use,1),size(midpoints,2));       

                % Find indices for the non-corner boundary edges

                addRows = ismember(1:size(midpoints,1), corner_coords_insert);

                oldDataInd = (1:size(midpoints,1)) + cumsum([0, addRows(1:end-1)]);

                % Add in the non-corner boundary edges

                boundary_coords(oldDataInd,:) = midpoints(:,:);

                % Find indices for the corner boundary edges

                newDataInd = (1:length(corner_coords_insert)) + corner_coords_insert';

                % Add in the corner boundary edges

                boundary_coords(newDataInd,:) = corner_coords2use(:,:);

            else

                % If there are no corners, just get the midpoints and use them

                Boundary_verts = boundary_edges_vert_id(edgesINDS,:);

                boundary_coords = (vertices(Boundary_verts(:, 1), :) + vertices(Boundary_verts(:, 2), :))./2;  

            end

            % The edges define a set of coordinates, and the order of the
            % defines how these coordinates need to be connected to form
            % the boundary. However the boundary needs to end where it
            % start so we make the first coordinate the final one as well
            
            BOUNDARY{nBounds} = [boundary_coords; boundary_coords(1,:)];
            
            % We don't know before hand how many distinct ROI boundaries
            % there will be, so we update it as we go.
            
            BOUNDARY_ROI_ID(nBounds) = i;

            nBounds = nBounds + 1;
            
            end

            case 'centroid'
                
                % Find the pairs of faces making up the boundary

                Boundary_faces = boundary_edges_face_id(edgesINDS,:); 
                
                % Find the face that is not in the next pair of faces. This
                % is the face we need to calculate the coordinates at that
                % step

                Boundary_faces_next = [Boundary_faces(2:end,:); Boundary_faces(1,:)];

                Boundary_face_order_mask = any(Boundary_faces(:,1)==Boundary_faces_next,2);

                Boundary_face_order_masked = Boundary_faces;
                Boundary_face_order_masked([Boundary_face_order_mask ~Boundary_face_order_mask]) = 0;

                % Get the order of faces
                
                Boundary_face_order = max(Boundary_face_order_masked,[],2);
                
                % Get the vertex ids of those order faces
                
                Boundary_face_order_verts = faces(Boundary_face_order,:);
                
                % Get the face centroid to use as the coordinates for the
                % boundary

                Boundary_coords = (vertices(Boundary_face_order_verts(:, 3), :) + vertices(Boundary_face_order_verts(:, 2), :) + vertices(Boundary_face_order_verts(:, 1), :))./3;

                BOUNDARY{nBounds} = [Boundary_coords; Boundary_coords(1,:)];

                BOUNDARY_ROI_ID(nBounds) = i;

                nBounds = nBounds + 1;   

            end

        end

    end
    
    case 'edge_vertices'

            % Get the number of unique ROIs
            
            roi_ids = unique(vertex_id);
            nrois = length(roi_ids);
            
            ROI_COMPONENTS = zeros(nrois,1);
            
            nBounds = 1;
            
            % Loop over each ROI
            
            for i = 1:nrois
                
                % Find the edges that connect two vertices of the same
                % ROI and is also a boundary edge
                
                index = sum([sum(unique_edges_metadata(:,3:4) == roi_ids(i),2)==2 unique_edges_metadata(:,8)==1],2) == 2;
                
                % Make an edge list of unique boundary edges for the
                % desired ROI
                
                edgeList_half = unique_edges_metadata(index,1:2);
                
                edgeList_orig = [edgeList_half; fliplr(edgeList_half)];
                
                unique_verts = unique(edgeList_orig(:));
                
                % Get the number of vertices that make up the boundary.
                
                nverts = length(unique_verts);
                
                % To make the next bit easier, relabel the vertex ids in 
                % edgeList_orig to be from 1:nverts
                
                new_vert_id = 1:nverts;
                old_vert_id = unique_verts;
                
                edgeList = edgeList_orig;

                for k = 1:numel(new_vert_id)
                    edgeList(edgeList_orig == old_vert_id(k)) = new_vert_id(k);
                end
                
                % Make a sparse adjacency matrix with the relabelled
                % vertices. if you graph this it should be a ring. This
                % ring indicates the order in which the vertices are
                % connected to define the boundary of the ROI. Basically
                % the nodes of this network represent vertices in the
                % original surface mesh
                                
                adj = sparse(edgeList(:, 1), edgeList(:, 2), 1, nverts, nverts);
                
                % Find the components of the adjacency matrix. Each component
                % represents a spatially continuous area for the current roi

                [comp,ROI_COMPONENTS(i)] = graphComponents(adj);

                % For the next part to work correctly we need to select a node
                % with a degree of two.

                NodesWithDeg2 = sum(full(adj))==2;

                % Loop over each component

                for j = 1:ROI_COMPONENTS(i)

                    % Find the "nodes" making up the component

                    Nodes2Use = find((comp == j).*NodesWithDeg2);

                    % Define a start "node"

                    N1 = Nodes2Use(1);

                    % Find the neighbours of the "node"

                    NodeNeighbors = find(adj( N1,:));

                    % Pick one of those neighbours to be the end "node"

                    N2 = NodeNeighbors(1);

                    % Find the connection linking the start and end "nodes". Once
                    % found,then remove it from the adjacency matrix 

                    adj(N1,N2) = 0;

                    adj(N2,N1) = 0;

                    % The boundary should be a cycle in the graph. If you remove
                    % one connection between "node" N1 and N2, if you find the 
                    % shortest path between N1 and N2, you get the ordering of
                    % vertices which define the boundary

                    G = graph(adj);

                    TR = shortestpathtree(G,N1,N2,'OutputForm','cell');
                    %save('error.mat')
                    cycle = [TR{1} TR{1}(1)];

                    % Insert the original vertex IDs, then extract the coordinates
                    % of those IDs. This is the boundary of this component of
                    % the ROI

                    new_vert_id = 1:nverts;
                    old_vert_id = unique_verts;

                    boundary_verts = cycle;

                    for k = 1:numel(new_vert_id)
                        boundary_verts(cycle == new_vert_id(k)) = old_vert_id(k);
                    end

                    BOUNDARY{nBounds} = vertices(boundary_verts,:);

                    BOUNDARY_ROI_ID(nBounds) = i;

                    nBounds = nBounds + 1;   
                end
                  
            end

            case 'edge_faces'
            % Honestly this is so janky if it gives weird results use
            % any of the other options

            % edges which belong to faces with different colours must be on
            % the boundary, extract only those
            boundary_edges_ind = unique_edges_metadata(:,5)==1;
            
            % This is the same as EDGES_unique_sameColor_faceID but for
            % boundary edges only
            boundary_edges_metadata = unique_edges_metadata(boundary_edges_ind,:);         

            % Clear up more large matrices/arrays that are no longer
            % needed
                        
            boundary_edges_facecolors = face_color(boundary_edges_metadata(:,6:7));

            % Get the number of unique ROIs
            roi_ids = unique(vertex_id);
            nrois = length(roi_ids);
            
            ROI_COMPONENTS = zeros(nrois,1);
            
            nBounds = 1;
            
            % Loop over each ROI
            
            for i = 1:nrois
                
                % Find the edges for who at least on of their associated
                % faces has a face colour matching the ROI of interest
                
                roiID = roi_ids(i);
                
                index = sum(boundary_edges_facecolors == roiID,2)>=1;
                
                % Make an edge list of unique boundary edges for the
                % desired ROI
                
                edgeList_half = boundary_edges_metadata(index,1:2);
                
                edgeList_orig = [edgeList_half; fliplr(edgeList_half)];
                
                unique_verts = unique(edgeList_orig(:));
                
                % Get the number of vertices that make up the boundary.
                
                nverts = length(unique_verts);
                
                % To make the next bit easier, relabel the vertex ids in 
                % edgeList_orig to be from 1:nverts
                
                new_vert_id = 1:nverts;
                old_vert_id = unique_verts;
                
                vert_new2old = [(1:nverts)' unique_verts];
                
                edgeList = edgeList_orig;

                for k = 1:numel(new_vert_id)
                    edgeList(edgeList_orig == old_vert_id(k)) = new_vert_id(k);
                end
                
                % Make a sparse adjacency matrix with the relabelled
                % vertices. if you graph this it should be a ring. This
                % ring indicates the order in which the vertices are
                % connected to define the boundary of the ROI. Basically
                % the nodes of this network represent vertices in the
                % original surface mesh
                                
                adj = sparse(edgeList(:, 1), edgeList(:, 2), 1, nverts, nverts);
                                
                % Find the components of the adjacency matrix. Each component
                % represents a spatially continuous area for the current roi

                [comp,ROI_COMPONENTS(i)] = graphComponents(adj);

                % For the next part to work correctly we need to select a node
                % with a degree of two.

                NodesWithDeg2 = sum(full(adj))==2;

                % Loop over each component

                for j = 1:ROI_COMPONENTS(i)

                    % Find the "nodes" making up the component. Nodes in this
                    % case are the vertices of the surface
                    NodesInComp = find(comp == j);
                    Nodes2Use = find((comp == j).*NodesWithDeg2);

                    % Define a start "node"

                    N1 = Nodes2Use(1);

                    % Find the neighbours of the "node"

                    NodeNeighbors = find(adj( N1,:));

                    % Pick one of those neighbours to be the end "node"

                    N2 = NodeNeighbors(1);

                    % Find the connection linking the start and end "nodes". Once
                    % found,then remove it from the adjacency matrix 

                    adj(N1,N2) = 0;

                    adj(N2,N1) = 0;

                    % The boundary should be a cycle in the graph. If you remove
                    % one connection between "node" N1 and N2, if you find the 
                    % shortest path between N1 and N2, you get the ordering of
                    % vertices which define the boundary

                    G = graph(adj);

                    TR = shortestpathtree(G,N1,N2,'OutputForm','cell');
                    attempts = 0;

                    % Remember when I said it gets janky? So it begins....

                    % Because the boundary when definied by faces can be so
                    % irregular, it throws up all kinds of problems with
                    % finding a path the connects all the edges which define
                    % said boundary (this is actually a somewhat constrained 
                    % version of the Chinese Postman Problem. Here is my best
                    % attempt

                    % First check if all the "nodes" are part of the path, if 
                    % not try some maths n' stuff 
                    while min(ismember(NodesInComp,TR{1})) == 0

                        % lets try again
                        adjFULL = full(adj);

                        % First we need to find any nodes which have a "cycle"
                        % attached to them. A cycle in this case refers to a
                        % sub-group of nodes which are connected (typically in
                        % a loop) to an origin node. Because they form a loop,
                        % when finding the shortest path as per above, all the
                        % nodes in the cycle are not found. A node has a cycle
                        % attached to it if it has a degree > 2
                        NodesWithCycles = find(sum(adjFULL)>=3);

                        % To make matters more complicated, if two nodes which
                        % have cycles are themselves connected, oh boy it gets
                        % fun
                        CheckNodesWithCyclesConnected = adjFULL(NodesWithCycles,NodesWithCycles);

                        % Check if nodes with cycles are connected
                        if max(CheckNodesWithCyclesConnected(:))~=0
                            % Nodes with cycles are connected. Yayyyyy.....

                            % find the nodes with cycles that are connected 
                            [Ci,Cj] = ind2sub(size(CheckNodesWithCyclesConnected),find(triu(CheckNodesWithCyclesConnected)==1));

                            % Loop over the pairs of nodes that are connected
                            for x = 1:length(Ci)

                                % Get the node ID of the connected pair
                                Chki = NodesWithCycles(Ci(x));
                                Chkj = NodesWithCycles(Cj(x));

                                % Get the original vertex ids
                                Chk_orig_vert = [vert_new2old(Chki,2) vert_new2old(Chkj,2)];

                                % Find the faces these belong to.
                                Chk_faces = boundary_edges_metadata((sum(ismember(boundary_edges_metadata(:,1:2),Chk_orig_vert),2)==2),4:5);      

                                % Find the ROI ID of those faces
                                Chk_faces_ROI = face_color(Chk_faces);

                                % We only want the face whose ROI ID
                                % corresponds to the ROI we are currently
                                % drawing the boundary of
                                FACE2Fix = Chk_faces((Chk_faces_ROI==roiID));

                                % Get the vertices which make up the face
                                FACE2Fix_verts = faces(FACE2Fix,:);

                                % Map from the original vertex IDs to the node
                                % IDs being used here
                                Verts2Fix = find(ismember(vert_new2old(:,2),FACE2Fix_verts));

                                % Because this algorithm is iterative and works
                                % by adding new nodes and assigning the
                                % original vertex IDs, to prevent anything
                                % becoming overly messed up, only use
                                % node/vertex IDs that were defined originally
                                Verts2Fix=Verts2Fix((Verts2Fix<=nverts));

                                % We now have found the node which connects the
                                % two nodes which form cycles. Chki,Chkj and 
                                % Chkk form a fully connected triad
                                Chkk = Verts2Fix(~ismember(Verts2Fix,[Chki Chkj]));

                                % Get the new nodes ID
                                NewNode = length(adjFULL)+1;

                                % Disconnect one of the two cycle node from the
                                % triad
                                adjFULL(Chki,Chkj) = 0;
                                adjFULL(Chkj,Chki) = 0;
                                adjFULL(Chkk,Chkj) = 0;
                                adjFULL(Chkj,Chkk) = 0;

                                % Reconnect the remaining nodes in the triad to
                                % a new node. This means the loop is now
                                % broken, and all nodes should be connected by
                                % a single path. This new node will correspond
                                % to the one we just disconnected above
                                adjFULL(NewNode,Chki) = 1;
                                adjFULL(Chki,NewNode) = 1;
                                adjFULL(NewNode,Chkk) = 1;
                                adjFULL(Chkk,NewNode) = 1;

                                adj = sparse(adjFULL);

                                % Sometimes, the two nodes with cycles are
                                % connected such that doing the above will
                                % result in one cycle becoming completely
                                % disconnected from the network. It will never
                                % be found on the path, the poor thing. We can
                                % fix this however.
                                [~,ChkCOmps] = graphComponents(adj);   
                                if ChkCOmps~=1

                                % Disconnected the newly formed node from Chkk,
                                % reconnedt Chkk and Chkj, and connect the
                                % newly formed node to Chkj. The cycle should
                                % now be broken and the graph is fully
                                % connected once again
                                adjFULL(NewNode,Chkj) = 1;
                                adjFULL(Chkj,NewNode) = 1;                                
                                adjFULL(NewNode,Chkk) = 0;
                                adjFULL(Chkk,NewNode) = 0;    
                                adjFULL(Chkk,Chkj) = 1;
                                adjFULL(Chkj,Chkk) = 1;    

                                end

                                % For the newly formed node, give it the vertex id 
                                % of the original node it corresponds to
                                vert_new2old(NewNode,:) = [NewNode vert_new2old(Chkj,2)];                            
                            end

                        else

                            % Find nodes which form cycles but aren't connected
                            % to other nodes which form cycles (these should
                            % have been delt with above)
                        for NodeWithCycleIND = 1:length(NodesWithCycles)

                            % Get the ID of the node with a cycle
                            CycleNodeID = NodesWithCycles(NodeWithCycleIND);

                            % Find the neighours of the node (i.e., the nodes
                            % it is connected to)
                        CycleNodeIDNeighbors = find(adj( CycleNodeID,:));

                            % The nodes in the cycle should not have been found
                            % on the shortest path originally defined, so we
                            % can use this information to extract just those
                            % nodes
                        CycleNodeIDNeighborsCycleIn = ismember(CycleNodeIDNeighbors,TR{1});

                            % If somehow the shortestpath already has explored
                            % the cycle, great! We don't need to worry.
                            % Otherwise we do the following:
                        if max(CycleNodeIDNeighborsCycleIn) == 1
                            % Find the neighbour nodes in the cycle
                            NeighborsInCycle = CycleNodeIDNeighbors(CycleNodeIDNeighborsCycleIn==1);
                            % Find those outside the cycle
                            NeighborsOutCycle = CycleNodeIDNeighbors(CycleNodeIDNeighborsCycleIn==0);
                            % Make our new node which will correspond to the
                            % original node which is the origin of the cycle
                            NewNode = length(adjFULL)+1;

                            % Disconnect one of the neighbour cycle nodes from
                            % the cycle node, also disconnect one of the
                            % non-cycle neighbour nodes as well
                            adjFULL(CycleNodeID,NeighborsInCycle(1)) = 0;
                            adjFULL(NeighborsInCycle(1),CycleNodeID) = 0;
                            adjFULL(CycleNodeID,NeighborsOutCycle(1)) = 0;
                            adjFULL(NeighborsOutCycle(1),CycleNodeID) = 0;
                            % Reconnect those two disconnected nodes to the new
                            % node, the cycle is now broken (hopefully)!!
                            adjFULL(NewNode,NeighborsInCycle(1)) = 1;
                            adjFULL(NeighborsInCycle(1),NewNode) = 1;
                            adjFULL(NewNode,NeighborsOutCycle(1)) = 1;
                            adjFULL(NeighborsOutCycle(1),NewNode) = 1;         

                            % For the newly formed node, give it the vertex id 
                            % of the original node it corresponds to
                            vert_new2old(NewNode,:) = [NewNode vert_new2old(CycleNodeID,2)];
                        end
                        end

                        end

                        % Recompute the shortest path
                        adj = sparse(adjFULL);

                        G = graph(adj);

                        TR = shortestpathtree(G,N1,N2,'OutputForm','cell');

                        % I can't promise to account for all the weird kinds of
                        % ways the boundary can be definied/connected under
                        % this approach. The following should catch any
                        % instances where things are so weird I cannot resolve
                        % it so it will just spit out its best effort.
                        attempts = attempts + 1;
                        if attempts >= 10
                            warning(['The boundary for ROI ',num2str(roiID),' has defeated my algorithm. Shame shame shame. The boundary is likely too complex in shape! Giving up and just using the final solution found'])
                            break
                        end
                    end

                    % Map the node IDs back to the original vertex IDs and use
                    % those to define the coordinates making up the boundary
                    cycle = [TR{1} TR{1}(1)];

                    boundary_verts = cycle;

                    for k = 1:size(vert_new2old,1)
                        boundary_verts(cycle == vert_new2old(k,1)) = vert_new2old(k,2);
                    end

                    BOUNDARY{nBounds} = vertices(boundary_verts,:);

                    BOUNDARY_ROI_ID(nBounds) = roiID;

                    nBounds = nBounds + 1;   
                end
            end

    end
otherwise
   error('Unrecognised ''boundary_method'' option')
end
