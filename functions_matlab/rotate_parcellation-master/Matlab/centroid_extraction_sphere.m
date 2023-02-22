% Function to extract centroids of a cortical parcellation on the
% Freesurfer sphere, for subsequent input to code for the performance of a
% spherical permutation. Runs on individual hemispheres.
%
% Requires read_surf.m and read_annotation.m functions, distributed with Freesurfer (Applications > freesurfer > matlab).
%
% Inputs:
% path_sphere   path to the Fressurfer .annot file describing coordinates
%               of each vertex on the sphere
%               example: [path of FreeSurfer parcellation '/surf/rh.sphere']
% path_annot    path to the Freesurfer .annot file describing parcel membership of
%               each vertex
%               example: [path of FreeSurfer parcellation '/label/rh.HCP.fsaverage.aparc.annot']
%
% Output:
% centroid      coordinates of the centroid of each region on the sphere
%
% Rafael Romero-Garcia, 2017

function centroid=centroid_extraction_sphere(path_sphere,path_annot)
[vertex_coords_sp, faces] = read_surf(path_sphere);                 % read in spherical vertex coordinates
[vertices, label_annot, colortable] = read_annotation(path_annot);  % read in parcel membership of vertices

ind=0;          % iteration counter
centroid=[];    % initialisation of centroid array
for ic=1:colortable.numEntries      % loop over parcellated structures
    if isempty(strfind(colortable.struct_names{ic},'unknown')) && isempty(strfind(colortable.struct_names{ic},'corpus')) % exclude "unknown" structures and corpus callosum from the parcellation 
        ind=ind+1;                                                  % increment counter for every valid region
        label=colortable.table(ic,5);                               % ID of current parcel
        vertices_reg=find(label_annot==label);                      % vertices with given ID
        centroid(ind,:)=mean(vertex_coords_sp(vertices_reg,:));     % average coordinates of all vertices within the current parcel to generate the centroid
    end
end
