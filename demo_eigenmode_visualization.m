%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% demo_eigenmode_visualization.m
%%%
%%% MATLAB script to demonstrate how to visualize surface and volume 
%%% eigenmodes
%%%
%%% NOTE 1: The script can also be used to visualize other types of surface 
%%%         eigenmodes (e.g., connectome eigenmodes). Just change the 
%%%         eigenmodes variable below. However, make sure that the 
%%%         variable is an array of size
%%%         [number of vertices x number of modes]. 
%%% NOTE 2: The script can also be used to visualize other types of volume 
%%%         eigenmodes (e.g., functional gradients). Just change the 
%%%         V_eigenmodes variable below. However, make sure that the 
%%%         variable is a 4D array of size
%%%         [number of voxels x number of voxels x number of voxels x number of modes].
%%%
%%% Original: James Pang, Monash University, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load relevant repository MATLAB functions

addpath(genpath('functions_matlab'));

%% Visualize surface eigenmodes

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';

% =========================================================================
%                       Load surface and eigenmodes                        
% =========================================================================

% Load surface file
[vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

% Load cortex mask
cortex = dlmread(sprintf('data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);

% Load fsLR_32k template midthickness surface eigenmodes
num_modes = 200;
eigenmodes = dlmread(sprintf('data/template_eigenmodes/%s_midthickness-%s_emode_%i.txt', surface_interest, hemisphere, num_modes));

% =========================================================================
%                         Visualize one eigenmode                          
% =========================================================================

% In this example, the 2nd eigenmode
mode_interest = 2;
surface_to_plot = surface_midthickness;
data_to_plot = eigenmodes(:, mode_interest);
medial_wall = find(cortex==0);

% with medial wall view
with_medial = 1;
fig = draw_surface_bluewhitered_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = sprintf('Single surface eigenmode with medial wall view (mode = %i)', mode_interest);

% without medial wall view
with_medial = 0;
fig = draw_surface_bluewhitered_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = sprintf('Single surface eigenmode without medial wall view (mode = %i)', mode_interest);

% =========================================================================
%                      Visualize multiple eigenmodes                       
% =========================================================================

% In this example, the 1st to 5th eigenmodes
mode_interest = [1:5];
surface_to_plot = surface_midthickness;
data_to_plot = eigenmodes(:, mode_interest);
medial_wall = find(cortex==0);

% with medial wall view
with_medial = 1;
fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Multiple surface eigenmodes with medial wall view';

% without medial wall view
with_medial = 0;
fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Multiple surface eigenmodes without medial wall view';

%% Visualize volume eigenmodes (in this example, eigenmodes of the thalamus)
% NOTE: Subcortical volume masks and eigenmodes for striatum (striatum) and
%       hippocampus (hippo) are also provided in the data folder for you to
%       use. Just change the structure variable below.
    
hemisphere = 'lh';
structure = 'tha';

% =========================================================================
%                     Load volume mask and eigenmodes                      
% =========================================================================

% Load mask file
mask_filename = sprintf('data/template_surfaces_volumes/hcp_%s-%s_thr25.nii.gz', structure, hemisphere);
V_mask = niftiread(mask_filename);

% Load volume eigenmodes
eigenmodes_filename = sprintf('data/results/subcortical/hcp_%s-%s_emode_30_noconstant_zscore.nii.gz', structure, hemisphere);
V_eigenmodes = niftiread(eigenmodes_filename);

% =========================================================================
%                         Visualize one eigenmode                          
% =========================================================================

% In this example, the 1st eigenmode
mode_interest = 1;
camera_view = [-27.5, 40]; % good camera view for striatum is [-61.5, 27] and hippocampus is [41, 43] 
markersize = 80;

volume_to_plot = V_mask;
data_to_plot = V_eigenmodes(:,:,:,mode_interest);

fig = draw_volume_bluewhitered(volume_to_plot, data_to_plot, camera_view, markersize);
fig.Name = sprintf('Single volume eigenmode (mode = %i)', mode_interest);

% =========================================================================
%                      Visualize multiple eigenmodes                       
% =========================================================================

% In this example, the 1st to 5th eigenmodes
mode_interest = [1:5];
camera_view = [-27.5, 40]; % good camera view for striatum is [-61.5, 27] and hippocampus is [41, 43] 
markersize = 80;

volume_to_plot = V_mask;
data_to_plot = V_eigenmodes(:,:,:,mode_interest);

fig = draw_volume_bluewhitered_gallery(volume_to_plot, data_to_plot, camera_view, markersize);
fig.Name = 'Multiple volume eigenmodes';
