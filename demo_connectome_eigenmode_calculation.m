%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% demo_connectome_eigenmode_calculation.m
%%%
%%% MATLAB script to demonstrate how to calculate various high-resolution 
%%% connectome-based eigenmodes. In particular, the script demonstrates 
%%% how to
%%% (1) calculate connectome eigenmodes and
%%% (2) calculate exponential distance rule (EDR) connectome eigenmodes
%%%
%%% NOTE 1: The script calculates connectome-based eigenmodes of a
%%%         synthetic connectome and EDR connectome on the fsaverage5
%%%         surface template with 10k vertices. This is for the demo to 
%%%         perform fast calculations. However, other surface resolutions 
%%%         are provided in the data folder for you to use. Just change the
%%%         surface_interest variable below.
%%% NOTE 2: An empirical high-resolution connectome can also be used. Just
%%%         change the connectome variable below. However make sure that
%%%         the variable is an array of size
%%%         [number of vertices x number of vertices].
%%% NOTE 3: Current demo uses 50 modes. For a proper analysis, we advise 
%%%         using between 100 to 200 modes.
%%%
%%% Original: James Pang, Monash University, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load relevant repository MATLAB functions

addpath(genpath('functions_matlab'));

%% Load surface files for fsaverage5

surface_interest = 'fsaverage5_10k';
hemisphere = 'lh';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

% Load cortex mask
cortex = dlmread(sprintf('data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);

% Calculate number of vertices
num_vertices = length(cortex);

%% Calculate synthetic connectome eigenmodes

hemisphere = 'lh';
surface_to_analyze = surface_midthickness;

% =========================================================================
%                   Generate surface local connectivity                    
% =========================================================================

surface_connectivity = calc_surface_connectivity(surface_to_analyze);

% Remove vertices corresponding to the medial wall
surface_connectivity = surface_connectivity(cortex_ind, cortex_ind);

% =========================================================================
%                      Generate synthetic connectome                       
% =========================================================================

% Replace the connectome variable below with your empirical data
% However, make sure that you remove the vertices corresponding to the
% medial wall first
connectome = zeros(size(surface_connectivity));

% Generate pseudorandom integers from imin to imax for the upper
% triangle elements of the connectome
rng(1)
imin = 0;
imax = 1e4;
triu_ind = calc_triu_ind(connectome);
triu_val = randi([imin, imax], size(triu_ind,1), size(triu_ind,2));

connectome(triu_ind) = triu_val;

% Symmetrize connectome 
connectome = connectome + connectome';

% Threshold connectome
threshold = 0.01; % preserve 1% of connections; change accordingly
connectome_threshold = threshold_edges_proportional(connectome, threshold);

% =========================================================================
%            Combine surface local connectivity and connectome             
% =========================================================================

surface_with_connectome = surface_connectivity + (connectome_threshold>0);
surface_with_connectome(surface_with_connectome>0) = 1;

% =========================================================================
%                           Calculate the modes                            
% =========================================================================

% The speed of this part depends on the resolution of surface_with_connectome. 
% Higher resolution matrices = Slower
% eig_vec = eigenvectors (eigenmodes)
% eig_val = eigenvalues

% Because of high computational requirements, precalculated data are stored 
% in data/examples if you only want to visualize the modes.
% For full calculation, set is_demo_calculate = 1
is_demo_calculate = 0;
if is_demo_calculate
    num_modes = 50;
    [eig_vec_temp, eig_val] = calc_network_eigenmode(surface_with_connectome, num_modes);
    
    % Bring back medial wall vertices with zero values
    eig_vec = zeros(num_vertices, num_modes);
    eig_vec(cortex_ind,:) = eig_vec_temp(:,1:num_modes);
    save(sprintf('data/examples/synthetic_connectome_eigenmodes-%s_%i.mat', hemisphere, num_modes), 'eig_val', 'eig_vec', '-v7.3')
else
    num_modes = 50;
    load(sprintf('data/examples/synthetic_connectome_eigenmodes-%s_%i.mat', hemisphere, num_modes), 'eig_vec')
end

% =========================================================================
%                      Some visualizations of results                      
% =========================================================================

% 1st to 5th modes with medial wall view
mode_interest = [1:5];
surface_to_plot = surface_midthickness;
data_to_plot = eig_vec(:, mode_interest);
medial_wall = find(cortex==0);
with_medial = 1;

% NOTE: The resulting connectome eigenmodes from the above demo calculation
%       will not necessarily resemble those in the paper. This is because
%       the connectome used above was synthetically generated 
%       (i.e., not empirical).
fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Multiple connectome eigenmodes with medial wall view';

%% Calculate synthetic EDR connectome eigenmodes

hemisphere = 'lh';
surface_to_analyze = surface_midthickness;

% =========================================================================
% Calculate Euclidean distance of surface vertices without the medial wall                    
% =========================================================================

surface_dist = squareform(pdist(surface_to_analyze.vertices(cortex_ind,:)));

% =========================================================================
%                    Generate synthetic EDR connectome                     
% =========================================================================

% Probability function
Pspace_func = @(scale, distance) exp(-scale*distance);

% Generate pseudorandom numbers to compare with probability function
rng(1)
rand_prob = rand(size(surface_dist));
rand_prob = triu(rand_prob,1) + triu(rand_prob,1)';
rand_prob(1:1+size(rand_prob,1):end) = 1;

% Calculate probability
% Scale = 0.120 matches empirical structural connectivity data. But you can 
% change its value accordingly.
scale = 0.120;
Pspace = Pspace_func(scale, surface_dist);
Pspace(1:1+size(Pspace,1):end) = 0;
Pspace = Pspace/max(Pspace(:));
    
% Generate EDR connectome
connectome = double(rand_prob < Pspace);
connectome(1:1+size(connectome,1):end) = 0;
    
% =========================================================================
%                           Calculate the modes                            
% =========================================================================

% The speed of this part depends on the resolution of connectome. 
% Higher resolution matrices = Slower
% eig_vec = eigenvectors (eigenmodes)
% eig_val = eigenvalues

% Because of high computational requirements, precalculated data are stored 
% in data/examples if you only want to visualize the modes.
% For full calculation, set is_demo_calculate = 1
is_demo_calculate = 0;
if is_demo_calculate
    num_modes = 50;
    [eig_vec_temp, eig_val] = calc_network_eigenmode(connectome, num_modes);
    
    % Bring back medial wall vertices with zero values
    eig_vec = zeros(num_vertices, num_modes);
    eig_vec(cortex_ind,:) = eig_vec_temp(:,1:num_modes);
    save(sprintf('data/examples/synthetic_EDRconnectome_eigenmodes-%s_%i.mat', hemisphere, num_modes), 'eig_val', 'eig_vec', '-v7.3')
else
    num_modes = 50;
    load(sprintf('data/examples/synthetic_EDRconnectome_eigenmodes-%s_%i.mat', hemisphere, num_modes), 'eig_vec')
end

% =========================================================================
%                      Some visualizations of results                      
% =========================================================================

% 1st to 5th modes with medial wall view
mode_interest = [1:5];
surface_to_plot = surface_midthickness;
data_to_plot = eig_vec(:, mode_interest);
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Multiple EDR connectome eigenmodes with medial wall view';
