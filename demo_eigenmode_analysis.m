%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% demo_eigenmode_analysis.m
%%%
%%% MATLAB script to demonstrate how to use surface eigenmodes to analyze 
%%% fMRI data. In particular, the script demonstrates how to
%%% (1) reconstruct a task fMRI spatial map,
%%% (2) reconstruct a resting-state fMRI spatiotemporal map and functional
%%%     connectivity (FC) matrix, and
%%% (3) calculate the eigenmode-based power spectral content of a spatial map
%%%
%%% NOTE 1: The script can also be used to analyze fMRI data using other
%%%         types of surface eigenmodes (e.g., connectome eigenmodes). 
%%%         Just change the eigenmodes variable below. However, make sure
%%%         that the variable is an array of size
%%%         [number of vertices x number of modes]. 
%%% NOTE 2: Current demo uses 50 modes. For a proper analysis, we advise 
%%%         using between 100 to 200 modes. 200 template geometric 
%%%         eigenmodes are provided in data/template_eigenmodes.
%%%
%%% Original: James Pang, Monash University, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load relevant repository MATLAB functions

addpath(genpath('functions_matlab'));

%% Load surface files for visualization

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

% Load cortex mask
cortex = dlmread(sprintf('data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);

%% Reconstruct a single-subject task fMRI spatial map

hemisphere = 'lh';
num_modes = 50;

% =========================================================================
%                    Load eigenmodes and empirical data                    
% =========================================================================

% Load 50 fsLR_32k template midthickness surface eigenmodes
eigenmodes = dlmread(sprintf('data/examples/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));

% Replace above line with the one below and make num_modes = 200 if using the 200 modes provided at data/template_eigenmodes
% eigenmodes = dlmread(sprintf('data/template_eigenmodes/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));

% Load example single-subject tfMRI z-stat data
data = load(sprintf('data/examples/subject_tfMRI_zstat-%s.mat', hemisphere));
data_to_reconstruct = data.zstat;

% =========================================================================
% Calculate reconstruction beta coefficients using 1 to num_modes eigenmodes
% =========================================================================

recon_beta = zeros(num_modes, num_modes);
for mode = 1:num_modes
    basis = eigenmodes(cortex_ind, 1:mode);
    
    recon_beta(1:mode,mode) = calc_eigendecomposition(data_to_reconstruct(cortex_ind), basis, 'matrix');
end

% =========================================================================
%     Calculate reconstruction accuracy using 1 to num_modes eigenmodes    
% =========================================================================

% reconstruction accuracy = correlation of empirical and reconstructed data

% At vertex level
recon_corr_vertex = zeros(1, num_modes);               
for mode = 1:num_modes
    recon_temp = eigenmodes(cortex_ind, 1:mode)*recon_beta(1:mode,mode);

    recon_corr_vertex(mode) = corr(data_to_reconstruct(cortex_ind), recon_temp);
end

% At parcellated level
parc_name = 'Glasser360';
parc = dlmread(sprintf('data/parcellations/fsLR_32k_%s-%s.txt', parc_name, hemisphere));

recon_corr_parc = zeros(1, num_modes);               
for mode = 1:num_modes
    recon_temp = eigenmodes(:, 1:mode)*recon_beta(1:mode,mode);

    recon_corr_parc(mode) = corr(calc_parcellate(parc, data_to_reconstruct), calc_parcellate(parc, recon_temp));
end

% =========================================================================
%                      Some visualizations of results                      
% =========================================================================

% Reconstruction accuracy vs number of modes at vertex and parcellated levels
figure('Name', 'tfMRI reconstruction - accuracy');
hold on;
plot(1:num_modes, recon_corr_vertex, 'k-', 'linewidth', 2, 'displayname', 'vertex')
plot(1:num_modes, recon_corr_parc, 'b-', 'linewidth', 2, 'displayname', 'parcellated')
hold off;
leg = legend('fontsize', 12, 'location', 'southeast', 'box', 'off');
set(gca, 'fontsize', 10, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1])
xlabel('number of modes', 'fontsize', 12)
ylabel('reconstruction accuracy', 'fontsize', 12)

% Reconstructed spatial map using N = num_modes modes
N = num_modes;
surface_to_plot = surface_midthickness;
data_to_plot = eigenmodes(:, 1:N)*recon_beta(1:N,N);
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = sprintf('tfMRI reconstruction - surface map using %i modes', N);

%% Reconstruct a single-subject resting-state fMRI spatiotemporal map and FC matrix

hemisphere = 'lh';
num_modes = 50;

% =========================================================================
%                    Load eigenmodes and empirical data                    
% =========================================================================

% Load 50 fsLR_32k template midthickness surface eigenmodes
eigenmodes = dlmread(sprintf('data/examples/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));

% Replace above line with the one below and make num_modes = 200 if using the 200 modes provided at data/template_eigenmodes
% eigenmodes = dlmread(sprintf('data/template_eigenmodes/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));

% Load example single-subject rfMRI time series data
data = load(sprintf('data/examples/subject_rfMRI_timeseries-%s.mat', hemisphere));
data_to_reconstruct = data.timeseries;
T = size(data_to_reconstruct, 2);

% =========================================================================
% Calculate reconstruction beta coefficients using 1 to num_modes eigenmodes
% =========================================================================

recon_beta = zeros(num_modes, T, num_modes);
for mode = 1:num_modes
    basis = eigenmodes(cortex_ind, 1:mode);
    
    recon_beta(1:mode,:,mode) = calc_eigendecomposition(data_to_reconstruct(cortex_ind,:), basis, 'matrix');
end

% =========================================================================
%     Calculate reconstruction accuracy using 1 to num_modes eigenmodes    
% =========================================================================

% reconstruction accuracy = correlation of empirical and reconstructed data

% At parcellated level
parc_name = 'Glasser360';
parc = dlmread(sprintf('data/parcellations/fsLR_32k_%s-%s.txt', parc_name, hemisphere));
num_parcels = length(unique(parc(parc>0)));

% Extract upper triangle indices
triu_ind = calc_triu_ind(zeros(num_parcels, num_parcels));

% Calculate empirical FC
data_parc_emp = calc_parcellate(parc, data_to_reconstruct);
data_parc_emp = calc_normalize_timeseries(data_parc_emp');
data_parc_emp(isnan(data_parc_emp)) = 0;

FC_emp = data_parc_emp'*data_parc_emp;
FC_emp = FC_emp/T;
FCvec_emp = FC_emp(triu_ind);

% Calculate reconstructed FC and accuracy (slow to run with more modes)
FCvec_recon = zeros(length(triu_ind), num_modes);
recon_corr_parc = zeros(1, num_modes);               
for mode = 1:num_modes
    recon_temp = eigenmodes(:, 1:mode)*squeeze(recon_beta(1:mode,:,mode));
 
    data_parc_recon = calc_parcellate(parc, recon_temp);
    data_parc_recon = calc_normalize_timeseries(data_parc_recon');
    data_parc_recon(isnan(data_parc_recon)) = 0;

    FC_recon_temp = data_parc_recon'*data_parc_recon;
    FC_recon_temp = FC_recon_temp/T;

    FCvec_recon(:,mode) = FC_recon_temp(triu_ind);
                    
    recon_corr_parc(mode) = corr(FCvec_emp, FCvec_recon(:,mode));
end

% =========================================================================
%                      Some visualizations of results                      
% =========================================================================

% Reconstruction accuracy vs number of modes at parcellated level
figure('Name', 'rfMRI reconstruction - accuracy');
hold on;
plot(1:num_modes, recon_corr_parc, 'b-', 'linewidth', 2)
hold off;
set(gca, 'fontsize', 10, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1])
xlabel('number of modes', 'fontsize', 12)
ylabel('reconstruction accuracy', 'fontsize', 12)

% Reconstructed FC using N = num_modes modes
N = num_modes;
FC_recon = zeros(num_parcels, num_parcels);
FC_recon(triu_ind) = FCvec_recon(:,N);
FC_recon = FC_recon + FC_recon';
FC_recon(1:(num_parcels+1):num_parcels^2) = 1;

figure('Name', sprintf('rfMRI reconstruction - FC matrix using %i modes', N));
imagesc(FC_recon)
caxis([-1 1])
colormap(bluewhitered)
cbar = colorbar;
set(gca, 'fontsize', 10, 'ticklength', [0.02 0.02])
xlabel('region', 'fontsize', 12)
ylabel('region', 'fontsize', 12)
ylabel(cbar, 'FC', 'fontsize', 12)
axis image

%% Calculate modal power spectral content of spatial maps

hemisphere = 'lh';
num_modes = 50;

% =========================================================================
%                   Load eigenmodes and empirical data
% =========================================================================

% Load 50 fsLR_32k template midthickness surface eigenmodes
eigenmodes = dlmread(sprintf('data/examples/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));

% Replace above line with the one below and make num_modes = 200 if using the 200 modes provided at data/template_eigenmodes
% eigenmodes = dlmread(sprintf('data/template_eigenmodes/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));

% Load example neurovault spatial map
if strcmpi(hemisphere, 'lh')
    data = gifti('data/examples/neurovault_map_100259.L.func.gii');
elseif strcmpi(hemisphere, 'rh')
    data = gifti('data/examples/neurovault_map_100259.R.func.gii');
end
data_to_reconstruct = data.cdata;

% =========================================================================
%                Calculate reconstruction beta coefficients                
% =========================================================================

basis = eigenmodes(cortex_ind, 1:num_modes);    
recon_beta = calc_eigendecomposition(data_to_reconstruct(cortex_ind), basis, 'matrix');

% =========================================================================
%                      Calculate modal power spectrum                      
% =========================================================================

[~, spectrum_norm] = calc_power_spectrum(recon_beta);

% =========================================================================
%                      Some visualizations of results                      
% =========================================================================

% Normalized power spectrum
figure('Name', 'rfMRI reconstruction - accuracy');
bar(1:num_modes, spectrum_norm)
set(gca, 'fontsize', 10, 'ticklength', [0.02 0.02], 'xlim', [2 num_modes], 'yscale', 'log')
xlabel('mode', 'fontsize', 12)
ylabel('normalized power (log scale)', 'fontsize', 12);
