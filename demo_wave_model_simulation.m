%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% demo_wave_model_simulation.m
%%%
%%% MATLAB script to demonstrate how to simulate the wave model on a
%%% cortical surface. In particular, the script demonstrates how to
%%% (1) simulate resting-state activity and
%%% (2) simulate evoked activity (targeting V1)
%%%
%%% NOTE 1: Current demo uses 50 modes. For a proper analysis, we advise 
%%%         using between 200 to 1000 modes. 200 template geometric 
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

%% Simulate resting-state activity

hemisphere = 'lh';
num_modes = 50;

% =========================================================================
%                    Load eigenmodes and eigenvalues data                  
% =========================================================================

% Load 50 fsLR_32k template midthickness surface eigenmodes
eigenmodes = dlmread(sprintf('data/examples/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
eigenvalues = dlmread(sprintf('data/examples/fsLR_32k_midthickness-%s_eval_%i.txt', hemisphere, num_modes));

% Replace above line with the one below and make num_modes = 200 if using the 200 modes provided at data/template_eigenmodes
% eigenmodes = dlmread(sprintf('data/template_eigenmodes/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
% eigenvalues = dlmread(sprintf('data/template_eigenmodes/fsLR_32k_midthickness-%s_eval_%i.txt', hemisphere, num_modes));

% =========================================================================
%         Simulate resting-state activity using Gaussian white noise       
% =========================================================================

param = loadParameters_wave_func;
param.tstep = 0.1; % in ms
param.tmax = 100;  % in ms
param.tspan = [0, param.tmax];
param.T = 0:param.tstep:param.tmax;

% Change value of param.is_time_ms to 1 because the time defined above is
% in ms. This is necessary as param.gamma_s needs to match the scale.
param.is_time_ms = 1;

% Method for solving the wave model (either 'ODE' or 'Fourier')
% 'Fourier' is faster for long time series
method = 'ODE';
% method = 'Fourier';

param.r_s = 30;      % (default) in mm
param.gamma_s = 116; % (default) in s^-1
if param.is_time_ms==1
    param.gamma_s = 116*1e-3;
end

% Create Gaussian white noise external input
% random number is set for now for reproducibility
rng(1)
ext_input = randn(size(eigenmodes,1), length(param.T));

[mode_activity_rest, simulated_activity_rest] = model_neural_waves(eigenmodes, eigenvalues, ext_input, param, method);

% =========================================================================
%                      Some visualizations of results                      
% =========================================================================

% Snapshot of activity snapshot every 10 ms
t_interest = [0:10:100];
t_interest_ind = dsearchn(param.T', t_interest');
surface_to_plot = surface_midthickness;
data_to_plot = simulated_activity_rest(:, t_interest_ind);
medial_wall = find(cortex==0);
with_medial = 0;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Activity snapshots without medial wall view';
colormap(fig, parula)

% Video of activity every 0.5 ms (increase this to better see the waves)
surface_to_plot = surface_midthickness;
data_to_plot = simulated_activity_rest;
t = param.T;
t_interest = [0:0.5:100];
is_time_ms = 1;
medial_wall = find(cortex==0);
cmap = parula;
output_filename = 'rest_waves';
save_video = 0;

fig = video_surface_activity(surface_to_plot, data_to_plot, hemisphere, t, t_interest, is_time_ms, medial_wall, cmap, output_filename, save_video);

%% Simulate evoked activity (stimulating V1)

hemisphere = 'lh';
num_modes = 50;

% =========================================================================
%                    Load eigenmodes and eigenvalues data                  
% =========================================================================

% Load 50 fsLR_32k template midthickness surface eigenmodes
eigenmodes = dlmread(sprintf('data/examples/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
eigenvalues = dlmread(sprintf('data/examples/fsLR_32k_midthickness-%s_eval_%i.txt', hemisphere, num_modes));

% Replace above line with the one below and make num_modes = 200 if using the 200 modes provided at data/template_eigenmodes
% eigenmodes = dlmread(sprintf('data/template_eigenmodes/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
% eigenvalues = dlmread(sprintf('data/template_eigenmodes/fsLR_32k_midthickness-%s_eval_%i.txt', hemisphere, num_modes));

% =========================================================================
%                    Load the HCPMMP1 (Glasser360) atlas                   
% =========================================================================

if strcmpi(hemisphere, 'lh')
    parc = gifti('data/parcellations/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.label.gii');
elseif strcmpi(hemisphere, 'rh')
    parc = gifti('data/parcellations/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.label.gii');
end

% =========================================================================
%               Find parcel number and index of vertices of V1             
% =========================================================================

if strcmpi(hemisphere, 'lh')
    ROI = 'L_V1_ROI';
elseif strcmpi(hemisphere, 'rh')
    ROI = 'R_V1_ROI';
end

parcel_number = parc.labels.key(strcmpi(parc.labels.name, ROI));
parcel_ind = find(parc.cdata==parcel_number);

% =========================================================================
%      Simulate evoked activity stimulating all vertices within V1 ROI     
% =========================================================================

param = loadParameters_wave_func;
param.tstep = 0.1; % in ms
param.tmax = 100;  % in ms
param.tspan = [0, param.tmax];
param.T = 0:param.tstep:param.tmax;

% Change value of param.is_time_ms to 1 because the time defined above is
% in ms. This is necessary as param.gamma_s needs to match the scale.
param.is_time_ms = 1;

% Method for solving the wave model (either 'ODE' or 'Fourier')
% 'Fourier' is faster for long time series
method = 'ODE';
% method = 'Fourier';

param.r_s = 30;      % (default) in mm
param.gamma_s = 116; % (default) in s^-1
if param.is_time_ms==1
    param.gamma_s = 116*1e-3;
end

% Create a 10 ms external input with amplitude = 20 to V1 ROI
% results are robust to amplitude
ext_input_time_range = [10, 20];
ext_input_time_range_ind = dsearchn(param.T', ext_input_time_range');
ext_input_amplitude = 20; 

ext_input = zeros(size(eigenmodes,1), length(param.T));
ext_input(parcel_ind, ext_input_time_range_ind(1):ext_input_time_range_ind(2)) = ext_input_amplitude; 

[mode_activity_evoke, simulated_activity_evoke] = model_neural_waves(eigenmodes, eigenvalues, ext_input, param, method);

% =========================================================================
%                      Some visualizations of results                      
% =========================================================================

% Snapshot of activity snapshot every 10 ms
t_interest = [10:10:100];
t_interest_ind = dsearchn(param.T', t_interest');
surface_to_plot = surface_midthickness;
data_to_plot = simulated_activity_evoke(:, t_interest_ind);
medial_wall = find(cortex==0);
with_medial = 0;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Activity snapshots without medial wall view';
colormap(fig, parula)

% Video of activity every 0.5 ms (increase this to better see the waves)
surface_to_plot = surface_midthickness;
data_to_plot = simulated_activity_evoke;
t = param.T;
t_interest = [10:0.5:100];
is_time_ms = 1;
medial_wall = find(cortex==0);
cmap = parula;
output_filename = 'evoke_waves';
save_video = 0;

fig = video_surface_activity(surface_to_plot, data_to_plot, hemisphere, t, t_interest, is_time_ms, medial_wall, cmap, output_filename, save_video);
