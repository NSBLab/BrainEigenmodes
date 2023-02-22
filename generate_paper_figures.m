%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate_paper_figures.m
%%%
%%% MATLAB script to generate the main figures of the paper
%%%
%%% NOTE : The configuration of your computer (e.g., screen resolution)
%%%        affects how the figures created will look. Hence, they will not 
%%%        100% visually match the figures in the paper, but the scientific 
%%%        contents are replicated.
%%%
%%% Original: James Pang, Monash University, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD REPO MATLAB FUNCTIONS

addpath(genpath('functions_matlab'))

%% MISCELLANEOUS VARIABLES

% NOTE: data provided is only for the below parameters, so don't change them
hemisphere = 'lh'; 
num_modes = 200;
parc_name = 'Glasser360';

data_empirical_folder = 'data/empirical';
data_results_folder = 'data/results';
data_figures_folder = 'data/figures';
data_template_surfaces_folder = 'data/template_surfaces_volumes';
data_template_eigenmodes_folder = 'data/template_eigenmodes';

%% LOAD SURFACES

mesh_interest = 'inflated';
[vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/fsLR_32k_%s-%s.vtk', mesh_interest, hemisphere));
surface_inflated.vertices = vertices';
surface_inflated.faces = faces';

mesh_interest = 'midthickness';
[vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/fsLR_32k_%s-%s.vtk', mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

%% LOAD EMPIRICAL AND RESULTS DATA

% HCP subject list
subjectfile = sprintf('%s/subject_list_HCP.txt', data_empirical_folder);
subject_list = dlmread(subjectfile);
num_subjects = length(subject_list);

% cortex and medial wall mask
cortex = dlmread(sprintf('%s/fsLR_32k_cortex-%s_mask.txt', data_template_surfaces_folder, hemisphere));
cortex_ind = find(cortex);
num_vertices = length(cortex);

% subjects in training and test sets
subject_list_training = 1:200;
subject_list_test = 201:255;
subject_list_all = 1:255;
num_subjects_training = length(subject_list_training);
num_subjects_test = length(subject_list_test);
num_subjects_all = length(subject_list_all);

% parcellation
filename_common_parcellation = @(parc_name, hemisphere) sprintf('data/parcellations/fsLR_32k_%s-%s.txt', parc_name, hemisphere);

% EMPIRICAL: raw task-activation maps
data_task = load(sprintf('%s/S255_tfMRI_ALLTASKS_raw_%s.mat', data_empirical_folder, hemisphere));

% EMPIRICAL: sample resting-state fMRI
data_dtseries = ft_read_cifti(sprintf('%s/subject_rfMRI_REST.dtseries.nii', data_empirical_folder));
if strcmpi(hemisphere, 'lh')
    data_resting = data_dtseries.dtseries(data_dtseries.brainstructure==1,:);
elseif strcmpi(hemisphere, 'rh')
    data_resting = data_dtseries.dtseries(data_dtseries.brainstructure==2,:);
end
T = size(data_resting,2);
clear data_dtseries

% EMPIRICAL: sample neurovault maps
temp = matfile(sprintf('%s/fsLR_32k_neurovault_maps_N=10000_%s.mat', data_empirical_folder, hemisphere));
data_neurovault = temp.neurovault_maps(:,1:10);

% EMPIRICAL: myelin map in HCPMMP1 (Glasser360 in our naming convention) atlas 
data_myelin = load(sprintf('%s/myelin_%s.mat', data_empirical_folder, parc_name));

% RESULTS: basis sets     
basis_geometric = dlmread(sprintf('%s/basis_geometric_midthickness-%s_evec_%i.txt', data_results_folder, hemisphere, num_modes));
basis_connectome = load(sprintf('%s/basis_connectome_midthickness-%s_evec_%i.mat', data_results_folder, hemisphere, num_modes));
basis_connectome = basis_connectome.eig_vec;
basis_connectome_density_matched = load(sprintf('%s/basis_connectome_density_matched_midthickness-%s_evec_%i.mat', data_results_folder, hemisphere, num_modes));
basis_connectome_density_matched = basis_connectome_density_matched.eig_vec;
basis_EDR = load(sprintf('%s/basis_connectome_EDR_midthickness-%s_evec_%i.mat', data_results_folder, hemisphere, num_modes));
basis_EDR = basis_EDR.eig_vec;
basis_PCA_task = load(sprintf('%s/basis_PCA_task_S255_training-%s_evec_%i.mat', data_results_folder, hemisphere, num_modes));
basis_PCA_task = basis_PCA_task.eig_vec;
basis_PCA_resting = load(sprintf('%s/basis_PCA_resting_S255_training-%s_evec_%i.mat', data_results_folder, hemisphere, num_modes));
basis_PCA_resting = basis_PCA_resting.eig_vec;
basis_Fourier = load(sprintf('%s/basis_Fourier_sphere-%s_evec_%i.mat', data_results_folder, hemisphere, num_modes));
basis_Fourier = basis_Fourier.eig_vec;

% RESULTS: task recon accuracy
% eigenmodes
load(sprintf('%s/task_reconstruction_corr_modes=%i_%s_%s.mat', data_results_folder, num_modes, parc_name, hemisphere), ...
             'task_recon_corr_geometric', 'task_recon_corr_connectome', ...
             'task_recon_corr_connectome_density_matched', 'task_recon_corr_EDR')
% removing low frequency
load(sprintf('%s/task_reconstruction_corr_remove_lowf_modes=%i_%s_%s.mat', data_results_folder, num_modes, parc_name, hemisphere), ...
             'task_recon_corr_remove_lowf_geometric', 'task_recon_corr_remove_lowf_connectome', ...
             'task_recon_corr_remove_lowf_connectome_density_matched', 'task_recon_corr_remove_lowf_EDR')    
% removing high frequency
load(sprintf('%s/task_reconstruction_corr_remove_highf_modes=%i_%s_%s.mat', data_results_folder, num_modes, parc_name, hemisphere), ...
             'task_recon_corr_remove_highf_geometric', 'task_recon_corr_remove_highf_connectome', ...
             'task_recon_corr_remove_highf_connectome_density_matched', 'task_recon_corr_remove_highf_EDR')  
% individual surfaces
load(sprintf('%s/task_reconstruction_corr_individual_modes=%i_%s_%s.mat', data_results_folder, num_modes, parc_name, hemisphere), ...
             'task_recon_corr_individual_geometric')

% RESULTS: resting recon accuracy
% eigenmodes
load(sprintf('%s/resting_reconstruction_corr_modes=%i_%s_%s.mat', data_results_folder, num_modes, parc_name, hemisphere), ...
             'resting_recon_corr_geometric', 'resting_recon_corr_connectome', ...
             'resting_recon_corr_connectome_density_matched', 'resting_recon_corr_EDR')
% individual surfaces
load(sprintf('%s/resting_reconstruction_corr_individual_modes=%i_%s_%s.mat', data_results_folder, num_modes, parc_name, hemisphere), ...
             'resting_recon_corr_individual_geometric')

% RESULTS: power spectrum (HCP, NeuroVault, noise)
load(sprintf('%s/spectrum_modes=%i_%s.mat', data_results_folder, num_modes, hemisphere), ...
             'spectrum_HCP', 'spectrum_neurovault', 'spectrum_noise')

% RESULTS: parcellated empirical resting-state time series 
load(sprintf('%s/S255_resting_empirical_%s.mat', data_results_folder, parc_name), 'resting_emp')

% RESULTS: model
model_fit = load(sprintf('%s/model_results_Glasser360_%s.mat', data_results_folder, hemisphere), ...
                         'FC_emp', 'FCD_emp', 'FC_model_wave', 'FCD_model_wave', ...
                         'FC_model_mass', 'FCD_model_mass', 'KS');
         
% RESULTS: wave model visual simulation
model_wave_visual = load(sprintf('%s/model_wave_neural_visual_%s.mat', data_results_folder, hemisphere), ...
                                 'tspan', 'simulated_neural_visual_vertex', 'simulated_neural_visual_parcel', ...
                                 'ROIS', 'ROI_names');

% RESULTS: wave model optimization
model_wave_optim = load(sprintf('%s/model_wave_rest_optimization_%s.mat', data_results_folder, hemisphere), ...
                                'rs_vec', 'optim_edge_FC', 'optim_node_FC', 'optim_FCD');

% RESULTS: spin rotations
load(sprintf('%s/spin_perm_id_%s_10000_%s.mat', data_results_folder, parc_name, hemisphere), 'perm_id')
   
%% DEFINE HCP TASK CONTRAST NAMES

tasks = {'SOCIAL', 'MOTOR', 'GAMBLING', 'WM', 'LANGUAGE', 'EMOTION', 'RELATIONAL'};

task_contrasts = get_HCP_task_contrasts;

contrasts = {};
counter = 0;
for ii = 1:length(tasks)
    task = tasks{ii};
    
    for jj=1:length(task_contrasts.(task))
        counter = counter+1;
        
        contrasts{counter,1} = task_contrasts.(task){jj};
    end
end

task_contrasts_ind = struct();
counter = 0;
for ii = 1:length(tasks)
    task = tasks{ii};
    
    task_contrasts_ind.(task) = counter + (1:length(task_contrasts.(task)));
    counter = counter + length(task_contrasts.(task));
end

% representative contrasts per task
representative_contrasts = {'social_tom_random'; ...
                            'motor_cue_avg'; ...
                            'gambling_punish_reward'; ...
                            'wm_2bk_0bk'; ...
                            'language_math_story'; ...
                            'emotion_faces_shapes'; ...
                            'relational_match_rel'};

%% DEFINE FIGURE-RELATED PROPERTIES

fontsize_axis = 10;
fontsize_label = 12;

%% FIGURE 1

% load all data relevant to Figure 1
data_Figure1 = load(sprintf('%s/Figure1.mat', data_figures_folder));

surface_to_plot = surface_midthickness;
num_modes = 200;
parc_name = 'Glasser360';
N_interest = [10, 100, 200];

cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];

% load parcellation
parc = dlmread(filename_common_parcellation(parc_name, hemisphere));
num_parcels = length(unique(parc(parc>0)));

fig = figure('Position', [200 200 1000 800]);

% =========================================================================
% A left: downsampled surface
% =========================================================================
mode_list = [1,2,3,4,20];
num_modes_to_plot = length(mode_list);

factor_x = 1.02;
factor_y = 1.08;
init_x = 0.02;
init_y = 0.52;
length_x = (0.77 - 1.*init_x)/(factor_x*(num_modes_to_plot-1) + 1);
length_y = (0.98 - init_y)/(factor_y*(2-1) + 1);

num_downsampled_vertices = 500;
surface_downsampled = gifti(sprintf('%s/fsLR_%i_midthickness-%s.surf.gii', data_template_surfaces_folder, num_downsampled_vertices, hemisphere));
data_to_plot = ones(size(surface_downsampled.vertices,1),1);

ax1 = axes('Position', [init_x init_y+factor_y*length_y*(2-1) length_x length_y]);    
patch(ax1, 'Vertices', surface_downsampled.vertices, 'Faces', surface_downsampled.faces, 'FaceVertexCData', data_to_plot, ...
           'EdgeColor', 'k', 'FaceColor', 'interp', 'FaceLighting', 'gouraud', 'linewidth', 0.5);
if strcmpi(hemisphere, 'lh')
    view([-90 0]);
elseif strcmpi(hemisphere, 'rh')
    view([90 0]);
end
colormap(ax1, 0.7*ones(1,3))
axis off
axis image
annotation(fig, 'textbox', [ax1.Position(1), ax1.Position(2)+ax1.Position(4)*0.75, ax1.Position(3), 0.05], 'string', {'surface mesh'}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
annotation(fig, 'arrow', [ax1.Position(1)+ax1.Position(3), init_x+factor_x*length_x*(1+1.3-1)], ...
                         (ax1.Position(2)+ax1.Position(4)*0.5)*ones(1,2), 'Color', 'k')  
annotation(fig, 'textbox', [ax1.Position(1)+ax1.Position(3), ax1.Position(2)+ax1.Position(4)*0.55, 0.05, 0.05], 'string', '$\Delta \psi = -\lambda \psi$', 'edgecolor', 'none', ...
                'fontsize', fontsize_axis+2, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex')

% =========================================================================
% A right: eigenmodes
% =========================================================================
for mode_ind=1:num_modes_to_plot
    mode = mode_list(mode_ind);

    if mode_ind==num_modes_to_plot
        ax1 = axes('Position', [0.03+init_x+factor_x*length_x*(mode_ind+1.3-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax1 = axes('Position', [init_x+factor_x*length_x*(mode_ind+1.3-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    end
    patch(ax1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', basis_geometric(:,mode), ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    camlight('headlight')
    material dull
    colormap(ax1,bluewhitered)
    axis off
    axis image
    if mode_ind==1
        annotation(fig, 'textbox', [ax1.Position(1), ax1.Position(2)+ax1.Position(4)*0.8, ax1.Position(3), 0.05], 'string', {'mode 1'; '\psi_1'}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
    elseif mode_ind==2
        annotation(fig, 'textbox', [ax1.Position(1), ax1.Position(2)+ax1.Position(4)*0.8, ax1.Position(3), 0.05], 'string', {'mode 2'; '\psi_2'}, 'edgecolor', 'none', ...
            'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle') 
    elseif mode_ind==3
        annotation(fig, 'textbox', [ax1.Position(1), ax1.Position(2)+ax1.Position(4)*0.8, ax1.Position(3), 0.05], 'string', {'mode 3'; '\psi_3'}, 'edgecolor', 'none', ...
            'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle') 
    elseif mode_ind==4
        annotation(fig, 'textbox', [ax1.Position(1), ax1.Position(2)+ax1.Position(4)*0.8, ax1.Position(3), 0.05], 'string', {'mode 4'; '\psi_4'}, 'edgecolor', 'none', ...
            'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle') 
    else
        annotation(fig, 'textbox', [ax1.Position(1), ax1.Position(2)+ax1.Position(4)*0.8, ax1.Position(3), 0.05], 'string', {'mode N'; '\psi_N'}, 'edgecolor', 'none', ...
            'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle') 
    end
end
    
annotation(fig, 'textbox', [ax1.Position(1)-0.056, ax1.Position(2)+ax1.Position(4)*0.45, 0.08, 0.02], 'string', '...', 'edgecolor', 'none', ...
            'fontsize', fontsize_axis*2, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
annotation(fig, 'arrow', [0.07+init_x+factor_x*length_x*(1.5+1.3-1), 0.85], ...
                         (ax1.Position(2)+0.025)*ones(1,2), 'Color', 'k')
annotation(fig, 'textbox', [0.0+init_x+factor_x*length_x*(1.06+1.3-1), ax1.Position(2)+0.01, 0.13, 0.02], 'string', {'low frequency /'; 'long wavelength'}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
annotation(fig, 'textbox', [0.855, ax1.Position(2)+0.01, 0.13, 0.02], 'string', {'high frequency /'; 'short wavelength'}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')

% =========================================================================
% B: mode decomposition
% =========================================================================
ax2 = axes('Position', [init_x+factor_x*length_x*(1-1) init_y+factor_y*length_y*(1-1) length_x*3.5 length_y*2]);
xstart = 0.07;
ystart = 0.6;
t1 = text(xstart+0.35, ystart-0.25, '$\displaystyle y(\mathbf{r},t) = \sum_{j=1}^N a_j(t)\psi_j(\mathbf{r})$', ...
    'fontsize', fontsize_label*1.3, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex');
t2 = text(xstart, ystart-0.35, '$\displaystyle y$', ...
    'fontsize', fontsize_label*1.2, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex');
text(xstart+0.22, t2.Position(2), '$\displaystyle \psi_1$', ...
    'fontsize', fontsize_label*1.2, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex');
text(xstart+0.22*2, t2.Position(2), '$\displaystyle \psi_2$', ...
    'fontsize', fontsize_label*1.2, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex');
text(xstart+0.22*3, t2.Position(2), '$\displaystyle \psi_3$', ...
    'fontsize', fontsize_label*1.2, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex');
t3 = text(xstart+0.11, ystart-0.45, '$\displaystyle = a_1 \times$', ...
    'fontsize', fontsize_label*1.2, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex');
text(xstart+0.11*3, t3.Position(2), '$\displaystyle + a_2 \times$', ...
    'fontsize', fontsize_label*1.2, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex');
text(xstart+0.11*5, t3.Position(2), '$\displaystyle + a_3 \times$', ...
    'fontsize', fontsize_label*1.2, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex');
text(xstart+0.11*7, t3.Position(2), '$\displaystyle + \ldots$', ...
    'fontsize', fontsize_label*1.2, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex');
axis off

for mode_ind=1:4
    if mode_ind>1
        mode = mode_list(mode_ind-1);
    end

    ax2_1 = axes('Position', [0.005+init_x+factor_x*length_x*(mode_ind-1)*0.76 init_y+factor_y*length_y*(0.8-1) length_x*0.4 length_y]);
    if mode_ind==1
        patch(ax2_1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_Figure1.task_map_emp.('social_tom_random'), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
    else
        patch(ax2_1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', basis_geometric(:,mode), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
    end
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    camlight('headlight')
    material dull
    colormap(ax2_1,bluewhitered)
    axis off
    axis image
end

% =========================================================================
% C: data
% =========================================================================
ax3 = axes('Position', [init_x+factor_x*length_x*(4.4-1) init_y+factor_y*length_y*(0.9-1) length_x*3 length_y*0.9]);
axis off
data_to_plot = data_resting;
data_to_plot = calc_normalize_timeseries(data_to_plot');
data_to_plot(isnan(data_to_plot)) = 0;
data_to_plot = data_to_plot';
for ii=1:3
    ax3_1 = axes('Position', [ax3.Position(1)+ax3.Position(3)*(0.15+0.05*(ii-1)) ax3.Position(2)+ax3.Position(4)*(0.3-0.23*(ii-1)) length_x*0.4 length_y]);
    if ii==1
        patch(ax3_1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_Figure1.task_map_emp.('motor_cue_avg'), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
    elseif ii==2
        patch(ax3_1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_Figure1.task_map_emp.('emotion_faces_shapes'), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
    elseif ii==3
        patch(ax3_1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_Figure1.task_map_emp.('wm_2bk_0bk'), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
    end
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    camlight('headlight')
    material dull
    colormap(ax3_1,bluewhitered)
    axis off
    axis image
    
    ax3_2 = axes('Position', [ax3_1.Position(1)+0.15 ax3_1.Position(2) ax3_1.Position(3) ax3_1.Position(4)]);
    patch(ax3_2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,ii), ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    camlight('headlight')
    material dull
    colormap(ax3_2,bluewhitered)
    axis off
    axis image
    
    annotation(fig, 'textbox', [ax3_1.Position(1)-0.065, ax3_1.Position(2)+ax3_1.Position(4)*0.4, 0.08, 0.02], 'string', sprintf('task %i', ii), 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
    if ii==1
        annotation(fig, 'arrow', [ax3.Position(1)+ax3.Position(3)*(0.45+0.05*(1-1)), ax3.Position(1)+ax3.Position(3)*(0.47+0.05*(3-1))], ...
                                 [length_y*0.5+ax3.Position(2)+ax3.Position(4)*(0.3-0.22*(1-1)), length_y*0.4+ax3.Position(2)+ax3.Position(4)*(0.3-0.22*(3-1))], ...
                                         'HeadLength', 6, 'HeadWidth', 6, 'Color', 'k')   
    end
    if ii==2
        annotation(fig, 'textbox', [ax3_2.Position(1)-0.068, ax3_2.Position(2)+ax3_2.Position(4)*0.35, 0.08, 0.02], 'string', 'time', 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
        annotation(fig, 'textbox', [ax3_1.Position(1)-0.05, ax3.Position(2)+ax3.Position(4)*1.05, ax3_1.Position(3)+0.03, 0.02], 'string', 'task-evoked', 'edgecolor', 'none', ...
                'fontsize', fontsize_label, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle') 
        annotation(fig, 'textbox', [ax3_2.Position(1)+0.04, ax3.Position(2)+ax3.Position(4)*1.05, ax3_2.Position(3)+0.05, 0.02], 'string', 'spontaneous', 'edgecolor', 'none', ...
                'fontsize', fontsize_label, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')     
    end
end

ax3_3 = axes('Position', [ax3.Position(1)+ax3.Position(3)*(0.62+0.05*(ii-1)) ax3.Position(2)+ax3.Position(4)*(1.4-1) length_x length_y*0.43]);
parc = dlmread(filename_common_parcellation(parc_name, hemisphere));
data_parc_emp = calc_parcellate(parc, data_resting);
data_parc_emp = calc_normalize_timeseries(data_parc_emp');
data_parc_emp(isnan(data_parc_emp)) = 0;

imagesc(data_parc_emp'*data_parc_emp/T)
caxis([-1 1])
colormap(ax3_3, bluewhitered)
% colorbar
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xtick', [], 'ytick', [])
axis image
xlabel('region', 'fontsize', fontsize_axis)
ylabel('region', 'fontsize', fontsize_axis)
title('FC', 'fontsize', fontsize_axis, 'fontweight', 'normal')

% =========================================================================
% D: reconstruction accuracy
% =========================================================================
ax5 = axes('Position', [0.05 0.1 0.35 0.35]);
hold on;
for jj=1:length(N_interest)
    xline(N_interest(jj), 'k:', 'linewidth', 1.5, 'HandleVisibility', 'off');
end
for ii=1:length(representative_contrasts)
    contrast = representative_contrasts{ii};
    plot(nanmean(task_recon_corr_geometric.(contrast),1), 'color', colors(ii,:), 'linewidth', 2, 'displayname', lower(tasks{ii}))
end
plot(nanmean(resting_recon_corr_geometric,1), 'color', colors(ii+1,:), 'linewidth', 2, 'displayname', lower('REST'))
hold off;
leg = legend('fontsize', fontsize_axis, 'location', 'southeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/2;
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1])
xlabel('number of modes', 'fontsize', fontsize_label)
ylabel('reconstruction accuracy', 'fontsize', fontsize_label)

% =========================================================================
% E: empirical and reconstructed task maps
% =========================================================================
factor_x = 1;
factor_y = 1.02;
init_x = ax5.Position(1)+ax5.Position(3)*1.28;
init_y = ax5.Position(2)*0.5;
length_x = (ax5.Position(3)*0.9)/(factor_x*(length(N_interest)+1-1) + 1);
length_y = (ax5.Position(2)+ax5.Position(4)-init_y)/(factor_y*(length(representative_contrasts)-1) + 1);
% empirical
for ii = 1:length(representative_contrasts)
    contrast = representative_contrasts{ii};
    task_name = lower(tasks{ii});
    
    data_to_plot = data_Figure1.task_map_emp.(contrast);
    
    ax6 = axes('Position', [init_x+factor_x*length_x*(1-1) init_y+factor_y*length_y*(length(representative_contrasts)-ii) length_x length_y]);
    patch(ax6, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    camlight('headlight')
    material dull
    colormap(ax6, bluewhitered)
    axis off
    axis image
    
    annotation(fig, 'textbox', [ax6.Position(1)-0.035, ax6.Position(2)+ax6.Position(4)*0.4, 0.01, 0.01], 'string', task_name, 'edgecolor', 'none', ...
            'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'middle')
    if ii==1
        annotation(fig, 'textbox', [ax6.Position(1)+ax6.Position(3)*0.44, ax6.Position(2)+ax6.Position(4)*0.93, 0.01, 0.01], 'string', 'data', 'edgecolor', 'none', ...
            'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
    end
end

% reconstruction
for ii = 1:length(representative_contrasts)
    contrast = representative_contrasts{ii};
    task_name = lower(tasks{ii});

    for jj = 1:length(N_interest)
        N = N_interest(jj);
        
        data_to_plot = data_Figure1.task_map_recon.(contrast)(:,jj);

        ax6 = axes('Position', [init_x+factor_x*length_x*(jj) init_y+factor_y*length_y*(length(representative_contrasts)-ii) length_x length_y]);
        patch(ax6, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
        if strcmpi(hemisphere, 'lh')
            view([-90 0]);
        elseif strcmpi(hemisphere, 'rh')
            view([90 0]);
        end
        camlight('headlight')
        material dull
        colormap(ax6, bluewhitered)
        axis off
        axis image
        
        if ii==1
            annotation(fig, 'textbox', [ax6.Position(1)+ax6.Position(3)*0.44-0.035, ax6.Position(2)+ax6.Position(4)*0.93, 0.08, 0.01], 'string', {'recon'; sprintf('(N = %i)', N)}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
        end
        
        if ii==1 && jj>=2
            annotation(fig, 'arrow', [ax6.Position(1)+ax6.Position(3)*0.6, ax6.Position(1)+ax6.Position(3)*0.57], ...
                                     [ax6.Position(2)+ax6.Position(4)*0.1, ax6.Position(2)+ax6.Position(4)*0.18], ...
                                     'LineStyle', 'none', 'HeadLength', 6, 'HeadWidth', 6, 'Color', 'k')
            annotation(fig, 'arrow', [ax6.Position(1)+ax6.Position(3)*0.87, ax6.Position(1)+ax6.Position(3)*0.84], ...
                                     [ax6.Position(2)+ax6.Position(4)*0.24, ax6.Position(2)+ax6.Position(4)*0.32], ...
                                     'LineStyle', 'none', 'HeadLength', 6, 'HeadWidth', 6, 'Color', 'k')                     
        elseif ii==2 && jj>=2
            annotation(fig, 'arrow', [ax6.Position(1)+ax6.Position(3)*0.82, ax6.Position(1)+ax6.Position(3)*0.79], ...
                                     [ax6.Position(2)+ax6.Position(4)*0.8, ax6.Position(2)+ax6.Position(4)*0.72], ...
                                     'LineStyle', 'none', 'HeadLength', 6, 'HeadWidth', 6, 'Color', 'k') 
        end
    end
end

% arrows from D to E and surface reconstructions
yloc_min = 0.1;
yloc_space = 0.2;
for jj = 1:length(N_interest)
    N = N_interest(jj);
    annotation(fig, 'line', (ax5.Position(1)+ax5.Position(3)*(N/num_modes))*ones(1,2), ...
                             [init_y*0.75, init_y*(yloc_min+yloc_space*(jj-1))], ...
                             'Color', 'k', 'LineStyle', ':', 'linewidth', 1.5) 
    annotation(fig, 'line', [ax5.Position(1)+ax5.Position(3)*(N/num_modes), init_x+factor_x*length_x*(jj)+length_x*0.5], ...
                             (init_y*(yloc_min+yloc_space*(jj-1)))*ones(1,2), ...
                             'Color', 'k', 'LineStyle', ':', 'linewidth', 1.5) 
    annotation(fig, 'arrow', init_x+factor_x*length_x*(jj)+length_x*0.5*ones(1,2), ...
                             [init_y*(yloc_min+yloc_space*(jj-1)), init_y*0.99], ...
                             'HeadLength', 6, 'HeadWidth', 6, 'Color', 'k', 'LineStyle', ':', 'linewidth', 1.5) 
    
    basis = basis_geometric(:,1:N);
    
    coeffs = calc_eigendecomposition(surface_to_plot.vertices, basis);
    new_vertices = basis*coeffs;
    surface_reconstruction = surface_to_plot;
    surface_reconstruction.vertices = new_vertices;
    
    inset_length = 0.05;
    ax5_inset = axes('Position', [ax5.Position(1)+ax5.Position(3)*(N/num_modes)-inset_length/2 ax5.Position(2)+ax5.Position(4) inset_length inset_length]);
    patch(ax5_inset, 'Vertices', surface_reconstruction.vertices, 'Faces', surface_reconstruction.faces, 'FaceVertexCData', zeros(size(surface_reconstruction.vertices,1),1), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    camlight('headlight')
    material dull
    colormap(ax5_inset, 0.7*ones(1,3))
    axis off
    axis image                     
end

% =========================================================================
% F: FC
% =========================================================================
extra_diagonal = 7;
factor_y = 1.6;
init_x = ax6.Position(1)+ax6.Position(3)*0.45;
init_y = init_y*0.95;
length_x = ax5.Position(3)*0.9;
length_y = (ax5.Position(2)+ax5.Position(4)-init_y)/(factor_y*(length(N_interest)-1) + 1);
for jj = 1:length(N_interest)
    N = N_interest(jj);
        
    point_start = 1-0.5; point_end = num_parcels+extra_diagonal+0.5;
    
    ax7 = axes('Position', [init_x init_y+factor_y*length_y*(length(N_interest)-jj) length_x length_y]);
    imagesc(data_Figure1.FC_combined(:,:,jj))
    hold on;
    plot([point_start+1+extra_diagonal point_end], point_start*ones(1,2), 'k-', 'linewidth', 1)
    plot(point_end*ones(1,2), [point_start point_end-1-extra_diagonal], 'k-', 'linewidth', 1)
    plot([point_start point_end-1-extra_diagonal], point_end*ones(1,2), 'k-', 'linewidth', 1)
    plot(point_start*ones(1,2), [point_start+1+extra_diagonal point_end], 'k-', 'linewidth', 1)
    for ii=1:num_parcels-1
        plot((point_start+1+extra_diagonal+(ii-1))*ones(1,2), point_start+[0 1]+(ii-1), 'k-', 'linewidth', 1)
        plot(point_start+1+extra_diagonal+[0 1]+(ii-1), (point_start+1+(ii-1))*ones(1,2), 'k-', 'linewidth', 1)

        plot((point_start+1+(ii-1))*ones(1,2), point_start+1+extra_diagonal+[0 1]+(ii-1), 'k-', 'linewidth', 1)
        plot(point_start+[0 1]+(ii-1), (point_start+1+extra_diagonal+(ii-1))*ones(1,2), 'k-', 'linewidth', 1)
    end
    hold off;
    caxis([-1 1])
    colormap(ax7, bluewhitered)
    cbar = colorbar;
    set(ax7, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02])
    ylabel(cbar, 'FC', 'fontsize', fontsize_label)
    axis image
    axis off
    
    text(-5, (num_parcels+extra_diagonal)/2, 'data', 'rotation', 90, 'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
    annotation(fig, 'textbox', [ax7.Position(1)+ax7.Position(3)*0.35, ax7.Position(2)+ax7.Position(4)*0.95, 0.08, 0.01], 'string', {'recon'; sprintf('(N = %i)', N)}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
end

%%% panel letters
annotation(fig, 'textbox', [0.01, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.01, 0.73, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.52, 0.73, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.01, 0.5, 0.01, 0.01], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.47, 0.5, 0.01, 0.01], 'string', 'E', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.84, 0.5, 0.01, 0.01], 'string', 'F', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% FIGURE 2

% load all data relevant to Figure 2
data_Figure2 = load(sprintf('%s/Figure2.mat', data_figures_folder));

surface_to_plot = surface_midthickness;
points_ind = data_Figure2.points_ind;
parcel_points = data_Figure2.parcel_points;

colors_schematic = lines(3);

mode_list = [1,2,3,4,20];
num_modes_to_plot = length(mode_list);

fig = figure('Position', [200 200 1000 600]);

% =========================================================================
% A: schematic
% =========================================================================
factor_x = 1.02;
factor_y = 0.6;
init_x = 0.02;
init_y = 0.57;
length_x = (0.89 - 1.*init_x)/(factor_x*(num_modes_to_plot+2-1) + 1);
length_y = (0.98 - init_y)/(factor_y*(2-1) + 1);

ax1_1 = axes('Position', [0.03+init_x+factor_x*length_x*(1-1) init_y+factor_y*length_y*(1.4-1) length_x*1.2 length_y*1.2]);
data_to_plot = zeros(size(surface_to_plot.vertices,1),1);    
patch(ax1_1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
       'EdgeColor', 'none', 'FaceColor', 'interp');
if strcmpi(hemisphere, 'lh')
    view([-90 0]);
elseif strcmpi(hemisphere, 'rh')
    view([90 0]);
end
camlight('headlight')
material dull
colormap(ax1_1, 0.8*ones(1,3))
axis off
axis image

for basis_ind=1:4
    if basis_ind==1
        connectome = data_Figure2.network_geometric;
        
        basis_name = {'geometric'};
        
        ax1_2_positions = [0.04 ax1_1.Position(2)+ax1_1.Position(4)*0.81 length_x*0.55 length_y*0.3];
        box1_positions = [ax1_1.Position(1)+ax1_1.Position(3)*0.25 ax1_1.Position(2)+ax1_1.Position(4)*0.55 ax1_1.Position(3)*0.1 ax1_1.Position(3)*0.15];
    elseif basis_ind==2
        connectome = data_Figure2.network_EDR;
        
        basis_name = {'EDR'};
        
        ax1_2_positions = [0.15 ax1_1.Position(2)+ax1_1.Position(4)*0.81 ax1_2_positions(3) ax1_2_positions(4)];
        box1_positions = [ax1_1.Position(1)+ax1_1.Position(3)*0.65 ax1_1.Position(2)+ax1_1.Position(4)*0.6 ax1_1.Position(3)*0.1 ax1_1.Position(3)*0.15];
    elseif basis_ind==3
        connectome = data_Figure2.network_connectome_1;
        
        basis_name = 'connectome';
        
        ax1_2_positions = [0.057 ax1_1.Position(2)-ax1_1.Position(4)*0.1 ax1_2_positions(3) ax1_2_positions(4)];
        box1_positions = [ax1_1.Position(1)+ax1_1.Position(3)*0.35 ax1_1.Position(2)+ax1_1.Position(4)*0.3 ax1_1.Position(3)*0.1 ax1_1.Position(3)*0.15];
    elseif basis_ind==4
        connectome = data_Figure2.network_connectome_2;
        
        ax1_2_positions = [0.13 ax1_1.Position(2)-ax1_1.Position(4)*0.1 ax1_2_positions(3) ax1_2_positions(4)];
        box1_positions = [ax1_1.Position(1)+ax1_1.Position(3)*0.65 ax1_1.Position(2)+ax1_1.Position(4)*0.4 ax1_1.Position(3)*0.1 ax1_1.Position(3)*0.15];    
    end
    
    [ind_y, ind_x] = find(triu(connectome>0));
    
    ax1_2 = axes('Position', ax1_2_positions);
    hold on;
    for ii=1:size(ind_x,1)
        if connectome(ind_y(ii), ind_x(ii))==1
            what_color = colors_schematic(1,:);
            linewidth = 1.5;
        elseif connectome(ind_y(ii), ind_x(ii))==2
            what_color = colors_schematic(2,:);
            linewidth = 1.5;
        elseif connectome(ind_y(ii), ind_x(ii))==3
            what_color = 'm';%colors_schematic(3,:);
            linewidth = 2.5;
        end
        plot3(parcel_points([ind_y(ii), ind_x(ii)],1), parcel_points([ind_y(ii), ind_x(ii)],2), parcel_points([ind_y(ii), ind_x(ii)],3), '-', 'linewidth', linewidth, 'color', what_color)
        
    end
    for ii=1:size(parcel_points,1)
        plot3(parcel_points(ii,1), parcel_points(ii,2), parcel_points(ii,3), 'ko', 'markersize', 2, 'markerfacecolor', 'k')
    end
    hold off;
    view([-90 0])
    set(ax1_2, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'ylim', [-68 -62], 'zlim', [22 28], 'xtick', [], 'ytick', [], 'ztick', [], 'linewidth', 2)
    axis off
    axis square    
    
    box1 = annotation('rectangle', box1_positions , ...
            'color', 'k', 'linewidth', 2);
    box2 = annotation('rectangle', [ax1_2_positions(1)+ax1_2_positions(3)*0.08 ax1_2_positions(2)-ax1_2_positions(4)*0.08 ax1_2_positions(3)*0.84 ax1_2_positions(4)+ax1_2_positions(4)*0.16], ...
                'color', 'k', 'linewidth', 2);
    if basis_ind==1 || basis_ind==2
        annotation('line', [box1.Position(1) box2.Position(1)], [box1.Position(2)+box1.Position(4) box2.Position(2)], ...  
                   'color', 'k', 'linewidth', 2);
        annotation('line', [box1.Position(1)+box1.Position(3) box2.Position(1)+box2.Position(3)], [box1.Position(2)+box1.Position(4) box2.Position(2)], ...  
                   'color', 'k', 'linewidth', 2); 
        annotation('textbox', [box2.Position(1), box2.Position(2)+box2.Position(4)*0.95, box2.Position(3), 0.02], ...  
                   'string', basis_name, 'edgecolor', 'none', 'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'bottom');
    elseif basis_ind==3
        annotation('line', [box1.Position(1) box2.Position(1)], [box1.Position(2) box2.Position(2)+box2.Position(4)], ...  
                   'color', 'k', 'linewidth', 2);
        annotation('line', [box1.Position(1)+box1.Position(3) box2.Position(1)+box2.Position(3)], [box1.Position(2) box2.Position(2)+box2.Position(4)], ...  
                   'color', 'k', 'linewidth', 2); 
        annotation('textbox', [box2.Position(1), box2.Position(2)-box2.Position(4)*0.15, box2.Position(3)*2.3, 0.02], ...  
                   'string', basis_name, 'edgecolor', 'none', 'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'top');
        annotation('line', [box2.Position(1)+box2.Position(3)*0.9 box2.Position(1)+box2.Position(3)*1.38], [box2.Position(2)+box2.Position(4)*0.75 box2.Position(2)+box2.Position(4)*0.65], ...  
                   'color', 'm', 'linewidth', 3);
        annotation('line', [box2.Position(1)+box2.Position(3)*0.88 box2.Position(1)+box2.Position(3)*1.38], [box2.Position(2)+box2.Position(4)*0.45 box2.Position(2)+box2.Position(4)*0.5], ...  
                   'color', 'm', 'linewidth', 3); 
        annotation('line', [box2.Position(1)+box2.Position(3)*0.95 box2.Position(1)+box2.Position(3)*1.38], [box2.Position(2)+box2.Position(4)*0.12 box2.Position(2)+box2.Position(4)*0.2], ...  
                   'color', 'm', 'linewidth', 3);       
    elseif basis_ind==4
        annotation('line', [box1.Position(1) box2.Position(1)], [box1.Position(2) box2.Position(2)+box2.Position(4)], ...  
                   'color', 'k', 'linewidth', 2);
        annotation('line', [box1.Position(1)+box1.Position(3) box2.Position(1)+box2.Position(3)], [box1.Position(2) box2.Position(2)+box2.Position(4)], ...  
                   'color', 'k', 'linewidth', 2);            
    end
end

% =========================================================================
% B: eigenmodes
% =========================================================================
for basis_ind=1:2
    
    if basis_ind==2
        data_to_plot_evec = basis_connectome;
        basis_name = {'connectome'};
    elseif basis_ind==1
        data_to_plot_evec = basis_EDR;
        basis_name = {'EDR'};
    end
    
    for mode_ind=1:num_modes_to_plot
        mode = mode_list(mode_ind);
        
        if mode_ind==num_modes_to_plot
            ax2 = axes('Position', [0.07+0.03+init_x+factor_x*length_x*(mode_ind+2-1) init_y+factor_y*length_y*(basis_ind-1) length_x length_y]);
        else
            ax2 = axes('Position', [0.07+init_x+factor_x*length_x*(mode_ind+2-1) init_y+factor_y*length_y*(basis_ind-1) length_x length_y]);
        end
        patch(ax2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot_evec(:,mode), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
        if strcmpi(hemisphere, 'lh')
            view([-90 0]);
        elseif strcmpi(hemisphere, 'rh')
            view([90 0]);
        end
        camlight('headlight')
        material dull
        colormap(ax2,bluewhitered)
        axis off
        axis image
        if basis_ind==1
            if mode_ind==1
                annotation(fig, 'textbox', [ax2.Position(1), ax2.Position(2)+ax2.Position(4)*1.4, ax2.Position(3), 0.05], 'string', {'mode 1'; '\psi_1'}, 'edgecolor', 'none', ...
                        'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
            elseif mode_ind==2
                annotation(fig, 'textbox', [ax2.Position(1), ax2.Position(2)+ax2.Position(4)*1.4, ax2.Position(3), 0.05], 'string', {'mode 2'; '\psi_2'}, 'edgecolor', 'none', ...
                        'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
            elseif mode_ind==3
                annotation(fig, 'textbox', [ax2.Position(1), ax2.Position(2)+ax2.Position(4)*1.4, ax2.Position(3), 0.05], 'string', {'mode 3'; '\psi_3'}, 'edgecolor', 'none', ...
                    'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle') 
            elseif mode_ind==4
                annotation(fig, 'textbox', [ax2.Position(1), ax2.Position(2)+ax2.Position(4)*1.4, ax2.Position(3), 0.05], 'string', {'mode 4'; '\psi_4'}, 'edgecolor', 'none', ...
                    'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle') 
            else
                annotation(fig, 'textbox', [ax2.Position(1), ax2.Position(2)+ax2.Position(4)*1.4, ax2.Position(3), 0.05], 'string', {'mode N'; '\psi_N'}, 'edgecolor', 'none', ...
                    'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle') 
            end
        end
        
        if mode_ind==1
            annotation(fig, 'textbox', [ax2.Position(1)-0.08, ax2.Position(2)+ax2.Position(4)*0.45, 0.08, 0.02], 'string', basis_name, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
        end
    end
    
    annotation(fig, 'textbox', [ax2.Position(1)-0.056, ax2.Position(2)+ax2.Position(4)*0.45, 0.08, 0.02], 'string', '...', 'edgecolor', 'none', ...
                'fontsize', fontsize_axis*2, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
end

% =========================================================================
% C: comparison of resting reconstruction accuracy
% =========================================================================
cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];

ax3 = axes('Position', [0.05 0.077 0.35 0.45]);
hold on;
plot(nanmean(resting_recon_corr_geometric,1), 'color', colors2(1,:), 'linewidth', 2, 'displayname', 'geometric')
plot(nanmean(resting_recon_corr_connectome,1), 'color', colors2(2,:), 'linewidth', 2, 'displayname', 'connectome')
plot(nanmean(resting_recon_corr_connectome_density_matched,1), 'color', colors2(3,:), 'linewidth', 2, 'displayname', 'connectome (density matched)')
plot(nanmean(resting_recon_corr_EDR,1), 'color', colors2(4,:), 'linewidth', 2, 'displayname', 'EDR')
hold off;
leg = legend('fontsize', fontsize_axis, 'location', 'southeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/2;
set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1])
xlabel('number of modes', 'fontsize', fontsize_label)
ylabel('reconstruction accuracy', 'fontsize', fontsize_label)

% =========================================================================
% D: comparison of task reconstruction accuracy
% =========================================================================
factor_x = 1.15;
init_x = ax3.Position(1)+ax3.Position(3)*0.85;
init_y = ax3.Position(2);
length_x = (0.9 - init_x)/(factor_x*(3+1-1) + 1);
length_y = ax3.Position(4);
yticks = [];
for ii=1:length(tasks)
    yticks = [yticks, mean(task_contrasts_ind.(tasks{ii}))];
end

data_1 = task_recon_corr_geometric;
recon_corr_combined = zeros(length(contrasts), num_modes, 3);
clims_combined = zeros(3,2);
for kk=1:3
    if kk==1
        data_2 = task_recon_corr_connectome;
    elseif kk==2
        data_2 = task_recon_corr_connectome_density_matched;
    elseif kk==3
        data_2 = task_recon_corr_EDR;
    end

    recon_corr = [];
    for ii=1:length(contrasts)
        contrast = contrasts{ii};

        temp_1 = nanmean(data_1.(contrast),1);
        temp_2 = nanmean(data_2.(contrast),1);
        recon_corr = cat(1, recon_corr, temp_1-temp_2);
    end
    recon_corr_combined(:,:,kk) = recon_corr;
    
    clim_max = max(recon_corr(:,2:end),[],'all');
    clim_min = min(recon_corr(:,2:end),[],'all');
    clims_combined(kk,:) = [clim_min, clim_max];
end

for kk=1:3
    if kk==1
        title_text = {'connectome'};
    elseif kk==2
        title_text = {'connectome'; '(density matched)'};
    elseif kk==3
        title_text = {'EDR'};
    end

    data_to_plot = recon_corr_combined(:,:,kk);
%     clim = clims_combined(kk,:);
    clim = [min(clims_combined(:,1)), max(clims_combined(:,2))]; % make absolute colorbar for all three eigenmodes

    ax4 = axes('Position', [init_x+factor_x*length_x*(kk) init_y length_x length_y]);
    imagesc(data_to_plot)
    for ii=1:length(contrasts)-1
        yline(ii+0.5, 'k-');
    end
    for ii=1:length(tasks)-1
        yline(task_contrasts_ind.(tasks{ii})(end)+0.5, 'k-', 'linewidth', 2);
    end
    caxis(clim)
    colormap(ax4, bluewhitered)
    % cbar = colorbar;
    set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ytick', yticks, 'yticklabel', lower(tasks), 'ticklabelinterpreter', 'none')
    if kk~=1
        set(ax4, 'yticklabel', {})
    end
    if kk==2
        xlabel('number of modes', 'fontsize', fontsize_label)
    end
    title(title_text, 'fontweight', 'normal', 'fontsize', fontsize_axis)
end

cbar = colorbar;
set(cbar, 'ticklength', 0.02, 'Position', [ax4.Position(1)+ax4.Position(3)*1.1 ax4.Position(2) 0.01 ax4.Position(4)], ...
    'fontsize', fontsize_axis-2);
ylabel(cbar, {'reconstruction accuracy difference'; ['(geometric ' char(8211) ' others)']}, 'fontsize', fontsize_label)

%%% panel letters
annotation(fig, 'textbox', [0.01, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.3, 0.99, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.01, 0.57, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.45, 0.57, 0.01, 0.01], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
    
%% FIGURE 3

% load all data relevant to Figure 3
data_Figure3 = load(sprintf('%s/Figure3.mat', data_figures_folder));

surface_to_plot = surface_midthickness;

cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];

fig = figure('Position', [200 200 1000 800]);

% =========================================================================
% A: power spectrum
% =========================================================================    

xlims = [2 num_modes+2];
ylims = [7.5e-4 1e-1];
N_interest = [50, 100, 200];
wavelengths = [60, 40, 30];

factor_x = 1.2;
init_x = 0.07;
init_y = 0.65;
length_x = (0.98 - init_x)/(factor_x*(2-1) + 1);
length_y = (0.94 - init_y)/(factor_y*(1-1) + 1);

% =========================================================================
% A left: power spectrum (HCP)
% =========================================================================   
data_to_plot = nanmean(spectrum_HCP,2);
ax1_1 = axes('Position', [init_x+factor_x*length_x*(1-1) init_y length_x length_y]);
hold on;
bar(1:num_modes, data_to_plot, 'displayname', 'HCP')
for jj=1:length(N_interest)
    if jj==1
        xline(N_interest(jj), 'k:', 'linewidth', 2, 'label', {'wavelength'; sprintf('~%i mm', wavelengths(jj))}, 'labelhorizontalalignment', 'left', 'HandleVisibility', 'off');
    else
        xline(N_interest(jj), 'k:', 'linewidth', 2, 'label', sprintf('~%i mm', wavelengths(jj)), 'labelhorizontalalignment', 'left', 'HandleVisibility', 'off');
    end
end
hold off;
leg = legend('fontsize', fontsize_label, 'location', 'east', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/2;
set(ax1_1, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', xlims, 'ylim', ylims, 'yscale', 'log')
xlabel('mode', 'fontsize', fontsize_label)
ylabel_h = ylabel('normalized power (log scale)', 'fontsize', fontsize_label);

% surface reconstructions
for jj = 1:length(N_interest)
    N = N_interest(jj);
    
    basis = basis_geometric(:,1:N);
    
    coeffs = calc_eigendecomposition(surface_to_plot.vertices, basis);
    new_vertices = basis*coeffs;
    surface_reconstruction = surface_to_plot;
    surface_reconstruction.vertices = new_vertices;
    
    inset_length = 0.05;
    ax1_1_inset = axes('Position', [ax1_1.Position(1)+ax1_1.Position(3)*(N/(xlims(2)+3))-inset_length/2 ax1_1.Position(2)+ax1_1.Position(4) inset_length inset_length]);
    patch(ax1_1_inset, 'Vertices', surface_reconstruction.vertices, 'Faces', surface_reconstruction.faces, 'FaceVertexCData', zeros(size(surface_reconstruction.vertices,1),1), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    camlight('headlight')
    material dull
    colormap(ax1_1_inset, 0.7*ones(1,3))
    axis off
    axis image                     
end

% =========================================================================
% A right: power spectrum (NeuroVault)
% =========================================================================   
data_to_plot = nanmean(spectrum_neurovault,2);
ax1_2 = axes('Position', [init_x+factor_x*length_x*(2-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
hold on;
bar(1:num_modes, data_to_plot, 'displayname', 'NeuroVault')
for jj=1:length(N_interest)
    if jj==1
        xline(N_interest(jj), 'k:', 'linewidth', 2, 'label', {'wavelength'; sprintf('~%i mm', wavelengths(jj))}, 'labelhorizontalalignment', 'left', 'HandleVisibility', 'off');
    else
        xline(N_interest(jj), 'k:', 'linewidth', 2, 'label', sprintf('~%i mm', wavelengths(jj)), 'labelhorizontalalignment', 'left', 'HandleVisibility', 'off');
    end
end
hold off;
% title('HCP', 'fontsize', fontsize_axis, 'fontweight', 'normal')
leg = legend('fontsize', fontsize_label, 'location', 'east', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/2;
% set(leg, 'visible', 'off')
set(ax1_2, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', xlims, 'ylim', ylims, 'yscale', 'log')
xlabel('mode', 'fontsize', fontsize_label)
% ylabel_h = ylabel('normalized power (log scale)', 'fontsize', fontsize_label);
% set(ylabel_h, 'Position', ylabel_h.Position+[0 0.4 0])

% surface reconstructions
for jj = 1:length(N_interest)
    N = N_interest(jj);
    
    basis = basis_geometric(:,1:N);
    
    coeffs = calc_eigendecomposition(surface_to_plot.vertices, basis);
    new_vertices = basis*coeffs;
    surface_reconstruction = surface_to_plot;
    surface_reconstruction.vertices = new_vertices;
    
    inset_length = 0.05;
    ax1_2_inset = axes('Position', [ax1_2.Position(1)+ax1_2.Position(3)*(N/(xlims(2)+3))-inset_length/2 ax1_2.Position(2)+ax1_2.Position(4) inset_length inset_length]);
    patch(ax1_2_inset, 'Vertices', surface_reconstruction.vertices, 'Faces', surface_reconstruction.faces, 'FaceVertexCData', zeros(size(surface_reconstruction.vertices,1),1), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    camlight('headlight')
    material dull
    colormap(ax1_2_inset, 0.7*ones(1,3))
    axis off
    axis image                     
end

% =========================================================================
% B: removing modes reconstruction accuracy
% =========================================================================  
factor_x = 1.3;
factor_y = 1.3;
init_y_1 = init_y-0.1;
init_y = 0.06;
length_x = (0.98 - init_x)/(factor_x*(4-1) + 1);
length_y = (init_y_1 - init_y)/(factor_y*(2-1) + 1);

N_interest = 51;
N_percent = ((N_interest-1)/num_modes)*100;

for ii=1:length(representative_contrasts)
    contrast = representative_contrasts{ii};
    task = lower(tasks{ii});
    
    data_to_plot_remove_highf = nanmean(task_recon_corr_remove_highf_geometric.(contrast),1);
    data_to_plot_remove_lowf = nanmean(task_recon_corr_remove_lowf_geometric.(contrast),1);
    
    if ii<=4
        ax4 = axes('Position', [init_x+factor_x*length_x*(ii-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax4 = axes('Position', [init_x+factor_x*length_x*(ii-4-1+0.25) init_y+factor_y*length_y*(1-1) length_x length_y]);
    end
    hold on;
    xline(N_percent, 'k:', 'HandleVisibility', 'off');
    if ii==length(representative_contrasts)
        plot(100*(0:(num_modes-1))/num_modes, data_to_plot_remove_lowf, '-', 'color', 'k', 'linewidth', 2, 'displayname', ['removing top long-' newline 'wavelength modes'])
        plot(100*(0:(num_modes-1))/num_modes, data_to_plot_remove_highf, '--', 'color', 'k', 'linewidth', 2, 'displayname', ['removing top short-' newline 'wavelength modes'])
    end
    plot(100*(0:(num_modes-1))/num_modes, data_to_plot_remove_lowf, '-', 'color', colors(ii,:), 'linewidth', 2, 'HandleVisibility', 'off')
    plot(100*(0:(num_modes-1))/num_modes, data_to_plot_remove_highf, '--', 'color', colors(ii,:), 'linewidth', 2, 'HandleVisibility', 'off')
    plot(N_percent, data_to_plot_remove_lowf(N_interest), 'k.', 'markersize', 25, 'HandleVisibility', 'off')
    plot(N_percent, data_to_plot_remove_highf(N_interest), 'k.', 'markersize', 25, 'HandleVisibility', 'off')
    hold off;
    if ii==length(representative_contrasts)
        leg = legend('fontsize', fontsize_axis, 'interpreter', 'none', 'box', 'off', 'numcolumns', 1, ...
                    'Position', [ax4.Position(1)+ax4.Position(3)*1.3 ax4.Position(2)+ax4.Position(4)*0.25 ax4.Position(3)*0.5 ax4.Position(3)*0.5]);
        leg.ItemTokenSize = leg.ItemTokenSize/1.4;
    end
    set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [0 100], 'ylim', [0 1], 'xtick', 0:25:100)
    if ii>4
        xlabel('percent of modes removed', 'fontsize', fontsize_label)
    end
    if ii==1 ||ii==5
        ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
    end
    title(task, 'fontsize', fontsize_axis, 'fontweight', 'normal', 'interpreter', 'none')
    
    surface_length = 0.05;
    %%% reconstructed map (removing high frequency modes)
    data_to_plot = data_Figure3.task_map_recon_remove_highf.(contrast);

    ax4_2 = axes('Position', [ax4.Position(1)+ax4.Position(3)*(N_percent/100 + 0.05) ax4.Position(2)+ax4.Position(4)*(data_to_plot_remove_highf(N_interest)-0.35) surface_length surface_length]);
    patch(ax4_2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    camlight('headlight')
    material dull
    colormap(ax4_2, bluewhitered)
    axis off
    axis image
    
    annotation(fig, 'arrow', [ax4.Position(1)+ax4.Position(3)*N_percent/100, ax4_2.Position(1)+ax4_2.Position(3)*0.2], ...
                             [ax4.Position(2)+ax4.Position(4)*(data_to_plot_remove_highf(N_interest)), ax4_2.Position(2)+ax4_2.Position(4)*0.8], ...
                             'LineStyle', '-', 'HeadLength', 6, 'HeadWidth', 6, 'Color', 'k') 
    
    %%% reconstructed map (removing low frequency modes)
    data_to_plot = data_Figure3.task_map_recon_remove_lowf.(contrast);

    ax4_3 = axes('Position', [ax4.Position(1)+ax4.Position(3)*0.05 ax4.Position(2)+ax4.Position(4)*(data_to_plot_remove_lowf(N_interest)-0.4) surface_length surface_length]);
    patch(ax4_3, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    camlight('headlight')
    material dull
    colormap(ax4_3, bluewhitered)
    axis off
    axis image
    
    annotation(fig, 'arrow', [ax4.Position(1)+ax4.Position(3)*N_percent/100, ax4_3.Position(1)+ax4_3.Position(3)*0.5], ...
                             [ax4.Position(2)+ax4.Position(4)*(data_to_plot_remove_lowf(N_interest)), ax4_3.Position(2)+ax4_3.Position(4)*0.9], ...
                             'LineStyle', '-', 'HeadLength', 6, 'HeadWidth', 6, 'Color', 'k')
                         
    %%% empirical map
    data_to_plot = data_Figure3.task_map_emp.(contrast);

    ax4_4 = axes('Position', [ax4.Position(1)+ax4.Position(3)-surface_length ax4.Position(2)+ax4.Position(4)-surface_length*0.8 surface_length surface_length]);
    patch(ax4_4, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    camlight('headlight')
    material dull
    colormap(ax4_4, bluewhitered)
    axis off
    axis image        
    title('data', 'fontsize', fontsize_axis-2, 'fontweight', 'normal', 'interpreter', 'none')
end

%%% panel letters
annotation(fig, 'textbox', [0.015, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.015, 0.6, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% FIGURE 4

% load all data relevant to Figure 4
% data_Figure4 = load(sprintf('%s/Figure4.mat', data_figures_folder));

surface_to_plot = surface_midthickness;

parc_name = 'Glasser360';
parc = dlmread(filename_common_parcellation(parc_name, hemisphere));
parcels = unique(parc(parc>0));
num_parcels = length(unique(parc(parc>0)));
triu_ind = calc_triu_ind(zeros(num_parcels, num_parcels));

cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];
colors_models = lines(2);

factor_x = 1.5;
factor_y = 1.4;
init_x = 0.06;
init_y = 0.07;
length_x = (1 - 1.*init_x)/(factor_x*(6-1) + 1);
length_y = (1 - 1.5*init_y)/(factor_y*(3-1) + 1);

fig = figure('Position', [200 200 1200 700]);

% =========================================================================
% A: model schematic
% =========================================================================
% wave model
ax1_1 = axes('Position', [init_x+factor_x*length_x*(1-1)+length_x*0.25 init_y+factor_y*length_y*(3.0-1)+length_y/2 length_x length_y/2]);
num_downsampled_vertices = 500;
surface_downsampled = gifti(sprintf('%s/fsLR_%i_midthickness-%s.surf.gii', data_template_surfaces_folder, num_downsampled_vertices, hemisphere));
data_to_plot = ones(size(surface_downsampled.vertices,1),1);
    
patch(ax1_1, 'Vertices', surface_downsampled.vertices, 'Faces', surface_downsampled.faces, 'FaceVertexCData', data_to_plot, ...
           'EdgeColor', 'k', 'FaceColor', 'interp', 'FaceLighting', 'gouraud', 'linewidth', 0.5);
if strcmpi(hemisphere, 'lh')
    view([-90 0]);
elseif strcmpi(hemisphere, 'rh')
    view([90 0]);
end
colormap(ax1_1, 0.7*ones(1,3))
axis off
axis image

annotation(fig, 'textbox', [ax1_1.Position(1)-0.06, ax1_1.Position(2)+ax1_1.Position(4)*0.45, 0.08, 0.02], 'string', {'wave'; 'model'}, 'edgecolor', 'none', ...
                'fontsize', fontsize_label, 'fontweight', 'bold', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'color', colors_models(1,:))
            
ax1_2 = axes('Position', [init_x+factor_x*length_x*(1.87-1) ax1_1.Position(2) 2.2*length_x length_y/2]);
ystart = 0.9;
yspacing = 0.25;
text(0, ystart, 'For each vertex at location $\mathbf{r}$:', ...
    'fontsize', fontsize_label, 'horizontalalignment', 'left', 'verticalalignment', 'middle', 'interpreter', 'latex')
text(0.5, ystart-yspacing*1.6, '$\displaystyle\left[\frac{1}{\gamma_s^2}\frac{\partial^2}{\partial t^2} + \frac{2}{\gamma_s}\frac{\partial}{\partial t} + 1 - r_s^2\nabla^2\right] \phi(\mathbf{r},t) = Q(\mathbf{r},t)$', ...
    'fontsize', fontsize_label, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex')
text(0.5, ystart-yspacing*1.6*2, '$M_{\rm fixed} = 1$; $M_{\rm free} = 1$', ...
    'fontsize', fontsize_label, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex')
axis off

% neural mass model
ax2_1 = axes('Position', [ax1_1.Position(1) init_y+factor_y*length_y*(2.93-1) length_x length_y/2]);
data_to_plot = ones(size(surface_to_plot.vertices,1),1);

parcel_points = [min(surface_to_plot.vertices(:,1)), -48, 50;
                 min(surface_to_plot.vertices(:,1)), -62, 10;
                 min(surface_to_plot.vertices(:,1)), -22, 19;
                 min(surface_to_plot.vertices(:,1)), 22, 43;
                 min(surface_to_plot.vertices(:,1)), 22, -5];
adjacency = zeros(size(parcel_points,1));
adjacency(1, [4]) = 1;
adjacency(2, 3) = 1;
adjacency(3, [4,5]) = 1;
adjacency = triu(adjacency,1)+triu(adjacency,1)';
% rng(1)
adjacency = rand(size(parcel_points,1)).*adjacency;
[ind_y, ind_x] = find(triu(adjacency>0));

patch(ax2_1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
           'EdgeColor', 'none', 'FaceColor', 'interp');   
hold on;
for ii=1:size(ind_x,1)
%     plot3(parcel_points([ind_y(ii), ind_x(ii)],1), parcel_points([ind_y(ii), ind_x(ii)],2), parcel_points([ind_y(ii), ind_x(ii)],3), 'k-', 'linewidth', 3*adjacency(ind_y(ii), ind_x(ii)))
    plot3(parcel_points([ind_y(ii), ind_x(ii)],1), parcel_points([ind_y(ii), ind_x(ii)],2), parcel_points([ind_y(ii), ind_x(ii)],3), 'k-', 'linewidth', 3*ii/size(ind_x,1))
end
for ii=1:size(parcel_points,1)
    plot3(parcel_points(ii,1), parcel_points(ii,2), parcel_points(ii,3), 'ko', 'markersize', 9, 'markerfacecolor', 'w')
    if ii==1
%         text(parcel_points(ii,1), parcel_points(ii,2), parcel_points(ii,3), 'i', 'fontsize', fontsize_axis, 'color', 'k', ...
%             'horizontalalignment', 'center', 'verticalalignment', 'middle')
    end
end
hold off;
if strcmpi(hemisphere, 'lh')
    view([-90 0]);
elseif strcmpi(hemisphere, 'rh')
    view([90 0]);
end
camlight('headlight')
material dull
colormap(ax2_1, 0.7*ones(1,3))
axis off
axis image

annotation(fig, 'textbox', [ax2_1.Position(1)-0.06, ax2_1.Position(2)+ax2_1.Position(4)*0.4, 0.08, 0.02], 'string', {'mass'; 'model'}, 'edgecolor', 'none', ...
                'fontsize', fontsize_label, 'fontweight', 'bold', 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'color', colors_models(2,:))

ax2_2 = axes('Position', [ax1_2.Position(1) ax2_1.Position(2) 2*length_x length_y/2]);
ystart = 0.9;
yspacing = 0.25;
text(0, ystart, 'For each region $i$:', ...
    'fontsize', fontsize_label, 'horizontalalignment', 'left', 'verticalalignment', 'middle', 'interpreter', 'latex')
text(0.5, ystart-yspacing*1.6, '$\displaystyle\frac{dS_i}{dt} = f(\mathbf{S}, \theta_i, C, G)$', ...
    'fontsize', fontsize_label, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex')
text(0.5, ystart-yspacing*1.6*2, '$M_{\rm fixed} = 15$; $M_{\rm free} = 4$', ...
    'fontsize', fontsize_label, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'interpreter', 'latex')
axis off

ax2_3 = axes('Position', [init_x+factor_x*length_x*(3.3-1) init_y+factor_y*length_y*(2.98-1) length_x*2.8 length_y]);
FC_example = corr(resting_emp.time_series(1:num_parcels,:,1)');

imagesc(FC_example)
caxis([-1 1])
colormap(ax2_3, bluewhitered)
cbar = colorbar('fontsize', fontsize_axis, 'ytick', [-1 0 1]);
set(ax2_3, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xtick', [], 'ytick', [])
xlabel('region', 'fontsize', fontsize_label)
ylabel('region', 'fontsize', fontsize_label)
ylabel(cbar, 'FC', 'fontsize', fontsize_label)
axis image
title('simulated FC', 'fontsize', fontsize_label, 'fontweight', 'normal')

annotation(fig, 'arrow', ax2_3.Position(1)+[0.02 0.04], ...
                         [ax1_2.Position(2)+ax1_2.Position(4)*0.5 ax2_3.Position(2)+ax2_3.Position(4)*0.6], 'Color', 'k')
annotation(fig, 'arrow', ax2_3.Position(1)+[0.02 0.04], ...
                         [ax2_2.Position(2)+ax2_2.Position(4)*0.5 ax2_3.Position(2)+ax2_3.Position(4)*0.4], 'Color', 'k')
annotation(fig, 'line', [ax1_1.Position(1)-0.05, ax1_2.Position(1)+ax1_2.Position(3)*1.05], ...
                        (ax1_1.Position(2)-0.011)*ones(1,2))
              
% =========================================================================
% B: data vs model
% =========================================================================
for model_ind = 1:2
    data_emp_FC = model_fit.FC_emp;
    data_emp_FCD = model_fit.FCD_emp;
    
    if model_ind==1
        data_model_FC = model_fit.FC_model_wave;
        data_model_FCD = model_fit.FCD_model_wave;
    elseif model_ind==2
        data_model_FC = model_fit.FC_model_mass;
        data_model_FCD = model_fit.FCD_model_mass;
    end

    % FC
    ax2_1 = axes('Position', [init_x+factor_x*length_x*(1-1) init_y+factor_y*length_y*(-0.98.*[model_ind-2]) length_x length_y*0.8]);
    extra_diagonal = 7;
    FC_combined_model = zeros(num_parcels+extra_diagonal);
    FC_combined_model(find(triu(ones(num_parcels+extra_diagonal),1+extra_diagonal))) = data_model_FC(triu_ind);
    FC_combined_emp = zeros(num_parcels+extra_diagonal);
    FC_combined_emp(find(triu(ones(num_parcels+extra_diagonal),1+extra_diagonal))) = data_emp_FC(triu_ind);

    FC_combined = FC_combined_model + FC_combined_emp';

    point_start = 1-0.5; point_end = num_parcels+extra_diagonal+0.5;

    imagesc(FC_combined)
    hold on;
    plot([point_start+1+extra_diagonal point_end], point_start*ones(1,2), 'k-', 'linewidth', 1)
    plot(point_end*ones(1,2), [point_start point_end-1-extra_diagonal], 'k-', 'linewidth', 1)
    plot([point_start point_end-1-extra_diagonal], point_end*ones(1,2), 'k-', 'linewidth', 1)
    plot(point_start*ones(1,2), [point_start+1+extra_diagonal point_end], 'k-', 'linewidth', 1)
    for ii=1:num_parcels-1
        plot((point_start+1+extra_diagonal+(ii-1))*ones(1,2), point_start+[0 1]+(ii-1), 'k-', 'linewidth', 1)
        plot(point_start+1+extra_diagonal+[0 1]+(ii-1), (point_start+1+(ii-1))*ones(1,2), 'k-', 'linewidth', 1)

        plot((point_start+1+(ii-1))*ones(1,2), point_start+1+extra_diagonal+[0 1]+(ii-1), 'k-', 'linewidth', 1)
        plot(point_start+[0 1]+(ii-1), (point_start+1+extra_diagonal+(ii-1))*ones(1,2), 'k-', 'linewidth', 1)
    end
    hold off;
    caxis([-1 1])
    colormap(ax2_1, bluewhitered)
    set(ax2_1, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02])
    axis off
    text(-5, (num_parcels+extra_diagonal)/2, 'data FC', 'rotation', 90, 'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
    annotation(fig, 'textbox', [ax2_1.Position(1)+ax2_1.Position(3)*0.12, ax2_1.Position(2)+ax2_1.Position(4)*0.96, 0.08, 0.02], 'string', 'model FC', 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
    
    if model_ind==1
        t = text(25, -45, 'wave model', 'fontsize', fontsize_label, 'fontweight', 'bold', 'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
            'color', colors_models(1,:));
        text(t.Position(1)+135, t.Position(2), {'($M_{\rm free} = 1$)'}, 'fontsize', fontsize_label, 'fontweight', 'bold', 'horizontalalignment', 'center', ...
            'verticalalignment', 'middle', 'interpreter', 'latex', 'color', colors_models(1,:))
    elseif model_ind==2
        t = text(25, -45, 'mass model', 'fontsize', fontsize_label, 'fontweight', 'bold', 'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
            'color', colors_models(2,:));
        text(t.Position(1)+135, t.Position(2), {'($M_{\rm free} = 4$)'}, 'fontsize', fontsize_label, 'fontweight', 'bold', 'horizontalalignment', 'center', ...
            'verticalalignment', 'middle', 'interpreter', 'latex', 'color', colors_models(2,:))
    end
    
    % edge FC
    ax2_2 = axes('Position', [init_x+factor_x*length_x*(2-1) ax2_1.Position(2) length_x ax2_1.Position(4)]);
    data_to_plot_x = real(atanh(data_model_FC(triu_ind)));
    data_to_plot_y = real(atanh(data_emp_FC(triu_ind)));
    hold on;
    plot(data_to_plot_x, data_to_plot_y, 'k.', 'markersize', 4)
    plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
        'r-', 'linewidth', 2);
    hold off;
    set(ax2_2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02])
    xlabel('model', 'fontsize', fontsize_axis)
    ylabel('data', 'fontsize', fontsize_axis)
    title('edge FC', 'fontsize', fontsize_label, 'fontweight', 'normal')
    [rho,pval] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
    text(max(get(gca,'xlim')), min(get(gca,'ylim'))+0.1*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], sprintf('r = %.2f', rho), 'color', 'r', ...
        'fontsize', fontsize_axis, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'right');
    box off;
    
    % node FC
    ax2_3 = axes('Position', [init_x+factor_x*length_x*(3-1) ax2_1.Position(2) length_x ax2_1.Position(4)]);
    data_to_plot_x = nanmean(data_model_FC,2);
    data_to_plot_y = nanmean(data_emp_FC,2);
    hold on;
    plot(data_to_plot_x, data_to_plot_y, 'k.', 'markersize', 10)
    plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
        'r-', 'linewidth', 2);
    hold off;
    set(ax2_3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02])
    xlabel('model', 'fontsize', fontsize_axis)
    ylabel('data', 'fontsize', fontsize_axis)
    title('node FC', 'fontsize', fontsize_label, 'fontweight', 'normal')
    [rho,pval] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
    text(max(get(gca,'xlim')), min(get(gca,'ylim'))+0.1*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], sprintf('r = %.2f', rho), 'color', 'r', ...
        'fontsize', fontsize_axis, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'right');
    box off
    
    cb = cbrewer('qual', 'Set3', 12, 'pchip');
    % phase FCD
    ax2_4 = axes('Position', [init_x+factor_x*length_x*(4-1) ax2_1.Position(2) length_x ax2_1.Position(4)]);
    hold on;
    histogram(data_emp_FCD(:), 'normalization', 'pdf', 'facecolor', cb(1,:), 'edgecolor', 'none', 'displayname', 'data', 'facealpha', 0.6)
    histogram(data_model_FCD(:), 'normalization', 'pdf', 'facecolor', cb(4,:), 'edgecolor', 'none', 'displayname', 'model', 'facealpha', 0.6)
    hold off;
    if model_ind==2
        objects = findobj(gca, 'type', 'histogram');
        set(ax2_4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 1], 'children', flipud(get(gca,'children')))
        leg = legend([objects(2) objects(1)], {'data', 'model'}, 'fontsize', fontsize_axis, 'location', 'northeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
    else
        objects = findobj(gca, 'type', 'histogram');
        set(ax2_4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 1])
        leg = legend('fontsize', fontsize_axis, 'location', 'northeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
    end
    leg.ItemTokenSize = leg.ItemTokenSize/2;
    leg_names = get(leg, 'string');
    xlabel('synchrony', 'fontsize', fontsize_axis)
    ylabel('pdf', 'fontsize', fontsize_axis)
    title('FCD', 'fontsize', fontsize_label, 'fontweight', 'normal')
    text(max(get(gca,'xlim')), min(get(gca,'ylim'))+0.1*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], sprintf('KS = %.2f', model_fit.KS(model_ind)), 'color', 'r', ...
        'fontsize', fontsize_axis, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'right');
    box off
    
    if model_ind==1
        annotation(fig, 'line', [ax2_1.Position(1)-0.03, ax2_4.Position(1)+ax2_4.Position(3)], (ax2_1.Position(2)-0.065)*ones(1,2))        
    end
end

% =========================================================================
% C: wave model visual stimulation time series
% =========================================================================
t_interest = [6, 16, 36, 42];
t_interest_ind = dsearchn(model_wave_visual.tspan', t_interest');
data_to_plot = model_wave_visual.simulated_neural_visual_vertex(:,t_interest_ind);
data_to_plot = data_to_plot./repmat(max(data_to_plot,[],1), size(data_to_plot,1), 1);
data_to_plot(isnan(data_to_plot)) = 0;
data_orig = data_to_plot;

surf_factor_x = 1.05;
surf_init_x = init_x+factor_x*length_x*(5-1);
surf_init_y = init_y+factor_y*length_y*(2.97-1);
surf_length_x = 2.5*length_x/(surf_factor_x*(length(t_interest_ind)-1) + 1);
surf_length_y = length_y/1.5;

for t_ind=1:length(t_interest_ind)
    t = t_interest(t_ind);
    
    nan_val = min(data_orig(find(cortex),t_ind));
    data_to_plot(isnan(data_to_plot(:,t_ind)),t_ind) = nan_val;
    data_to_plot(find(~cortex),t_ind) = nan_val-1e-2;
    clims = [min(data_to_plot(:,t_ind)), max(data_to_plot(:,t_ind))];
    if clims(2)<=0
        clims(2) = 0.01;
    end
    
    ax4_1 = axes('Position', [surf_init_x+surf_factor_x*surf_length_x*(t_ind-1) surf_init_y+length_y/2.3 surf_length_x surf_length_y]);
    patch(ax4_1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,t_ind), ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([-90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([90 0]);
    end
    caxis(clims)
    camlight('headlight')
    material dull
    colormap(ax4_1,[0.5 0.5 0.5; parula])
    axis off
    axis image
    obj_title = title(sprintf('t = %i ms', (t)), 'fontsize', fontsize_axis, 'fontweight', 'normal');
    obj_title.Position = obj_title.Position + [0 0 18];
    
    ax4_2 = axes('Position', [surf_init_x+surf_factor_x*surf_length_x*(t_ind-1) surf_init_y surf_length_x surf_length_y]);
    patch(ax4_2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,t_ind), ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    if strcmpi(hemisphere, 'lh')
        view([90 0]);
    elseif strcmpi(hemisphere, 'rh')
        view([-90 0]);
    end
    caxis(clims)
    camlight('headlight')
    material dull
    colormap(ax4_2,[0.5 0.5 0.5; parula])
    axis off
    axis image
    
    if t_ind==1
        annotation(fig, 'arrow', ax4_2.Position(1)+ax4_2.Position(3).*[0.09, 0.01], ...
                                 ax4_2.Position(2)+ax4_2.Position(4).*[0.5, 0.75], ...
                                 'LineStyle', '-', 'HeadLength', 8, 'HeadWidth', 8, 'Color', 'r')
        annotation(fig, 'arrow', ax4_2.Position(1)+ax4_2.Position(3).*[0.1, 0.19], ...
                                 ax4_2.Position(2)+ax4_2.Position(4).*[0.4, 0.15], ...
                                 'LineStyle', '-', 'HeadLength', 8, 'HeadWidth', 8, 'Color', 'r')
    elseif t_ind==2
        annotation(fig, 'arrow', ax4_2.Position(1)+ax4_2.Position(3).*[0.09, 0.01], ...
                                 ax4_2.Position(2)+ax4_2.Position(4).*[0.58, 0.8], ...
                                 'LineStyle', '-', 'HeadLength', 8, 'HeadWidth', 8, 'Color', 'r')
        annotation(fig, 'arrow', ax4_2.Position(1)+ax4_2.Position(3).*[0.15, 0.21], ...
                                 ax4_2.Position(2)+ax4_2.Position(4).*[0.36, 0.14], ...
                                 'LineStyle', '-', 'HeadLength', 8, 'HeadWidth', 8, 'Color', 'r')                     
    elseif t_ind==3
        annotation(fig, 'arrow', ax4_1.Position(1)+ax4_1.Position(3).*[0.91, 0.68], ...
                                 ax4_1.Position(2)+ax4_1.Position(4).*[0.66, 0.81], ...
                                 'LineStyle', '-', 'HeadLength', 8, 'HeadWidth', 8, 'Color', 'r')
        annotation(fig, 'arrow', ax4_1.Position(1)+ax4_1.Position(3).*[0.85, 0.55], ...
                                 ax4_1.Position(2)+ax4_1.Position(4).*[0.32, 0.23], ...
                                 'LineStyle', '-', 'HeadLength', 8, 'HeadWidth', 8, 'Color', 'r')                     
    elseif t_ind==length(t_interest_ind)  
        annotation(fig, 'arrow', ax4_1.Position(1)+ax4_1.Position(3).*[0.8, 0.55], ...
                                 ax4_1.Position(2)+ax4_1.Position(4).*[0.73, 0.85], ...
                                 'LineStyle', '-', 'HeadLength', 8, 'HeadWidth', 8, 'Color', 'r')
        annotation(fig, 'arrow', ax4_1.Position(1)+ax4_1.Position(3).*[0.7, 0.4], ...
                                 ax4_1.Position(2)+ax4_1.Position(4).*[0.28, 0.18], ...
                                 'LineStyle', '-', 'HeadLength', 8, 'HeadWidth', 8, 'Color', 'r') 
                             
        cbar = colorbar(ax4_2,'southoutside');
        set(cbar, 'fontsize', fontsize_axis-4, 'ticklength', 0.02, 'ytick', [min(data_to_plot(:,t_ind)), max(data_to_plot(:,t_ind))], 'yticklabel', {}, ...
            'position', [ax4_2.Position(1)-ax4_2.Position(3)*0.2, ax4_2.Position(2)+0.01, ax4_2.Position(3)*0.6, 0.01])
        ylabel(cbar, 'amplitude', 'fontsize', fontsize_axis)
        annotation(fig, 'textbox', [cbar.Position(1)-0.03, cbar.Position(2)*1, 0.04, 0.01], 'string', 'min', 'edgecolor', 'none', ...
               'fontsize', fontsize_axis-2, 'horizontalalignment', 'center', 'verticalalignment', 'middle')
        annotation(fig, 'textbox', [cbar.Position(1)+cbar.Position(3)-0.008, cbar.Position(2)*1, 0.04, 0.01], 'string', 'max', 'edgecolor', 'none', ...
                   'fontsize', fontsize_axis-2, 'horizontalalignment', 'center', 'verticalalignment', 'middle') 
    end
end

% =========================================================================
% D: hierarchy amplitude time series
% =========================================================================
ROIS = model_wave_visual.ROIS;
data_to_plot = model_wave_visual.simulated_neural_visual_parcel(ROIS,:);
[~, max_ind] = max(data_to_plot, [], 2);
tpeak = model_wave_visual.tspan(max_ind);
[~,tpeak_ind] = sort(tpeak,'ascend');
data_to_plot = data_to_plot + abs(min(data_to_plot(tpeak_ind,:),[],'all'));
data_to_plot = data_to_plot + min(data_to_plot(tpeak_ind,:),[],'all');

ROIS_ordered = ROIS(tpeak_ind);

parc_tpeak = zeros(size(parc,1),1);
parcels = unique(parc(parc>0));
for ii=1:length(ROIS)
    parc_tpeak(parc==parcels(ROIS_ordered(ii))) = ii;
end

colors3 = turbo(length(ROIS)+2);
colors3 = colors3(3:end,:);

ax5 = axes('Position', [surf_init_x+factor_x*length_x*0.25 init_y+factor_y*length_y*(2-1) length_x*2 length_y]);
hold on;
for jj=1:length(ROIS)
    ROI_name = model_wave_visual.ROI_names{ROIS_ordered(jj)};
    ROI_name = ROI_name(3:end);
    ROI_name = ROI_name(1:(end-4));
    
	plot(model_wave_visual.tspan, data_to_plot(tpeak_ind(jj),:), 'Color', colors3(jj,:), 'linewidth', 2, 'displayname', ROI_name)
end
hold off;
set(ax5, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 50])
xlabel('time (ms)', 'fontsize', fontsize_label)
ylabel('amplitude', 'fontsize', fontsize_label)

data_to_plot = parc_tpeak;
data_to_plot(find(~cortex)) = min(data_to_plot)-1;
clims = [min(data_to_plot), max(data_to_plot)];

% surface inset
boundary_method = 'midpoint';
BOUNDARY = findROIboundaries(surface_to_plot.vertices,surface_to_plot.faces,parc_tpeak,boundary_method);
    
ax5_1 = axes('Position', [ax5.Position(1)+ax5.Position(3)*0.5 ax5.Position(2)+ax5.Position(4)*0.55 ax5.Position(3)*0.55 ax5.Position(4)*0.45]);
patch(ax5_1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
           'EdgeColor', 'none', 'FaceColor', 'interp');
hold on;
for ii = 1:length(BOUNDARY)
    plot3(BOUNDARY{ii}(:,1), BOUNDARY{ii}(:,2), BOUNDARY{ii}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
end
hold off;
if strcmpi(hemisphere, 'lh')
    view([-90 0]);
elseif strcmpi(hemisphere, 'rh')
    view([90 0]);
end
caxis(clims)
camlight('headlight')
material dull
colormap(ax5_1,[0.5 0.5 0.5; 0.8 0.8 0.8; colors3])
axis off
axis image

ax5_2 = axes('Position', [ax5_1.Position(1) ax5_1.Position(2)-ax5_1.Position(4)*1.02 ax5_1.Position(3) ax5_1.Position(4)]);
patch(ax5_2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
           'EdgeColor', 'none', 'FaceColor', 'interp');
hold on;
for ii = 1:length(BOUNDARY)
    plot3(BOUNDARY{ii}(:,1), BOUNDARY{ii}(:,2), BOUNDARY{ii}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
end
hold off;
if strcmpi(hemisphere, 'lh')
    view([90 0]);
elseif strcmpi(hemisphere, 'rh')
    view([-90 0]);
end
caxis(clims)
camlight('headlight')
material dull
colormap(ax5_2,[0.5 0.5 0.5; 0.8 0.8 0.8; colors3])
axis off
axis image

% =========================================================================
% E: time to peak vs myelin
% =========================================================================
ax6 = axes('Position', [surf_init_x+factor_x*length_x*0.25 init_y+factor_y*length_y*(1-1) ax5.Position(3) ax5.Position(4)]);
data_to_plot_x = tiedrank(data_myelin.myelin(ROIS));
data_to_plot_y = tiedrank(tpeak');
hold on;
plot(data_to_plot_x, data_to_plot_y, 'k.', 'markersize', 25)
plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
    'r-', 'linewidth', 2);
hold off;
set(ax6, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [min(data_to_plot_x) max(data_to_plot_x)], 'ylim', [min(data_to_plot_y) max(data_to_plot_y)])
xlabel('T1w:T2w (rank)', 'fontsize', fontsize_label)
ylabel('time to peak (rank)', 'fontsize', fontsize_label)

stat = 'linear';
if strcmpi(stat, 'rank')
    [rho,~] = corr(data_to_plot_x, data_to_plot_y, 'type', 'spearman');
elseif strcmpi(stat, 'linear')
    [rho,~] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
end

temp = NaN(size(perm_id));
for ii=1:length(ROIS)
    ROI = ROIS(ii);
    
    ROI_ind = find(perm_id==ROI);
    temp(ROI_ind) = ii;
end

ROI_perm_id = zeros(length(ROIS), size(perm_id,2));
for jj=1:size(perm_id,2)
    temp2 = temp(:,jj);
    temp2(isnan(temp2)) = [];
    
    ROI_perm_id(:,jj) = temp2;
end

pval = perm_sphere_p(data_to_plot_x, data_to_plot_y, ROI_perm_id, 'pearson');

text(min(get(gca,'xlim'))+0.05*[max(get(gca,'xlim'))-min(get(gca,'xlim'))], min(get(gca,'ylim'))+0.25*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], ['r = ', sprintf('%.2f', rho)], 'color', 'r', ...
    'fontsize', fontsize_axis, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'left');
text(min(get(gca,'xlim'))+0.05*[max(get(gca,'xlim'))-min(get(gca,'xlim'))], min(get(gca,'ylim'))+0.12*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], extract_pvalue_text(pval,1,'spin'), 'color', 'r', ...
    'fontsize', fontsize_axis, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'left');
box off

%%% panel letters
a1 = annotation(fig, 'textbox', [0.01, 0.98, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center');
a2 = annotation(fig, 'textbox', [a1.Position(1), a1.Position(2)-0.32, a1.Position(3), a1.Position(4)], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center');
a3 = annotation(fig, 'textbox', [a1.Position(1)+0.7, a1.Position(2), a1.Position(3), a1.Position(4)], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center');
a4 = annotation(fig, 'textbox', [a3.Position(1), a2.Position(2), a1.Position(3), a1.Position(4)], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center');
a5 = annotation(fig, 'textbox', [a3.Position(1), a4.Position(2)-0.32, a1.Position(3), a1.Position(4)], 'string', 'E', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center');

%% FIGURE 5

% load all data relevant to Figure 5
data_Figure5 = load(sprintf('%s/Figure5.mat', data_figures_folder));

structures = {'tha', 'striatum', 'hippo'};
structure_names = {'thalamus', 'striatum', 'hippocampus'};
mode_interest = [1, 2, 3];

subcortical_gradients = data_Figure5.subcortical_gradients;
subcortical_emode = data_Figure5.subcortical_emode;
subcortical_corr = data_Figure5.subcortical_corr;
subcortical_locs = data_Figure5.subcortical_locs;
subcortical_gradients_variance = data_Figure5.subcortical_gradients_variance;

factor_x = 1.35;
factor_y = 1.25;
init_x = 0.06;
init_y = -0.05;
length_x = (1 - 1.5*init_x)/(factor_x*(3-1) + 1);
length_y = (1 - 1.5*init_y)/(factor_y*(2-1) + 1);
markersize = 10;
views = [-27.5 40;
         -61.5, 27;
         41, 43];

fig = figure('Position', [200 200 1000 800]);

for structure_ind = 1:length(structures)
    structure = structures{structure_ind};
    structure_name = structure_names{structure_ind};
    
    % =====================================================================
    % A-C: gradient vs eigenmode
    % =====================================================================
    lim_gradients = double([min(subcortical_gradients.(structure)(:,mode_interest), [], 'all'), max(subcortical_gradients.(structure)(:,mode_interest), [], 'all')]);
    for ii = 1:length(mode_interest)
        mode = mode_interest(ii);
        
        ax1 = axes('Position', [init_x+length_x*0.67+factor_x*length_x*(structure_ind-1) init_y+factor_y*length_y*(2-1)-0.05+length_y/(length(mode_interest)-0.)*(-(ii-length(mode_interest))) length_x*0.33 length_y/(length(mode_interest)+0.85)]);
        data_to_plot_x = subcortical_gradients.(structure)(:,mode);
        data_to_plot_y = subcortical_emode.(structure)(:,mode);
        [rho,~] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
        if rho<0
            data_to_plot_y = -data_to_plot_y;
        end
        hold on;
        plot(data_to_plot_x, data_to_plot_y, 'k.', 'markersize', 6)
        plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
            'r-', 'linewidth', 2);
        hold off;
        set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', lim_gradients, ...
            'ylim', double([min(data_to_plot_y), max(data_to_plot_y)])+[-0.7 0])
        xlabel(sprintf('gradient %i', mode), 'fontsize', fontsize_label-2)
        if ii~=length(mode_interest)
            set(ax1, 'xticklabel', {})
        end
        ylabel(sprintf('mode %i', mode), 'fontsize', fontsize_label-2)
        [rho,pval] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
        text(lim_gradients(2), min(get(gca,'ylim'))+0.09*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], sprintf('r = %.3f', rho), 'color', 'r', ...
            'fontsize', fontsize_axis, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'right');
        box off;
        
        %%% eigenmode 3D visualization
        data_to_plot = subcortical_emode.(structure)(:,mode);
        
        [rho,~] = corr(subcortical_gradients.(structure)(:,mode), data_to_plot, 'type', 'pearson');
        if rho<0
            data_to_plot = -data_to_plot;
        end
        
        fig_temp = figure('Visible', 'off');
        imagesc(data_to_plot);
        new_map = bluewhitered(size(data_to_plot,1));
        close(fig_temp)

        [~, sort_ind] = sort(data_to_plot, 'ascend');
        
        ax1_1 = axes('Position', [ax1.Position(1)-0.21 ax1.Position(2)+ax1.Position(4)*0 length_x*0.34 ax1.Position(4)]);
        scatter3(subcortical_locs.(structure)(sort_ind,1), subcortical_locs.(structure)(sort_ind,2), subcortical_locs.(structure)(sort_ind,3), markersize, new_map, 'filled');
        set(ax1_1, 'xlim', [min(subcortical_locs.(structure)(sort_ind,1)), max(subcortical_locs.(structure)(sort_ind,1))], ...
                 'ylim', [min(subcortical_locs.(structure)(sort_ind,2)), max(subcortical_locs.(structure)(sort_ind,2))], ...
                 'zlim', [min(subcortical_locs.(structure)(sort_ind,3)), max(subcortical_locs.(structure)(sort_ind,3))])
        axis square
        view([views(structure_ind,1) views(structure_ind,2)])
        axis off
        
        annotation(fig, 'textbox', [ax1_1.Position(1), ax1_1.Position(2)+ax1_1.Position(4)*0.88, ax1_1.Position(3), 0.01], 'string', sprintf('mode %i', mode), 'edgecolor', 'none', ...
        'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')   
        
        %%% gradient 3D visualization
        data_to_plot = subcortical_gradients.(structure)(:,mode);
        
        fig_temp = figure('Visible', 'off');
        imagesc(data_to_plot);
        new_map = bluewhitered(size(data_to_plot,1));
        close(fig_temp)

        [~, sort_ind] = sort(data_to_plot, 'ascend');
        
        ax1_2 = axes('Position', [ax1_1.Position(1)+ax1_1.Position(3)*1 ax1_1.Position(2) ax1_1.Position(3) ax1_1.Position(4)]);
        p = scatter3(subcortical_locs.(structure)(sort_ind,1), subcortical_locs.(structure)(sort_ind,2), subcortical_locs.(structure)(sort_ind,3), markersize, new_map, 'filled');
        set(ax1_2, 'xlim', [min(subcortical_locs.(structure)(sort_ind,1)), max(subcortical_locs.(structure)(sort_ind,1))], ...
                 'ylim', [min(subcortical_locs.(structure)(sort_ind,2)), max(subcortical_locs.(structure)(sort_ind,2))], ...
                 'zlim', [min(subcortical_locs.(structure)(sort_ind,3)), max(subcortical_locs.(structure)(sort_ind,3))])
        axis square
        view([views(structure_ind,1) views(structure_ind,2)])
        axis off
        
        annotation(fig, 'textbox', [ax1_2.Position(1), ax1_2.Position(2)+ax1_2.Position(4)*0.88, ax1_2.Position(3), 0.01], 'string', sprintf('gradient %i', mode), 'edgecolor', 'none', ...
        'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle') 
        
        if ii==length(mode_interest)
            % draw direction labels
            xfrac = 0.1; yfrac = 0.1; zfrac = 0.15;
            xstart = min(subcortical_locs.(structure)(:,1)); xend = min(subcortical_locs.(structure)(:,1)) + xfrac*(max(subcortical_locs.(structure)(:,1))-min(subcortical_locs.(structure)(:,1)));
            % ystart = min(y); yend = min(y) + yfrac*(max(y)-min(y));
            ystart = max(subcortical_locs.(structure)(:,2)); yend = max(subcortical_locs.(structure)(:,2)) - yfrac*(max(subcortical_locs.(structure)(:,2))-min(subcortical_locs.(structure)(:,2)));
            zstart = min(subcortical_locs.(structure)(:,3)); zend = min(subcortical_locs.(structure)(:,3)) + zfrac*(max(subcortical_locs.(structure)(:,3))-min(subcortical_locs.(structure)(:,3)));

            if strcmpi(structure, 'tha')
                ax1_1_dir = axes('Position', [ax1_1.Position(1)+0.01 ax1_1.Position(2)-0.12 0.3 0.3]);
            elseif strcmpi(structure, 'striatum')
                ax1_1_dir = axes('Position', [ax1_1.Position(1)+0.01 ax1_1.Position(2)-0.07 ax1_1_dir.Position(3) ax1_1_dir.Position(4)]);
            elseif strcmpi(structure, 'hippo')
                ax1_1_dir = axes('Position', [ax1_1.Position(1)-0.12 ax1_1.Position(2)-0.18 ax1_1_dir.Position(3) ax1_1_dir.Position(4)]);
            end
            hold on;
            plot3([xstart xend], ystart*ones(1,2), zstart*ones(1,2), 'k-')
            plot3(xstart*ones(1,2), [ystart yend], zstart*ones(1,2), 'k-')
            plot3(xstart*ones(1,2), ystart*ones(1,2), [zstart zend], 'k-')
            hold off;
            text(min(subcortical_locs.(structure)(:,1)) + (xfrac+0.05*xfrac)*(max(subcortical_locs.(structure)(:,1))-min(subcortical_locs.(structure)(:,1))), ystart, zstart, 'R') 
            text(xstart, max(subcortical_locs.(structure)(:,2)) - (yfrac+0.1*yfrac)*(max(subcortical_locs.(structure)(:,2))-min(subcortical_locs.(structure)(:,2))), zstart, 'P') 
            text(xstart, ystart, min(subcortical_locs.(structure)(:,3)) + (zfrac+0.35*zfrac)*(max(subcortical_locs.(structure)(:,3))-min(subcortical_locs.(structure)(:,3))), 'D') 
            xlim([min(subcortical_locs.(structure)(:,1)) max(subcortical_locs.(structure)(:,1))])
            ylim([min(subcortical_locs.(structure)(:,2)) max(subcortical_locs.(structure)(:,2))])
            zlim([min(subcortical_locs.(structure)(:,3)) max(subcortical_locs.(structure)(:,3))])
            [az, el] = view(ax1_1);
            view([az el])
            axis off
        end
    end
    
    % =====================================================================
    % D-F: correlation matrix and max absolute correlation
    % =====================================================================
    cmap = [1 1 1; cbrewer('seq', 'Reds', 100, 'pchip')];
    
    ax2 = axes('Position', [init_x+factor_x*length_x*(structure_ind-1) init_y+factor_y*length_y*(0.98-1) length_x length_y]);
    imagesc(abs(subcortical_corr.(structure)))
    caxis([0 1])
    colormap(cmap)
    cbar = colorbar('fontsize', fontsize_axis-2);
    set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xtick', 5:5:20, 'ytick', 5:5:20)
    xlabel('gradient', 'fontsize', fontsize_label)
    ylabel('mode', 'fontsize', fontsize_label)
    ylabel(cbar, 'absolute correlation, |r|', 'fontsize', fontsize_label-2)
    axis square
    
    ax2_2 = axes('Position', [ax2.Position(1) ax2.Position(2)+ax2.Position(4)*0.78 ax2.Position(3) ax2.Position(4)*0.18]);
    colororder({'k', 'b'})
    yyaxis(ax2_2, 'left')
    bar(1:size(subcortical_corr.(structure),1), max(abs((subcortical_corr.(structure))),[],1), 'facecolor', 0.8*ones(1,3))
    set(ax2_2, 'fontsize', fontsize_axis-2, 'ticklength', [0.02, 0.02], 'xtick', 5:5:20, 'xticklabel', {}, 'xlim', [1-0.5, size(subcortical_corr.(structure),1)+0.5], 'ylim', [0 1])
    ylabel('max |r|', 'fontsize', fontsize_label-2)
    box on
    
    yyaxis(ax2_2, 'right')
    plot(1:size(subcortical_corr.(structure),1), subcortical_gradients_variance.(structure)(1:size(subcortical_corr.(structure),1)), 'b.-', ...
        'markersize', 15, 'linewidth', 1.5)
    set(ax2_2, 'fontsize', fontsize_axis-2, 'ticklength', [0.02, 0.02])
    ylabel({'% variance'; 'explained'}, 'fontsize', fontsize_label-2)
    
    %%% structure names
    annotation(fig, 'textbox', [ax1.Position(1)-0.19, 0.97, length_x, 0.01], 'string', structure_name, 'edgecolor', 'none', ...
            'fontsize', fontsize_label, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')    
end

%%% panel letters
a1 = annotation(fig, 'textbox', [0.01, 0.98, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center');
a2 = annotation(fig, 'textbox', [a1.Position(1), 0.42, a1.Position(3), a1.Position(4)], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center');
a3 = annotation(fig, 'textbox', [a1.Position(1)+0.33, a1.Position(2), a1.Position(3), a1.Position(4)], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center');
annotation(fig, 'textbox', [a3.Position(1), a2.Position(2), a1.Position(3), a1.Position(4)], 'string', 'E', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
a4 = annotation(fig, 'textbox', [a1.Position(1)+0.66, a1.Position(2), a1.Position(3), a1.Position(4)], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center');
annotation(fig, 'textbox', [a4.Position(1), a2.Position(2), a1.Position(3), a1.Position(4)], 'string', 'F', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
