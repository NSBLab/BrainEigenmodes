%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate_paper_figures_extended_Nature.m
%%%
%%% MATLAB script to generate the final Extended Data figures of the paper
%%%
%%% NOTE : The configuration of your computer (e.g., screen resolution)
%%%        affects how the figures created will look. Hence, they will not 
%%%        100% visually match the figures in the paper, but the scientific 
%%%        contents are replicated.
%%%
%%% Original: James Pang, Monash University, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load relevant repository MATLAB functions

addpath(genpath('functions_matlab'));

%% MISCELLANEOUS VARIABLES

% NOTE: data provided is only for the below parameters, so please don't change them
hemisphere = 'lh'; 
num_modes = 200;
parc_name = 'Glasser360';

data_empirical_folder = 'data/empirical';
data_results_folder = 'data/results';
data_figures_folder = 'data/figures_Nature';
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

% RESULTS: resting recon accuracy
% eigenmodes
load(sprintf('%s/resting_reconstruction_corr_modes=%i_%s_%s.mat', data_results_folder, num_modes, parc_name, hemisphere), ...
             'resting_recon_corr_geometric', 'resting_recon_corr_connectome', ...
             'resting_recon_corr_connectome_density_matched', 'resting_recon_corr_EDR')

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

%% EXTENDED DATA FIGURE 1
% bioRxiv Supplementary Figure 1

surface_to_plot = surface_midthickness;

medial_wall = find(~cortex);
mode_list = [1:5, 10,25,50,100,200];
num_modes_to_plot = length(mode_list);

factor_x_small = 1.05;
factor_x_big = 2.3;
factor_y = 1.1;
init_x = 0.05;
init_y = 0.01;
length_x = (0.99-init_x)/(factor_x_small + factor_x_big*(4-1) + 1);
length_y = (0.95-init_y)/(factor_y*(num_modes_to_plot-1) + 1);

fig = figure('Position', [200 200 1000 900]);
for basis_ind=1:4
    if basis_ind==1
        data_to_plot = basis_geometric;
        basis_name = 'geometric';
    elseif basis_ind==2
        data_to_plot = basis_connectome;
        basis_name = 'connectome';
    elseif basis_ind==3
        data_to_plot = basis_connectome_density_matched;
        basis_name = {'connectome'; '(density matched)'};
    elseif basis_ind==4
        data_to_plot = basis_EDR;
        basis_name = 'EDR';
    end

    for mode_ind=1:num_modes_to_plot
        mode = mode_list(mode_ind);
        
        data_to_plot(medial_wall,mode) = min(data_to_plot(:,mode))*1.1;
        clims = [min(data_to_plot(:,mode)), max(data_to_plot(:,mode))];
        if clims(2)<=0
            clims(2) = 0.01;
        end

        ax1 = axes('Position', [init_x+factor_x_big*length_x*(basis_ind-1) init_y+factor_y*length_y*(-[mode_ind-num_modes_to_plot]) length_x length_y]);
        obj1 = patch(ax1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,mode), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
        if strcmpi(hemisphere, 'lh')
            view([-90 0]);
        elseif strcmpi(hemisphere, 'rh')
            view([90 0]);
        end
        caxis(clims)
        camlight('headlight')
        material dull
        colormap(ax1,[0.5,0.5,0.5; bluewhitered])
        axis off
        axis image
        
        if basis_ind==1
            annotation(fig, 'textbox', [ax1.Position(1)-0.03, ax1.Position(2)+ax1.Position(4)*0.45, 0.01, 0.01], 'string', {'mode'; sprintf('%i',mode)}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
        end
        
        ax2 = axes('Position', [init_x+factor_x_small*length_x+factor_x_big*length_x*(basis_ind-1) init_y+factor_y*length_y*(-[mode_ind-num_modes_to_plot]) length_x length_y]);
        obj2 = patch(ax2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,mode), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
        if strcmpi(hemisphere, 'lh')
            view([90 0]);
        elseif strcmpi(hemisphere, 'rh')
            view([-90 0]);
        end
        caxis(clims)
        camlight('headlight')
        material dull
        colormap(ax2,[0.5,0.5,0.5; bluewhitered])
        axis off
        axis image
        
        if mode_ind==1
            annotation(fig, 'textbox', [ax1.Position(1), ax1.Position(2)+ax1.Position(4)*1, ax2.Position(1)-ax1.Position(1)+ax2.Position(3), 0.01], 'string', basis_name, 'edgecolor', 'none', ...
                'fontsize', fontsize_label, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
        end
    end
end

%% EXTENDED DATA FIGURE 2
% bioRxiv Supplementary Figure 7

% load all data relevant to Extended Data Figure 2
data_ExtFigure2 = load(sprintf('%s/ExtFigure2.mat', data_figures_folder));

parc_name_list = {'Schaefer100', 'Schaefer200', 'Glasser360', 'Schaefer400', 'Schaefer600', 'Schaefer800', 'Schaefer1000'};

cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];

factor_x = 1.3;
factor_y = 1.3;
init_x = 0.05;
init_y = 0.07;
length_x = (1 - 1.5*init_x)/(factor_x*(4-1) + 1);
length_y = (1 - 1.8*init_y)/(factor_y*(2-1) + 1);

fig = figure('Position', [200 200 1000 600]);
for parc_name_ind = 1:length(parc_name_list)
    parc_name_interest = parc_name_list{parc_name_ind};
    
    if parc_name_ind<=4
        ax1 = axes('Position', [init_x+factor_x*length_x*(parc_name_ind-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax1 = axes('Position', [init_x+factor_x*length_x*(parc_name_ind-4-1)+0.12 init_y+factor_y*length_y*(1-1) length_x length_y]);
    end
    hold on;
    plot(data_ExtFigure2.resting_recon_corr_parcellations_geometric.(parc_name_interest), 'color', colors2(1,:), 'linewidth', 2, 'displayname', 'geometric')
    plot(data_ExtFigure2.resting_recon_corr_parcellations_connectome.(parc_name_interest), 'color', colors2(2,:), 'linewidth', 2, 'displayname', 'connectome')
    plot(data_ExtFigure2.resting_recon_corr_parcellations_connectome_density_matched.(parc_name_interest), 'color', colors2(3,:), 'linewidth', 2, 'displayname', 'connectome (density matched)')
    plot(data_ExtFigure2.resting_recon_corr_parcellations_EDR.(parc_name_interest), 'color', colors2(4,:), 'linewidth', 2, 'displayname', 'EDR')
    hold off;
    if parc_name_ind==length(parc_name_list)
        leg = legend('fontsize', fontsize_axis, 'interpreter', 'none', 'box', 'off', 'numcolumns', 1, ...
                    'Position', [ax1.Position(1)+ax1.Position(3)*0.7 ax1.Position(2)+ax1.Position(4)*0.1 ax1.Position(3)*0.5 ax1.Position(3)*0.5]);
        leg.ItemTokenSize = leg.ItemTokenSize/2;
    end
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1])
    if parc_name_ind>4
        xlabel('number of modes', 'fontsize', fontsize_label)
    end
    if parc_name_ind==1 || parc_name_ind==5
        ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
    end
    title(parc_name_interest, 'fontsize', fontsize_axis, 'fontweight', 'normal')
end

%% EXTENDED DATA FIGURE 3
% bioRxiv Supplementary Figure 8

% load all data relevant to Extended Data Figure 3
data_ExtFigure3 = load(sprintf('%s/ExtFigure3.mat', data_figures_folder));

parc_name_list = {'Schaefer100', 'Schaefer200', 'Glasser360', 'Schaefer400', 'Schaefer600', 'Schaefer800', 'Schaefer1000'};

factor_x = 1.3;
factor_y = 1.3;
init_x = 0.08;
init_y = 0.07;
length_x = (1 - 1.5*init_x)/(factor_x*(4-1) + 1);
length_y = (1 - 1.8*init_y)/(factor_y*(2-1) + 1);

yticks = [];
for ii=1:length(tasks)
    yticks = [yticks, mean(task_contrasts_ind.(tasks{ii}))];
end

recon_corr_combined = zeros(length(contrasts), num_modes, length(parc_name_list));
clims_combined = zeros(length(parc_name_list),2);
for parc_name_ind = 1:length(parc_name_list)
    parc_name_interest = parc_name_list{parc_name_ind};
    
    data_1 = data_ExtFigure3.task_recon_corr_parcellations_geometric.(parc_name_interest);
    
    % connectome eigenmodes
    data_2 = data_ExtFigure3.task_recon_corr_parcellations_connectome.(parc_name_interest);
    
    recon_corr = [];
    for ii=1:length(contrasts)
        contrast = contrasts{ii};

        temp_1 = nanmean(data_1.(contrast),1);
        temp_2 = nanmean(data_2.(contrast),1);
        recon_corr = cat(1, recon_corr, temp_1-temp_2);
    end
    recon_corr_combined(:,:,parc_name_ind) = recon_corr;
    
    clim_max = max(recon_corr(:,2:end),[],'all');
    clim_min = min(recon_corr(:,2:end),[],'all');
    clims_combined(parc_name_ind,:) = [clim_min, clim_max];
end

fig = figure('Position', [200 200 1000 600]);
for parc_name_ind = 1:length(parc_name_list)
    parc_name_interest = parc_name_list{parc_name_ind};
    
    data_to_plot = recon_corr_combined(:,:,parc_name_ind);
    clim = [min(clims_combined(:,1)), max(clims_combined(:,2))]; % make absolute colorbar for all parcellations
    
    if parc_name_ind<=4
        ax1 = axes('Position', [init_x+factor_x*length_x*(parc_name_ind-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax1 = axes('Position', [init_x+factor_x*length_x*(parc_name_ind-4-1)+0.12 init_y+factor_y*length_y*(1-1) length_x length_y]);
    end
    imagesc(data_to_plot)
    for ii=1:length(contrasts)-1
        yline(ii+0.5, 'k-');
    end
    for ii=1:length(tasks)-1
        yline(task_contrasts_ind.(tasks{ii})(end)+0.5, 'k-', 'linewidth', 2);
    end
    caxis(clim)
    colormap(ax1, bluewhitered)
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ytick', yticks, 'yticklabel', {}, 'ticklabelinterpreter', 'none')
    if parc_name_ind>4
        xlabel('number of modes', 'fontsize', fontsize_label)
    end
    if parc_name_ind==1 || parc_name_ind==5
        set(ax1, 'yticklabel', lower(tasks))
    end
    title(parc_name_interest, 'fontsize', fontsize_axis, 'fontweight', 'normal')
end

cbar = colorbar;
set(cbar, 'ticklength', 0.02, 'Position', [ax1.Position(1)+ax1.Position(3)*1.1 ax1.Position(2) 0.01 ax1.Position(4)], ...
    'fontsize', fontsize_axis-2);
ylabel(cbar, {'reconstruction accuracy difference'; ['(geometric ' char(8211) ' connectome)']}, 'fontsize', fontsize_label)

%% EXTENDED DATA FIGURE 4
% bioRxiv Supplementary Figure 13

% load all data relevant to Extended Data Figure 4
data_ExtFigure4 = load(sprintf('%s/ExtFigure4.mat', data_figures_folder));

cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];

fig = figure('Position', [200 200 1000 600]);

% =========================================================================
% A: reconstrution accuracy (resting)
% =========================================================================
ax1 = axes('Position', [0.05 0.35 0.35 0.45]);
hold on;
plot(data_ExtFigure4.resting_recon_corr_geometric, 'color', colors2(1,:), 'linewidth', 2, 'displayname', 'geometric eigenmodes')
plot(data_ExtFigure4.resting_recon_corr_PCA, 'color', colors2(5,:), 'linewidth', 2, 'displayname', 'PCA (rest)')
hold off;
leg = legend('fontsize', fontsize_axis, 'location', 'southeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/2;
set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1])
xlabel('number of modes/components', 'fontsize', fontsize_label)
ylabel('reconstruction accuracy', 'fontsize', fontsize_label)

% =========================================================================
% B: reconstrution accuracy (task)
% =========================================================================
factor_x = 1.15;
init_x = ax1.Position(1)+ax1.Position(3)*0.9;
init_y = 0.54;
length_x = (0.98 - init_x)/(factor_x*(4+1-1) + 1);
length_y = 0.36;
yticks = [];
for ii=1:length(tasks)
    yticks = [yticks, mean(task_contrasts_ind.(tasks{ii}))];
end

data_1 = data_ExtFigure4.task_recon_corr_geometric;
recon_corr_combined = zeros(length(contrasts), num_modes, 7);
clims_combined = zeros(7,2);
for kk=1:7
    contrast_to_ind = representative_contrasts{kk};
    
    data_2 = data_ExtFigure4.task_recon_corr_PCA.(contrast_to_ind);

    recon_corr = [];
    for ii=1:length(contrasts)
        contrast = contrasts{ii};

        temp_1 = data_1.(contrast);
        temp_2 = data_2.(contrast);
        recon_corr = cat(1, recon_corr, temp_1-temp_2);
    end
    recon_corr_combined(:,:,kk) = recon_corr;
    
    clim_max = max(recon_corr(:,4:end),[],'all');
    clim_min = min(recon_corr(:,4:end),[],'all');
    clims_combined(kk,:) = [clim_min, clim_max];
end

for kk=1:4
    contrast_to_ind = representative_contrasts{kk};
    task_name = {'PCA'; sprintf('(%s)', lower(tasks{kk}))};
    
    data_to_plot = recon_corr_combined(:,:,kk);
    clim = [min(clims_combined(:,1)), max(clims_combined(:,2))]; % make absolute colorbar for all PCAs

    ax2 = axes('Position', [init_x+factor_x*length_x*(kk) init_y  length_x length_y]);
    imagesc(data_to_plot)
    for ii=1:length(contrasts)-1
        yline(ii+0.5, 'k-');
    end
    for ii=1:length(tasks)-1
        yline(task_contrasts_ind.(tasks{ii})(end)+0.5, 'k-', 'linewidth', 2);
    end
    text(num_modes+2, find(strcmpi(contrasts, contrast_to_ind))-0.5, '\ast', 'fontsize', 10, 'fontweight', 'b', 'verticalalignment', 'middle')
    caxis(clim)
    colormap(ax2, bluewhitered)
    set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ytick', yticks, 'yticklabel', lower(tasks), 'ticklabelinterpreter', 'none')
    if kk~=1
        set(ax2, 'yticklabel', {})
    end
    title(task_name, 'fontweight', 'normal', 'fontsize', fontsize_axis)
end

for kk=1:3
    contrast_to_ind = representative_contrasts{4+kk};
    task_name = {'PCA'; sprintf('(%s)', lower(tasks{4+kk}))};
    
    data_to_plot = recon_corr_combined(:,:,4+kk);
    clim = [min(clims_combined(:,1)), max(clims_combined(:,2))]; % make absolute colorbar for all PCAs

    ax2 = axes('Position', [init_x+factor_x*length_x*(kk)+0.05 init_y-length_y*1.3 length_x length_y]);
    imagesc(data_to_plot)
    for ii=1:length(contrasts)-1
        yline(ii+0.5, 'k-');
    end
    for ii=1:length(tasks)-1
        yline(task_contrasts_ind.(tasks{ii})(end)+0.5, 'k-', 'linewidth', 2);
    end
    text(num_modes+2, find(strcmpi(contrasts, contrast_to_ind))-0.5, '\ast', 'fontsize', 10, 'fontweight', 'b', 'verticalalignment', 'middle')
    caxis(clim)
    colormap(ax2, bluewhitered)
    set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ytick', yticks, 'yticklabel', lower(tasks), 'ticklabelinterpreter', 'none')
    if kk~=1
        set(ax2, 'yticklabel', {})
    end
    if kk==2
        xlabel('number of modes/components', 'fontsize', fontsize_label)
    end
    title(task_name, 'fontweight', 'normal', 'fontsize', fontsize_axis)
end

cbar = colorbar;
set(cbar, 'ticklength', 0.02, 'Position', [ax2.Position(1)+ax2.Position(3)*1.1 ax2.Position(2) 0.01 ax2.Position(4)], ...
    'fontsize', fontsize_axis-2);
ylabel(cbar, {'reconstruction accuracy difference'; ['(geometric eigenmodes ' char(8211) ' PCA)']}, 'fontsize', fontsize_label)

%%% panel letters
annotation(fig, 'textbox', [0.01, 0.86, 0.01, 0.01], 'string', 'a', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.43, 0.98, 0.01, 0.01], 'string', 'b', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% EXTENDED DATA FIGURE 5
% bioRxiv Supplementary Figure 14

% load all data relevant to Extended Data Figure 5
data_ExtFigure5 = load(sprintf('%s/ExtFigure5.mat', data_figures_folder));

surface_to_plot = surface_midthickness;

medial_wall = find(~cortex);
mode_list = [1:4, 10,50,100];
num_modes_to_plot = length(mode_list);

basis_name_list = {'geometric', 'F_1 (reg)', 'F_2 (reg)', 'F_3 (reg)', ...
                   'F_1 (irreg)', 'F_2 (irreg)', 'F_3 (irreg)'};
basis_interest = {};
for basis_ind = 1:6
    basis_interest{basis_ind} = nansum(basis_Fourier{basis_ind},3);
end

fig = figure('Position', [200 200 1000 1000]);

init_x_1 = 0.04;
init_y_1 = 0.58;
init_x_2 = 0.05;
init_y_2 = 0.05;

% =========================================================================
% A: basis sets
% =========================================================================
factor_x_small = 0.52; 
factor_x_big = 1.1;
factor_y = 1.1;
init_x = init_x_1;
init_y = init_y_1;
length_x = (0.98)/(factor_x_small + factor_x_big*(6-1) + 1);
length_y = (0.96-init_y)/(factor_y*(num_modes_to_plot-1) + 1);

for basis_ind=1:length(basis_name_list)-1
    data_to_plot = basis_interest{basis_ind};
    basis_name = basis_name_list{basis_ind+1};

    for mode_ind=1:1:num_modes_to_plot
        mode = mode_list(mode_ind);
        
        if min(data_to_plot(:,mode))==0
            data_to_plot(medial_wall,mode) = -0.005;
        else
            data_to_plot(medial_wall,mode) = min(data_to_plot(:,mode))*1.1;
        end
        clims = [min(data_to_plot(:,mode)), max(data_to_plot(:,mode))];
        if clims(2)<=0
            clims(2) = 0.01;
        end

        ax1 = axes('Position', [init_x+factor_x_big*length_x*(basis_ind-1) init_y+factor_y*length_y*(-[mode_ind-num_modes_to_plot]) length_x length_y]);
        obj1 = patch(ax1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,mode), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
        if strcmpi(hemisphere, 'lh')
            view([-90 0]);
        elseif strcmpi(hemisphere, 'rh')
            view([90 0]);
        end
        caxis(clims)
        camlight('headlight')
        material dull
        colormap(ax1,[0.5,0.5,0.5; bluewhitered])
        axis off
        axis image
        
        if basis_ind==1
            annotation(fig, 'textbox', [ax1.Position(1)-0.0, ax1.Position(2)+ax1.Position(4)*0.4, 0.01, 0.01], 'string', {'mode'; sprintf('%i',mode)}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
        end
        
        ax2 = axes('Position', [init_x+factor_x_small*length_x+factor_x_big*length_x*(basis_ind-1) init_y+factor_y*length_y*(-[mode_ind-num_modes_to_plot]) length_x length_y]);
        obj2 = patch(ax2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,mode), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
        if strcmpi(hemisphere, 'lh')
            view([90 0]);
        elseif strcmpi(hemisphere, 'rh')
            view([-90 0]);
        end
        caxis(clims)
        camlight('headlight')
        material dull
        colormap(ax2,[0.5,0.5,0.5; bluewhitered])
        axis off
        axis image
        
        if mode_ind==1
            annotation(fig, 'textbox', [ax1.Position(1), ax1.Position(2)+ax1.Position(4)*1, ax2.Position(1)-ax1.Position(1)+ax2.Position(3), 0.01], 'string', basis_name, 'edgecolor', 'none', ...
                'fontsize', fontsize_label, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'bottom', 'interpreter', 'none')
        end
    end
end

% =========================================================================
% B: reconstruction accuracy
% =========================================================================
cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];

factor_x = 1.3;
factor_y = 1.3;
init_x = init_x_2;
init_y = init_y_2;
length_x = (1 - 1.3*init_x)/(factor_x*(4-1) + 1);
length_y = (init_y_1*0.95 - 1.8*init_y)/(factor_y*(2-1) + 1);

for ii=1:length(representative_contrasts)
    contrast = representative_contrasts{ii};
    task_name = lower(tasks{ii});
    
    color_shade_cmap = colors2(1:length(basis_name_list)-3,:);
    
    if ii<=4
        ax3 = axes('Position', [init_x+factor_x*length_x*(ii-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax3 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
    end
    hold on;
    for basis_ind = 1:length(basis_name_list)
        basis_name = basis_name_list{basis_ind};
        
        if basis_ind==1
            data_to_plot = data_ExtFigure5.task_recon_corr_geometric.(contrast);
        else
            data_to_plot = data_ExtFigure5.task_recon_corr_Fourier{basis_ind-1}.(contrast);
        end
        
        if basis_ind<=4
            plot(data_to_plot, '-', 'color', color_shade_cmap(basis_ind,:), 'linewidth', 2, 'displayname', basis_name)
        else
            plot(data_to_plot, '--', 'color', color_shade_cmap(basis_ind-3,:), 'linewidth', 2, 'displayname', basis_name)
        end
    end
    hold off;
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1])
    if ii>4
        xlabel('number of modes', 'fontsize', fontsize_label)
    end
    if ii==1 || ii==5
        ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
    end
    title(task_name, 'fontsize', fontsize_axis, 'fontweight', 'normal')
end

% rest
ii = 8;
color_shade_cmap = colors2(1:length(basis_name_list)-3,:);

ax4 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
hold on;
for basis_ind = 1:length(basis_name_list)
    basis_name = basis_name_list{basis_ind};

    if basis_ind==1
        data_to_plot = data_ExtFigure5.resting_recon_corr_geometric;
    else
        data_to_plot = data_ExtFigure5.resting_recon_corr_Fourier{basis_ind-1};
    end
    
    if basis_ind<=4
        plot(data_to_plot, '-', 'color', color_shade_cmap(basis_ind,:), 'linewidth', 2, 'displayname', basis_name)
    else
        plot(data_to_plot, '--', 'color', color_shade_cmap(basis_ind-3,:), 'linewidth', 2, 'displayname', basis_name)
    end
end
hold off;
leg = legend('fontsize', fontsize_axis, 'position', [init_x init_y+factor_y*length_y*(2-1)+length_y*1.15 1-init_x 0.01], 'interpreter', 'none', 'box', 'off', 'numcolumns', 7);
leg.ItemTokenSize = leg.ItemTokenSize/1;
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1])
xlabel('number of modes', 'fontsize', fontsize_label)
if ii==1
    ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
end
title('rest', 'fontsize', fontsize_axis, 'fontweight', 'normal')

%%% panel letters
annotation(fig, 'textbox', [0.015, 0.99, 0.01, 0.01], 'string', 'a', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.015, 0.56, 0.01, 0.01], 'string', 'b', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% EXTENDED DATA FIGURE 6
% bioRxiv Supplementary Figure 15

surface_to_plot = surface_midthickness;

tmap = data_neurovault(:,1);
threshold_1 = 1.6;
threshold_2 = 2.3;
tmap_threshold_1 = zeros(size(tmap));
tmap_threshold_1(abs(tmap)>threshold_1) = 1;
tmap_threshold_2 = zeros(size(tmap));
tmap_threshold_2(abs(tmap)>threshold_2) = 1;
 
fig = figure('Position', [200 200 700 400]);

% =========================================================================
% A: line plot
% =========================================================================
x = 0:0.01:8;
data_to_plot = 0.8*(sin(x) + 0.8*sin(2*x) + 2.5*sin(3*x));

ax1 = axes('Position', [0.1 0.08 0.5 0.84]);
plot(x, data_to_plot, 'k', 'linewidth', 2)
yline(threshold_1, 'b--', 'linewidth', 2);
yline(-threshold_1, 'b--', 'linewidth', 2);
yline(threshold_2, 'm--', 'linewidth', 2);
yline(-threshold_2, 'm--', 'linewidth', 2);
set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xtick', [], 'ytick', [-max(abs(data_to_plot)), 0, max(abs(data_to_plot))], ...
            'yticklabel', {'min', '0', 'max'}, 'ylim', [-max(abs(data_to_plot)), max(abs(data_to_plot))], 'xlim', [0 7.5])
text(max(get(gca,'xlim'))*1.02, threshold_1, 'threshold_1', 'fontsize', fontsize_axis, 'color', 'b')
text(max(get(gca,'xlim'))*1.02, -threshold_1, 'threshold_1', 'fontsize', fontsize_axis, 'color', 'b')
text(max(get(gca,'xlim'))*1.02, threshold_2, 'threshold_2', 'fontsize', fontsize_axis, 'color', 'm')
text(max(get(gca,'xlim'))*1.02, -threshold_2, 'threshold_2', 'fontsize', fontsize_axis, 'color', 'm')
xlabel('vertex', 'fontsize', fontsize_label)
ylabel('statistic', 'fontsize', fontsize_label)
box off

% =========================================================================
% B: unthresholded and thresholded maps
% =========================================================================
factor_x = 1;
factor_y = 1.5;
init_x = ax1.Position(1)+ax1.Position(3)*1.08;
init_y = ax1.Position(2)*0.5;
length_x = 0.95-ax1.Position(3);
length_y = (ax1.Position(2)+ax1.Position(4)-init_y)/(factor_y*(3-1) + 1);

% unthresholded
data_to_plot = tmap;

ax2 = axes('Position', [init_x init_y+factor_y*length_y*(3-1) length_x length_y]);
patch(ax2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');       
if strcmpi(hemisphere, 'lh')
    view([-90 0]);
elseif strcmpi(hemisphere, 'rh')
    view([90 0]);
end
camlight('headlight')
material dull
caxis([-max(abs(data_to_plot)), max(abs(data_to_plot))])
colormap(ax2, bluewhitered)
axis off
axis image
title('unthresholded', 'fontsize', fontsize_axis, 'fontweight', 'normal')

cbar = colorbar(ax2,'southoutside');
set(cbar, 'fontsize', fontsize_axis-2, 'ticklength', 0.02, 'ytick', [-max(abs(data_to_plot)), 0, max(abs(data_to_plot))], 'yticklabel', {'min', '0', 'max'}, ...
    'position', [ax2.Position(1)+ax2.Position(3)*0.6, ax2.Position(2)+0.01, ax2.Position(3)*0.12, 0.015])

% threshold_1
data_to_plot = tmap_threshold_1;

ax3 = axes('Position', [init_x init_y+factor_y*length_y*(2-1) length_x length_y]);
patch(ax3, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
if strcmpi(hemisphere, 'lh')
    view([-90 0]);
elseif strcmpi(hemisphere, 'rh')
    view([90 0]);
end
camlight('headlight')
material dull
colormap(ax3, [1 1 1; 0 0 0])
axis off
axis image
title('threshold_1', 'fontsize', fontsize_axis, 'fontweight', 'normal', 'color', 'b')

% threshold_2
data_to_plot = tmap_threshold_2;

ax4 = axes('Position', [init_x init_y+factor_y*length_y*(1-1) length_x length_y]);
patch(ax4, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
if strcmpi(hemisphere, 'lh')
    view([-90 0]);
elseif strcmpi(hemisphere, 'rh')
    view([90 0]);
end
camlight('headlight')
material dull
colormap(ax4, [1 1 1; 0 0 0])
axis off
axis image
title('threshold_2', 'fontsize', fontsize_axis, 'fontweight', 'normal', 'color', 'm')

%%% panel letters
annotation(fig, 'textbox', [0.02, 0.98, 0.01, 0.01], 'string', 'a', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.75, 0.98, 0.01, 0.01], 'string', 'b', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% EXTENDED DATA FIGURE 7
% bioRxiv Supplementary Figure 16

cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];

xlims = [2 num_modes];
ylims = [2e-3 1];

factor_x = 1.3;
factor_y = 1.3;
init_x = 0.07;
init_y = 0.07;
length_x = (1 - 1.5*init_x)/(factor_x*(4-1) + 1);
length_y = (1 - 1.8*init_y)/(factor_y*(2-1) + 1);

fig = figure('Position', [200 200 1000 600]);
for ii=1:length(representative_contrasts)
    contrast = representative_contrasts{ii};
    task_name = lower(tasks{ii});
    
    data_to_plot = spectrum_HCP(:,strcmpi(contrasts, contrast));
    data_to_plot = data_to_plot(2:end);
    
    if ii<=4
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-4-1)+0.12 init_y+factor_y*length_y*(1-1) length_x length_y]);
    end
    hold on;
    bar(2:num_modes, data_to_plot/max(data_to_plot), 'facecolor', colors(ii,:), 'displayname', 'data')
    hold off;
    title(task_name, 'fontsize', fontsize_axis, 'fontweight', 'normal')
    if ii==length(representative_contrasts)
        leg = legend('fontsize', fontsize_axis, 'interpreter', 'none', 'box', 'off', 'numcolumns', 1, ...
                    'Position', [ax1.Position(1)+ax1.Position(3)*1.2 ax1.Position(2)+ax1.Position(4)*0.4 ax1.Position(3)*0.5 ax1.Position(3)*0.5]);
        leg.ItemTokenSize = leg.ItemTokenSize/2;
        set(leg, 'visible', 'off')
    end
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', xlims, 'ylim', ylims, 'yscale', 'log')
    if ii>4
        xlabel('mode', 'fontsize', fontsize_label)
    end
    if ii==1 || ii==5
        ylabel('normalized power (log scale)', 'fontsize', fontsize_label)
    end
end

%% EXTENDED DATA FIGURE 8
% bioRxiv Supplementary Figure 17

% load all data relevant to Extended Data Figure 8
data_ExtFigure8 = load(sprintf('%s/ExtFigure8.mat', data_figures_folder));

cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];

fwhm_list = 0:2:50;

xlims = [2 num_modes];
ylims = [5e-5 1e-1];

fig = figure('Position', [200 200 1000 500]);

factor_x = 1.1;
factor_y = 1.35;
init_x = 0.07;
init_y = 0.09;
length_x = 0.37;
length_y = (0.95 - init_y)/(factor_y*(2-1) + 1);

% =========================================================================
% A top: power spectrum (HCP)
% =========================================================================   
fwhm_interest_list = [40,30,20,10,0];
color_shade_cmap = cbrewer('seq', 'YlGnBu', length(fwhm_interest_list), 'pchip');
color_shade_cmap = color_shade_cmap(1:end,:);

data_to_plot = nanmean(spectrum_HCP,2);
ax1_1 = axes('Position', [init_x+factor_x*length_x*(1-1) init_y+factor_y*length_y*(2-1) length_x*1.08 length_y]);
hold on;
bar(1:num_modes, data_to_plot, 'facecolor', 0.2*ones(1,3), 'displayname', 'data')
for jj=1:length(fwhm_interest_list)
    fwhm_interest = fwhm_interest_list(jj);
    fwhm_interest_index = find(fwhm_list==fwhm_interest);
    
    if fwhm_interest==20
        plot(1:num_modes, spectrum_noise{fwhm_interest_index}, 'color', color_shade_cmap(jj,:), 'linewidth', 2, 'displayname', sprintf('FWHM = %i mm (best fit)', fwhm_interest))
        text(num_modes+3, spectrum_noise{fwhm_interest_index}(end-3), sprintf('FWHM = %i mm (best fit)', fwhm_interest), 'fontsize', fontsize_axis-2, 'fontweight', 'normal', ...
            'horizontalalignment', 'left', 'verticalalignment', 'middle')
    else
        plot(1:num_modes, spectrum_noise{fwhm_interest_index}, 'color', color_shade_cmap(jj,:), 'linewidth', 2, 'displayname', sprintf('FWHM = %i mm', fwhm_interest))
        text(num_modes+3, spectrum_noise{fwhm_interest_index}(end-3), sprintf('FWHM = %i mm', fwhm_interest), 'fontsize', fontsize_axis-2, 'fontweight', 'normal', ...
            'horizontalalignment', 'left', 'verticalalignment', 'middle')
    end
end
hold off;
title('HCP', 'fontsize', fontsize_axis, 'fontweight', 'normal')
leg = legend('fontsize', fontsize_axis, 'location', 'eastoutside', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/2;
set(leg, 'visible', 'off')
set(ax1_1, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', xlims, 'ylim', ylims, 'yscale', 'log')
% xlabel('mode', 'fontsize', fontsize_label)
% ylabel_h = ylabel('normalized power (log scale)', 'fontsize', fontsize_label);

% =========================================================================
% A bottom: power spectrum (NeuroVault)
% =========================================================================   
fwhm_interest_list = [40,30,20,10,0]; %[0,10,20,30,40];
color_shade_cmap = cbrewer('seq', 'YlGnBu', length(fwhm_interest_list), 'pchip');
color_shade_cmap = color_shade_cmap(1:end,:);

data_to_plot = nanmean(spectrum_neurovault,2);
ax1_2 = axes('Position', [init_x+factor_x*length_x*(1-1) init_y+factor_y*length_y*(1-1) length_x*1.08 length_y]);
% ax1_2 = axes('Position', [ax0_positions(1) ax0_positions(2) ax0_positions(3) ax0_positions(4)*0.5]);
hold on;
bar(1:num_modes, data_to_plot, 'facecolor', 0.2*ones(1,3), 'displayname', 'data')
for jj=1:length(fwhm_interest_list)
    fwhm_interest = fwhm_interest_list(jj);
    fwhm_interest_index = find(fwhm_list==fwhm_interest);

    if fwhm_interest==20
        plot(1:num_modes, spectrum_noise{fwhm_interest_index}, 'color', color_shade_cmap(jj,:), 'linewidth', 2, 'displayname', sprintf('FWHM = %i mm (best fit)', fwhm_interest))
        text(num_modes+3, spectrum_noise{fwhm_interest_index}(end-3), sprintf('FWHM = %i mm (best fit)', fwhm_interest), 'fontsize', fontsize_axis-2, 'fontweight', 'normal', ...
            'horizontalalignment', 'left', 'verticalalignment', 'middle')
    else
        plot(1:num_modes, spectrum_noise{fwhm_interest_index}, 'color', color_shade_cmap(jj,:), 'linewidth', 2, 'displayname', sprintf('FWHM = %i mm', fwhm_interest))
        text(num_modes+3, spectrum_noise{fwhm_interest_index}(end-3), sprintf('FWHM = %i mm', fwhm_interest), 'fontsize', fontsize_axis-2, 'fontweight', 'normal', ...
            'horizontalalignment', 'left', 'verticalalignment', 'middle')
    end
end
hold off;
title('NeuroVault', 'fontsize', fontsize_axis, 'fontweight', 'normal')
leg = legend('fontsize', fontsize_axis, 'location', 'eastoutside', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/2;
set(leg, 'visible', 'off')
set(ax1_2, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', xlims, 'ylim', ylims, 'yscale', 'log')
xlabel('mode', 'fontsize', fontsize_label)
ylabel_h = ylabel('normalized power (log scale)', 'fontsize', fontsize_label);
set(ylabel_h, 'Position', ylabel_h.Position+[0 0.4 0])

% =========================================================================
% B: average MSLE
% =========================================================================   
ax2 = axes('Position', [init_x+factor_x*length_x*(1.95-1) init_y+factor_y*length_y*(1-1) length_x*0.67 0.96-init_y]);
hold on;
plot(fwhm_list, data_ExtFigure8.MSLE_HCP_average_map, 'ko-', 'linewidth', 2, 'markersize', 8, 'displayname', 'HCP');
plot(fwhm_list, data_ExtFigure8.MSLE_neurovault_average_map, 'k^-', 'linewidth', 2, 'markersize', 8, 'displayname', 'NeuroVault');
hold off;
leg = legend('fontsize', fontsize_axis, 'location', 'northwest', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/1;
set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'yscale', 'linear', 'xlim', [0 50])
xlabel('FWHM (mm)', 'fontsize', fontsize_label)
ylabel('mean square logarithmic error (MSLE)', 'fontsize', fontsize_label)

% =========================================================================
% C top: MSLE (HCP)
% =========================================================================   
cmap = cbrewer('seq', 'YlOrBr', 1000, 'pchip');

data_to_plot = data_ExtFigure8.MSLE_HCP_all_maps;
[~, min_ind] = min(data_to_plot,[],2);
[~, min_sort_ind] = sort(min_ind, 'ascend');
clims = prctile(data_to_plot(:), [2 98]);

ax3_1 = axes('Position', [init_x+factor_x*length_x*(2.7-1) init_y+factor_y*length_y*(2-1) length_x*0.45 length_y]);
imagesc(fwhm_list, 1:length(min_ind), data_to_plot(min_sort_ind, :))
hold on;
plot([fwhm_list(min_ind(min_sort_ind(1))), fwhm_list(min_ind(min_sort_ind)), fwhm_list(min_ind(min_sort_ind(end)))], ...
     [0.5, 1:length(min_ind), length(min_ind)+0.5], 'k-', 'linewidth', 2)
hold off;
caxis(clims)
colormap(ax3_1, cmap)
cbar = colorbar;
title('HCP', 'fontsize', fontsize_axis, 'fontweight', 'normal')
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'ytick', [], 'ticklabelinterpreter', 'none')
ylabel('map', 'fontsize', fontsize_label)
set(cbar, 'fontsize', fontsize_axis, 'ticklength', 0.02, ...
    'position', [ax3_1.Position(1)+ax3_1.Position(3)*1.5, ax3_1.Position(2), ax3_1.Position(3)*0.1, ax3_1.Position(4)])

% =========================================================================
% C bottom: MSLE (NeuroVault)
% =========================================================================  
data_to_plot = data_ExtFigure8.MSLE_neurovault_all_maps;
data_to_plot = data_to_plot(~isnan(sum(data_to_plot,2)), :);
[~, min_ind] = min(data_to_plot,[],2);
[~, min_sort_ind] = sort(min_ind, 'ascend');
clims = prctile(data_to_plot(:), [2 98]);

ax3_2 = axes('Position', [init_x+factor_x*length_x*(2.7-1) init_y+factor_y*length_y*(1-1) length_x*0.45 length_y]);
imagesc(fwhm_list, 1:length(min_ind), data_to_plot(min_sort_ind, :))
hold on;
plot([fwhm_list(min_ind(min_sort_ind(1))), fwhm_list(min_ind(min_sort_ind)), fwhm_list(min_ind(min_sort_ind(end)))], ...
     [0.5, 1:length(min_ind), length(min_ind)+0.5], 'k-', 'linewidth', 2)
hold off;
caxis(clims)
colormap(ax3_2, cmap)
cbar = colorbar;
title('NeuroVault', 'fontsize', fontsize_axis, 'fontweight', 'normal')
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'ytick', [], 'ticklabelinterpreter', 'none')
xlabel('FWHM (mm)', 'fontsize', fontsize_label)
ylabel('map', 'fontsize', fontsize_label)
ylabel_h = ylabel(cbar, 'mean square logarithmic error (MSLE)', 'fontsize', fontsize_label);
set(ylabel_h, 'Position', ylabel_h.Position+[0 0.9 0])
set(cbar, 'fontsize', fontsize_axis, 'ticklength', 0.02, ...
    'position', [ax3_2.Position(1)+ax3_2.Position(3)*1.5, ax3_2.Position(2), ax3_2.Position(3)*0.1, ax3_2.Position(4)])

%%% panel letters
annotation(fig, 'textbox', [0.015, 0.98, 0.01, 0.01], 'string', 'a', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.40, 0.98, 0.01, 0.01], 'string', 'b', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.735, 0.98, 0.01, 0.01], 'string', 'c', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% EXTENDED DATA FIGURE 9
% bioRxiv Supplementary Figure 19

% load all data relevant to Extended Data Figure 9
data_ExtFigure9 = load(sprintf('%s/ExtFigure9.mat', data_figures_folder));

parc_name = 'Glasser360';
parc = dlmread(filename_common_parcellation(parc_name, hemisphere));
parcels = unique(parc(parc>0));
num_parcels = length(unique(parc(parc>0)));

surface_to_plot = surface_midthickness;

medial_wall = find(~cortex);

type_name_list = {'data', 'wave model', 'mass model'};

fig = figure('Position', [200 200 800 1000]);

diff_y = 0.27;
init_x_1 = 0.02;
init_y_1 = 0.82;
init_x_2 = 0.13;
init_y_2 = init_y_1 - diff_y;
init_x_3 = init_x_2;
init_y_3 = init_y_2 - diff_y;
init_x_4 = init_x_2;
init_y_4 = init_y_3 - diff_y;

% =========================================================================
% A: lag matrices
% =========================================================================  
factor_x = 0.8;
init_x = init_x_1;
init_y = init_y_1+0.015;
length_x = (1 - 2.*init_x)/(factor_x*(3-1) + 1);
length_y = length_x*0.38;

for type_ind = 1:length(type_name_list)
    ax1 = axes('Position', [init_x+factor_x*length_x*(type_ind-1) init_y length_x length_y]);
    imagesc(data_ExtFigure9.peak_lags{type_ind})
    colormap(ax1, bluewhitered)
    cbar = colorbar('fontsize', fontsize_axis);
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xtick', [], 'ytick', [])
    xlabel('region', 'fontsize', fontsize_label)
    ylabel('region', 'fontsize', fontsize_label)
    if type_ind==3
        ylabel(cbar, 'time delay (s)', 'fontsize', fontsize_label)
    end
    title(type_name_list{type_ind}, 'fontsize', fontsize_label, 'fontweight', 'normal')
    axis image
end

% =========================================================================
% B left: mean lags on surface
% ========================================================================= 
data_to_plot_temp = {};
for type_ind = 1:length(type_name_list)
    data_to_plot_temp{type_ind} = data_ExtFigure9.peak_lags_mean{type_ind}';
end
factor_x_small = 1.05;
factor_x_big = 1.5;
factor_y_big = 1.02;
init_x = init_x_2;
init_y = init_y_2;
length_x = (0.4-init_x)/(factor_x_small + factor_x_big*(1-1) + 1);
length_y = (init_y_1 - 0.03 - init_y)/(factor_y_big*(length(type_name_list)-1) + 1);

for type_ind = 1:length(type_name_list)
    data_orig = zeros(num_vertices, 1);
    for ii=1:num_parcels
        ind_parcel = find(parc==parcels(ii));
        data_orig(ind_parcel,:) = repmat(data_to_plot_temp{type_ind}(ii), length(ind_parcel), 1);
    end

    boundary_method = 'midpoint';
    BOUNDARY = findROIboundaries(surface_to_plot.vertices,surface_to_plot.faces,parc,boundary_method);

    data_to_plot = data_orig;
    data_to_plot(medial_wall) = min(data_to_plot)*1.1;
    clims = [min(data_to_plot), max(data_to_plot)];
    if clims(2)<=0
        clims(2) = 0.01;
    end

    ax2_1 = axes('Position', [init_x+factor_x_big*length_x*(1-1) init_y+factor_y_big*length_y*(-[type_ind-length(type_name_list)]) length_x length_y]);
    obj1 = patch(ax2_1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
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
    colormap(ax2_1,[0.5,0.5,0.5; bluewhitered])
    axis off
    axis image
    
    annotation(fig, 'textbox', [0, ax2_1.Position(2), 0.132, length_y], 'string', {'mean lag'; sprintf('(%s)',type_name_list{type_ind})}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')

    ax2_2 = axes('Position', [init_x+factor_x_small*length_x+factor_x_big*length_x*(1-1) init_y+factor_y_big*length_y*(-[type_ind-length(type_name_list)]) length_x length_y]);
    obj2 = patch(ax2_2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
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
    colormap(ax2_2,[0.5,0.5,0.5; bluewhitered])
    axis off
    axis image
    
    if type_ind==3
        % colorbar
        ax2_temp = axes('position', [ax2_2.Position(1)-ax2_2.Position(3)*0.28, ax2_2.Position(2)+0.006, ax2_2.Position(3)*0.5, 0.006]);
        imagesc(linspace(-1,1,1000), ones(1,1000), linspace(-1,1,1000))
        set(ax2_temp, 'fontsize', fontsize_axis-2, 'ytick', [], 'xtick', [-1 0 1], 'xticklabel', {'neg', '0', 'pos'});
        colormap(ax2_temp,bluewhitered)
    end
end

% =========================================================================
% B right: mean lags correlation 
% =========================================================================
factor_x = 1.4;
init_x = length_x*(factor_x_small + factor_x_big*(1-1) + 1) + init_x_2 + 0.08;
init_y = init_y_2 + 0.05;
length_x = (0.98 - init_x)/(factor_x*(2-1) + 1);
length_y = length_y*(factor_y_big*(length(type_name_list)-1) + 1) - (init_y - init_y_2) - 0.01;

for type_ind = 2:length(type_name_list)
    ax3 = axes('Position', [init_x+factor_x*length_x*(type_ind-2) init_y length_x length_y]);
    data_to_plot_x = data_to_plot_temp{type_ind};
    data_to_plot_y = data_to_plot_temp{1};
    hold on;
    plot(data_to_plot_x, data_to_plot_y, 'k.', 'markersize', 15)
    plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
        'r-', 'linewidth', 2);
    hold off;
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', 1.2*[min(data_to_plot_x) max(data_to_plot_x)], 'ylim', [1.5 1.2].*[min(data_to_plot_y) max(data_to_plot_y)])
    xlabel(type_name_list{type_ind}, 'fontsize', fontsize_label)
    ylabel(type_name_list{1}, 'fontsize', fontsize_label)

    stat = 'linear';
    % stat = 'rank';
    if strcmpi(stat, 'rank')
        [rho,~] = corr(data_to_plot_x, data_to_plot_y, 'type', 'spearman');
    elseif strcmpi(stat, 'linear')
        [rho,~] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
    end
    pval = perm_sphere_p(data_to_plot_x, data_to_plot_y, perm_id, 'pearson');
    text(min(get(gca,'xlim'))+1*[max(get(gca,'xlim'))-min(get(gca,'xlim'))], min(get(gca,'ylim'))+0.2*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], ['r = ', sprintf('%.2f', rho)], 'color', 'r', ...
        'fontsize', fontsize_axis, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'right');
    text(min(get(gca,'xlim'))+1*[max(get(gca,'xlim'))-min(get(gca,'xlim'))], min(get(gca,'ylim'))+0.1*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], extract_pvalue_text(pval,1,'spin'), 'color', 'r', ...
        'fontsize', fontsize_axis, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'right');
    box off
end

% =========================================================================
% C left: PC1 lags on surface
% ========================================================================= 
component = 1;
data_to_plot_temp = {};
for type_ind = 1:length(type_name_list)
    data_to_plot_temp{type_ind} = data_ExtFigure9.score{type_ind}(:,component);
end
init_x = init_x_3;
init_y = init_y_3;
length_x = (0.4-init_x)/(factor_x_small + factor_x_big*(1-1) + 1);
length_y = (init_y_2 - 0.03 - init_y)/(factor_y_big*(length(type_name_list)-1) + 1);

for type_ind = 1:length(type_name_list)
    data_orig = zeros(num_vertices, 1);
    for ii=1:num_parcels
        ind_parcel = find(parc==parcels(ii));
        data_orig(ind_parcel,:) = repmat(data_to_plot_temp{type_ind}(ii), length(ind_parcel), 1);
    end

    boundary_method = 'midpoint';
    BOUNDARY = findROIboundaries(surface_to_plot.vertices,surface_to_plot.faces,parc,boundary_method);

    data_to_plot = data_orig;
    data_to_plot(medial_wall) = min(data_to_plot)*1.1;
    clims = [min(data_to_plot), max(data_to_plot)];
    if clims(2)<=0
        clims(2) = 0.01;
    end

    ax2_1 = axes('Position', [init_x+factor_x_big*length_x*(1-1) init_y+factor_y_big*length_y*(-[type_ind-length(type_name_list)]) length_x length_y]);
    obj1 = patch(ax2_1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
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
    colormap(ax2_1,[0.5,0.5,0.5; bluewhitered])
    axis off
    axis image
    
    annotation(fig, 'textbox', [0, ax2_1.Position(2), 0.132, length_y], 'string', {'PC1 lag'; sprintf('(%s)',type_name_list{type_ind})}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')

    ax2_2 = axes('Position', [init_x+factor_x_small*length_x+factor_x_big*length_x*(1-1) init_y+factor_y_big*length_y*(-[type_ind-length(type_name_list)]) length_x length_y]);
    obj2 = patch(ax2_2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
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
    colormap(ax2_2,[0.5,0.5,0.5; bluewhitered])
    axis off
    axis image
    
    if type_ind==3
        % colorbar
        ax2_temp = axes('position', [ax2_2.Position(1)-ax2_2.Position(3)*0.28, ax2_2.Position(2)+0.006, ax2_2.Position(3)*0.5, 0.006]);
        imagesc(linspace(-1,1,1000), ones(1,1000), linspace(-1,1,1000))
        set(ax2_temp, 'fontsize', fontsize_axis-2, 'ytick', [], 'xtick', [-1 0 1], 'xticklabel', {'neg', '0', 'pos'});
        colormap(ax2_temp,bluewhitered)
    end
    
    annotation(fig, 'textbox', [ax2_1.Position(1), ax2_1.Position(2)+ax2_1.Position(4)*0.85, ax2_2.Position(1)-ax2_1.Position(1)+ax2_1.Position(3), 0.01], 'string', sprintf('var = %.1f%%', data_ExtFigure9.explained{type_ind}(component)), 'edgecolor', 'none', ...
        'fontsize', fontsize_axis-2, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'bottom', 'interpreter', 'none')
end

% =========================================================================
% C right: PC1 lags correlation 
% =========================================================================
factor_x = 1.4;
init_x = length_x*(factor_x_small + factor_x_big*(1-1) + 1) + init_x_3 + 0.08;
init_y = init_y_3 + 0.05;
length_x = (0.98 - init_x)/(factor_x*(2-1) + 1);
length_y = length_y*(factor_y_big*(length(type_name_list)-1) + 1) - (init_y - init_y_3) - 0.01;

for type_ind = 2:length(type_name_list)
    ax3 = axes('Position', [init_x+factor_x*length_x*(type_ind-2) init_y length_x length_y]);
    data_to_plot_x = data_to_plot_temp{type_ind};
    data_to_plot_y = data_to_plot_temp{1};
    hold on;
    plot(data_to_plot_x, data_to_plot_y, 'k.', 'markersize', 15)
    plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
        'r-', 'linewidth', 2);
    hold off;
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', 1.2*[min(data_to_plot_x) max(data_to_plot_x)], 'ylim', [1.5 1.2].*[min(data_to_plot_y) max(data_to_plot_y)])
    xlabel(type_name_list{type_ind}, 'fontsize', fontsize_label)
    ylabel(type_name_list{1}, 'fontsize', fontsize_label)

    stat = 'linear';
    % stat = 'rank';
    if strcmpi(stat, 'rank')
        [rho,~] = corr(data_to_plot_x, data_to_plot_y, 'type', 'spearman');
    elseif strcmpi(stat, 'linear')
        [rho,~] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
    end
    pval = perm_sphere_p(data_to_plot_x, data_to_plot_y, perm_id, 'pearson');
    text(min(get(gca,'xlim'))+1*[max(get(gca,'xlim'))-min(get(gca,'xlim'))], min(get(gca,'ylim'))+0.2*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], ['r = ', sprintf('%.2f', rho)], 'color', 'r', ...
        'fontsize', fontsize_axis, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'right');
    text(min(get(gca,'xlim'))+1*[max(get(gca,'xlim'))-min(get(gca,'xlim'))], min(get(gca,'ylim'))+0.1*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], extract_pvalue_text(pval,1,'spin'), 'color', 'r', ...
        'fontsize', fontsize_axis, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'right');
    box off
end

% =========================================================================
% D left: PC2 lags on surface 
% =========================================================================
component = 2;
data_to_plot_temp = {};
for type_ind = 1:length(type_name_list)
    data_to_plot_temp{type_ind} = data_ExtFigure9.score{type_ind}(:,component);
end
init_x = init_x_4;
init_y = init_y_4;
length_x = (0.4-init_x)/(factor_x_small + factor_x_big*(1-1) + 1);
length_y = (init_y_3 - 0.03 - init_y)/(factor_y_big*(length(type_name_list)-1) + 1);

for type_ind = 1:length(type_name_list)
    data_orig = zeros(num_vertices, 1);
    for ii=1:num_parcels
        ind_parcel = find(parc==parcels(ii));
        data_orig(ind_parcel,:) = repmat(data_to_plot_temp{type_ind}(ii), length(ind_parcel), 1);
    end

    boundary_method = 'midpoint';
    BOUNDARY = findROIboundaries(surface_to_plot.vertices,surface_to_plot.faces,parc,boundary_method);

    data_to_plot = data_orig;
    data_to_plot(medial_wall) = min(data_to_plot)*1.1;
    clims = [min(data_to_plot), max(data_to_plot)];
    if clims(2)<=0
        clims(2) = 0.01;
    end

    ax2_1 = axes('Position', [init_x+factor_x_big*length_x*(1-1) init_y+factor_y_big*length_y*(-[type_ind-length(type_name_list)]) length_x length_y]);
    obj1 = patch(ax2_1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
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
    colormap(ax2_1,[0.5,0.5,0.5; bluewhitered])
    axis off
    axis image
    
    annotation(fig, 'textbox', [0, ax2_1.Position(2), 0.132, length_y], 'string', {'PC2 lag'; sprintf('(%s)',type_name_list{type_ind})}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')

    ax2_2 = axes('Position', [init_x+factor_x_small*length_x+factor_x_big*length_x*(1-1) init_y+factor_y_big*length_y*(-[type_ind-length(type_name_list)]) length_x length_y]);
    obj2 = patch(ax2_2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot, ...
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
    colormap(ax2_2,[0.5,0.5,0.5; bluewhitered])
    axis off
    axis image
    
    if type_ind==3
        % colorbar
        ax2_temp = axes('position', [ax2_2.Position(1)-ax2_2.Position(3)*0.28, ax2_2.Position(2)+0.006, ax2_2.Position(3)*0.5, 0.006]);
        imagesc(linspace(-1,1,1000), ones(1,1000), linspace(-1,1,1000))
        set(ax2_temp, 'fontsize', fontsize_axis-2, 'ytick', [], 'xtick', [-1 0 1], 'xticklabel', {'neg', '0', 'pos'});
        colormap(ax2_temp,bluewhitered)
    end
    
    annotation(fig, 'textbox', [ax2_1.Position(1), ax2_1.Position(2)+ax2_1.Position(4)*0.85, ax2_2.Position(1)-ax2_1.Position(1)+ax2_1.Position(3), 0.01], 'string', sprintf('var = %.1f%%', data_ExtFigure9.explained{type_ind}(component)), 'edgecolor', 'none', ...
        'fontsize', fontsize_axis-2, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'bottom', 'interpreter', 'none')
end

% =========================================================================
% D right: PC2 lags correlation 
% =========================================================================
factor_x = 1.4;
init_x = length_x*(factor_x_small + factor_x_big*(1-1) + 1) + init_x_4 + 0.08;
init_y = init_y_4 + 0.05;
length_x = (0.98 - init_x)/(factor_x*(2-1) + 1);
length_y = length_y*(factor_y_big*(length(type_name_list)-1) + 1) - (init_y - init_y_4) - 0.01;

for type_ind = 2:length(type_name_list)
    ax3 = axes('Position', [init_x+factor_x*length_x*(type_ind-2) init_y length_x length_y]);
    data_to_plot_x = data_to_plot_temp{type_ind};
    data_to_plot_y = data_to_plot_temp{1};
    hold on;
    plot(data_to_plot_x, data_to_plot_y, 'k.', 'markersize', 15)
    plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
        'r-', 'linewidth', 2);
    hold off;
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', 1.2*[min(data_to_plot_x) max(data_to_plot_x)], 'ylim', [1.5 1.2].*[min(data_to_plot_y) max(data_to_plot_y)])
    xlabel(type_name_list{type_ind}, 'fontsize', fontsize_label)
    ylabel(type_name_list{1}, 'fontsize', fontsize_label)

    stat = 'linear';
    % stat = 'rank';
    if strcmpi(stat, 'rank')
        [rho,~] = corr(data_to_plot_x, data_to_plot_y, 'type', 'spearman');
    elseif strcmpi(stat, 'linear')
        [rho,~] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
    end
    pval = perm_sphere_p(data_to_plot_x, data_to_plot_y, perm_id, 'pearson');
    text(min(get(gca,'xlim'))+1*[max(get(gca,'xlim'))-min(get(gca,'xlim'))], min(get(gca,'ylim'))+0.2*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], ['r = ', sprintf('%.2f', rho)], 'color', 'r', ...
        'fontsize', fontsize_axis, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'right');
    text(min(get(gca,'xlim'))+1*[max(get(gca,'xlim'))-min(get(gca,'xlim'))], min(get(gca,'ylim'))+0.1*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], extract_pvalue_text(pval,1,'spin'), 'color', 'r', ...
        'fontsize', fontsize_axis, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'right');
    box off
end

%%% panel letters
annotation(fig, 'textbox', [0.015, 0.99, 0.01, 0.01], 'string', 'a', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.015, 0.8, 0.01, 0.01], 'string', 'b', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.015, 0.53, 0.01, 0.01], 'string', 'c', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.015, 0.26, 0.01, 0.01], 'string', 'd', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% EXTENDED DATA FIGURE 10
% bioRxiv Supplementary Figure 20

cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];

fig = figure;
hold on;
plot(model_wave_optim.rs_vec, model_wave_optim.optim_edge_FC, 'o-', 'color', colors2(1,:), 'linewidth', 2, 'markersize', 8, 'displayname', 'edge FC correlation')
plot(model_wave_optim.rs_vec, model_wave_optim.optim_node_FC, 'o-', 'color', colors2(2,:), 'linewidth', 2, 'markersize', 8, 'displayname', 'node FC correlation')
plot(model_wave_optim.rs_vec, model_wave_optim.optim_FCD, 'o-', 'color', colors2(4,:), 'linewidth', 2, 'markersize', 8, 'displayname', 'FCD KS statistic')
hold off;
leg = legend('fontsize', fontsize_label, 'location', 'southeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/1;
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'ylim', [0 1], 'xlim', [min(model_wave_optim.rs_vec), max(model_wave_optim.rs_vec)])
xlabel('$r_s$', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('optimization metric', 'fontsize', fontsize_label)
