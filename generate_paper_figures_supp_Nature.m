%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate_paper_figures_supp_Nature.m
%%%
%%% MATLAB script to generate the final Supplementary figures of the paper
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

%% SUPPLEMENTARY FIGURE 1
% bioRxiv Supplementary Figure 2

cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];

fig = figure;
hold on;
for ii=1:length(tasks)
    task = tasks{ii};
    
    num_contrasts = length(task_contrasts_ind.(task));
    
    for jj=1:num_contrasts
        contrast = contrasts{task_contrasts_ind.(task)(jj)};
        if jj==1
            plot(nanmean(task_recon_corr_geometric.(contrast),1), 'color', colors(ii,:), 'linewidth', 2, 'displayname', lower(tasks{ii}))
        else
            plot(nanmean(task_recon_corr_geometric.(contrast),1), 'color', colors(ii,:), 'linewidth', 2, 'HandleVisibility', 'off')
        end
    end
end
hold off;
leg = legend('fontsize', fontsize_axis, 'location', 'southeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/2;
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1], 'xtick', 50:50:200)
xlabel('number of modes', 'fontsize', fontsize_label)
ylabel('reconstruction accuracy', 'fontsize', fontsize_label)

%% SUPPLEMENTARY FIGURE 2
% bioRxiv Supplementary Figure 3

% load all data relevant to Supplementary Figure 2
data_SuppFigure2 = load(sprintf('%s/SuppFigure2.mat', data_figures_folder));

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
for ii=1:length(representative_contrasts)
    contrast = representative_contrasts{ii};
    task_name = lower(tasks{ii});
    
    color_shade_cmap = repmat(colors(ii,:),length(parc_name_list),1);
    brighten_vec = linspace(0,0.85,length(parc_name_list));
    for parc_name_ind=1:length(parc_name_list)
        color_shade_cmap(parc_name_ind,:) = brighten(color_shade_cmap(parc_name_ind,:), brighten_vec(parc_name_ind));
    end
    color_shade_cmap = flipud(color_shade_cmap);
    
    if ii<=4
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
    end
    hold on;
    for parc_name_ind = 1:length(parc_name_list)
        parc_name_interest = parc_name_list{parc_name_ind};
        
        plot(data_SuppFigure2.task_recon_corr_parcellations.(parc_name_interest).(contrast), 'color', color_shade_cmap(parc_name_ind,:), 'linewidth', 2, 'displayname', parc_name_interest)
    end
    hold off;
    leg = legend('fontsize', fontsize_axis, 'location', 'southeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
    leg.ItemTokenSize = leg.ItemTokenSize/2;
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1])
    if ii>4
        xlabel('number of modes', 'fontsize', fontsize_label)
    end
    if ii==1 || ii==5
        ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
    end
    title(task_name, 'fontsize', fontsize_axis, 'fontweight', 'normal')
    
    if ii==1
        annotation('textarrow', ax1.Position(1)+ax1.Position(3).*[0.7, 0.05], ...
                                ax1.Position(2)+ax1.Position(4).*[0.72, 0.92], 'String', {'decreasing'; 'resolution'}, ...
                                'fontsize', fontsize_axis, 'horizontalalignment', 'center');
    else
        annotation('textarrow', ax1.Position(1)+ax1.Position(3).*[0.7, 0.05], ...
                                ax1.Position(2)+ax1.Position(4).*[0.72, 0.92], 'String', {''}, ...
                                'fontsize', fontsize_axis, 'horizontalalignment', 'center');
    end
end

% rest
ii = 8;
brighten_vec = linspace(0,0.8,length(parc_name_list));
color_shade_cmap = repmat(brighten_vec', 1, 3); 
color_shade_cmap = flipud(color_shade_cmap);

ax1 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
hold on;
for parc_name_ind = 1:length(parc_name_list)
    parc_name_interest = parc_name_list{parc_name_ind};
    
    plot(data_SuppFigure2.resting_recon_corr_parcellations.(parc_name_interest), 'color', color_shade_cmap(parc_name_ind,:), 'linewidth', 2, 'displayname', parc_name_interest)
end
hold off;
leg = legend('fontsize', fontsize_axis, 'location', 'southeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/2;
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1])
xlabel('number of modes', 'fontsize', fontsize_label)
if ii==1
    ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
end
title('rest', 'fontsize', fontsize_axis, 'fontweight', 'normal')

annotation('textarrow', ax1.Position(1)+ax1.Position(3).*[0.7, 0.05], ...
                                ax1.Position(2)+ax1.Position(4).*[0.735, 0.9], 'String', {''}, ...
                                'fontsize', fontsize_axis, 'horizontalalignment', 'center');
                            
%% SUPPLEMENTARY FIGURE 3
% bioRxiv Supplementary Figure 4

% load all data relevant to Supplementary Figure 3
data_SuppFigure3 = load(sprintf('%s/SuppFigure3.mat', data_figures_folder));

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
for ii=1:(length(representative_contrasts)+1)
    if ii>length(representative_contrasts)
        titlename = 'rest';
        data_to_plot_template = nanmean(resting_recon_corr_geometric,1);
        data_to_plot_individual = data_SuppFigure3.resting_recon_corr_individual;
    else
        contrast = representative_contrasts{ii};
        task = lower(tasks{ii});
        titlename = task;
        data_to_plot_template = nanmean(task_recon_corr_geometric.(contrast),1);
        data_to_plot_individual = data_SuppFigure3.task_recon_corr_individual.(contrast);
    end
    
    if ii<=4
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
    end
    hold on;
    plot(1:num_modes, data_to_plot_template, 'color', colors(ii,:), 'linewidth', 2, 'displayname', 'template')
    plot(1:num_modes, data_to_plot_individual, '--', 'color', colors(ii,:), 'linewidth', 2, 'displayname', 'individual-specific')
    hold off;
    if ii==length(representative_contrasts)+1
        leg = legend('fontsize', fontsize_axis-1, 'interpreter', 'none', 'box', 'off', 'numcolumns', 1, ...
                    'Position', [ax1.Position(1)+ax1.Position(3)*0.5 ax1.Position(2)+ax1.Position(4)*0.61 ax1.Position(3)*0.5 ax1.Position(3)*0.2]);
        leg.ItemTokenSize = leg.ItemTokenSize/1.3;
    end
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1], 'xtick', 50:50:200)
    if ii>4
        xlabel('number of modes', 'fontsize', fontsize_label)
    end
    if ii==1 ||ii==5
        ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
    end
    title(titlename, 'fontsize', fontsize_axis, 'fontweight', 'normal', 'interpreter', 'none')
    
    ax2 = axes('Position', [ax1.Position(1)+ax1.Position(3)*0.42 ax1.Position(2)+ax1.Position(4)*0.18 ax1.Position(3)*0.5 ax1.Position(4)*0.3]);
    plot(1:num_modes, data_to_plot_template-data_to_plot_individual, 'k-', 'linewidth', 2)
    yline(0, 'k:');
    set(ax2, 'fontsize', fontsize_axis-2, 'ticklength', [0.02 0.02], 'xlim', [2 num_modes], 'ylim', 0.05*[-1 1], 'xtick', 100:100:200)
    xlabel('number of modes', 'fontsize', fontsize_label-2)
    ylabel('difference', 'fontsize', fontsize_label-2)
end

%% SUPPLEMENTARY FIGURE 4
% bioRxiv Supplementary Figure 5

% load all data relevant to Supplementary Figure 4
data_SuppFigure4 = load(sprintf('%s/SuppFigure4.mat', data_figures_folder));

n = 200;

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
for ii=1:(length(representative_contrasts)+1)
    if ii>length(representative_contrasts)
        titlename = 'rest';
        data_to_plot_template = resting_recon_corr_geometric(:,n);
        data_to_plot_individual = data_SuppFigure4.resting_recon_corr_individual;
    else
        contrast = representative_contrasts{ii};
        task = lower(tasks{ii});
        titlename = task;
        data_to_plot_template = task_recon_corr_geometric.(contrast)(:,n);
        data_to_plot_individual = data_SuppFigure4.task_recon_corr_individual.(contrast);
    end
    
    if ii<=4
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
    end
    hold on;
    plot([0 1], [0 1], 'k:')
    plot(data_to_plot_individual, data_to_plot_template, '.', 'color', colors(ii,:), 'markersize', 10)
    hold off;
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [0.8 1], 'ylim', [0.8 1], 'xtick', 0.8:0.1:1, 'ytick', 0.8:0.1:1)
    if ii>4
        xlabel('individual-specific accuracy', 'fontsize', fontsize_label)
    end
    if ii==1 ||ii==5
        ylabel('template accuracy', 'fontsize', fontsize_label)
    end
    title(titlename, 'fontsize', fontsize_axis, 'fontweight', 'normal', 'interpreter', 'none')

    ax2 = axes('Position', [ax1.Position(1)+ax1.Position(3)*0.7 ax1.Position(2)+ax1.Position(4)*0.15 ax1.Position(3)*0.3 ax1.Position(4)*0.22]);
    data_to_plot_diff = data_to_plot_template-data_to_plot_individual;
    limits = [min(data_to_plot_diff), max(data_to_plot_diff)];
    histogram(data_to_plot_diff, 12, 'facecolor', 0.7*ones(1,3))
    set(ax2, 'fontsize', fontsize_axis-2, 'ticklength', [0.02 0.02], 'xlim', max(abs(limits))*[-1 1], ...
        'xtick', [-max(abs(limits)), 0, max(abs(limits))], 'xticklabel', {sprintf('%.2f', -max(abs(limits))), '0', sprintf('%.2f', max(abs(limits)))})
    xlabel('difference', 'fontsize', fontsize_label-2)
    ylabel('count', 'fontsize', fontsize_label-2)
    box on
end

%% SUPPLEMENTARY FIGURE 5
% bioRxiv Supplementary Figure 6

% load all data relevant to Supplementary Figure 5
data_SuppFigure5 = load(sprintf('%s/SuppFigure5.mat', data_figures_folder));

num_modes_temp = 500;
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
for ii=1:(length(representative_contrasts)+1)
    if ii>length(representative_contrasts)
        titlename = 'rest';
        data_to_plot_template = data_SuppFigure5.resting_recon_corr_geometric_500;
        data_to_plot_individual = data_SuppFigure5.resting_recon_corr_individual_500;
    else
        contrast = representative_contrasts{ii};
        task = lower(tasks{ii});
        titlename = task;
        data_to_plot_template = data_SuppFigure5.task_recon_corr_geometric_500.(contrast);
        data_to_plot_individual = data_SuppFigure5.task_recon_corr_individual_500.(contrast);
    end
    
    if ii<=4
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
    end
    hold on;
    plot(1:num_modes_temp, data_to_plot_template, 'color', colors(ii,:), 'linewidth', 2, 'displayname', 'template')
    plot(1:num_modes_temp, data_to_plot_individual, '--', 'color', colors(ii,:), 'linewidth', 2, 'displayname', 'individual-specific')
    hold off;
    if ii==length(representative_contrasts)+1
        leg = legend('fontsize', fontsize_axis-1, 'interpreter', 'none', 'box', 'off', 'numcolumns', 1, ...
                    'Position', [ax1.Position(1)+ax1.Position(3)*0.5 ax1.Position(2)+ax1.Position(4)*0.61 ax1.Position(3)*0.5 ax1.Position(3)*0.2]);
        leg.ItemTokenSize = leg.ItemTokenSize/1.35;
    end
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [200 num_modes_temp], 'ylim', [0.8 1], 'xtick', 200:100:500, 'ytick', 0.8:0.1:1)
    if ii>4
        xlabel('number of modes', 'fontsize', fontsize_label)
    end
    if ii==1 ||ii==5
        ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
    end
    title(titlename, 'fontsize', fontsize_axis, 'fontweight', 'normal', 'interpreter', 'none')
    
    ax2 = axes('Position', [ax1.Position(1)+ax1.Position(3)*0.42 ax1.Position(2)+ax1.Position(4)*0.18 ax1.Position(3)*0.5 ax1.Position(4)*0.3]);
    plot(1:num_modes_temp, data_to_plot_template-data_to_plot_individual, 'k-', 'linewidth', 2)
    yline(0, 'k:');
    set(ax2, 'fontsize', fontsize_axis-2, 'ticklength', [0.02 0.02], 'xlim', [200 num_modes_temp], 'ylim', 0.025*[-1 1], 'xtick', 200:100:500)
    xlabel('number of modes', 'fontsize', fontsize_label-2)
    ylabel('difference', 'fontsize', fontsize_label-2)
end

%% SUPPLEMENTARY FIGURE 6
% bioRxiv Supplementary Figure 9

% load all data relevant to Supplementary Figure 6
data_SuppFigure6 = load(sprintf('%s/SuppFigure6.mat', data_figures_folder));

density_list = 1:13;
density_interest_list = 1:2:13;
density_interest_list_ind = dsearchn(density_list', density_interest_list');

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
for ii=1:length(representative_contrasts)
    contrast = representative_contrasts{ii};
    task_name = lower(tasks{ii});
    
    color_shade_cmap = repmat(colors(ii,:),length(density_interest_list),1);
    brighten_vec = linspace(0,0.85,length(density_interest_list));
    for density_ind=1:length(density_interest_list)
        color_shade_cmap(density_ind,:) = brighten(color_shade_cmap(density_ind,:), brighten_vec(density_ind));
    end
    color_shade_cmap = flipud(color_shade_cmap);
    
    if ii<=4
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
    end
    hold on;
    plot(data_SuppFigure6.task_recon_corr_geometric.(contrast), '--', 'color', 'k', 'linewidth', 2, 'displayname', 'geometric')
   
    for density_ind=1:length(density_interest_list)
        density = density_interest_list(density_ind);
        density_num = density_interest_list_ind(density_ind);
        
        plot(data_SuppFigure6.task_recon_corr_connectome_density.(contrast)(density_num,:), 'color', color_shade_cmap(density_ind,:), 'linewidth', 2, 'displayname', sprintf('density = %.1f%%', density))
    end
    hold off;
    leg = legend('fontsize', fontsize_axis-1, 'location', 'southeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
    leg.ItemTokenSize = leg.ItemTokenSize/1.4;
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
brighten_vec = linspace(0,0.8,length(density_interest_list));
color_shade_cmap = repmat(brighten_vec', 1, 3); 
color_shade_cmap = flipud(color_shade_cmap);

ax1 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
hold on;
plot(data_SuppFigure6.resting_recon_corr_geometric, '--', 'color', 'k', 'linewidth', 2, 'displayname', 'geometric')

for density_ind=1:length(density_interest_list)
    density = density_interest_list(density_ind);
    density_num = density_interest_list_ind(density_ind);

    plot(data_SuppFigure6.resting_recon_corr_connectome_density(density_num,:), 'color', color_shade_cmap(density_ind,:), 'linewidth', 2, 'displayname', sprintf('density = %.1f%%', density))
end
hold off;
leg = legend('fontsize', fontsize_axis-1, 'location', 'southeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/1.4;
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1])
if ii>4
    xlabel('number of modes', 'fontsize', fontsize_label)
end
if ii==1 || ii==5
    ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
end
title('rest', 'fontsize', fontsize_axis, 'fontweight', 'normal')

%% SUPPLEMENTARY FIGURE 7
% bioRxiv Supplementary Figure 10

% load all data relevant to Supplementary Figure 7
data_SuppFigure7 = load(sprintf('%s/SuppFigure7.mat', data_figures_folder));

density_list = 0.1:0.1:5;
density_interest_list = [0.1, 1.55];
density_interest_ind = dsearchn(density_list', density_interest_list');

cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];

factor_x = 1.3;
factor_y = 1.3;
init_x = 0.06;
init_y = 0.07;
length_x = (1 - 1.4*init_x)/(factor_x*(4-1) + 1);
length_y = (1 - 1.8*init_y)/(factor_y*(2-1) + 1);

fig = figure('Position', [200 200 1000 600]);
for ii=1:length(representative_contrasts)
    contrast = representative_contrasts{ii};
    task_name = lower(tasks{ii});
    
    if ii<=4
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
    end
    
    data_to_plot = data_SuppFigure7.task_recon_corr_connectome_density.(contrast);
    
    hold on;
    plot(density_interest_list(1)*ones(1,2), [min(data_to_plot), data_to_plot(density_interest_ind(1))-0.001], 'k:', 'linewidth', 1.5);
    plot(density_interest_list(2)*ones(1,2), [min(data_to_plot), data_to_plot(density_interest_ind(2))-0.001], 'k:', 'linewidth', 1.5);
    plot(density_list, data_to_plot, '.-', 'color', colors(ii,:), 'markersize', 25)
    yline(data_SuppFigure7.task_recon_corr_geometric.(contrast), 'k-', 'geometric', 'linewidth', 1.5, 'fontsize', fontsize_axis, 'labelverticalalignment', 'bottom');
    hold off;
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [-0.2, 5], 'xtick', [0,1:5], ...
        'ylim', [min(data_to_plot), max(get(gca,'ylim'))])
    if ii>4
        xlabel('density (%)', 'fontsize', fontsize_label)
    end
    if ii==1 || ii==5
        ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
    end
    title(task_name, 'fontsize', fontsize_axis, 'fontweight', 'normal')
end

% rest
ii = 8;
data_to_plot = data_SuppFigure7.resting_recon_corr_connectome_density;

ax1 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
hold on;
plot(density_interest_list(1)*ones(1,2), [min(data_to_plot), data_to_plot(density_interest_ind(1))-0.001], 'k:', 'linewidth', 1.5);
plot(density_interest_list(2)*ones(1,2), [min(data_to_plot), data_to_plot(density_interest_ind(2))-0.001], 'k:', 'linewidth', 1.5);
plot(density_list, data_to_plot, '.-', 'color', colors(ii,:), 'markersize', 25)
yline(data_SuppFigure7.resting_recon_corr_geometric, 'k-', 'geometric', 'linewidth', 1.5, 'fontsize', fontsize_axis, 'labelverticalalignment', 'bottom');
hold off;
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [-0.2, 5], 'xtick', [0,1:5], ...
        'ylim', [min(data_to_plot), max(get(gca,'ylim'))])
if ii>4
    xlabel('density (%)', 'fontsize', fontsize_label)
end
if ii==1 || ii==5
    ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
end
title('rest', 'fontsize', fontsize_axis, 'fontweight', 'normal')

%% SUPPLEMENTARY FIGURE 8
% bioRxiv Supplementary Figure 11

% load all data relevant to Supplementary Figure 8
data_SuppFigure8 = load(sprintf('%s/SuppFigure8.mat', data_figures_folder));

parc_name_list = {'Schaefer400', 'Schaefer600', 'Schaefer800', 'Schaefer1000'};

cmap1 = lines(7);
cmap2 = cbrewer('qual', 'Set1', 8, 'pchip');
colors = [cmap1(4,:); cmap2(3,:); cmap1(6,:); cmap2(1,:); cmap1(3,:); cmap2(7,:); cmap2(8,:); 0 0 0];
colors2 = [0 0 0; cmap1([4,5,1:3,6:7],:)];

factor_x = 1.3;
factor_y = 1.3;
init_x = 0.05;
init_y = 0.08;
length_x = (1 - 1.5*init_x)/(factor_x*(4-1) + 1);
length_y = (1 - 1.8*init_y)/(factor_y*(2-1) + 1);

fig = figure('Position', [200 200 1000 600]);
for ii=1:length(representative_contrasts)
    contrast = representative_contrasts{ii};
    task_name = lower(tasks{ii});

    color_shade_cmap = repmat(colors(ii,:),length(parc_name_list),1);
    brighten_vec = linspace(0,0.85,length(parc_name_list));
    for parc_name_ind = 1:length(parc_name_list)
        color_shade_cmap(parc_name_ind,:) = brighten(color_shade_cmap(parc_name_ind,:), brighten_vec(parc_name_ind));
    end
    color_shade_cmap = flipud(color_shade_cmap);

    if ii<=4
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-1) init_y+factor_y*length_y*(2-1) length_x length_y]);
    else
        ax1 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
    end
    hold on;
    data_to_plot = data_SuppFigure8.task_recon_corr_geometric.(contrast);
    data_to_plot_x = 100*(1:num_modes)/sum(cortex);
    plot(data_to_plot_x, data_to_plot, '--', 'color', 'k', 'linewidth', 2, 'displayname', 'geometric')

    for parc_name_ind = 1:length(parc_name_list)
        parc_name_interest = parc_name_list{parc_name_ind};

        parc = dlmread(filename_common_parcellation(parc_name_interest, hemisphere));
        parcels = unique(parc(parc>0));
        num_parcels = length(unique(parc(parc>0)));

        if num_parcels<200
            num_modes = num_parcels;
        else
            num_modes = 200;
        end  

        data_to_plot = data_SuppFigure8.task_recon_corr_discrete_connectome.(parc_name_interest).(contrast);
        data_to_plot_x = 100*(1:num_modes)/num_parcels;

        plot(data_to_plot_x, data_to_plot, 'color', color_shade_cmap(parc_name_ind,:), 'linewidth', 2, 'displayname', parc_name_interest)
    end
    hold off;
    leg = legend('fontsize', fontsize_axis-1, 'location', 'southeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
    leg.ItemTokenSize = leg.ItemTokenSize/1.4;
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02 0.02], 'xlim', [1e-4, 1e2], 'ylim', [0 1], 'xscale', 'linear')
    if ii>4
        xlabel('% of available modes', 'fontsize', fontsize_label)
    end
    if ii==1 || ii==5
        ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
    end
    title(task_name, 'fontsize', fontsize_axis, 'fontweight', 'normal')

    ax1_inset = axes('Position', [ax1.Position(1)+ax1.Position(3)*0.5 ax1.Position(2)+ax1.Position(4)*0.55 length_x*0.48 length_y*0.15]);
    data_to_plot = data_SuppFigure8.task_recon_corr_geometric.(contrast);
    data_to_plot_x = 100*(1:num_modes)/sum(cortex);
    plot(data_to_plot_x, data_to_plot, '--', 'color', 'k', 'linewidth', 2, 'displayname', 'geometric')
    set(gca, 'fontsize', fontsize_axis-2, 'ticklength', [0.02 0.02], 'xtick', [0.2, 0.4, 0.6], 'ylim', [0 1], 'xscale', 'linear')
    xlabel('%', 'fontsize', fontsize_label-2)
    box on
end

% rest
ii = 8;
brighten_vec = linspace(0,0.8,length(parc_name_list));
color_shade_cmap = repmat(brighten_vec', 1, 3); 
color_shade_cmap = flipud(color_shade_cmap);

ax1 = axes('Position', [init_x+factor_x*length_x*(ii-4-1) init_y+factor_y*length_y*(1-1) length_x length_y]);
hold on;
data_to_plot = data_SuppFigure8.resting_recon_corr_geometric;
data_to_plot_x = 100*(1:num_modes)/sum(cortex);
plot(data_to_plot_x, data_to_plot, '--', 'color', 'k', 'linewidth', 2, 'displayname', 'geometric')

for parc_name_ind=1:length(parc_name_list)
    parc_name_interest = parc_name_list{parc_name_ind};

    parc = dlmread(filename_common_parcellation(parc_name_interest, hemisphere));
    parcels = unique(parc(parc>0));
    num_parcels = length(unique(parc(parc>0)));

    if num_parcels<200
        num_modes = num_parcels;
    else
        num_modes = 200;
    end  

    data_to_plot = data_SuppFigure8.resting_recon_corr_discrete_connectome.(parc_name_interest);
    data_to_plot_x = 100*(1:num_modes)/num_parcels;

    plot(data_to_plot_x, data_to_plot, 'color', color_shade_cmap(parc_name_ind,:), 'linewidth', 2, 'displayname', parc_name_interest)
end
hold off;
leg = legend('fontsize', fontsize_axis-1, 'location', 'southeast', 'interpreter', 'none', 'box', 'off', 'numcolumns', 1);
leg.ItemTokenSize = leg.ItemTokenSize/1.4;
if ii>4
    xlabel('% of available modes', 'fontsize', fontsize_label)
end
if ii==1 || ii==5
    ylabel('reconstruction accuracy', 'fontsize', fontsize_label)
end
title('rest', 'fontsize', fontsize_axis, 'fontweight', 'normal')

ax1_inset = axes('Position', [ax1.Position(1)+ax1.Position(3)*0.5 ax1.Position(2)+ax1.Position(4)*0.55 length_x*0.48 length_y*0.15]);
data_to_plot = data_SuppFigure8.resting_recon_corr_geometric;
data_to_plot_x = 100*(1:num_modes)/sum(cortex);
plot(data_to_plot_x, data_to_plot, '--', 'color', 'k', 'linewidth', 2, 'displayname', 'geometric')
set(gca, 'fontsize', fontsize_axis-2, 'ticklength', [0.02 0.02], 'xtick', [0.2, 0.4, 0.6], 'ylim', [0 1], 'xscale', 'linear')
xlabel('%', 'fontsize', fontsize_label-2)
box on

%% SUPPLEMENTARY FIGURE 9
% bioRxiv Supplementary Figure 12

surface_to_plot = surface_midthickness;

medial_wall = find(~cortex);
mode_list = [1:5, 10,25,50,100,200];
num_modes_to_plot = length(mode_list);

factor_x_small = 0.6;
factor_x_big = 1.25;
factor_y = 1.1;

fig = figure('Position', [200 200 600 1000]);

init_x = 0.11;
init_y = 0.505;
length_x = (1-init_x)/(factor_x_small+factor_x_big*(4-1) + 1);
length_y = (0.95-init_y)/(factor_y*(num_modes_to_plot-1) + 1);

for basis_ind=1:4
    if basis_ind==8
        data_to_plot = basis_PCA_resting;
        basis_name = {'PCA'; '(rest)'};
    else
        contrast = representative_contrasts{basis_ind};
        task = lower(tasks{basis_ind});
        data_to_plot = basis_PCA_task.(contrast);
        basis_name = {'PCA'; sprintf('(%s)', task)};
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
            annotation(fig, 'textbox', [ax1.Position(1)-0.04, ax1.Position(2)+ax1.Position(4)*0.4, 0.01, 0.01], 'string', {'component'; sprintf('%i',mode)}, 'edgecolor', 'none', ...
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

init_y = 0.005;

for basis_ind=5:8
    if basis_ind==8
        data_to_plot = basis_PCA_resting;
        basis_name = {'PCA'; '(rest)'};
    else
        contrast = representative_contrasts{basis_ind};
        task = lower(tasks{basis_ind});
        data_to_plot = basis_PCA_task.(contrast);
        basis_name = {'PCA'; sprintf('(%s)', task)};
    end

    for mode_ind=1:num_modes_to_plot
        mode = mode_list(mode_ind);
        
        data_to_plot(medial_wall,mode) = min(data_to_plot(:,mode))*1.1;
        clims = [min(data_to_plot(:,mode)), max(data_to_plot(:,mode))];
        if clims(2)<=0
            clims(2) = 0.01;
        end

        ax1 = axes('Position', [init_x+factor_x_big*length_x*(basis_ind-5) init_y+factor_y*length_y*(-[mode_ind-num_modes_to_plot]) length_x length_y]);
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
        
        if basis_ind==5
            annotation(fig, 'textbox', [ax1.Position(1)-0.04, ax1.Position(2)+ax1.Position(4)*0.4, 0.01, 0.01], 'string', {'component'; sprintf('%i',mode)}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
        end
        
        ax2 = axes('Position', [init_x+factor_x_small*length_x+factor_x_big*length_x*(basis_ind-5) init_y+factor_y*length_y*(-[mode_ind-num_modes_to_plot]) length_x length_y]);
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

%% SUPPLEMENTARY FIGURE 10
% bioRxiv Supplementary Figure 18

circle_locs_x = [0.2, 0.8];
circle_locs_y = 0.7*ones(1,2);
curve_line_x = linspace(circle_locs_x(1)-0.02, circle_locs_x(2)+0.02, 40);
curve_line_y = linspace(circle_locs_y(1), circle_locs_y(1)-0.5, length(curve_line_x)/2);
curve_line_y = cat(2, curve_line_y, linspace(curve_line_y(end), circle_locs_y(2), length(curve_line_x)/2));
rng(1)
curve_line_y = curve_line_y + [0, 0, 0, 0.05*randn(1,length(curve_line_y)-6), 0, 0, 0];
p = polyfit(curve_line_x, curve_line_y, 5);
v = polyval(p, curve_line_x);

fig = figure('Position', [200 200 560 200]);
ax = axes('Position', [0.01 0.01 0.98 0.98]);
hold on
plot([0,1], circle_locs_y, 'k-', 'linewidth', 2);
plot(curve_line_x, v, 'r-', 'linewidth', 2);
plot(circle_locs_x(1), circle_locs_y(1), 'ko', 'markersize', 20, 'markerfacecolor', 'b')
plot(circle_locs_x(2), circle_locs_y(2), 'ko', 'markersize', 20, 'markerfacecolor', 'b')
hold off;
set(gca, 'xtick', [], 'ytick', [], 'xlim', [0 1], 'ylim', [0 1])

% texts
text(0, circle_locs_y(1)-0.04, 'cortex', 'horizontalalignment', 'left', 'fontsize', fontsize_label)
text(circle_locs_x(1), circle_locs_y(1)-0.12, '$\mathbf{r''}$', 'horizontalalignment', 'center', 'fontsize', fontsize_label+2, 'interpreter', 'latex')
text(circle_locs_x(2), circle_locs_y(2)-0.12, '$\mathbf{r}$', 'horizontalalignment', 'center', 'fontsize', fontsize_label+2, 'interpreter', 'latex')
text(circle_locs_x(1), circle_locs_y(1)+0.14, '$Q(\mathbf{r''},t'')$', 'horizontalalignment', 'center', 'fontsize', fontsize_label+2, 'interpreter', 'latex')
text(circle_locs_x(2), circle_locs_y(2)+0.14, '$\phi(\mathbf{r},t)$', 'horizontalalignment', 'center', 'fontsize', fontsize_label+2, 'interpreter', 'latex')
text(mean(circle_locs_x), min(v)-0.07, '$W(\mathbf{r},t;\mathbf{r''},t'') := W(\mathbf{r}-\mathbf{r''},t-t'')$', 'horizontalalignment', 'center', 'fontsize', fontsize_label+2, 'interpreter', 'latex')
axis off

%% SUPPLEMENTARY FIGURE 11
% bioRxiv Supplementary Figure 21

parc_name = 'Glasser360';
parc = dlmread(filename_common_parcellation(parc_name, hemisphere));
num_parcels = length(unique(parc(parc>0)));

fig = figure;
[~, max_ind] = max(model_wave_visual.simulated_neural_visual_parcel, [], 2);
data_to_plot_x = tiedrank(data_myelin.myelin(1:num_parcels));
data_to_plot_y = tiedrank(model_wave_visual.tspan(max_ind)');
hold on;
plot(data_to_plot_x, data_to_plot_y, 'k.', 'markersize', 25)
plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
    'r-', 'linewidth', 2);
hold off;
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [min(data_to_plot_x) max(data_to_plot_x)], 'ylim', [min(data_to_plot_y) max(data_to_plot_y)])
xlabel('T1w:T2w (rank)', 'fontsize', fontsize_label)
ylabel('time to peak (rank)', 'fontsize', fontsize_label)

stat = 'linear';
if strcmpi(stat, 'rank')
    [rho,~] = corr(data_to_plot_x, data_to_plot_y, 'type', 'spearman');
elseif strcmpi(stat, 'linear')
    [rho,~] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
end
pval = perm_sphere_p(data_to_plot_x, data_to_plot_y, perm_id, 'pearson');
text(min(get(gca,'xlim'))+1*[max(get(gca,'xlim'))-min(get(gca,'xlim'))], min(get(gca,'ylim'))+0.95*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], ['r = ', sprintf('%.2f', rho)], 'color', 'r', ...
    'fontsize', fontsize_label, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'right');
text(min(get(gca,'xlim'))+1*[max(get(gca,'xlim'))-min(get(gca,'xlim'))], min(get(gca,'ylim'))+0.87*[max(get(gca,'ylim'))-min(get(gca,'ylim'))], extract_pvalue_text(pval,1,'spin'), 'color', 'r', ...
    'fontsize', fontsize_label, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'right');
box off

%% SUPPLEMENTARY VIDEO 1

surface_to_plot = surface_midthickness;
data_to_plot = model_wave_visual.simulated_neural_visual_vertex;
t = model_wave_visual.tspan;
t_interest = 1:0.1:90;
is_time_ms = 1;
medial_wall = find(~cortex);
with_medial = 1;
cmap = parula;
show_colorbar = 1;
output_filename = 'SuppVideo1';
save_video = 0;

fig = video_surface_activity(surface_to_plot, data_to_plot, hemisphere, t, ...
                             t_interest, is_time_ms, medial_wall, with_medial, ...
                             cmap, show_colorbar, output_filename, save_video);
                         