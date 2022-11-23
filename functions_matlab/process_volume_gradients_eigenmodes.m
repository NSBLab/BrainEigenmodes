%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% process_volume_gradients_eigenmodes.m
%%%
%%% MATLAB script to process volume gradients and eigenmodes from HCP data.
%%% In particular, the script demonstrates how to 
%%% (1) calculate the z-scores of the gradients,
%%% (2) calculate the z-scores of the eigenmodes and removing the first
%%%     (constant) eigenmode, and
%%% (3) combine the gradients and eigenmodes into a single file and
%%%     calculate correlations between them
%%%
%%% Original: James Pang, Monash University, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load preliminary variables

hemispheres = {'lh', 'rh'};
structures = {'tha', 'striatum', 'hippo'};

data_dir = 'data';

%% Calculate z-scores of the gradients

for structure_ind = 1:length(structures)
    structure = structures{structure_ind};
    
    for hemisphere_ind = 1:length(hemispheres)
        hemisphere = hemispheres{hemisphere_ind};
        
        input_filename = sprintf('%s/results/subcortical/hcp_%s-%s_gradient_20.nii.gz', data_dir, structure, hemisphere);
        output_filename = sprintf('%s/results/subcortical/hcp_%s-%s_gradient_20_zscore.nii', data_dir, structure, hemisphere);
        mask_filename = sprintf('%s/template_surfaces_volumes/hcp_%s-%s_thr25.nii.gz', data_dir, structure, hemisphere);

        V = niftiread(input_filename);
        V_info = niftiinfo(input_filename);
        
        V_mask = niftiread(mask_filename);
        ind = find(V_mask~=0);
        
        num_gradients = size(V,4);

        V_zscore = zeros(size(V));
        for gradient_ind = 1:num_gradients
            temp = V(:,:,:,gradient_ind);
            temp2 = NaN(size(temp));
            temp2(ind) = zscore(temp(ind));
            
            V_zscore(:,:,:,gradient_ind) = temp2;
        end
        
        V_zscore_info = V_info;
        V_zscore_info.Filename = [output_filename, '.gz'];
        niftiwrite(single(V_zscore), output_filename, V_zscore_info, 'Compressed', true);
    end
end

%% Calculate z-scores of the eigenmodes and removing the first one 
% we need to remove the first one as it represents a constant mode and
% to match with the functional gradients, which do not have a constant
% gradient

for structure_ind = 1:length(structures)
    structure = structures{structure_ind};
    
    for hemisphere_ind = 1:length(hemispheres)
        hemisphere = hemispheres{hemisphere_ind};
        
        input_filename = sprintf('%s/template_eigenmodes/hcp_%s-%s_emode_31.nii.gz', data_dir, structure, hemisphere);
        output_filename = sprintf('%s/results/subcortical/hcp_%s-%s_emode_30_noconstant_zscore.nii', data_dir, structure, hemisphere);
        mask_filename = sprintf('%s/template_surfaces_volumes/hcp_%s-%s_thr25.nii.gz', data_dir, structure, hemisphere);

        V = niftiread(input_filename);
        V_info = niftiinfo(input_filename);
        
        V_mask = niftiread(mask_filename);
        ind = find(V_mask~=0);
        
        num_modes = size(V,4);

        V_zscore = zeros(size(V));
        for mode_ind = 1:num_modes
            temp = V(:,:,:,mode_ind);
            temp2 = NaN(size(temp));
            temp2(ind) = zscore(temp(ind));
            
            V_zscore(:,:,:,mode_ind) = temp2;
        end
        
        V_noconstant_zscore = V_zscore(:,:,:,2:end); 
        V_noconstant_zscore_info = V_info;
        V_noconstant_zscore_info.Filename = [output_filename, '.gz'];
        V_noconstant_zscore_info.ImageSize = [V_info.ImageSize(1), V_info.ImageSize(2), V_info.ImageSize(3), V_info.ImageSize(4)-1];
        niftiwrite(single(V_noconstant_zscore), output_filename, V_noconstant_zscore_info, 'Compressed', true);
    end
end

%% Combine gradients and eigenmodes into a single file and calculate correlations

gradients = struct();
eigenmodes = struct();
correlation = struct();

for structure_ind = 1:length(structures)
    structure = structures{structure_ind};
    
    for hemisphere_ind = 1:length(hemispheres)
        hemisphere = hemispheres{hemisphere_ind};
        
        gradients_filename = sprintf('%s/results/subcortical/hcp_%s-%s_gradient_20_zscore.nii.gz', data_dir, structure, hemisphere);
        eigenmodes_filename = sprintf('%s/results/subcortical/hcp_%s-%s_emode_30_noconstant_zscore.nii.gz', data_dir, structure, hemisphere);
        mask_filename = sprintf('%s/template_surfaces_volumes/hcp_%s-%s_thr25.nii.gz', data_dir, structure, hemisphere);

        V_gradients = niftiread(gradients_filename);
        V_eigenmodes = niftiread(eigenmodes_filename);
        V_mask = niftiread(mask_filename);
        ind = find(V_mask~=0);
        
        num_modes = size(V_gradients,4);
        
        for mode_ind = 1:num_modes
            temp1 = V_gradients(:,:,:,mode_ind);
            temp2 = V_eigenmodes(:,:,:,mode_ind);

            gradients.(structure).(hemisphere)(:,mode_ind) = temp1(ind);
            eigenmodes.(structure).(hemisphere)(:,mode_ind) = temp2(ind);
        end
        
        correlation.(structure).(hemisphere) = corr(eigenmodes.(structure).(hemisphere), gradients.(structure).(hemisphere));
    end
end

save(sprintf('%s/results/subcortical/data_gradients_eigenmodes_20.mat', data_dir), 'gradients', 'eigenmodes', 'correlation', '-v7.3')

