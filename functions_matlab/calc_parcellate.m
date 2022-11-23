function data_parcellated = calc_parcellate(parc, data_input)
% calc_parcellate.m
%
% Parcellate data
%
% Inputs: parc          : parcellation [Nx1]
%         data_input    : data to parcellate [NxM]
%
% Output: data_parcellated : parcellated data [num_parcelsxM]
%
% Original: James Pang, Monash University, 2021

%%

num_vertices = size(parc,1);
parcels = unique(parc(parc>0));
num_parcels = length(parcels);

if size(data_input,1) ~= num_vertices
    data_input = data_input';
end

data_parcellated = zeros(num_parcels, size(data_input,2));

for parcel_ind = 1:num_parcels
    parcel_interest = parcels(parcel_ind);

    ind_parcel = find(parc==parcel_interest);
    
    data_parcellated(parcel_ind,:) = nanmean(data_input(ind_parcel,:));
end
