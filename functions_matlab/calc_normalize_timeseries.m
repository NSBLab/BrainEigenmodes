function data_normalized = calc_normalize_timeseries(data)
% calc_normalize_timeseries.m
%
% Normalize timeseries with respect to mean and std
%
% Input: data             : fMRI data [TxN]
%                           T = length of time
%                           N = number of points
%
% Output: data_normalized : normalized fMRI data [TxN]
%
% Original: James Pang, Monash University, 2021

%%

T = size(data,1);
data_normalized = detrend(data, 'constant');
data_normalized = data_normalized./repmat(std(data_normalized),T,1);