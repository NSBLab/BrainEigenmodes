function [peak_lags, peak_corr, zero_lag_corr, lag_proj_unweighted, lag_proj_weighted] = calc_lag_projection(data, TR, lag_lim, max_TR_shift, is_data_normalized)

% data = [T x V x num_subjects]
% Adopted from:

if nargin<5
    is_data_normalized = 0;
end

num_timepoints = size(data, 1);
num_points = size(data, 2);
num_subjects = size(data, 3);

% Setup
% max_TR_shift = round(lag_lim/TR);
lags = -max_TR_shift:max_TR_shift;

if num_subjects==1
    signal = data;
    if ~is_data_normalized
        signal = calc_normalize_timeseries(signal);
    end

    [peak_lags, peak_corr, zero_lag_corr] = single_subject_analysis(signal, TR, lags, lag_lim);
    zero_lag_corr = tanh(zero_lag_corr); % inverse Fisher z transform
else
    % Initialize matrices
    peak_lags = zeros(num_points);
    peak_corr = zeros(num_points);
    zero_lag_corr = zeros(num_points);
    
    peak_lags_nan = zeros(num_points);
    peak_corr_nan = zeros(num_points);
    zero_lag_corr_nan = zeros(num_points);
    
    for sub_ind = 1:num_subjects
        signal = data(:,:,sub_ind);
        if ~is_data_normalized
            signal = calc_normalize_timeseries(signal);
        end
        
        [peak_lags_temp, peak_corr_temp, zero_lag_corr_temp] = single_subject_analysis(signal, TR, lags, lag_lim);
        
        peak_lags = nansum(cat(3, peak_lags, peak_lags_temp), 3);
        peak_corr = nansum(cat(3, peak_corr, peak_corr_temp), 3);
        zero_lag_corr = nansum(cat(3, zero_lag_corr, zero_lag_corr_temp), 3);
        
        peak_lags_nan = peak_lags_nan + isnan(peak_lags_temp);
        peak_corr_nan = peak_corr_nan + isnan(peak_corr_temp);
        zero_lag_corr_nan = zero_lag_corr_nan + isnan(zero_lag_corr_temp);
    end
    
    % Calculate group averages
    peak_lags = peak_lags./(num_subjects - peak_lags_nan);
    peak_corr = peak_corr./(num_subjects - peak_corr_nan);
    zero_lag_corr = zero_lag_corr./(num_subjects - zero_lag_corr_nan);
    zero_lag_corr = tanh(zero_lag_corr); % inverse Fisher z transform 
end

% Calculate unweighted lag projection
lag_proj_unweighted = nanmean(peak_lags);

% Calculate Weighted lag projection 
% (inversely weight lags by correlation magnitude to reduce sampling error)
lag_weights = tan((pi/2)*(1-abs(zero_lag_corr))).^(-2);    % weighted by 1/f^2(r); f(r) = tan[(pi/2)(1-|r|)]
lag_weights(logical(eye(size(lag_weights)))) = 0;          % zero-out diagonal weights
lag_weights(isnan(peak_lags)) = nan;
lag_proj_weighted = peak_lags.*lag_weights;
lag_proj_weighted = nansum(lag_proj_weighted)./nansum(lag_weights);

end

function [peak_lags, peak_corr, zero_lag_corr] = single_subject_analysis(signal, TR, lags, lag_lim)

num_timepoints = size(signal, 1);

% Calculate covariance matrix
Cov = lagged_cov(signal, signal, max(lags));
Cov = Cov/num_timepoints;

% Calculate peak lag and correlation matrices using parabolic interpolation
[peak_lags, peak_corr] = parabolic_interp(Cov, TR);
peak_lags(abs(peak_lags) > lag_lim) = nan;  % Exclude long lags (generally occur when CCF is flat)

% Calculate zero-lag correlation
temp = Cov(:,:,lags==0);    % zero-lag correlation
d = zeros(size(temp));
d(logical(eye(length(temp)))) = sqrt(diag(temp));
temp = d^(-1)*temp/d;
temp = atanh(temp);         % Fisher z transform
temp(isnan(peak_lags)) = nan;
zero_lag_corr = temp;

end

function r = lagged_cov(Avg1,Avg2,L)
% Computes (unnormalized) cross-covariance function out to +/- L lags 
% (TR shifts) between each column of Avg1 and Avg2

% Avg1 and Avg2 are time x region matrices/vectors

% See tdmx_template.m for appropriate normalization of lagged_cov output

L1 = size(Avg1,2);
L2 = size(Avg2,2);
r = single(zeros(L1,L2,2*L+1));

k = 1;
for i = -L:L
    tau = abs(i);

    if i >=0
        Avg1_lagged = Avg1(1:end-tau,:);
        Avg2_lagged = Avg2(1+tau:end,:);
    else
        Avg1_lagged = Avg1(1+tau:end,:);
        Avg2_lagged = Avg2(1:end-tau,:);
    end    

    r(:,:,k) = Avg1_lagged'*Avg2_lagged;
    k = k+1;
end

end

function [peak_lag,peak_cov] = parabolic_interp(lcc,tr)
%This function uses parabolic interpolation to find the lag using the
%extremum of the lagged cross correlation.

%The function first uses the
%sign of the cross correlation at zero to decide whether to find a maximum
%or minmum. Next, we look for the global max/min.

%lcc is the empirical lagged covariance curve, lags is a vector with the timepoints
%in each temporal direction (e.g. -8:2:8 for +/- 8 seconds with a 2 second TR). 

s = size(lcc);
peak_lag = nan([1,s(1)*s(2)]);
peak_cov = peak_lag;

% linearize
lcc = reshape(lcc,[s(1)*s(2),s(3)])';

% find index of extremum (max or min determined by sign at zero-lag)
[~,I]= max(bsxfun(@times,lcc,sign(lcc((s(3)+1)/2,:))),[],1);

% ensure extremum is not at an endpoint (this would preclude parabolic interpolation)
use = I>1 & I<s(3);
lcc = lcc(:,use);

% place peaks at center
x0 = I(use) - (s(3)+1)/2;

% set up three-point ccf for interpolation (y1,y2,y3)
i = sub2ind([size(lcc),sum(use)],I(use),1:sum(use));
lcc = [lcc(i-1);lcc(i);lcc(i+1)];

% fit parabola: tau = TR * (y1-y3) / (2*(y1-2y2+y3))
b = (lcc(3,:) - lcc(1,:))/2;
a = (lcc(1,:) + lcc(3,:) - 2*lcc(2,:))/2;
peak_lag(use) =  (-b./(2*a));

% construct parabola to get covariance (y = ax^2 + bx + c)
peak_cov(use) = a.*(peak_lag(use).^2) + b.*peak_lag(use) + lcc(2,:);

% put back TR information
peak_lag(use) = (peak_lag(use) + x0)*tr;

peak_lag = reshape(peak_lag,[s(1) s(2)]);
peak_cov = reshape(peak_cov,[s(1) s(2)]);
	
end


