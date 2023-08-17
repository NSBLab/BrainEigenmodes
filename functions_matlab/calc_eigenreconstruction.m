function [corrCoeffs, recon, betaCoeffs] = calc_eigenreconstruction(data, eigenvectors, method)
% calc_eigenreconstruction.m
%
% Reconstruct data using eigenvectors and calculate the coefficient of 
% contribution of each vector
%
% Inputs: data         : data [MxP]
%                        M = number of points, P = number of independent data
%         eigenvectors : eigenvectors [MxN]
%                        M = number of points, N = number of eigenvectors
%         method       : type of calculation
%                        'matrix', 'matrix_separate', 'regression'
% Output: corrCoeffs   : correlation coefficient values [NxP]
%         recon        : reconstructions using 1..P modes [MxNxP]
%         betaCoeffs   : coefficient values {Nx1}
%
% Original: Mehul Gajwani, Monash University, 2023

%%

if nargin < 3; method = 'matrix'; end

[M,P] = size(data); N = size(eigenvectors, 2);
assert(size(data,1) == size(eigenvectors, 1));

recon = nan(M,N,P);
corrCoeffs = nan(N,P);
betaCoeffs = cell(N, 1);

for n = 1:N
    betaCoeffs{n} = calc_eigendecomposition(data, eigenvectors(:,1:n), method); %#ok<*AGROW>
    recon(:,n,:) = eigenvectors(:,1:n)*betaCoeffs{n};
    corrCoeffs(n,:) = diag(corr(data, squeeze(recon(:,n,:))));
end

% if the first column of `eigenvectors` is already constant, then the first
% row of `corrCoeffs` may be NaNs
assert(~any(isnan(corrCoeffs(2:end,:)), 'all')); 
assert(~any(isnan(recon), 'all'));

end

