function [corrCoeffs,recon,betaCoeffs,constOffset,fcCorr,fcRecon] = ...
    calc_eigenreconstruction(data,eigenvectors,method,modesq)
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
%         modesq       : modes to query to get reconstruction [1xN]
% Output: corrCoeffs   : correlation coefficient values [NxP]
%         recon        : reconstructions using 1..N modes [MxNxP]
%         betaCoeffs   : coefficient values {Nx1}
%         constOffset  : offset if eigenvectors do not contain constant column {Nx1} 
%         fcCorr       : treating `data` as spanning P timepoints, reconstruct FC and calculate corr coeff between original and new FC matrices [Nx1] 
%         fcRecon      : reconstructions of FC using 1..N modes [MxMxN] 
%
% Original: Mehul Gajwani, Monash University, 2024

%%

% Prelims
[M,P] = size(data);
assert(size(data,1) == size(eigenvectors, 1));

if nargin < 3; method = 'matrix'; end
if nargin < 4; modesq = 1:size(eigenvectors, 2); end
nq = length(modesq);

recon = nan(M,nq,P);
corrCoeffs = nan(nq,P);
betaCoeffs = cell(nq,1);
constOffset = cell(nq,1);

if nargout > 4
    tril2vec = @(A,k) reshape(A(tril(true(size(A)), k)), [], 1);
    fcRecon = nan(M,M,nq);
    fcCorr = nan(nq,1);
    fcOrig = tril2vec(corr(data.'),-1);
end

% Calculations
for n = 1:nq
    [betaCoeffs{n}, constOffset{n}] = calc_eigendecomposition(data, eigenvectors(:,1:modesq(n)), method);
    recon(:,n,:) = eigenvectors(:,1:modesq(n))*betaCoeffs{n} + constOffset{n};
    corrCoeffs(n,:) = diag(corr(data, squeeze(recon(:,n,:))));
    if nargout > 4
        fcRecon(:,:,n) = corr(squeeze(recon(:,n,:)).');
        fcCorr(n) = corr( fcOrig , tril2vec(corr(fcRecon(:,:,n).'),-1) );
    end
end

% if the first column of `eigenvectors` is already constant, then the first
% row of `corrCoeffs` may be NaNs
assert(~any(isnan(corrCoeffs(2:end,:)), 'all')); 
assert(~any(isnan(recon), 'all'));

end




