function [coeffs, constOffset] = calc_eigendecomposition(data, eigenvectors, method)
% calc_eigendecomposition.m
%
% Decompose data using eigenvectors and calculate the coefficient of 
% contribution of each vector
%
% Inputs: data         : data [MxP]
%                        M = number of points, P = number of independent data
%         eigenvectors : eigenvectors [MxN]
%                        M = number of points, N = number of eigenvectors
%         method       : type of calculation
%                        'matrix', 'matrix_separate', 'regression'
% Output: coeffs       : coefficient values [NxP]
%         constOffest  : offset needed for reconstruction if `eigenvectors` does not include constant column [1xP]
%
% Original: James Pang, Monash University, 2022
% Revised: Mehul Gajwani, Monash University, 2024

%%

[M,P] = size(data);
% [~,N] = size(eigenvectors);

if nargin<3
    method = 'matrix';
end

%%% ENSURE THAT `eigenvectors` CONTAINS A COLUMN OF CONSTANTS %%%
eigenvectors = [eigenvectors, ones(M, 1)];
% if eigenvectors already had a constant column, then the rank will not
% change (note that N <= M by force)
hasConstantCol = rank(eigenvectors) == size(eigenvectors, 2)-1;
% if it had a constant column, remove the new addition
if hasConstantCol; eigenvectors = eigenvectors(:,1:(end-1)); end
constOffset = zeros(P,1);

switch method
    case 'matrix'
        coeffs = (eigenvectors.'*eigenvectors)\(eigenvectors.'*data);
    case 'matrix_separate'
        coeffs = zeros(size(eigenvectors, 2),P);
        
        for p = 1:P
            coeffs(:,p) = (eigenvectors.'*eigenvectors)\(eigenvectors.'*data(:,p));
        end
    case 'regression'
        coeffs = zeros(size(eigenvectors, 2),P);
        
        for p = 1:P
            coeffs(:,p) = regress(data(:,p), eigenvectors);
        end
end

% if it didn't have a constant column, remove the extra coefficient here
if ~hasConstantCol; constOffset = coeffs(end,:); coeffs = coeffs(1:(end-1),:); end
    
end