function triu_ind = calc_triu_ind(matrix)
% calc_triu_ind.m
%
% Calculate upper triangle indices of matrix
%
% Input: matric : data matrix (array)
%
% Output: triu_ind : indices of upper triangle (vector)
%
% Original: James Pang, Monash University, 2022

%%

[M, N] = size(matrix);

triu_ind = find(triu(ones(M,N),1));

