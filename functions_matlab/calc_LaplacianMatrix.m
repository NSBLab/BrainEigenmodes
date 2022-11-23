function [L, Lnorm] = calc_LaplacianMatrix(W)
% calc_LaplacianMatrix.m
%
% Calculate Laplacian and normalized Laplacian
%
% Input: W       : connectivity matrix [NxN]
%
% Outputs: L     : Laplacian matrix [NxN]
%          Lnorm : normalized Laplacian matrix [NxN]
%
% Original: James Pang, Monash University, 2021

%%

d = sum(W,2);   % degree or strength of each node
D = diag(d);

L = D - W;

Dhalf = diag(1./sqrt(d));
Lnorm = Dhalf*(L*Dhalf);

