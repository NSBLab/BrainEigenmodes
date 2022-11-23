function [eig_vec, eig_val] = calc_network_eigenmode(network, num_modes)
% calc_network_eigenmode.m
%
% Calculate the eigenvectors (eigenmodes) and eigenvalues of the input 
% network
%
% Inputs: network   : symmetric connectivity matrix [N x N]
%                     N = number of row/columns
%         num_modes : number of modes to return (int)
%
% Outputs: eig_vec  : eigenvectors (eigenmodes) [N x num_modes]
%          eig_val  : eigenvalues [num_modes x 1]
%
% Original: James Pang, Monash University, 2022

%%
if nargin<2
    num_modes = size(network,1);
end

% Remove self-connections
network(1:(1+size(network,1)):end) = 0;

% Calculate normalized Laplacian matrix
[~, Lnorm] = calc_LaplacianMatrix(network);
Lnorm(isnan(Lnorm)) = 0;

% Calculate the eigenvectors and eigenvalues
[eig_vec, eig_val] = eig(Lnorm);

% Extract first num_modes modes
eig_vec = eig_vec(:, 1:num_modes);
eig_val = diag(eig_val(1:num_modes, 1:num_modes));