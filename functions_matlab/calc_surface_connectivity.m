function surface_connectivity = calc_surface_connectivity(surface)
% calc_surface_connectivity.m
%
% Calculate surface-based connectivity matrix
%
% Input: surface               : surface structure with fields
%                                vertices and faces
%
% Output: surface_connectivity : surface connectivity matrix [NxN]
%
% Original: James Pang, Monash University, 2021

%%

num_vertices = size(surface.vertices, 1);
f = surface.faces;

temp = ...
    sparse( ...
        reshape(f           ,  [], 1), ...
        reshape(f(:, [2 3 1]), [], 1), ...
        true, num_vertices, num_vertices ...
    );

surface_connectivity = full(temp | temp'); % can replace with +full if desired