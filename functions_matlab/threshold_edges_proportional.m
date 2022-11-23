function W = threshold_edges_proportional(W, p, weight_to_remove)
% threshold_edges_proportional.m
%
% This function "thresholds" the connectivity matrix by preserving a
% proportion p (0<p<1) of the strongest or weakest weights. All other weights,
% and all weights on the main diagonal (self-self connections) are set to 0.
%
% Based on Brain Connectivity Toolbox's threshold_proportional
%
% Inputs: W                : connectivity matrix [NxN]
%         p                : proportion of weights to preserve (float)
%                            p=1 (all weights preserved)
%                            p=0 (no weights preserved)
%         weight_to_remove : strongest or weakest (string)
%
% Output: W                : thresholded matrix [NxN]
%
% Original: James Pang, QIMR Berghofer, 2020

%%

if nargin<3
    weight_to_remove = 'weakest';
end

N = size(W,1);                                %number of nodes
W(1:N+1:end)=0;                             %clear diagonal

if max(max(abs(W-W.'))) < 1e-10             %if symmetric matrix
    W=triu(W);                              %ensure symmetry is preserved
    ud=2;                                   %halve number of removed links
else
    ud=1;
end

ind=find(W);                                %find all links
E=sortrows([ind W(ind)], -2);               %sort by magnitude
en=round(length(ind)*p);                     %number of links to be preserved

if strcmpi(weight_to_remove, 'weakest')
    W(E(en+1:end,1))=0;                         %apply threshold
elseif strcmpi(weight_to_remove, 'strongest')
    W(E(1:(length(ind)-en),1))=0;
end

if ud==2                                    %if symmetric matrix
    W=W+W.';                                %reconstruct symmetry
end

