function out = coord2freq_1D(T, w, padN)
%% coord2freq_1D.m
%
% Calculates the 1D continuous Fourier transform of T.
%
% Inputs: T      : array of 1D transfer function or signal in coordinate space 
%                  size(T) = [1, length(w)]
%         w      : vector of temporal frequencies
%         padN   : number extending the length of T to N by padding zeros
%                  (not a required input)
%
% Output: out    : array of transfer function or signal in frequency space 
% 
% Example:
% >> w = linspace(-1,1,100); 
% >> T = rand(1,length(w));
% >> out = coord2freq_1D(T, w) % gives out  Fourier transform of T
%
% Original: Kevin Aquino, University of Sydney, 2014
%           James Pang, University of Sydney, 2016
% Version 1.2: James Pang, University of Sydney, Jan 2018
% Version 2: James Pang, Monash University, 2022

%%

% calculating the -1 vectors needed for the Fourier transform
wM = (-1).^(1:length(w));

% performing the Fourier transform 
if (nargin > 2)
    out = wM.*ifft(wM.*T, padN, 2);
else
    out = wM.*ifft(wM.*T, [], 2);
end
