function [ft, w] = ctfft(f, t)
% ctfft.m
%
% Calculate continuous fourier transform using FFT
%
% Inputs: f   : function to take fourier transform [1xT]
%               T = length of time
%               NOTE: f must be centered around 0
%         t   : time vector [1xT]
%
% Outputs: ft : Fourier transform of function f [1xT]
%          w  : frequency vector [1xT]
%
% Original: Kevin Aquino, Monash University, 2021
% Edited: James Pang, Monash University, 2022

%%

Nt = length(t);

a = t(end)*2;
beta = a/Nt;
w = (t/beta)*2*pi/(Nt*beta);

ft = ((-1).^[1:Nt]).*beta.*ifft((-1).^[1:Nt].*f,[],2);
return