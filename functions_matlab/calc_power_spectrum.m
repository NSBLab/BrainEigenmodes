function [power_spectrum, power_spectrum_norm] = calc_power_spectrum(data)
% calc_power_spectrum.m
%
% Calculate spatial power spectrum from eigenmode analysis
%
% Inputs: data                 : data [NxP]
%                                N = number of modes, P = number of independent data
%
% Outputs: power_spectrum      : power spectrum [NxP]
%          power_spectrum_norm : normalized power spectrum [NxP]
%
% Original: James Pang, Monash University, 2022

%%

[N,P] = size(data);

power_spectrum = abs(data).^2;
power_spectrum_norm = power_spectrum./repmat(nansum(power_spectrum,1), N, 1);


end