function sim_activity = model_neural_waves(eig_vec, eig_val, ext_input, param, method)
% model_neural_waves.m
%
% Simulate neural waves on a surface using eigenmodes.
%
% Inputs: eig_vec      : eigenvectors (eigenmodes) [V x num_modes]
%                        V = number of vertices
%                        num_modes = number of modes
%         eig_val      : eigenvalues [num_modes x 1]
%         ext_input    : spatiotemporal external input [V x T]
%                        T = number of time points
%         param        : model parameters (struct)
%                        Create instance using loadParameters_wave_func.m
%         method       : method for calculating the activity (string)
%                        ODE = via solving ODEs
%                        Fourier = via solving Fourier transform
%
% Output: sim_activity : simulated wave activity [V x T]
%
% Original: James Pang, Monash University, 2022

%%

if nargin<5
    method = 'ODE';
end

% if time is in ms, param.gamma_s must be in ms^-1
% hence uncomment below to rescale to ms^-1
% param.gamma_s = param.gamma_s*1e-3; 

num_modes = size(eig_vec,2);

switch method
    case 'ODE'
    
        % mode decomposition of external input
        ext_input_coeffs = calc_eigendecomposition(ext_input, eig_vec, 'matrix');

        % initialize simulated activity vector
        sim_activity = zeros(num_modes, size(ext_input_coeffs,2));

        for mode_ind = 1:num_modes
            mode_coeff = ext_input_coeffs(mode_ind,:);
            lambda = eig_val(mode_ind);

            [tout, yout] = ode45(@(t,y) wave_ODE(t, y, mode_coeff, lambda, param), param.T, [mode_coeff(1); 0]);

            sim_activity(mode_ind,:) = yout(:,1);
        end
        
    case 'Fourier'
        
        % append time vector with negative values to have a zero center
        T_append = [-param.tmax:param.tstep:param.tmax];
        Nt = length(T_append);
        t0_ind = dsearchn(T_append', 0);
        
        % mode decomposition of external input
        ext_input_coeffs_temp = calc_eigendecomposition(ext_input, eig_vec, 'matrix');

        % append external input coefficients for negative time values
        ext_input_coeffs = zeros(num_modes, Nt);
        ext_input_coeffs(:,t0_ind:end) = ext_input_coeffs_temp;
        
        % initialize simulated activity vector
        sim_activity = zeros(num_modes, size(ext_input_coeffs,2));

        for mode_ind = 1:num_modes
            mode_coeff = ext_input_coeffs(mode_ind,:);
            lambda = eig_val(mode_ind);

            yout = wave_Fourier(mode_coeff, lambda, T_append, param);

            sim_activity(mode_ind,:) = yout;
        end
        
        sim_activity = sim_activity(:,t0_ind:end);
end

% combine mode time series with mode spatial map
sim_activity = eig_vec*sim_activity;

end

function out = wave_ODE(t, y, mode_coeff, lambda, param)
% wave_ODE.m
%
% Calculate the temporal activity of one mode via solving an ODE.
%
% Inputs: t          : time variable 
%         y          : activity variable
%         mode_coeff : coefficient of the mode [1 x T]
%         lambda     : eigenvalue of the mode (float)
%         param      : model parameters (struct)
%
% Output: out        : activity and its first-order derivative [2 x 1]
%
% Original: James Pang, Monash University, 2022


out = zeros(2,1);

% t_ind = dsearchn(param.T', t);
coef_interp = interp1(param.T, mode_coeff, t); % interpolated coefficient at t

out(1) = y(2);
out(2) = param.gamma_s^2*(coef_interp - (2/param.gamma_s)*y(2) - y(1)*(1 + param.r_s^2*lambda));

end

function out = wave_Fourier(mode_coeff, lambda, T, param)
% wave_Fourier.m
%
% Calculate the temporal activity of one mode via Fourier transform.
%
% Inputs: mode_coeff : coefficient of the mode [1 x T]
%         lambda     : eigenvalue of the mode (float)
%         T          : time vector with zero center [1 x T]
%         param      : model parameters (struct)
%
% Output: out        : activity [1 x T]
%
% Original: James Pang, Monash University, 2022

Nt = length(T);
Nw = Nt;
wsamp = 1/mean(param.tstep)*2*pi;
wMat = (-1).^(1:Nw);
jvec = 0:Nw-1;
w = (wsamp)*1/Nw*(jvec - Nw/2);

% mode_coeff_fft = ctfft(mode_coeff, param.T);	
mode_coeff_fft = coord2freq_1D(mode_coeff, w);	

out_fft = param.gamma_s^2*mode_coeff_fft./(-w.^2 - 2*1i*w*param.gamma_s + param.gamma_s^2*(1 + param.r_s^2*lambda));

% calculate inverse Fourier transform
out = real(freq2coord_1D(out_fft, w));

end



