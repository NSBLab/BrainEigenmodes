function [mode_activity, sim_activity] = model_BOLD_balloon(eig_vec, ext_input, param, method)
% model_BOLD_balloon.m
%
% Simulate BOLD balloon model on a surface using eigenmodes.
%
% Inputs: eig_vec       : eigenvectors (eigenmodes) [V x num_modes]
%                         V = number of vertices
%                         num_modes = number of modes
%         ext_input     : spatiotemporal external input [V x T]
%                         T = number of time points
%         param         : model parameters (struct)
%                         Create instance using loadParameters_balloon_func.m
%         method        : method for calculating the activity (string)
%                         ODE = via solving ODEs
%                         Fourier = via solving Fourier transform
%
% Outputs: mode_activity : simulated mode activity [num_modes x T]
%          sim_activity  : simulated BOLD activity [V x T]
%          
% Original: James Pang, Monash University, 2022

%%

if nargin<4
    method = 'ODE';
end

num_modes = size(eig_vec,2);

switch method
    case 'ODE'
    
        % mode decomposition of external input
        ext_input_coeffs = calc_eigendecomposition(ext_input, eig_vec, 'matrix');

%         % initialize simulated activity vector
%         sim_activity = zeros(num_modes, size(ext_input_coeffs,2));
% 
%         for mode_ind = 1:num_modes
%             mode_coeff = ext_input_coeffs(mode_ind,:);
% 
%             [tout, yout] = ode45(@(t,y) balloon_ODE(t, y, mode_coeff, param), param.T, mode_coeff(1)*ones(5,1));
% 
%             sim_activity(mode_ind,:) = yout(:,5);
%         end

        
        % initialize activity vectors
        sol.z = zeros(num_modes, length(param.T));
        sol.f = zeros(num_modes, length(param.T));
        sol.v = zeros(num_modes, length(param.T));
        sol.q = zeros(num_modes, length(param.T));
        sol.BOLD = zeros(num_modes, length(param.T));

        % setting initial condition
%         F0 = repmat(ext_input_coeffs(:,1), 1, 4);
        F0 = repmat(0.001*ones(num_modes,1), 1, 4);
        
        F = F0;
        sol.z(:,1) = F(:,1);
        sol.f(:,1) = F(:,2);
        sol.v(:,1) = F(:,3);
        sol.q(:,1) = F(:,4);
        
        for k = 2:length(param.T)
            dF = balloon_ODE(ext_input_coeffs(:,k-1), F, param);
            F = F + dF*param.tstep;
            sol.z(:,k) = F(:,1);
            sol.f(:,k) = F(:,2);
            sol.v(:,k) = F(:,3);
            sol.q(:,k) = F(:,4);
        end

        sol.BOLD = 100*param.V0*(param.k1*(1 - sol.q) + param.k2*(1 - sol.q./sol.v) + param.k3*(1 - sol.v));

        sim_activity = sol.BOLD;


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

            yout = balloon_Fourier(mode_coeff, T_append, param);

            sim_activity(mode_ind,:) = yout;
        end
        
        sim_activity = sim_activity(:,t0_ind:end);
end

mode_activity = sim_activity;

% combine mode time series with mode spatial map
sim_activity = eig_vec*sim_activity;


end

% function out = balloon_ODE(t, y, mode_coeff, param)
% % balloon_ODE.m
% %
% % Calculate the temporal activity of one mode via solving an ODE.
% %
% % Inputs: t          : time variable 
% %         y          : activity variable
% %         mode_coeff : coefficient of the mode [1 x T]
% %         param      : model parameters (struct)
% %
% % Output: out        : activity variables [5 x 1]
% %
% % Original: James Pang, Monash University, 2022
% 
% 
% out = zeros(5,1);
% 
% % t_ind = dsearchn(param.T', t);
% coef_interp = interp1(param.T, mode_coeff, t); % interpolated coefficient at t
% 
% S = coef_interp;
% z = y(1);
% f = y(2);
% v = y(3);
% q = y(4);
% Y = y(5);
% 
% out(1) = S - param.kappa*z - param.gamma*(f - 1);
% out(2) = z;
% out(3) = (1/param.tau)*(f - v.^(1/param.alpha));
% out(4) = (1/param.tau)*((f/param.rho).*(1 - (1 - param.rho).^(1./f)) - q.*v.^(1/param.alpha - 1));
% out(5) = 100*param.V0*(param.k1*(1 - q) + param.k2*(1 - q./v) + param.k3*(1 -  v));
% 
% end

function dF = balloon_ODE(S, F, param)
% balloon_ODE.m
%
% Calculate the temporal activity by solving an ODE.
%
% Inputs: S          : spatiotemporal external input [V x 1]
%         F          : solutions at one time point
%         param      : model parameters (struct)
%
% Output: dF         : time derivative of variables [4 x 1]
%
% Original: James Pang, Monash University, 2022

z = F(:,1);
f = F(:,2);
v = F(:,3);
q = F(:,4);

dF(:,1) = S - param.kappa*z - param.gamma*(f - 1);
dF(:,2) = z;
dF(:,3) = (1/param.tau)*(f - v.^(1/param.alpha));
dF(:,4) = (1/param.tau)*((f/param.rho).*(1 - (1 - param.rho).^(1./f)) - q.*v.^(1/param.alpha - 1));

end

function out = balloon_Fourier(mode_coeff, T, param)
% balloon_Fourier.m
%
% Calculate the temporal activity of one mode via Fourier transform.
%
% Inputs: mode_coeff : coefficient of the mode [1 x T]
%         T          : time vector with zero center [1 x T]
%         param      : model parameters (struct)
%
% Output: out        : activity [1 x T]
%
% Original: James Pang, Monash University, 2022

Nt = length(T);
Nw = Nt;
wsamp = 1/mean(param.tstep)*2*pi;
jvec = 0:Nw-1;
w = (wsamp)*1/Nw*(jvec - Nw/2);

% mode_coeff_fft = ctfft(mode_coeff, param.T);	
mode_coeff_fft = coord2freq_1D(mode_coeff, w);
	
% Transfer functions from Robinson et al. 2006, Aquino et al. 2012, 2014,
% Pang et al. 2016, 2018
T_Fz = 1 ./ (-(w + 1i*0.5*param.kappa).^2 + param.w_f^2);
T_yF = param.V_0 * (param.alpha*(param.k2 + param.k3)*(1 - 1i*param.tau*w) - (param.k1 + param.k2)*(param.alpha + param.beta - 1 - 1i*param.tau*param.alpha*param.beta*w))./((1 - 1i*param.tau*w).*(1 - 1i*param.tau*param.alpha*w));
T_yz = T_yF.*T_Fz;

out_fft = T_yz.*mode_coeff_fft;

% calculate inverse Fourier transform
out = real(freq2coord_1D(out_fft, w));

end



