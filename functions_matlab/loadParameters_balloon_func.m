function param = loadParameters_balloon_func
%% loadParameters_balloon_func.m     
%
% Contains all the parameters of the model and the 
% computational parameters for the simulations/calculations. It is
% necessary to make an instance of the class before you can use the 
% parameters.
%
% Example:
% >> param = loadParameters_balloon_func;
% 
%
% Some important notes:
% 1. All the parameters can be changed by overwriting the existing instance 
%   (example: param.kappa = 0.5 if you want to change the signal decay rate).
%
% Original: James Pang, Monash University, 2022
%
%%
    % =====================================================================
    %               DEFAULT INDEPENDENT MODEL PARAMETERS
    % ===================================================================== 
        
    param.kappa    = 0.65;       % signal decay rate [s^-1]
    param.gamma    = 0.41;       % rate of elimination [s^-1]
    param.tau      = 0.98;       % hemodynamic transit time [s]
    param.alpha    = 0.32;       % Grubb's exponent [unitless]
    param.rho      = 0.34;       % resting oxygen extraction fraction [unitless]
    param.V0       = 0.02;       % resting blood volume fraction [unitless]

%     % scanner parameters
%     param.k1       = 3.72;       % 
%     param.k2       = 0.53;       % 
%     param.k3       = 0.53;       %

    param.w_f      = 0.56;
    param.Q0       = 1;
    param.rho_f    = 1000;
    param.eta      = 0.3;
    param.Xi_0     = 1;
    param.beta     = 3;
    param.V_0      = 0.02;
%     param.k1       = 7*param.rho;
%     param.k2       = 2;
%     param.k3       = 2*param.rho - 0.2;
    param.k1       = 3.72;
    param.k2       = 0.527;
    param.k3       = 0.48;
    param.beta     = (param.rho + (1 - param.rho)*log(1 - param.rho))/param.rho;

    % =====================================================================
    %                    COMPUTATIONAL PARAMETERS
    % ===================================================================== 

    param.tstep    = 0.01;       % time step
    param.tmax     = 100;        % maximum time
        
    % =====================================================================
    %                     DEPENDENT PARAMETERS
    % ===================================================================== 
    
    param.tspan    = [0, param.tmax];             % time period limits
    param.T        = 0:param.tstep:param.tmax;    % time vector
        
end
