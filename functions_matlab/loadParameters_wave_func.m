function param = loadParameters_wave_func
%% loadParameters_wave_func.m     
%
% Contains all the parameters of the wave model and the 
% computational parameters for the simulations/calculations. It is
% necessary to make an instance of the function before you can use the 
% parameters.
%
% Example:
% >> param = loadParameters_wave_func;
% 
%
% Some important notes:
% 1. All the parameters can be changed by overwriting the existing instance 
%    (example: param.r_s = 10 if you want to change the length scale).
%
% 2. The dependent parameters need to be manually updated when other independent 
%    parameters are changed.
%
% Original: James Pang, Monash University, 2022

%%
    % =====================================================================
    %               DEFAULT INDEPENDENT MODEL PARAMETERS
    % ===================================================================== 
        
    param.r_s        = 30;   % length scale [mm]
    param.gamma_s    = 116;  % damping rate [s^-1]
        
    % =====================================================================
    %                    COMPUTATIONAL PARAMETERS
    % ===================================================================== 
    
    param.is_time_ms = 0;    % 1 = if time is in ms; 0 = if time is in s
    param.tstep      = 0.01; % time step [s]
    param.tmax       = 1;    % maximum time [s]
        
    % =====================================================================
    %                     DEPENDENT PARAMETERS
    % ===================================================================== 
        
    param.tspan    = [0, param.tmax];             % time period limits
    param.T        = 0:param.tstep:param.tmax;    % time vector
    
end
