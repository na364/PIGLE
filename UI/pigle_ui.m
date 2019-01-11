

%% parameters for config_model.m

z_enabled = 0;
dKz_include_in_isf = 0;
theta_enabled = 0;
zero_p_init = 0; % set initial momentum be set to zero? (if set to 0, p_init will correspond to thermal distribution)
interactions_active = 1;
N_runs = 1;
run_parallel = 1;

% Specify dK as a 2D vector, 3rd dim is azimuths.
dK = [0.05 0.1 0.15 0.2:0.1:5];
azim_1 = [1 0];
azim_2 = [cosd(30) sind(30)];

% specify beam parameters and geometrical parameters for scatering calculations
theta_tot = 44.4; % Degrees
beam_ki = 3.3977; % ANgstrom ^{-1} 

% Specify simulation time parameters
% (those will be adjusted by the program, see below if interested)
sample_time = 2e-3;
isf_sample_time = 5e-2;
thermalizing_time = 30;
stop_time = 1024*0.1;

% N_steps and N_ISF_steps are calculated after PIGLE adjusts the requested time parameters
max_N_steps = 1e9;
max_N_ISF_steps = 6e5;

