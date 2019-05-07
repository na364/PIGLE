

%% parameters for config_model.m

z_enabled = 0; % Enable/Disable (1/0) motion prependicular to the surface
dKz_include_in_isf = 0; % Enable/Disable (1/0) calculation of the ISF for prependicular momentum transfer
theta_enabled = 1; % Enable/Disable (1/0) rotation of rigid body (rotation axis prependicular to the surface)
zero_p_init = 0; % set initial momentum be set to zero? (if set to 0, p_init will correspond to thermal distribution)
interactions_active = 1; % Enable/Disable (1/0)
N_runs = 2; % Number of runs
run_parallel = 1; % Enable/Disable (1/0) of parallel computing. If parallel computing toolbox is not installed - set this option to zero.

% Specify momentum transfer (dK) for calculating the ISF at.
% The program will calculate the ISF for the defined dK, on the defined
% azimuths (two are allowed). If interactions are enabled, the dK values at each azimuth
% are adjusted to conform with the 'in-phase' scattering requirement, see
% the Computer Physics Communications paper by Avidor et. al.
dk = [0.05 0.1 0.15 0.2:0.5:2.6]; % magnitude
azim_1 = [1 0];
azim_2 = [cosd(30) sind(30)];

% Specify simulation time parameters
% (those will be adjusted by the program, see below if interested)
sample_time = 1e-3;
isf_sample_time = 5e-2; % time step for which the position will be recorded, for later use in the ISF calculations
thermalizing_time = 200; % during this initial period, the trajectory won't be recorded.
stop_time = 1024*0.1; % Total simulation time.

% N_steps and N_ISF_steps are calculated after PIGLE adjusts the requested time parameters
max_N_steps = 1e9; % The simulation won't start if the calculated number of steps is greater than max_N_steps. Prevents too-long runs.
max_N_ISF_steps = 6e5; % The simulation won't start if the calculated number of steps is greater than max_N_ISF_steps. Prevents insensible memory explotation

