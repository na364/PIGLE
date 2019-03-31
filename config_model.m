% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

% Config_model.m allows one to both configure the simulation and the post analysis
% routines (if enabled). For the simulation, here you can define the dimentions of
% the simulation (2D-4D), enable inter-adsorbate interactions (and define), enable
% parallelization and time constants (step, decimation/isf, total, thermalization). For the post-analysis (calculation of the ISF), define whether to
% include the prependicular direction in calculatin gthe ISF, the number of repetitions
% (N_run) and the momentum transfer in which the ISF is to be calculated. 

%% Fundamental Constants (Do not touch)
k1=9.6485335; % meV/Angstrm to a.m.u * Angstrm/(ps)^2 ... 1.6021766208e-19 / 1000 * 1/1.660539e-27 * 1e20 / 1e24
k2=1.2023514; % 1/k_B
params.k_B = 0.8317035; % Boltzmann constant in Aˆ2 amu psˆ-2 Kˆ-1 .... 8.6173303e-5 eV/K * k1

%% Sim general params (user defined)

pigle_ui

params.z_enabled = z_enabled;
params.dKz_include_in_isf = dKz_include_in_isf;
params.theta_enabled = theta_enabled;
params.zero_p_init = zero_p_init; % set initial momentum be set to zero? (if set to 0, p_init will correspond to thermal distribution)
params.interactions.active = interactions_active;
params.N_runs = N_runs;
params.run_parallel = run_parallel;

% Specify simulation time parameters
% (those will be adjusted by the program, see below if interested)
params.sample_time = sample_time;
params.sample_time_clist = sample_time_clist;
params.nclist            = floor(sample_time_clist/sample_time);
if params.nclist < 1, error('sample_time > sample_time_clist, can not proceed'); end
params.isf_sample_time = isf_sample_time;
params.thermalizing_time = thermalizing_time;
params.stop_time = stop_time;

%% Parameters for scattering calculations

% Specify dK as a 2D vector, 3rd dim is azimuths.
params.dK(:,:,1) = azim_1.*dK';
params.dK(:,:,2) = azim_2.*dK';

% beam props and total scattering angle
params.beam_ki                = beam_ki;
params.theta_tot                = theta_tot;

%% Config dimentions and parallelization

if ~params.z_enabled, params.dKz_include_in_isf = 0; end

params.model_dim = 4-(~params.z_enabled)-(~params.theta_enabled);

if params.run_parallel ~= 0
    c = parcluster('local'); % build the 'local' cluster object
    jobStorageLocation = [pwd '/tmp_cluster_data'];
    if ~exist(jobStorageLocation,'dir'), mkdir(jobStorageLocation); end
    c.JobStorageLocation = jobStorageLocation;
    params.n_parall_runs = c.NumWorkers; % get the number of workers
else
    params.n_parall_runs = 1;
end
if params.n_parall_runs > params.N_runs
    params.n_parall_runs = params.N_runs;
end
if params.n_parall_runs > 1 && params.N_runs > 1 && isempty(gcp('nocreate'))
   parpool(c, params.n_parall_runs)
end

%% Adjust timing

% adjust time parameters to be
% (a) Power of 2 - for FFT efficiency
% (b) Decimation compatible, etc ...
ntscat=2^nextpow2(params.stop_time/params.isf_sample_time);
params.stop_time=(ntscat-1)*params.isf_sample_time;
params.decimation=round(params.isf_sample_time/params.sample_time);
params.sample_time = params.isf_sample_time/params.decimation;
params.t=linspace(0,params.stop_time,ntscat);

% Add thermalization time, and recalculate steps
params.stop_time = params.stop_time + params.thermalizing_time;
params.N_steps = params.stop_time/params.sample_time;
params.N_ISF_steps = params.stop_time/params.isf_sample_time;
disp(['N_steps = ' num2str(params.N_steps) ' .... N_ISF_steps = ' num2str(params.N_ISF_steps)])

if params.N_steps > max_N_steps, error('N_steps is too high'); end
if params.N_ISF_steps > max_N_ISF_steps, error('N_ISF_steps is too high'); end

%% Add surface and Adsorbates

% Load surface parameters
surface_params

% Calculate the rest of the simulation parameters based on these.
params = calculate_sim_params(params, A_strct, A_theta_strct, r_conf);

%% Apply In-Phase Scattering Condition
% If interactions are enabled:
% phase from shadow scatterers (outside of unitcell) must be equal
% exp(-i*((R+a)*dK)) = exp(-i*((R)*dK)) ==> 2*pi*n = dK*a where a is of the
% supercell
if params.interactions.active
    disp('#############################################')
    disp('# dK is changed to comply with interactions #')
    disp('#############################################')
    G_supercell = 2*pi./params.supercell.celldim(1:2);
    for i=1:size(params.dK,3)
        params.dK(:,:,i) = G_supercell.*round(params.dK(:,:,i)./G_supercell);
    end
end

% Calculate particle configuration (inc. form factors)
for i=1:length(params.prtcl)
    params.prtcl(i).conf = prepare_configuration(r_conf(i),'CoM_form_factor_conf',CoM_form_factor_conf(i),'form_factor_conf',form_factor_conf(i),'dK',params.dK,'beam_ki',params.beam_ki,'theta_tot',params.theta_tot);
end

%% Calculate initial conditions:
params = prepFuncs.r_init(params);
if params.zero_p_init == 1, disp('initial momentum is set to 0'); end
params = prepFuncs.p_init(params);

%%
if exist('./pigle_data.mat','file'), disp('pigle_data.mat exists. DELETING'); delete pigle_data.mat; end
save('pigle_data.mat', 'params');
assignin('base', 'params', params);
