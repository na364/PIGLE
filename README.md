% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

PIGLE - Particles Interacting in Generalized Langevin Equation simulator
Release: 1.0rc1 (a - alpha , b - beta, rc - release candidate)
Please cite PIGLE software using doi:// xxxx

PIGLE simulator solves the Generalized Langevin Equation for interacting particles, in a 4D potential energy surface.
The potential has up to 3 spatial dimensions, and a rigid-body like rotation. Particles doesn't neceseraly share the
same properties (such as mass or friction). When the third spatial dimension ('z') is enabled, the user can choose 
either constant pressure (where desorbing particles re-appear at the system), or in simple mode where desorbing particles
are being frozen for later filtration. However, currently there is not much physical use for that feature, and it might
be removed in future releases. Following the calculation of the particles trajectories, each particle (center of mass)
can be replaced with a configuration of particles, for example an eight memeber ring of equally contributing scattering
centers. Such an approach allows the calculations of the intermediate scattering function with a less simplified
representation of the molecular form factor, and also allows to visualize rotations.

####################
# Getting started: #
####################

Installation:
-------------
- Install Matlab (tested with 2017b). If working under Linux, make sure to obtain a supported version of a supported comipiler.
- Edit prep_environment.m

Execution of jobs:
------------------
- Create a potential of your choise and save it as 'mat' file, or create a function which generates the potential.
  One can use the functions and scripts in the subfolder generatePES.
- Configure the parameters in m-files under the subfolder UI.
- Run 'run_pigle.m'

Visualization:
--------------
The 'make_movie' file will help you to visuallize the dynamics.

Parallel Computing:
-------------------
PIGLE support parallel computing (via matlab/Simulink support).


####################
# List of Files:   #
####################

PIGLE:
aux_files        config_model.m   generatePES  make_movie.m     pigle_sim               prep_environment.m  run_pigle.m      sweepParams  UI
benchmark  f_interaction.m  LICENSE.txt  pigle_data.mat  prepare_configuration.m  README.txt          surface_params.m  TODO.txt

./aux_files:
calc_kinetic_energy.m    calculate_sim_params.m                        generate_A_from_function.m  make_data_path.m  resample_data.m  timeinput.m
calculate_average_isf.m  generate_A_from_frequencies_multiple_gamma.m  hlp_f.m                     prepFuncs.m       sim_gle_nd.m

./benchmark:
analytic_gle.m  analytic_le.m  benchmark_biexp_gle.m  biexp_gle_isf.m  incoh_brownian_rot.m  incoh_langevin_rot.m  plot_diffusion_models.m

./generatePES:
loadPES.m  params_for_function_prepare_potential.m  PES_library  prepare_potential.m

./generatePES/PES_library:
hexagonal6interp.m  hexagonal.m

./pigle_sim:
create_model.m  delete_unconnected_lines.m  sl_interactions.slx  sl_pigle_Population.slx

./sweepParams:
config_job.m  config_job_params.m  pigle_run_single_task.m  pigle_shell_params.m  run_job.m

./UI:
pigle_ui.m  pigle_ui_surface_params.m  pigle_wrapper_params.m



