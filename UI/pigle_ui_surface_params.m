
%% params for surface_params.m
T=140;
Nprtcl_total = 3;
mass_list = [16 16];
radius = [0.2 1.5];
number_density = 0.01*[2 1];
eta =4 ; eta2=4;
eta_theta = 6; eta_theta2 = 6;
tau = [1];

a1=3.6147/sqrt(2);                          % Copper 111 lattice constant in Angstrom
%a1=2.71;                          % Ru 0001 lattice constant in Angstrom
x0 = 0; nx = 30; xdim = a1;                 % x dimention params of the unitcell/PES
y0 = 0; ny = 50; ydim = a1*sqrt(3.0);       % y dimention params of the unitcell/PES
z0 = 0; nz = 20; zdim = 10;                 % z dimention params of the unitcell/PES
theta0 = 0; ntheta = 20; thetadim = 2*pi/6; % theta dimention params of the unitcell/PES
numOfPrmtvCells = [1 2]; % How many primitive cells exist in the XY potential

unitcell = prepFuncs.make_unitcell([nx xdim x0],[ny ydim y0],'z',[nz zdim z0 z_enabled],'theta',[ntheta thetadim theta0 theta_enabled],'numOfPrmtvCells',numOfPrmtvCells);

angular_mass_list = mass_list.*radius.^2;

Nmass = length(mass_list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adsorbate configuration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See prepare_configuration.m.
% If new adsorbate_conf_case_num (in addition to '1') are implemented in
% prepare_configuration.m, surface_params.m needs to be updated.
% Case '1' is for top symmetric molecule
% Each parameter needs to be stored either as a single cell (which will be distributed to all
% populations), or as a cell-array (with the i'th element distributed to the i'th
% population)
%
r_conf_case_num = {1};
r_conf_radius   = {0};
r_conf_Natoms   = {1};
CoM_form_factor_case_num = {0}; % The configuration case refer to form_factor.m.
CoM_form_factor_hemisphere_radius = {2.5};
form_factor_case_num = {repmat(0,r_conf_Natoms{1},2)}; % The configuration case refer to form_factor.m.
form_factor_hemisphere_radius = {repmat(2.5,r_conf_Natoms{1},2)};

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Translational Friction %
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Depending on the case (see calc_A function in calculate_sim_params.m),
% different fields of the structure A_struct are expected. For new cases in
% calculate_sim_params, ammend also surface_params.m
% Each parameter needs to be stored either as a single cell (which will be distributed to all
% populations), or as a cell-array (with the i'th element distributed to the i'th
% population)
%
A_case = {1};
A_w0   = {eta};
A_dw   = {1./tau};
A_eta  = {eta, eta2};
A_tau  = {tau};
A_spatial_depended_friction = {0}; % enable/disable spatial depended friction

% To define spatial dependent friction, we provide function+variables to
% generate a PES-like multi-dimensional array, which will be used to scale
% each element of the friction matrix. For example, if the PES-like has 4
% dimensions (ie xyz+theta), and the 'A' matrix of the friction (for a
% specific point in space) is 2D, the final 'A' would be of 6 dimensions.

% Pointers to the population specific functions for generating the PES-like
% scaling factor for the friction
friction_func_list = {@prepare_potential};

% Define variables to hold the arguments which are stored in friction_arg_list (see below)
params_for_pos_depended_spatial_friction

% Define the arguments for friction generation.
friction_arg_list = {unitcell, friction_strct(1)};

%%%%%%%%%%%%%%%%%%%%
% Theta (Angular) Friction %
%%%%%%%%%%%%%%%%%%%%
% Depending on the case (see calc_A function in calculate_sim_params.m),
% different fields of the structure A_struct are expected.
% Each parameter needs to be stored either as a single cell (which will be distributed to all
% populations), or as a cell-array (with the i'th element distributed to the i'th
% population)
%
A_theta_case = {1};
A_theta_w0   = {eta_theta};
A_theta_dw   = {1};
A_theta_eta  = {eta_theta};
A_theta_tau  = {1};
A_spatial_depended_theta_friction = {0}; % enable/disable spatial depended angular friction

% To define spatial dependent friction, we provide function+variables to
% generate a PES-like multi-dimensional array, which will be used to scale
% each element of the friction matrix. For example, if the PES-like has 4
% dimensions (ie xyz+theta), and the 'A' matrix of the friction (for a
% specific point in space) is 2D, the final 'A' would be of 6 dimensions.

% Pointers to the population specific functions for generating the PES-like
% scaling factor for the friction
theta_friction_func_list = {@prepare_potential};

% Define variables to hold the arguments which are stored in friction_arg_list (see below)
params_for_pos_depended_spatial_theta_friction

% Define the arguments for friction generation.
theta_friction_arg_list = {unitcell, theta_friction_strct(1)};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface-Adsorbate Potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pointers to the population specific functions for generating the PES
PES_func_list = {@prepare_potential,@prepare_potential};
%PES_func_list = {@loadPES, @loadPES, @loadPES};

% Define variables to hold the arguments which are stored in PES_arg_list (see below)
params_for_function_prepare_potential

% Define the arguments for PES generation.
PES_arg_list = {unitcell, pot_strct(1);unitcell, pot_strct(2)};

%% Parameters for Interactions

prepare_params_for_interactions

% assign functions to species:
% For each pair in f_perm, a function case is defined in f_func. This
% function case is taken from f_interaction.m - and f_func_params contain
% the arguments for each function.
f_perm = [1 1;1 2;2 2];
f_func = [repmat(1,1,3)];
f_func_params = {[fparam1_12 13 fparam1_6 7],[fparam2_12 13 fparam2_6 7],[fparam3_12 13 fparam3_6 7]};

% Define the boundaries for interactions:
% out_cutoff_r - The supercell must be larger than that number (see calculate_sim_params.m).
%                TODO: include in connection lists, once implemented
% in_cutoff_r -  the force between particles will be calculated for r >= r_in
out_cutoff_r = norm(unitcell.celldim)*10*0+49;
in_cutoff_r = 0.5;

% x_interactions - the points in which the force is to be calculated
x_min = in_cutoff_r/10; % in Angstrom
x_max = out_cutoff_r + 1; % in Angstrom
numOfPoints_interactions = 500;
x_interactions = linspace(x_min, x_max, numOfPoints_interactions); %x_min/max in Angstrom

