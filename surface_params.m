% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

% surface_params.m holds all the info required to configure the latice
% properties, adsorbates, and the interactions (PES and interadsorbate
% forces). While some effort was made to place user interface variables in
% the beginning of the file, one needs to go through it all to configure
% all the relevant information.

%% Pre-defined variables

% override/add if supp_surface_params.m exists
% This option is used (for example) when generating multi-run job (for
% probing a parameter space).
if exist([pwd '/supp_surface_params.m'],'file')
    supp_surface_params
else
    disp(['No supp_surface_params.m at ' pwd])
end

pigle_ui_surface_params

params.T = T; % Temperature / K

%% Lattice Properties

% unitcel - a struct with the following fields:
%   numOfPrmtvCells - number of primitive cells in the unitcell (along the lateral dimentions)
%   x,y,z,theta - vectors containing the grid
%   celldim     - vector containing the size of the cell [xmax-xmin ymax-ymin .. etc]
params.unitcell = unitcell;

if ~params.z_enabled && isfield('z',params.unitcell)
    params.unitcell = rmfield(params.unitcell,'z');
    params.unitcell.celldim(3)=[];
end
if ~params.theta_enabled && isfield('theta',params.unitcell)
    params.unitcell = rmfield(params.unitcell,'theta');
    params.unitcell.celldim(end)=[];
end

%% Adsorbate(s) Properties

% Temporary total number of particles, actual number will be calculated besed on coverage
params.Nprtcl = Nprtcl_total;

% Mass for each population of particles [amu]
params.mass_list = mass_list;
Nmass = length(params.mass_list);

% Angular mass for each population of particles [amu * Angstrm^2]
params.angular_mass_list = angular_mass_list;

% Density of each population (requested)
% The program will find the size of supercell which can fulfill 
% the densities (based on the requested number of particles), then
% will adjust the number of particles in each population to match the density.
% The true densities will be slightly different than the ones requested,
% and can be calculated in post simulation time. The modified densities and
% number of particles are reported to the user.

params.number_density = [number_density];

if length(params.number_density) ~= length(params.mass_list) ...
        || length(params.angular_mass_list) ~= length(params.mass_list)
    error('numbe_density, mass_list, and angular_mass_list - must have the same length lengthes')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adsorbate configuration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The configuration case refer to prepare_configuration.m.
% For example, case '1' ia for top symmetric molecule, with two parameters,
% r0 and Natoms. The configuration will be a ring with the requested number
% of particles. Later on, after the simmulation had completed, the ISF will
% be calculated using these configurations.

[r_conf(1:Nmass).caseNum]                         = deal(r_conf_case_num{:});
[r_conf(1:Nmass).r0]                              = deal(r_conf_radius{:});
[r_conf(1:Nmass).Natoms]                          = deal(r_conf_Natoms{:});
[CoM_form_factor_conf(1:Nmass).caseNum]           = deal(CoM_form_factor_case_num{:}); % The configuration case refer to form_factor.m.
[CoM_form_factor_conf(1:Nmass).hemisphere_radius] = deal(CoM_form_factor_hemisphere_radius{:});
[form_factor_conf(1:Nmass).caseNum]               = deal(form_factor_case_num{:}); % The configuration case refer to form_factor.m.
[form_factor_conf(1:Nmass).hemisphere_radius]     = deal(form_factor_hemisphere_radius{:});

% TODO remove this awful cod ing practice and in general, transfer all the particle handling to
% proper OOP
for i=1:Nmass
    for j=1:length(form_factor_conf(i).caseNum)
        scatCntr(j).caseNum = form_factor_conf(i).caseNum(j);
        scatCntr(j).hemisphere_radius = form_factor_conf(i).hemisphere_radius(j);
    end
    form_factor_conf(i).scatCntr = scatCntr;
end
form_factor_conf = rmfield(form_factor_conf,'caseNum');
form_factor_conf = rmfield(form_factor_conf,'hemisphere_radius');

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Translational Friction %
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Depending on the case (see calc_A function in calculate_sim_params.m),
% different fields of the structure A_struct are expected. For new cases in
% calculate_sim_params, ammend also surface_params.m
[A_strct(1:Nmass).A_case] = deal(A_case{:});
[A_strct(1:Nmass).w0]   = deal(A_w0{:});
[A_strct(1:Nmass).dw]   = deal(A_dw{:});
[A_strct(1:Nmass).eta]  = deal(A_eta{:});
[A_strct(1:Nmass).tau]  = deal(A_tau{:});

[params.prtcl(1:Nmass).A_spatial_depended_friction]    = deal(A_spatial_depended_friction{:});

% Assign PESlike nD-array to each population, and scale friction
% Order of dimensions: 1st - 'y', 2nd - 'x'. If 'z' dimension is enabled,
% then the 'z' dimension is assumed to be 3rd, last dimension is angular (either 3rd or 4th).
% So for example, PotMatrix(i,j) is the surface coordinate (j,i).
for i=1:length(params.mass_list)
    if params.prtcl(i).A_spatial_depended_friction == 1
        [PESlike] = friction_func_list{i}(friction_arg_list{i,:});
        PESlike = prepFuncs.adjust_PES(params,PESlike); % resuce PotMatrix to match the enabled dimensions
        params.prtcl(i).friction.scaleMat = PESlike;
    else
        params.prtcl(i).friction.scaleMat = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%
% Angular Friction %
%%%%%%%%%%%%%%%%%%%%
% Depending on the case (see calc_A function in calculate_sim_params.m),
% different fields of the structure A_struct are expected.
[A_theta_strct(1:Nmass).A_case] = deal(A_theta_case{:});
[A_theta_strct(1:Nmass).w0]     = deal(A_theta_w0{:});
[A_theta_strct(1:Nmass).dw]     = deal(A_theta_dw{:});
[A_theta_strct(1:Nmass).eta]    = deal(A_theta_eta{:});
[A_theta_strct(1:Nmass).tau]    = deal(A_theta_tau{:});

[params.prtcl(1:Nmass).A_spatial_depended_theta_friction]    = deal(A_spatial_depended_theta_friction{:});

% Assign PESlike nD-array to each population, and scale friction
% Order of dimensions: 1st - 'y', 2nd - 'x'. If 'z' dimension is enabled,
% then the 'z' dimension is assumed to be 3rd, last dimension is angular (either 3rd or 4th).
% So for example, PotMatrix(i,j) is the surface coordinate (j,i).
for i=1:length(params.mass_list)
    if params.prtcl(i).A_spatial_depended_theta_friction == 1
        [PESlike] = theta_friction_func_list{i}(theta_friction_arg_list{i,:});
        PESlike = prepFuncs.adjust_PES(params,PESlike); % resuce PotMatrix to match the enabled dimensions
        params.prtcl(i).friction.theta_scaleMat = PESlike;
    else
        params.prtcl(i).friction.theta_scaleMat = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface-Adsorbate Potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The PES is constructed for each population, with few help functions.
%
%    PES_func_list - cell-array of pointers (one per population) to functions
%                 used to generate the PES. A list of functions can be
%                 found under <PIGLE root>/generatePES/PES_library.
%                 If the PES is saved as a mat file (and there is no need for
%                 generation function), use the function loadPES.
%
%    PES_arg_list  - Cell-array.  PES_arg_list{i} holds the parameters to be
%                    passed to the function PES_func_list{i}.
%
%    params_for_function_prepare_potential - a matlab file to hold variables
%                 to be placed in PES_arg_list. The file serves as a 'buffer',
%                 to avoid changing the values in this (surface_params.m) file.

% Pointers to the population specific functions for generating the PES
if ~exist('PES_func_list','var')
    disp('PES_func_list is not defined. Assuming PES_func_list={@prepare_potential,@prepare_potential}')
    PES_func_list = {@prepare_potential, @prepare_potential};
    %PES_func_list = {@loadPES, @loadPES, @loadPES};
end

% Define the arguments for PES generation.
% PES_arg_list(i) is the list of variables to pass to the function PES_func_list(i). 
% If PES_func_list(i)=@loadPES (PES saved as a mat file),
% use PES_arg_list(i)={'filename', 'var name'}, where 'filename' is the full path to
% the mat file, and 'var name' is the variable in the mat file which contains the PES.
if ~exist('PES_arg_list','var')
    disp('PES_arg_list is not defined. Assuming PES_arg_list={params.unitcell, pot_strct(1); params.unitcell, pot_strct(1)}')
    
    % Define variables to hold the arguments which are stored in PES_arg_list (see below)
    if exist('params_for_function_prepare_potential','file')
        params_for_function_prepare_potential;
    end
    PES_arg_list = {params.unitcell, pot_strct(1); ...
        params.unitcell, pot_strct(1)};
end

% Assign PES to each population, and workout forces.
% Order of dimensions: 1st - 'y', 2nd - 'x'. If 'z' dimension is enabled,
% then the 'z' dimension is assumed to be 3rd, last dimension is angular (either 3rd or 4th).
% So for example, PotMatrix(i,j) is the surface coordinate (j,i).
for i=1:length(params.mass_list)
    [PotMatrix] = PES_func_list{i}(PES_arg_list{i,:});
    PotMatrix = prepFuncs.adjust_PES(params,PotMatrix); % resuce PotMatrix to match the enabled dimensions
    params.prtcl(i).pes.PotMatrix = PotMatrix;
     
    % Calc force from the PES
    % gradient acts funny, gives the deriv on sec dimension first,
    %then deriv on 1st dimension, then the rest. But for spacing - acts normal (PES,spacing1, spacing2, ....)
    dx = diff(params.unitcell.x(1:2)); dy = diff(params.unitcell.y(1:2));
    if params.model_dim ~= length(size(PotMatrix))
        error('PotMatrix must have the dimentions of the model, even if dummy dimentions')
    elseif params.model_dim == 2
        [fx,fy] = gradient(PotMatrix,dy,dx);
    elseif params.model_dim == 3 && params.z_enabled
        dz = diff(params.unitcell.z(1:2));
        [fx,fy,fz] = gradient(PotMatrix,dy,dx,dz);
    elseif params.model_dim == 3 && params.theta_enabled
        dtheta = diff(params.unitcell.theta(1:2));
        [fx,fy,ftheta] = gradient(PotMatrix,dy,dx,dtheta);
    else % 4D
        dz = diff(params.unitcell.z(1:2)); dtheta = diff(params.unitcell.theta(1:2));
        [fx,fy,fz,ftheta] = gradient(PotMatrix,dy,dx,dz,dtheta);
    end
        
    params.prtcl(i).pes.fx = -fx*k1;
    params.prtcl(i).pes.fy = -fy*k1;
    if exist('fz','var'), params.prtcl(i).pes.fz = -fz*k1; end
    if exist('ftheta','var'), params.prtcl(i).pes.ftheta = -ftheta*k1; end
end

%% Parameters for Interactions

% To turn the interaction off/on, see config_model.m

% assign functions to species:
% For each pair in f_perm, a function case is defined in f_func. This
% function case is taken from f_interaction.m - and f_func_params contain
% the arguments for each function. The force dimensions are expected to be
% in meV/Angstrom
params.interactions.f_perm = f_perm;
params.interactions.f_func = f_func;
params.interactions.f_func_params = f_func_params;

% Define the boundaries for interactions:
% out_cutoff_r - The supercell must be larger than that number (see calculate_sim_params.m).
%                TODO: include in connection lists, once implemented
% in_cutoff_r -  the force between particles will be calculated for r >= r_in
params.interactions.out_cutoff_r = out_cutoff_r;
params.interactions.in_cutoff_r  = in_cutoff_r;

% Decimate the interadsorbate force
% x - the points in which the force is to be calculated
% Fint(i,x) = the decimated force for the ith permutation. There should be
% N*(N+1)/2 permutations in total, where N is the number of populations.

params.interactions.x = x_interactions; % in Angstrom

for i=1:length(params.interactions.f_func)
    % force is in meV/Angstrm
    params.interactions.Fint(:,i) = ...
        f_interaction(params.interactions.x,params.interactions.f_func(i),params.interactions.f_func_params{i})*k1;
end
