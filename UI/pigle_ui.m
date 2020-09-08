

%% parameters for config_model.m

z_enabled = 0;

dKz_include_in_isf = 0;

theta_enabled = 1;

thetatilt_enabled=1;%Please do not use tilt if theta is not enabled

% theta will be the Euler angle theta and tilt will be the Euler angle phi
%Currently only implemented for diatomic molecules propely
%Recommended if using angular simulations coupled to an external potential
%
%Currently implemented so that theta goes between 0 and 2 pi in this
%context, giving each oreintation of the molecule two possible
%descriptions: (theta,tilt) and (2pi-theta,tilt+pi). This is to cope with
%the singularity at the poles. Potentials will need to be prepared for
%this. Additionally, to avoid the tilt moment of inertia going to zero when
%theta =0,pi, the moment of inertia instead becomes the moment of inertia 
%of the smallest theta quantisation unit 

zero_p_init = 0; % set initial momentum be set to zero? (if set to 0, p_init will correspond to thermal distribution)
interactions_active = 0;%interabsorbate interactions
uniform_dist=1; % set initial x,y particle distribution to be uniform? If set to 0, the positions will be randomly distributed through the unit cell
N_runs = 5;
run_parallel = 0;

conf3D2D=0; %allows use of a 3D conformation in a 2D spatial sim. Only works for 2 angular dimensions

% Specify dK as a 2D vector, 3rd dim is azimuths.
%dK = [0.05 0.1 0.15 0.2:0.1:1 1.2:0.2:5];
dK=[0.1:0.1:3.6];%[0.05:0.05:3.6];
azim_1 = [1 0];
azim_2 = [cosd(45) sind(45)];

% specify beam parameters and geometrical parameters for scatering calculations
% specify beam parameters and geometrical parameters for scatering calculations
theta_tot = 44.4; % Degrees
beam_ek=9.5;%/meV
beam_ki =KE2k(beam_ek);%3.3977; % Angstrom ^{-1} 

% Specify simulation time parameters
% (those will be adjusted by the program, see below if interested)
sample_time = 1e-3;
sample_time_clist = 1e-2;
isf_sample_time = 5e-2;
thermalizing_time = 50;%50;
stop_time = 1024;

% N_steps and N_ISF_steps are calculated after PIGLE adjusts the requested time parameters
max_N_steps = 1e9;%1e9
max_N_ISF_steps = 6e5;%6e5

