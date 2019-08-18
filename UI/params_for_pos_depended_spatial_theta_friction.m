
% defind a structure to be passed as argument to the function which create
% the PES-like array to scale the friction

%% PES-like 1

theta_friction_strct(1).ref_De = [1600 1600 1600 1600 1600 1600]*0;

theta_friction_strct(1).V = [0.5 0.5 ; 0.5 0.5; 0.5 0.5; 1 1 ; 1 1 ; 1 1 ]*1; % top, slope1,slop2,bridge,hcp,fcc

theta_friction_strct(1).is_potval = [1 0 0 1 1 1];
theta_friction_strct(1).a = [1 NaN NaN 1 1 1]*1;
theta_friction_strct(1).r_e = [2 NaN NaN 2 2 2];
theta_friction_strct(1).f_2D = @hexagonal6interp;

%% PES-like 2

theta_friction_strct(2) = theta_friction_strct(1);
%theta_friction_strct(2).V = [1 1 ; 0.5 0.5; 0.5 0.5; 1 1 ; 1 1 ; 1 1 ]*1; % top, slope1,slop2,bridge,hcp,fcc
