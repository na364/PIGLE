
% defind a structure to be passed as argument to the function which create
% the PES-like array to scale the friction

%% PES-like 1

friction_strct(1).ref_De = [1600 1600 1600 1600 1600 1600]*0;

% friction_strct(1).V = [0.1 0.1 ; 0.1 0.1; 0.1 0.1; 1 1 ; 1 1 ; 1 1 ]*1; % top, slope1,slop2,bridge,hcp,fcc
friction_strct(1).V = [1 ; 0.01; 0.01; 0.2 ; 0.2 ; 1 ]*1; % top, slope1,slop2,bridge,hcp,fcc

friction_strct(1).is_potval = [1 0 0 1 1 1];
friction_strct(1).a = [1 NaN NaN 1 1 1]*1;
friction_strct(1).r_e = [2 NaN NaN 2 2 2];
friction_strct(1).f_2D = @hexagonal6interp;

%% PES-like 2

friction_strct(2) = friction_strct(1);
% friction_strct(2).V = [1 1 ; 0.8 0.8; 0.8 0.8; 1 1 ; 1 1 ; 1 1 ]*1; % top, slope1,slop2,bridge,hcp,fcc
