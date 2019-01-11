
% defind a structure to be passed as argument to the function which create the potential

%% Pot 1

pot_strct(1).ref_De = [1600 1600 1600 1600 1600 1600]*0;

%pot_strct(1).V = [192 168 ; 0.5 0.5; 0.5 0.5; 18 62 ; 112 104 ; 156 118 ]*1; % top, slope1,slop2,bridge,hcp,fcc
%pot_strct(1).V = [0 0 ; 0.5 0.5; 0.5 0.5; 80 80 ; 80 80 ; 80 80 ]*1; % top, slope1,slop2,bridge,hcp,fcc
pot_strct(1).V = [250 ; 0.5 ; 0.5 ; 110 ; 22 ; 0]; % top, slope1,slop2,bridge,hcp,fcc

pot_strct(1).is_potval = [1 0 0 1 1 1];
pot_strct(1).a = [1 NaN NaN 1 1 1]*1;
pot_strct(1).r_e = [2 NaN NaN 2 2 2];
pot_strct(1).f_2D = @hexagonal6interp;

%% Pot 2

pot_strct(2) = pot_strct(1);
pot_strct(2).V = [exp(15*number_density(2))-1 ; 0.5 ; 0.5 ; 240 ; 270 ; 280]; % top, slope1,slop2,bridge,hcp,fcc
