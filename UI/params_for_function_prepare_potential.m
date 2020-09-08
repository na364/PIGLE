
%Defined a structure to be passed as argument to the function which create the potential
%potential values in meV assumed

%If you wish to use a specific 5D potential, use @loadPES on line 168 of
%pigle_ui_surface_params.m, and create a file with a function that
%generates a 5D matrix of potential values. This should take the form
%PotMatrix(x,y,z,theta,tilt).

%% Pot 1


pot_strct(1).ref_De = [1600 1600 1600 1600 1600 1600]*0;
%Sets depth of Morse potential well in z direction, set to 0 if not using z in simulation


%pot_strct(i).V contains parameters for the surface-adsorbate potential for
%the ith species. What each parameter means depends on the potential
%function provided. Potential values are typically in meV.

% pot_strct(1).V = [192 168 ; 0.5 0.5; 0.5 0.5; 18 62 ; 112 104 ; 156 118 ]*1; % top, slope1,slop2,bridge,hcp,fcc
% pot_strct(1).V = [0 45 ; 0.5 0.5; 0.5 0.5; 80 100 ; 80 100 ; 80 100 ]*1; % top, slope1,slop2,bridge,hcp,fcc
% pot_strct(1).V = [0 ; 0.5 ;  0.5; 250 ; 210 ; 60 ]*1; % top, slope1,slop2,bridge,hcp,fcc
% pot_strct(1).V = [250 ; 0.5 ; 0.5 ; 110 ; 22 ; 0]; % top, slope1,slop2,bridge,hcp,fcc
%pot_strct(1).V = [200 ; 0.1 ; 0.1 ; 5   ; 5  ;50]; % top, slope1,slop2,bridge,hcp,fcc
%pot_strct(1).V=[200;50;100;0.5]; %max,min,saddle,slopef
% pot_strct(1).V=[135;0;115;0.5]%; %CO Cu-100 empirical
%pot_strct(1).V=[235;115;0];%Hollow, bridge, Tmode

%Use prepare_potential to see exact functional dependance
pot_strct(1).anglecase=3; %0=no angle dependence % 1=theta dependence %2=theta and tilt dependence 3 cu100

pot_strct(1).V=[115;235;0;32;300;26;235;135;12;16;0.2];%;115]; %diffusion barrier,
%hollow, top, top Tmode, top R mode, bridge T mode, bridge Rmode1 -transverse, bridge
%Rmode2, mass carbon, mass oxygen, decay length


 pot_strct(1).is_potval=[1 1 1 0 0 0 0 0 0 0 0 0 1];
 %Is the corresponding entry in pot_strct.V a potential value? Causes scaling of this value with z potential if set to 1
 
 
pot_strct(1).a=[1 1 1 NaN NaN NaN NaN NaN NaN NaN NaN 1];
%Sets width of Morse potential in z direction in angstroms


 pot_strct(1).r_e=[2 2 2 NaN NaN NaN NaN NaN NaN NaN NaN 2];
 %Sets equilibrium height in z direction in angstroms



%If a parameter is not a potential value, set minmax= 0 0
pot_strct(1).theta_minmax=[5 0; 5 0; 0 0];%amplitude /meV, angle of minimum energy in radians
pot_strct(1).tilt_minmax=[5 0;0  0;5 0;0 0];%amplitude /meV, angle of minimum energy in radians
%used for angular potential cases 1 & 2. This will modify each of the
%potential values in pot_strct.V by as
%V'=V+amplitude*(cos(angle-minimum_angle+pi/2))^2


% pot_strct(1).is_potval=[1 1 0];
% pot_strct(1).a=[1 1 NaN];
% pot_strct(1).r_e=[2 2 NaN];


% pot_strct(1).is_potval = [1 0 0 1 1 1];%Is each entry a potential or a scale factor?
% pot_strct(1).a = [1 NaN NaN 1 1 1]*1;%width of z pot well
% pot_strct(1).r_e = [2 NaN NaN 2 2 2]; % equilibrium height
% pot_strct(1).is_potval=[1 1 1 0];
% pot_strct(1).a=[1 1 1 NaN];
% pot_strct(1).r_e=[2 2 2 NaN];

pot_strct(1).f_2D = @cu100_2;%set to be a function on the path which creates the desired 2D potential.

% pot_strct(1).ref_De = [1600 1600 1600 1600 1600 1600]*0;
% pot_strct(1).V = [40 ; 1 ;  4 ]; % scalingF, A, p : for hexagonal.m
% pot_strct(1).is_potval = [0 0 0];
% pot_strct(1).a = [NaN NaN NaN]*1;
% pot_strct(1).r_e = [NaN NaN NaN];
% 
% pot_strct(1).f_2D = @hexagonal;


%% Pot 2

%PES for the second species of adsorbate

%pot_strct(2) = pot_strct(1);
%pot_strct(2).V = [exp(15*number_density(2))-1 ; 0.5 ; 0.5 ; 240 ; 270 ; 280]; % top, slope1,slop2,bridge,hcp,fcc
