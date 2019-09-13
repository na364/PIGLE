
% %a1 - unitcell dimension, should be defined already
unitcell_area = (sqrt(3)/2)*(a1^2);

%% dipole dipole
% mu0=3 ;
% n0= 1 * 1/unitcell_area; alpha=20;
% mu_eff=dipole_dipole_repulsion.effective_dipole_moment(mu0,alpha,n0,sum(number_density));
mu_eff=5;
[fparam1, ~] = dipole_dipole_repulsion.f_Kohn_Lau(mu_eff,1);


%% H & CO
% % mass1
% mu0=6 ; n0= 1 * 1/unitcell_area; alpha=20;
% mu_eff=dipole_dipole_repulsion.effective_dipole_moment(mu0,alpha,n0,sum(number_density));
% [fparam1, ~] = dipole_dipole_repulsion.f_Kohn_Lau(mu_eff,1);

% %mass28
% mu0=6 ; n0=1 * 1/unitcell_area; alpha=10;
% mu_eff=dipole_dipole_repulsion.effective_dipole_moment(mu0,alpha,n0,sum(number_density));
% [fparam2, ~] = dipole_dipole_repulsion.f_Kohn_Lau(mu_eff,1);
% 
% % between mass1 and mass28
% mu0=6 ; n0=1 * 1/unitcell_area; alpha=20;
% mu_eff=dipole_dipole_repulsion.effective_dipole_moment(mu0,alpha,n0,sum(number_density));
% [fparam3, ~] = dipole_dipole_repulsion.f_Kohn_Lau(mu_eff,1);


%% u=4*epsilon*(sigma^12.*r.^(-12) - sigma^6.*r.^(-6)); plot(r,u)
% 
% % H-H
% u0 = -500; r0=2;
% epsilon = u0/(r0^(-4) - 2*r0^(-2));
% sigma=(0.5*r0^4)^(1/6);
% 
% fparam1_12 = -12 * 4*epsilon*sigma^12; % r^-{13}
% fparam1_6  =  -6 * 4*epsilon*sigma^6; % r^-{7}
% 
% %H-O
% u0 = -500; r0=2;
% epsilon = u0/(r0^(-4) - 2*r0^(-2));
% sigma=(0.5*r0^4)^(1/6);
% 
% fparam2_12 = -12 * 4*epsilon*sigma^12; % r^-{13}
% fparam2_6  =  -6 * 4*epsilon*sigma^6; % r^-{7}
% 
% % O-O
% u0 = -500; r0=2;
% epsilon = u0/(r0^(-4) - 2*r0^(-2));
% sigma=(0.5*r0^4)^(1/6);
% 
% fparam3_12 = -12 * 4*epsilon*sigma^12; % r^-{13}
% fparam3_6  =  -6 * 4*epsilon*sigma^6; % r^-{7}

% O-O LJ from TIP4P/2005
sigma = 3.1589; %Angstrm
epsilon = 93.2*params.k_B; % K ... Boltzmann constant in Aˆ2 amu psˆ-2 Kˆ-1
fparam1_12 = -12 * 4*epsilon*sigma^12; % r^-{13}
fparam1_6  =  6 * 4*epsilon*sigma^6; % r^-{7}

fparam2_12 = -12 * 4*epsilon*sigma^12; % r^-{13}
fparam2_6  =  6 * 4*epsilon*sigma^6; % r^-{7}

fparam3_12 = -12 * 4*epsilon*sigma^12; % r^-{13}
fparam3_6  =  6 * 4*epsilon*sigma^6; % r^-{7}
