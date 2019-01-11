
%a1 - unitcell dimension, should be defined already
unitcell_area = (sqrt(3)/2)*(a1^2);
% mass1
mu0=3 ; n0= 1 * 1/unitcell_area; alpha=20;
mu_eff=dipole_dipole_repulsion.effective_dipole_moment(mu0,alpha,n0,sum(number_density));
[fparam1, ~] = dipole_dipole_repulsion.f_Kohn_Lau(mu_eff,1);

%mass28
mu0=1 ; n0=1 * 1/unitcell_area; alpha=10;
mu_eff=dipole_dipole_repulsion.effective_dipole_moment(mu0,alpha,n0,sum(number_density));
[fparam2, ~] = dipole_dipole_repulsion.f_Kohn_Lau(mu_eff,1);

% between mass1 and mass28
mu0=2 ; n0=1 * 1/unitcell_area; alpha=15;
mu_eff=dipole_dipole_repulsion.effective_dipole_moment(mu0,alpha,n0,sum(number_density));
[fparam3, ~] = dipole_dipole_repulsion.f_Kohn_Lau(mu_eff,1);