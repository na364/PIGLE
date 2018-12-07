% Script to evalute incoherent rotational diffusion lineshape and plot for
% a reasonable parameter set

%close all;
r=4; % radius of a molecule
I=104.15*4^2; % something like n*m*r^2 where n = num atoms, m = atomic mass [amu], r = molecular radius [angstroms]
kb=0.8314; % Boltzmann constant in [Angstrom]ˆ2 amu psˆ-2 Kˆ-1
T=298.92;
g=2.0; % rotational friction
N=3000; % number of bessel functions to include in sum
DK=1.0; % experimentally realistic dk
t=0:0.05:200; % short tSE as fast motion, no barriers

%% Plot a selection of models
figure;
isf=incoh_langevin_rot(r,I,kb,T,g,N,DK,t);
plot(t,isf,'b-','LineWidth',2);
ylim([0 1]);
hold on;
isf=incoh_brownian_rot(r,I,kb,T,g,N,DK,t);
plot(t,isf,'r--','LineWidth',2);
set(gca,'TickLabelInterpreter','latex','FontSize',20);
xlabel('$t\,$(ps)','interpreter','latex','FontSize',20);
ylabel('$I(\Delta K,t)$','interpreter','latex','FontSize',20);