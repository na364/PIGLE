% Copyright (c) 2018, Peter Townsend
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

% benchmark_biexp_gle.m

% Script to compare single and biexp kernel gle isf
% Motivation: use the biexp as benchmark for lampd

%close all;

m_amu=104.15; % flavour of the month
amu_cmu_conversion=0.10375; % one amu in cmu mass units
mass=m_amu*amu_cmu_conversion;
friction1=1;
omega1=1;
friction2=5; % equal strength of friction
omega2=0.2;
kb=0.08625;
T=150;

%% Make a nice plot of analytic ISFs
analyticFlatISFs=figure;
hold on;
box on;
setime=1:0.01:100;
DK=1.0;
colorstring='rgm';
annotation_horiz=[0.20 0.40 0.47 0.15];
annotation_heights=[0.20 0.24 0.42 0.90]; % heights at which to place text box for each curve

isf=analytic_gle(DK,setime,friction1,omega1,kb,T,mass); % do the first line for infinite wc outside the k-loop
plot(setime,isf,[colorstring(1),'-'],'LineWidth',2);
an=annotation('textbox',[annotation_horiz(1) annotation_heights(1) 0.1 0.1],'String',['$\omega_c=',num2str(omega1),'$'],'FitBoxToText','on');
set(an,'interpreter','latex','FontSize',14);

isf=analytic_gle(DK,setime,friction2,omega2,kb,T,mass); % plot one-exp gle result for comparison
plot(setime,isf,[colorstring(2),'-'],'LineWidth',2);
an=annotation('textbox',[annotation_horiz(2) annotation_heights(2) 0.1 0.1],'String',['$\omega_{c}=',num2str(omega2),'$']);
set(an,'interpreter','latex','FontSize',14)

isf=biexp_gle_isf(DK,setime,friction1,omega1,friction2,omega2,kb,T,mass); % plot biexp gle result with negligible second cpt
plot(setime,isf,[colorstring(3),'-'],'LineWidth',2);
an=annotation('textbox',[annotation_horiz(3) annotation_heights(3) 0.1 0.1],'String',['$\omega_{1}=',num2str(omega1),', \omega_{2}=',num2str(omega2),'$']);
set(an,'interpreter','latex','FontSize',14)

xlabel('$t$ (ps)','interpreter','latex','FontSize',20);
ylabel('$I(\Delta K,t)$','interpreter','latex','FontSize',20);
set(gca,'FontSize',20);
set(gca,'TickLabelInterpreter','latex');
set(gca,'XScale','log');
ax=gca;
ax.LineWidth=1;
