% Copyright (c) 2018, Peter Townsend.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function isf=analytic_gle(DK,time,gamma,wD,kb,T,m)

% Computes analytic ISF for particle obeying GLE on flat surface.

C=wD-gamma;
S=sqrt(wD)*(wD-3*gamma)/sqrt(wD-4*gamma);
omega=0.5*sqrt(wD^2-4*gamma*wD);
X=(kb*T/(m*gamma^2))*(gamma/wD-1+gamma*time+(exp(-wD*time/2)/wD).*(C*cosh(omega*time)+S*sinh(omega*time)));

exponent=-DK^2*X;
isf=exp(exponent);