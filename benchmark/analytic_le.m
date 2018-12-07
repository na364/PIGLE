% Copyright (c) 2018, Peter Townsend
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function isf=analytic_le(DK,time,gamma,kb,T,m)

exponent=-DK^2*(kb*T/(m*gamma^2))*(exp(-gamma*time)+gamma*time-1);
isf=exp(exponent);