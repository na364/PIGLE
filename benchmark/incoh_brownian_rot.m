% Copyright (c) 2018, Peter Townsend.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function isf=incoh_brownian_rot(r,I,kb,T,g,N,DK,t)

% function to return incoherent isf for a molecule of moment of inertia I,
% radius r as for incoh_langevin_rot but assuming perfect brownian motion,
% ie. no inertial effects

t=(t(:))'; % row vec
n=-N:N; % small but irrelevant loss of efficiency using double sided sum
n=n(:); % column
[n,t]=meshgrid(n,t); 
n=n'; t=t';% still don't understand meshgrid indexing
isf=besselj(n,DK*r).*besselj(n,DK*r).*exp(-kb*T*n.^2.*t./(I*g));
isf=sum(isf,1); % leave us with a row