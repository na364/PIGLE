% Copyright (c) 2018, Peter Townsend
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function isf=incoh_langevin_rot(r,I,kb,T,g,N,DK,t)

% function to return incoherent isf for a molecule of moment of inertia I,
% radius r

% the angular co-ordinate is assumed to obey the LE 
% I*theta'' = -friction*theta' + G(t)
% where G(t) is a random torque obeying FDT
% in analogy with translational LE and FDT

% function accepts a single DK for simplicity, but takes a vector of time
% inputs

% g = 'gamma', friction in units of per-ps

t=(t(:))'; % row vec
n=-N:N; % small but irrelevant loss of efficiency using double sided sum
n=n(:); % column
[n,t]=meshgrid(n,t); 
n=n'; t=t';% still don't understand meshgrid indexing
isf=besselj(n,DK*r).*besselj(n,DK*r).*exp(-kb*T*n.^2.*(exp(-g*t)+g*t-1)/(I*g^2));
isf=sum(isf,1); % leave us with a row


