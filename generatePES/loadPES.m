% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function [PES] = loadPES(filename,PESfield)
%LOADPES Summary of this function goes here
%   Detailed explanation goes here
tmp = load(filename,PESfield);
PES=tmp.(PESfield);
end

