% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function [ f ] = f_interaction( x,f_func,f_func_params )
%F_INTERACTION computes the inter-particle force
% The force dimensions are expected to be in meV/Angstrom
%   Inputs:
%       x             - a vector. Points along the distance between adsorbed particles at which the force is to be calculated
%       f_func        - Case number. The type of force to consider.
%       f_func_params - Parameters for the force equation.
%
%   Output:
%       f             - the inter-adsorbate force, f(x)
%

switch f_func
    case 1
        % Simple sum of -A_i*(x.^(-pwr_i)
        f = zeros(size(x));
        for i=1:2:length(f_func_params)
            A = f_func_params(i);
            if A < 0, warning('Prefactor of inter-adsorbate force is negative. When multiplied by (dx,dy)/|r| will result attraction'); end
            pwr = f_func_params(i+1);
            f = f + A*(x.^(-pwr));
        end
        
    case 2
        %No interactions:
        f = 0*x;
        
    otherwise
        error('No such force function defind');
        
end


end

