% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.


function [KE, Rot_KE] = calc_kinetic_energy(data_prtcl,params_prtcl,k_B,z_enabled,theta_enabled,time_steps)
% CALC_KINETIC_ENERGY both translational and rotational
% A function to calculate both rotational and translational KE, taking into
% account whetehr translation includes 2 or three dimensions.
%
%   Inputs:
%       data_prtcl    - Simulation result struct for one population of particles
%       params_prtcl  - Pre-simulation configuration structure for one population of particles
%       k_B           - Boltzman constant
%       z_enabled     - is the third spatial dimension enabled
%       theta_enabled - is the angular dimension enabled
%       time_steps    - time steps to consider. If not defined, the default is 1:size(data_prtcl.p,3)
%
%   Outputs:
%       KE - kinetic energy in kelvin
%       Rot_KE - rotational kinetic energy in kelvin
%
    if ~exist('time_steps','var') || (exist('time_steps','var') && isempty(time_steps))
        time_steps = 1:size(data_prtcl.p,3);
    end
    KE = 0.5*(data_prtcl.p(1,:,time_steps).^2 + ...
              data_prtcl.p(2,:,time_steps).^2)./params_prtcl.mass;
    if z_enabled > 0
          KE = KE + 0.5 * data_prtcl.p(3,:,time_steps).^2./params_prtcl.mass;
    end
    KE = squeeze(KE);
    if z_enabled > 0
        for j=1:length(data_prtcl.freeze)
            if data_prtcl.freeze(j) > 0
                KE(j,data_prtcl.freeze(j):end)=NaN;
            end
        end
    end
    KE = reshape(KE,1,[]);
    KE = KE(~isnan(KE));
    KE = mean(2*KE/k_B)/(2+int16(z_enabled > 0));
    
    if theta_enabled
        Rot_KE = 0.5*(data_prtcl.p(3+int16(z_enabled > 0),:,:).^2)./params_prtcl.angular_mass;
        Rot_KE = squeeze(Rot_KE);
        if z_enabled > 0
            for j=1:length(data_prtcl.freeze)
                if data_prtcl.freeze(j) > 0
                    Rot_KE(j,data_prtcl.freeze(j):end)=NaN;
                end
            end
        end
        Rot_KE = reshape(Rot_KE,1,[]); Rot_KE = Rot_KE(~isnan(Rot_KE));
        Rot_KE = mean(2*Rot_KE/k_B);
    else
        Rot_KE = NaN;
    end

end