% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

classdef hlp_f
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods (Static)
        
        function r1 = calc_new_r(r,r_conf,z_enabled,theta_enabled)
            % CALC_NEW_R spans center-of-mass trajectories into trajectories of a
            % configuration of scattering centers
            
            spatial_dim = 2+double(z_enabled>0);
            Nprtcl = size(r,2);
            Natoms = size(r_conf,1);
            tSteps = size(r,3);
            r1 = zeros(spatial_dim,Nprtcl*Natoms,tSteps);
            for i=1:tSteps
                for j=1:Nprtcl
                    r_CoM = r(1:spatial_dim,j,i); r_CoM1 = repmat(r_CoM,1,Natoms);
                    if theta_enabled
                        rotAngle = r(spatial_dim+1,j,i);
                    else
                        rotAngle = 0;
                    end
                    % For now, assume z==0, TODO rotate 3D according to
                    % given axis
                    r_tmp = [hlp_f.Rmat(rotAngle)*r_conf(:,1:2)'];
                    if z_enabled
                        r_conf_new = [r_tmp; zeros(1,size(r_tmp,2))];
                    else
                        r_conf_new = r_tmp;
                    end
                    r1(:,(j-1)*Natoms+1:(j)*Natoms,i) = r_conf_new + r_CoM1;
                end
            end
        end

        function r_supercell = calc_r_supercell(params,prtcl_num,r)
            % CALC_R_SUPERCELL works out the particle position in the
            % supercell for a given 'absolute' trajectory
            
            tmp = 1:2;            
            if params.theta_enabled, tmp = [tmp params.model_dim]; end
            
            supercell_dim = repmat(params.supercell.celldim(tmp)',1,params.prtcl(prtcl_num).Nprtcl,size(r,3));            
            r_supercell = mod(r(tmp,:,:),supercell_dim);
            
            if params.z_enabled
                if params.theta_enabled
                    r_supercell(4,:,:) = r_supercell(3,:,:);                    
                else
                    r_supercell(3,:,:) = r(3,:,:);
                end
            end
        end
        
        function m = Rmat(x)
            m = [cos(x) sin(x); -sin(x) cos(x)];
        end
    end
end

