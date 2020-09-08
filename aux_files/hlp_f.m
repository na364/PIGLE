% Copyright (c) 2018, Nadav Avidor.
% Copyright (c) 2020, Daniel Cropper.
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
        
        function r1 = calc_new_r(r,r_conf,z_enabled,theta_enabled,tilt_enabled,conf3D2D)
            % CALC_NEW_R spans center-of-mass trajectories into trajectories of a
            % configuration of scattering centers
            
            spatial_dim = 2+double(z_enabled>0);
            Nprtcl = size(r,2);
            Natoms = size(r_conf,1);
            tSteps = size(r,3);
            r1 = zeros(spatial_dim,Nprtcl*Natoms,tSteps);
            for i=1:tSteps
                for j=1:Nprtcl
                    if z_enabled
                    r_CoM = r(1:3,j,i);
                    elseif conf3D2D
                      r_CoM=r(1:2,j,i);
                      r_CoM(3)=0;
                    else
                         r_CoM=r(1:2,j,i);
                    end
                    r_CoM1 = repmat(r_CoM,1,Natoms);
                    if theta_enabled
                        tmp=3+z_enabled;
                        rotAngle = r(tmp,j,i);
                    else
                        rotAngle = 0;
                    end
                    if tilt_enabled
                        tmp=4+z_enabled;
                        tiltAngle=r(tmp,j,i);
                    else
                        tiltAngle=0;
                    end
                    % For now, assume z==0, TODO rotate 3D according to
                    % given axis
                    if tilt_enabled
                        r_tmp=[hlp_f.Rmat3(tiltAngle,rotAngle)*r_conf(:,1:3)'];
                        if ~conf3D2D
                            r_tmp=r_tmp(1:2,:);
                        end
                    else
                    r_tmp = [hlp_f.Rmat(rotAngle)*r_conf(:,1:2)'];
                    end
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
            if params.theta_enabled, tmp = [tmp 3+params.z_enabled]; end
            if params.thetatilt_enabled, tmp=[tmp 4+params.z_enabled]; end
            
            supercell_dim = repmat(params.supercell.celldim(tmp)',1,params.prtcl(prtcl_num).Nprtcl,size(r,3));            
            r_supercell = mod(r(tmp,:,:),supercell_dim);
            
            if params.z_enabled
                if params.theta_enabled && ~params.thetatilt_enabled
                    r_supercell(4,:,:) = r_supercell(3,:,:); 
                    r_supercell(3,:,:) = r(3,:,:);
                elseif params.theta_enabled && params.thetatilt_enabled
                    
                    r_supercell(5,:,:)=r_supercell(4,:,:);
                    r_supercell(4,:,:) = r_supercell(3,:,:);
                    r_supercell(3,:,:) = r(3,:,:);
                else
                    r_supercell(3,:,:) = r(3,:,:);
                end
            end
            
%             
           
        end
        
        function m = Rmat(x)
            m = [cos(x) sin(x); -sin(x) cos(x)];
        end
        function m=Rmat3(x,y)
            m1 =[cos(x) sin(x) 0; -sin(x) cos(x) 0; 0 0 1];
            m2 =[1 0 0; 0 cos(y) -sin(y); 0 sin(y) cos(y)];
            m=m2*m1;
        end
    end
end

