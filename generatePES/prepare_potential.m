% Copyright (c) 2018, Nadav Avidor.
% Copyright (c) 2020, Daniel Cropper.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function PotMatrix = prepare_potential(unitcell, pot_strct)
% PREPARE_POTENTIAL constructs a 5D PES for adsorbate.
% The function creates a 2D PES, then modify it for each theta, then modify
% for the 'z' dimension, assuming a morse potential.
%
% Alternatively a user defined function can be used to generate the potential. 
%PotMatrix should be of the form V(x,y,z,theta,tilt). See line 
%
% Inputs:
%       unitcell - struct with fields that holds the grid: x,y (mandatory),
%                  and z,theta (optional). Also, celldim holds the 'size'
%                  of each dimension.
%       pot_strct - struct with fields:
%                   V(i,:) - holds angular dependend values for the i'th variable
%                            which is passed to the function which creates the
%                            lateral 2D PES.
%                is_potval - defines if V(i,:) are values of the PES, or
%                            other values which needs to be passed to the function
%                            creating the 2D latteral PES.
%                   ref_De - depth of the 'z' potential at r_e
%                        a - controls the width of the well (within the morse potential)
%                      r_e - equilibrium distance for the morse potential
%                     f_2D - pointer to a function for creating 2D lateral PES
%

if ~isfield(unitcell,'theta')
    unitcell.theta=0;
end
if ~isfield(unitcell,'z')
    unitcell.z=0;
end
if ~isfield(unitcell, 'tilt')
    unitcell.tilt=0;
end
n_xyzthetatilt = [length(unitcell.x) length(unitcell.y) length(unitcell.z) length(unitcell.theta) length(unitcell.tilt)];

    
%Generates a PotMatrix by varying the parameters in pot_strct by angle and
%by z as desired



V_morse = @(De,a,r_e,r) De.*(exp(-2*a*(r-r_e))-2*exp(-a*(r-r_e)));
%Vcor_f = @(r) Vcor_/(z(end)-r_e)*(r-r_e)-Vcor_; % control decay in corrugation as function of 'z'
% Vcor = Vcor_f(z); 

angle1 = unitcell.theta;
angle2= unitcell.tilt;
z     = unitcell.z;
n_xyzthetatilt = [length(unitcell.x) length(unitcell.y) length(unitcell.z) length(unitcell.theta) length(unitcell.tilt)];
pot_vals=zeros(length(pot_strct.V),n_xyzthetatilt(3),n_xyzthetatilt(4),n_xyzthetatilt(5));
switch pot_strct.anglecase
    case 0
        %No angular potential, only z variation
        for i=1:length(pot_strct.V)
     v = pot_strct.V(i,:);

     if pot_strct.is_potval(i)
         tmp = V_morse(pot_strct.ref_De(i), pot_strct.a(i), pot_strct.r_e(i), z')+v;
     else
         tmp = repmat(v,length(z),1);
     end
     %tmp = repmat(tmp,1,length(angle)/length(v));
     pot_vals(i,:,:,:) = tmp;
        end
        
        
    case 1
        %Angular variation with only theta
        
        
for i=1:length(pot_strct.V)
     v=zeros(n_xyzthetatilt(3),n_xyzthetatilt(4),n_xyzthetatilt(5));
     for dum1=1:n_xyzthetatilt(3)
         for dum2=1:n_xyzthetatilt(4)
             for dum3 = 1:n_xyzthetatilt(5)
    v(dum1,dum2,dum3) = pot_strct.V(i,:);
             end
         end
     end
     amp=pot_strct.theta_minmax(i,:);
    
     t = [0:2*pi/n_xyzthetatilt(4):2*pi-2*pi/n_xyzthetatilt(4)];
          for dum1=1:n_xyzthetatilt(3)
         for dum2=1:n_xyzthetatilt(4)
             for dum3 = 1:n_xyzthetatilt(5)
     v(dum1,dum2,dum3)=amp(1)*(cos((t(dum2)-amp(2))+pi/2)).^2+v(dum1,dum2,dum3);
             end
         end
          end
     for dum1=1:n_xyzthetatilt(3)
         for dum2=1:n_xyzthetatilt(4)
             for dum3 = 1:n_xyzthetatilt(5)
     if pot_strct.is_potval(i)
         Vzpot=V_morse(pot_strct.ref_De(i), pot_strct.a(i), pot_strct.r_e(i), z');
         v(dum1,dum2,dum3) = Vzpot(dum3)+v(dum1,dum2,dum3);

     end
             end
         end
     end
     %tmp = repmat(tmp,1,length(angle)/length(v));
     pot_vals(i,:,:,:) = v;
end

        case 2
            %Angular variation with theta and tilt
             
        for i=1:length(pot_strct.V)
            v=zeros(n_xyzthetatilt(3),n_xyzthetatilt(4),n_xyzthetatilt(5));
            for dum1=1:n_xyzthetatilt(4)
                for dum2=1:n_xyzthetatilt(5)
                    for dum3=1:n_xyzthetatilt(3)
            v(dum3,dum1,dum2)=pot_strct.V(i);
                    end
                end
            end
            amp_theta=pot_strct.theta_minmax(i,:);
            amp_tilt=pot_strct.theta_minmax(i,:);
            
            
            t=[0:2*pi/n_xyzthetatilt(4):2*pi-2*pi/n_xyzthetatilt(4)];
            for dum1=1:n_xyzthetatilt(4)
                for dum2=1:n_xyzthetatilt(5)
                    for dum3=1:n_xyzthetatilt(3)
            v(dum3,dum1,dum2)=amp_theta(1)*(cos((t(dum1)-amp_theta(2))+pi/2)).^2+v(dum3,dum1,dum2);
                    end
                end
            end
            
            s=[0:2*pi/n_xyzthetatilt(5):2*pi-2*pi/n_xyzthetatilt(5)];
            for dum1=1:n_xyzthetatilt(4)
                for dum2=1:n_xyzthetatilt(5)
                    for dum3=1:n_xyzthetatilt(3)
            v(dum3,dum1,dum2)=amp_tilt(1)*(cos((s(dum2)-amp_tilt(2))+pi/2)).^2+v(dum3,dum1,dum2);
                    end
                end
            end
            
            if pot_strct.is_potval(i)
                morse=V_morse(pot_strct.ref_De(i), pot_strct.a(i), pot_strct.r_e(i), z');
                for dum1=1:n_xyzthetatilt(4)
                for dum2=1:n_xyzthetatilt(5)
                    for dum3=1:n_xyzthetatilt(3)
          v(dum3,dum1,dum2) = morse(dum3)+v(dum3,dum1,dum2);
          tmp=v;
                    end
                end
                end
            else
          tmp =v; %repmat(v,length(z),1);
            end
     %tmp = repmat(tmp,1,length(angle)/length(v));
     pot_vals(i,:,:,:) = tmp;
        end
            
end

    
% f_2D_args_=zeros(length(z),length(angle1), length(angle2));
for i=1:length(z)
    for j=1:length(angle1)
        for k=1:length(angle2)
      tmp = pot_vals(:,i,j,k)';
      f_2D_args_(i,j,k) = {num2cell(tmp)};
        end
    end
end



[PotMatrix] = stack_2Dpot_2_5Dpot(unitcell.celldim, n_xyzthetatilt, pot_strct.f_2D, f_2D_args_);

end

function [PotMatrix] = stack_2Dpot_2_5Dpot(celldim, n_xyzthetatilt, f_2D, f_2D_args_)
% STACK_2DPOT_2_4DPOT
%

nx=n_xyzthetatilt(1); ny=n_xyzthetatilt(2);
nz=n_xyzthetatilt(3); ntheta=n_xyzthetatilt(4);
ntilt=n_xyzthetatilt(5);

for i=1:nz
    for j=1:ntheta
        for k=1:ntilt
        f_2D_args = f_2D_args_{i,j,k};
        V = f_2D(celldim,nx,ny,f_2D_args{:});
        if size(V,1) == nx
            % make V(y,x), so surface(V) gives the correct plot
            V = V';
        end
        end
        PotMatrix(:,:,i,j,k)=real(V);
    end
end

end
