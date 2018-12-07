% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function PotMatrix = prepare_potential(unitcell, pot_strct)
% PREPARE_POTENTIAL constructs a 4D PES for adsorbate.
% The function creates a 2D PES, then modify it for each theta, then modify
% for the 'z' dimension, assuming a morse potential.
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

V_morse = @(De,a,r_e,r) De.*(exp(-2*a*(r-r_e))-2*exp(-a*(r-r_e)));
%Vcor_f = @(r) Vcor_/(z(end)-r_e)*(r-r_e)-Vcor_; % control decay in corrugation as function of 'z'
% Vcor = Vcor_f(z); 

angle = unitcell.theta;
z     = unitcell.z;
n_xyztheta = [length(unitcell.x) length(unitcell.y) length(unitcell.z) length(unitcell.theta)];

for i=1:length(pot_strct.V)
     v = pot_strct.V(i,:);
     
     % if there are 2 values for each 'v', assume that its for rotational PES
     if length(v) == 2
         t = [0:pi/n_xyztheta(4):pi-pi/n_xyztheta(4)];
         v = diff(v)*(cos(t+(v(1)<v(2))*pi/2)).^2+v(1);
     end
     
     if pot_strct.is_potval(i)
         tmp = V_morse(pot_strct.ref_De(i), pot_strct.a(i), pot_strct.r_e(i), z')+v;
     else
         tmp = repmat(v,length(z),1);
     end
     %tmp = repmat(tmp,1,length(angle)/length(v));
     pot_vals(i,:,:) = tmp;
end

for i=1:length(z)
    for j=1:length(angle)
      tmp = pot_vals(:,i,j)';
      f_2D_args_(i,j) = {num2cell(tmp)};
    end
end



[PotMatrix] = stack_2Dpot_2_4Dpot(unitcell.celldim, n_xyztheta, pot_strct.f_2D, f_2D_args_);

end

function [PotMatrix] = stack_2Dpot_2_4Dpot(celldim, n_xyztheta, f_2D, f_2D_args_)
% STACK_2DPOT_2_4DPOT
%

nx=n_xyztheta(1); ny=n_xyztheta(2);
nz=n_xyztheta(3); ntheta=n_xyztheta(4);

for i=1:nz
    for j=1:ntheta
        f_2D_args = f_2D_args_{i,j};
        V = f_2D(celldim,nx,ny,f_2D_args{:});
        if size(V,1) == nx
            % make V(y,x), so surface(V) gives the correct plot
            V = V';
        end
        PotMatrix(:,:,i,j)=real(V);
    end
end

end