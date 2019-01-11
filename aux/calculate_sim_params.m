% Copyright (c) 2018, Nadav Avidor.
% Copyright (c) 2016, Ryan Collingham.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the simulation parameters based on physical inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = calculate_sim_params(params, A_strct, A_theta_strct, r_conf)

%% Calculate supercell properties

% mass and number of particles
Nprtcl_list = ceil(params.Nprtcl .* (params.number_density)/sum(params.number_density));

% Supercell (for interactions) Properties
desired_prmtvCellsInSprCell = max(ceil(Nprtcl_list./params.number_density));
superCellDim = [0 0];
while 1
    if desired_prmtvCellsInSprCell < PrmtvCellsInSprCell(superCellDim+[1 0], params)
        break;
    end
    superCellDim = superCellDim + [1 0];
    if desired_prmtvCellsInSprCell < PrmtvCellsInSprCell(superCellDim+[0 1], params)
        break;
    end
        superCellDim = superCellDim + [0 1];
end
 
params.supercell.celldim = params.unitcell.celldim(1:2) .* superCellDim;
if params.interactions.out_cutoff_r > min(params.supercell.celldim) && params.interactions.active
    error('out_cutoff_r > params.supercell.celldim. Reduce out_cutoff_r or number_density, or, increase Nprtcl')
end

if params.model_dim > 2
    params.supercell.celldim(3) = params.unitcell.celldim(3);
end
if params.theta_enabled
    params.supercell.celldim(params.model_dim) = 2*pi;
end

% Now, after optimizing the supercell size to as close as possible to
% requirements, adjust number of particles to make coverage as accurate as
% possible.
prmtvCellsInSprCell = PrmtvCellsInSprCell(superCellDim, params);

%% Calculate prtcl properties

Nprtcl_list = ceil(params.number_density * prmtvCellsInSprCell);

if params.Nprtcl * 1.3 < sum(Nprtcl_list) && sum(Nprtcl_list) > 6
    error(['Nprtcl which satisfy the number density is calculated to be ' ...
        num2str(sum(Nprtcl_list)) ' (more than 1.3 * requested Nprtcl. Abborting']);
else
    params.Nprtcl = sum(Nprtcl_list);
    disp(['Nprtcl was changed to accommodate the number densities. New Nprtcl = ' num2str(params.Nprtcl)]);
    total_calculated_number_density = params.Nprtcl / ... 
        (params.supercell.celldim(1)*params.supercell.celldim(2) / (params.unitcell.celldim(1)*params.unitcell.celldim(2)) *  (params.unitcell.numOfPrmtvCells(1)*params.unitcell.numOfPrmtvCells(2)));
    disp(['Total number density is: ' num2str(total_calculated_number_density)]);
end

for i=1:length(params.mass_list)    
    params.prtcl(i).mass = params.mass_list(i);
    params.prtcl(i).angular_mass = params.angular_mass_list(i);
    params.prtcl(i).Nprtcl = Nprtcl_list(i);
    params.prtcl(i).A_strct = A_strct(i);
    params.prtcl(i).A = calc_A(A_strct(i));
    params.prtcl(i).A_theta_strct = A_theta_strct(i);
    params.prtcl(i).A_theta = calc_A(A_theta_strct(i));
    params.prtcl(i).momenta_dimension = length(params.prtcl(i).A);
    params.prtcl(i).momenta_dimension_theta = length(params.prtcl(i).A_theta);
    params.prtcl(i).B = calculate_B(params.prtcl(i).A, params.prtcl(i).mass, params.k_B, params.T);
    params.prtcl(i).B_theta = calculate_B(params.prtcl(i).A_theta, params.prtcl(i).angular_mass, params.k_B, params.T);
end

%% Set defaults for data output (output everything).
params.last_position = inf;
params.last_momenta = inf;
end

function A = calc_A(A_strct)
eta = A_strct.eta;
tau = A_strct.tau;
% A coefficient matrix:
switch A_strct.A_case
    case 1 % No filter
        A = eta;
    case 2 % Low pass
        s = sqrt(eta / tau);
        A = [0, -s ; ...
            s, 1/tau];
    case 3 % spike
        dw = A_strct.dw;
        w0 = A_strct.w0;
        % to see the filter:
        % w=0:0.01:50; tau=[6 6]; eta=[2 2]; dw=[0.1 0.1]; w0=[0.28 20];
        % y = zeros(size(w));
        % for i=1:length(eta)
        % y=y+2*sqrt(2/pi)*(eta(i)*dw(i)/tau(i))*(dw(i)^2+w0(i)^2+w.^2)./((dw(i)^2+(w-w0(i)).^2).*(dw(i)^2+(w+w0(i)).^2));
        % end
        % plot(w,y)
        A = generate_A_from_frequencies_multiple_gamma(w0, dw, eta, tau);
    case 4 % band pass
        
        % w=0:0.01:50; tau=[1 6]; eta=[2];
        % y = sqrt(2/pi)*eta*(1./(1+(tau(1)*w).^2)-1./(1+(tau(2)*w).^2))
        % plot(w,y)
        
        sL = sqrt(eta/tau(1));
        sH = sqrt(eta/tau(2));
        A = [0   sH        -sL; ...
            sH   1/tau(2)   0 ; ...
            sL   0          1/tau(1)];
        
    case 5 % bi low-pass
        
        s1 = sqrt(eta(1) / tau(1));
        s2 = sqrt(eta(2) / tau(2));
        A = [0,     -s1   ,   -s2   ; ...
            s1,   1/tau(1),    0    ; ...
            s2,      0    , 1/tau(2)];
        
end

end

function B = calculate_B(A, m, k_B, T)
% Calculate the B matrix coefficients.
B = real(sqrtm(m * k_B * T * (A + A.')));
end

function prmtvCellsInSprCell = PrmtvCellsInSprCell(superCellDim, params)
prmtvCellsInSprCell = superCellDim(1)*params.unitcell.numOfPrmtvCells(1) ...
    *superCellDim(2)*params.unitcell.numOfPrmtvCells(2);
end