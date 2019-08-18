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

% Currently, spatial dependent friction is supported only for single
% particle in each population. That's because the GLE subsystem in the sim model 
% isn't designed to accept friction which is particle specific - only population specific.
% Since there's not much sense in having multiple species with single particle in each, its banned as well.
% Furthermore, spatial dependent friction isn't supported with GLE at the
% moment, so only friction matrix of dimension 1x1 is supported. This is
% because the interpulation doesn't support extracting a sub-array (see the
% sim model).
if ~isempty(find([A_strct(:).A_case] > 1,1)) || ~isempty(find([A_theta_strct(:).A_case] > 1,1))
    if sum([params.prtcl.A_spatial_depended_friction])+sum([params.prtcl.A_spatial_depended_theta_friction]) > 0
        error('Spatial dependent friction is supported ONLY for non-GLE friction');
    end
end

if params.Nprtcl > 1
    if sum([params.prtcl.A_spatial_depended_friction])+sum([params.prtcl.A_spatial_depended_theta_friction]) > 0
        warning('Spatial dependent friction is supported ONLY for non-GLE friction');
    end
end

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
    warning('out_cutoff_r > params.supercell.celldim. Reduce out_cutoff_r or number_density, or, increase Nprtcl.')
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
    params.prtcl(i).A = calc_A(A_strct(i),params.prtcl(i).friction.scaleMat,params.prtcl(i).A_spatial_depended_friction);
    params.prtcl(i).A_theta_strct = A_theta_strct(i);
    params.prtcl(i).A_theta = calc_A(A_theta_strct(i),params.prtcl(i).friction.theta_scaleMat,params.prtcl(i).A_spatial_depended_theta_friction);
    params.prtcl(i).momenta_dimension = size(params.prtcl(i).A,length(size(params.prtcl(i).friction.scaleMat))*params.prtcl(i).A_spatial_depended_friction+1);
    params.prtcl(i).momenta_dimension_theta = size(params.prtcl(i).A_theta,length(size(params.prtcl(i).friction.theta_scaleMat))*params.prtcl(i).A_spatial_depended_theta_friction+1);
    params.prtcl(i).B = calculate_B(params.prtcl(i).A, params.prtcl(i).mass, params.k_B, params.T,params.prtcl(i).momenta_dimension>1,params.prtcl(i).A_spatial_depended_friction);
    params.prtcl(i).B_theta = calculate_B(params.prtcl(i).A_theta, params.prtcl(i).angular_mass, params.k_B, params.T,params.prtcl(i).momenta_dimension_theta>1,params.prtcl(i).A_spatial_depended_theta_friction);
end

%% Set defaults for data output (output everything).
params.last_position = inf;
params.last_momenta = inf;
end

function A = calc_A(A_strct,scaleMat,spatial_depended_friction_enabled)
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

if spatial_depended_friction_enabled == 0
    return;
end

% Now scale 'A' to get spatial dependent friction, using a PES-like array
% The following code assumes that scaleMat (the PES-like scaling array) can
% have different dimensions - that needs both better definision and also
% better coding!
Amat = zeros([size(scaleMat),size(A)]);
for i=1:size(A,1)
    for j=1:size(A,2)
        switch length(size(scaleMat))
            case 1 % no spatial depended friction
                Amat(i,j) = A(i,j);
            case 2
                Amat(:,:,i,j) = A(i,j) * scaleMat;
            case 3
                Amat(:,:,:,i,j) = A(i,j) * scaleMat;
            case 4 % up to 4 dimensions at the moment
                Amat(:,:,:,:,i,j) = A(i,j) * scaleMat;
        end
    end
end

A = Amat;

end

function B = calculate_B(Amat, m, k_B, T, isA2D, spatial_depended_friction_enabled)
% Calculate the B matrix coefficients.

% If spatial-dependent-friction is not enabled, do simple calc of B
if spatial_depended_friction_enabled == 0
    A = Amat;
    B = real(sqrtm(m * k_B * T * (A + A.')));
    return;
end


% Amat is expected to have dimensions as follows:
% dimensions of the PES + {0,2}, depends if the friction is GLE or not.
% So if the PES is 4D, Amat=Amat(x,y,z,theta) or Amat=Amat(x,y,z,theta,:,:)
% The B matrix needs to be calculated seperatly for each possible position
% in the unit-cell.

size_Amat = size(Amat);
Bmat = zeros(size_Amat);

nD_Amat = length(size_Amat);

% Span all possible positions at the unitcell
for i=1:(nD_Amat - isA2D*2)
    vecs4perm{i} = [1:size_Amat(i)];
end
paramIndMat = prepFuncs.allPerm(vecs4perm{:});


for i=1:size(paramIndMat,2)
    if ~isA2D
        A = Amat(paramIndMat{:,i});
        B = real(sqrtm(m * k_B * T * (A + A.')));
        Bmat(paramIndMat{:,i}) = B;
    else
        A = squeeze(Amat(paramIndMat{:,i},:,:));
        warning('off');
        B = real(sqrtm(m * k_B * T * (A + A.')));
        if B^2 ~= (m * k_B * T * (A + A.'))
            error('Calc B matrix faild')
        end
        warning('on');
        Bmat(paramIndMat{:,i},:,:) = B;
    end
end

B = Bmat;

end

function prmtvCellsInSprCell = PrmtvCellsInSprCell(superCellDim, params)
prmtvCellsInSprCell = superCellDim(1)*params.unitcell.numOfPrmtvCells(1) ...
    *superCellDim(2)*params.unitcell.numOfPrmtvCells(2);
end