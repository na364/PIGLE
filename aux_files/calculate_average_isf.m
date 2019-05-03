% Copyright (c) 2018, Nadav Avidor.
% Copyright (c) 2016, Ryan Collingham.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate ISFs for multiple runs (same parameters, different random seeds).
% Return the average ISF. The purpose of this is to remove noisy parts.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data_last, isf_c_CoM,     isf_inc_CoM,     isf_1_CoM,     isf_c,     isf_inc,     isf_1, ...
                 max_isf_c_CoM, max_isf_inc_CoM, max_isf_1_CoM, max_isf_c, max_isf_inc, max_isf_1] = calculate_average_isf(N, params, dK)

    % Reshape dK to matrix of (:,2)
    dK_reshaped = reshape(permute(dK,[1 3 2]),[],2);    

    % Run N simulations and sum the FFTs of the scattering functions.
    [data_last, isf_c_CoM ,isf_inc_CoM ,isf_1_CoM , isf_c,isf_inc,isf_1] = sum_isf(N, params, dK_reshaped);
    
        % Reshape isf to matrix with dimention (dK,length(isf),num_of_azimuth)
    if size(dK,3) > 1
        isf_c_CoM = permute(reshape(isf_c_CoM,[],2,size(isf_c_CoM,2)),[1 3 2]);
        isf_inc_CoM = permute(reshape(isf_inc_CoM,[],2,size(isf_inc_CoM,2)),[1 3 2]);
        isf_1_CoM = permute(reshape(isf_1_CoM,[],2,size(isf_1_CoM,2)),[1 3 2]);
        isf_c = permute(reshape(isf_c,[],2,size(isf_c,2)),[1 3 2]);
        isf_inc = permute(reshape(isf_inc,[],2,size(isf_inc,2)),[1 3 2]);
        isf_1 = permute(reshape(isf_1,[],2,size(isf_1,2)),[1 3 2]);
    end
    
    % Normalise the ISF to unity.
    [isf_c_CoM, max_isf_c_CoM] = normalise(isf_c_CoM);
    [isf_inc_CoM, max_isf_inc_CoM] = normalise(isf_inc_CoM);
    [isf_1_CoM, max_isf_1_CoM] = normalise(isf_1_CoM);
    [isf_c, max_isf_c] = normalise(isf_c);
    [isf_inc, max_isf_inc] = normalise(isf_inc);
    [isf_1, max_isf_1] = normalise(isf_1);
    
end

function [data, isf_c_CoM ,isf_inc_CoM ,isf_1_CoM , isf_c,isf_inc,isf_1] = sum_isf(N, params, dK)
    % Sum the FFTs of the scattering functions of N simulations.

    z_enabled = params.z_enabled;
    dKz_include_in_isf = params.dKz_include_in_isf;
    theta_enabled = params.theta_enabled;
    
    sum_Skw_c_CoM = 0;
    sum_Skw_inc_CoM = 0;
    sum_Skw_1_CoM = 0;
    sum_Skw_c = 0;
    sum_Skw_inc = 0;
    sum_Skw_1 = 0;

    % Run simulations in parallel to reduce execution time.
    N_itter = ceil(N/params.n_parall_runs);
    n_parall_runs = params.n_parall_runs;
    numOfActualSim = 0;
    for i = 1:N_itter
        disp(['Iter ' num2str(i)]); toc
        
        % Run simulation to get position data.
        if i*params.n_parall_runs > N
            n_parall_runs = N-(i-1)*n_parall_runs;
        end
        data = sim_position(params,n_parall_runs);
        numOfActualSim = numOfActualSim + length(data);
        
        % Attach conf to 'data' struct, which is used in the parfor.
        for j=1:length(data)
            for k=1:length(data(j).prtcl)
                data(j).prtcl(k).conf = params.prtcl(k).conf;
            end
        end
        
        % Calculate the scattered amplitudes for our value of dK.
        if length(data) > 1
            parfor j=1:length(data)
                [Skw_c_CoM,Skw_inc_CoM,Skw_1_CoM,Skw_c,Skw_inc,Skw_1] = scattered_amplitudes(dK, data(j),z_enabled,dKz_include_in_isf,theta_enabled);

                % Add to the running total Skw.
                sum_Skw_c_CoM = sum_Skw_c_CoM + Skw_c_CoM;
                sum_Skw_inc_CoM = sum_Skw_inc_CoM + Skw_inc_CoM;
                sum_Skw_1_CoM = sum_Skw_1_CoM + Skw_1_CoM;
                sum_Skw_c = sum_Skw_c + Skw_c;
                sum_Skw_inc = sum_Skw_inc + Skw_inc;
                sum_Skw_1 = sum_Skw_1 + Skw_1;
            end
        else % reduce the overhead of parfor for one itteration
            [Skw_c_CoM,Skw_inc_CoM,Skw_1_CoM,Skw_c,Skw_inc,Skw_1] = scattered_amplitudes(dK, data, z_enabled,dKz_include_in_isf,theta_enabled);
            
            % Add to the running total Skw.
            sum_Skw_c_CoM = sum_Skw_c_CoM + Skw_c_CoM;
            sum_Skw_inc_CoM = sum_Skw_inc_CoM + Skw_inc_CoM;
            sum_Skw_1_CoM = sum_Skw_1_CoM + Skw_1_CoM;
            sum_Skw_c = sum_Skw_c + Skw_c;
            sum_Skw_inc = sum_Skw_inc + Skw_inc;
            sum_Skw_1 = sum_Skw_1 + Skw_1;
        end
    end
    
    sum_Skw_c_CoM = sum_Skw_c_CoM / numOfActualSim;
    sum_Skw_inc_CoM = sum_Skw_inc_CoM / numOfActualSim;
    sum_Skw_1_CoM = sum_Skw_1_CoM / numOfActualSim;
    sum_Skw_c = sum_Skw_c / numOfActualSim;
    sum_Skw_inc = sum_Skw_inc / numOfActualSim;
    sum_Skw_1 = sum_Skw_1 / numOfActualSim;
    
    if length(data)>1, data(2:end)=[]; end
    
    % Work out the ISF
    isf_c_CoM =   real(ifft(sum_Skw_c_CoM,[],2));
    isf_inc_CoM = real(ifft(sum_Skw_inc_CoM,[],2));
    isf_1_CoM =   real(ifft(sum_Skw_1_CoM,[],2));
    isf_c =   real(ifft(sum_Skw_c,[],2));
    isf_inc = real(ifft(sum_Skw_inc,[],2));
    isf_1 =   real(ifft(sum_Skw_1,[],2));
    
    
    % Reduce the ISF to positive SE time, to save space
    indx = length(params.t)/2+1;
    
    isf_c_CoM = isf_c_CoM(:,1:indx,:);    
    isf_inc_CoM = isf_inc_CoM(:,1:indx,:);
    isf_1_CoM = isf_1_CoM(:,1:indx,:);

    if exist('isf_c','var') && ~isempty(isf_c)
        isf_c = isf_c(:,1:indx,:);
    end
    
    if exist('isf_inc','var') && ~isempty(isf_inc)
        isf_inc = isf_inc(:,1:indx,:);
    end
    
    if exist('isf_1','var') && ~isempty(isf_1)
        isf_1 = isf_1(:,1:indx,:);
    end
    
end

function data = sim_position(params,n_parall_runs)
% Simulate using params, return position data and the last position and
% momenta.
data = sim_gle_nd(params,n_parall_runs);
%data.prtcl = rmfield(data.prtcl,{'r_supercell','p'});

end

function [Skw_c_CoM,Skw_inc_CoM,Skw_1_CoM,Skw_c,Skw_inc,Skw_1] = scattered_amplitudes(delta_k, data,z_enabled,dKz_include_in_isf,theta_enabled)
% Calculate the scattering amplitude for a given delta k vector

%% Center of Mass
A_c_CoM = zeros(size(delta_k,1),size(data.prtcl(1).r,3));
Skw_inc_CoM = A_c_CoM;

if z_enabled > 0 && dKz_include_in_isf
    if size(delta_k,2) < 3, delta_k(:,3) = zeros(size(delta_k,1),1); end
elseif size(delta_k,2) == 3
    error('delta_k is 3D, but z dimension is disabled')
end    

for j=1:length(data.prtcl)
    r = data.prtcl(j).r;    
    FF_CoM = data.prtcl(j).conf.form_factor.FF_CoM;
    
    % Reshape F_CoM to a vector compatible with the reshaped dK
    FF_CoM = reshape(FF_CoM,[],1);    
    
    for i=1:size(r,2)
        r_i = squeeze(r(1:(2 + (z_enabled > 0 && dKz_include_in_isf)),i,:));
        A_i = FF_CoM .* exp(1i * (delta_k * r_i));
        A_c_CoM = A_c_CoM + A_i;
        Skw_1_CoM = conj(fft(A_i,[],2)).*fft(A_i,[],2);
        Skw_inc_CoM = Skw_inc_CoM + Skw_1_CoM;
    end
end

Skw_c_CoM = conj(fft(A_c_CoM,[],2)).*fft(A_c_CoM,[],2);

% Skw_inc_CoM = real(Skw_inc_CoM);
% Skw_c_CoM = real(Skw_c_CoM);
% Skw_1_CoM = real(Skw_1_CoM);

%% Full configuration

% is there any particle which is not a point particle?
not_a_point=0;
for j=1:length(data.prtcl)
    if size(data.prtcl(j).conf.r_conf,1) > 1
        not_a_point=1;
    end
end

if not_a_point == 0
    Skw_c=[]; Skw_inc=[]; Skw_1=[];
else

    A_c = zeros(size(delta_k,1),size(data.prtcl(1).r,3));
    Skw_inc = A_c;

    for j=1:length(data.prtcl)
        r = data.prtcl(j).r;
        r_conf = data.prtcl(j).conf.r_conf;
        r1 = hlp_f.calc_new_r(r,r_conf,z_enabled,theta_enabled);
        
        % Prep FF for the calculation.
        % TODO check if it works for Natoms>1 etc.
        FF_tmp = data.prtcl(j).conf.form_factor.FF;
        % Reshape F_tmp to a vector compatible with the reshaped dK
        FF_tmp = reshape(FF_tmp,[],size(FF_tmp,2));
        % Duplicate the FF_tmp to match r1 (with length of num of scattering centers * num of
        % CoM trajectories)
        FF = repmat(FF_tmp,1,size(r,2));

        for i=1:size(r1,2)
            r_i = squeeze(r1(1:(2 + (z_enabled > 0 && dKz_include_in_isf)),i,:));
            A_i = FF(:,i) .* exp(1i * (delta_k * r_i));
            A_c = A_c + A_i;
            Skw_1 = conj(fft(A_i,[],2)).*fft(A_i,[],2);
            Skw_inc = Skw_inc + Skw_1;
        end
    end

    Skw_c = conj(fft(A_c,[],2)).*fft(A_c,[],2);

    % Skw_inc = real(Skw_inc);
    % Skw_c = real(Skw_c);
    % Skw_1 = real(Skw_1);
end

end

function [I,maxI] = normalise(I)
% Normalise the ISF to unity
maxI=max(I,[],2);
I=I./maxI;
end

