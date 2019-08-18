% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface between Simulink model and Matlab code for n-D
% trajectories. Shuffles RNG and returns 'data' struct, which contains
% position, momentum and freeze vector for each particle in each population.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = sim_gle_nd(params,n_parall_runs)
% SIM_GLE_ND Simulate the GLE using calculated parameters, in parallel or series.

sim_name = 'sl_pigle_main_current'; % Name of the Simulink model on disk.
seed_limit = 100000; % Upper limit for integer seed values.

% Seed the rng.
rng('shuffle');

% set seed_values, seperate from params struct, to allow compilation(?)
for j=1:n_parall_runs
for i=1:length(params.prtcl)
    seed_values_master(j,i).val = ...
        randi(seed_limit, params.prtcl(i).momenta_dimension,params.prtcl(i).Nprtcl,2+params.z_enabled);
    theta_seed_values_master(j,i).val = ...
        randi(seed_limit, params.prtcl(i).momenta_dimension_theta,params.prtcl(i).Nprtcl,1);
end
end

%% for running interactively

if 0 % for debug
    % population 1
    seed_values_1 = seed_values_master(1,1).val;
    theta_seed_values_1 = theta_seed_values_master(1,1).val;
    assignin('base','seed_values_1',seed_values_1)
    assignin('base','theta_seed_values_1',theta_seed_values_1)
    save('pigle_data.mat', 'params','seed_values_1','theta_seed_values_1','-append');
    clear seed_values_1 theta_seed_values_1;
    
    % population 2
    seed_values_2 = seed_values_master(1,2).val;
    theta_seed_values_2 = theta_seed_values_master(1,2).val;
    assignin('base','seed_values_2',seed_values_2)
    assignin('base','theta_seed_values_2',theta_seed_values_2)
    save('pigle_data.mat', 'params','seed_values_2','theta_seed_values_2','-append');
    clear seed_values_2 theta_seed_values_2;
    
    % population 3
    seed_values_3 = seed_values_master(1,3).val;
    theta_seed_values_3 = theta_seed_values_master(1,3).val;
    assignin('base','seed_values_3',seed_values_3)
    assignin('base','theta_seed_values_3',theta_seed_values_3)
    save('pigle_data.mat', 'params','seed_values_3','theta_seed_values_3','-append');
    clear seed_values_3 theta_seed_values_3;
    
    assignin('base','params',params);

end
%%

pause(1); clear params; load('./pigle_data.mat','params'); % workaround timing issue

simIn(1:n_parall_runs) = Simulink.SimulationInput(sim_name);
for j=1:n_parall_runs
    for i=1:length(params.prtcl)
        simIn(j) = setVariable(simIn(j),['seed_values_' num2str(i)],seed_values_master(j,i).val);
        simIn(j) = setVariable(simIn(j),['theta_seed_values_' num2str(i)],theta_seed_values_master(j,i).val);
        simIn(j) = setModelParameter(simIn(j),'SimulationMode', 'rapid');
    end
end

if length(simIn)>1
    simOut = parsim(simIn, 'ShowProgress', 'on');
else
    simOut = sim(simIn);
end


%data.t = simOut.;

steps2thermalized = params.thermalizing_time/params.isf_sample_time;

error_indx = [];
for i=1:n_parall_runs
    
    if ~isempty(simOut(i).ErrorMessage)
        error_indx = [error_indx i];
        continue;
    end
        
    for j=1:length(params.prtcl)
        
        data(i).prtcl(j).r_supercell=simOut(i).(['pos_supercell_' num2str(j)]);
        data(i).prtcl(j).r=simOut(i).(['pos' num2str(j)]);
        data(i).prtcl(j).p=simOut(i).(['p' num2str(j)]);
        
        data(i).prtcl(j).r(:,:,1:steps2thermalized)=[];
        data(i).prtcl(j).r_supercell(:,:,1:steps2thermalized)=[];
        data(i).prtcl(j).p(:,:,1:steps2thermalized)=[];
        
        if params.z_enabled
            data(i).prtcl(j).freeze=squeeze(simOut(i).(['freeze_' num2str(j)]));
            freeze = zeros(1,params.prtcl(j).Nprtcl);
            for k=1:params.prtcl(j).Nprtcl
                tmp = find(data(i).prtcl(j).freeze(k,:),1,'first');
                if ~isempty(tmp), freeze(k) = tmp; end
            end
            data(i).prtcl(j).freeze = freeze;
        end
    end
end
            
if ~isempty(error_indx), disp(['##### Errors at sim no. :' num2str(error_indx) '   #####']); end
data(error_indx) = [];

% data.t(1:steps2thermalized)=[];
% data.t=data.t-data.t(1);

% % Re-seed RNG and apply initial conditions for y-direction
% params.seed_values = randi(seed_limit, params.momenta_dimension, 1);
% params.r_init = r_init(2);
% params.p_init = p_init(2, :);
% 
% % Simulate the gle for y-direction.
% assignin('base', 'params', params);
% sim(sim_name);
% r(:, 2) = position.Data;
% p(:, 2) = momentum.Data;

end
