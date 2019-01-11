% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

global data_file
%global data_path
%global sub_job_path
%global pigle_path
global proj_name
global isGraphicOn
%global data_path_script_name

if ~exist('proj_name','var') || ...
   (exist('proj_name','var') && isempty(proj_name))
        prep_environment;
end

%% Wrapper parameters
pigle_wrapper_params

if ~isSave
    if  ~isGraphicOn
        % there's no obvious motivation for using that mode
        isSave = 1;
        disp('Turning SAVE MODE ON')
    else
        isSave = str2num(timeinput(10,0,'Turn on Save Mode?'));
    end
end

if isGraphicOn
    reduceData = str2num(timeinput(10,reduceData,'Reduce data saved? 0 - dont. 1 - remove p and r_supercell. 2 - leave one population. 3 - clear the data structure'));
    if exist('params','var')
        clearParams = str2num(timeinput(10,1,'Clear variable param and reconfigure?'));
    end
end

%% Prepare save file if needed

if isSave
    
    if isGraphicOn
        ISF2save = str2num(timeinput(20,num2str(ISF2save),'Which ISF to save: 1-Inc. 2-C. 3-single.'));
    end
    
    prepFuncs.get_data_filename(proj_name)    
    disp(['data_file was set to ' data_file])
end

%% Config the model (create params variable)

if (exist('clearParams','var') && clearParams) || ~exist('clearParams','var')
    clear params clearParams
    if exist('./pigle_data.mat','file')
        delete 'pigle_data.mat'
    end
    config_model
end

%% Create the model
create_model

%% Run the model, potentially calculate ISF
if isISF
    % Calcualate the average ISF of N runs trajectories
    tic;
    [data, isf_c_CoM, isf_inc_CoM, isf_1_CoM, isf_c, isf_inc, isf_1, ...
     max_isf_c_CoM, max_isf_inc_CoM, max_isf_1_CoM, max_isf_c, max_isf_inc, max_isf_1] ...
     = calculate_average_isf(params.N_runs, params, params.dK);
    toc;
    
    params.t_isf = params.t(1:size(isf_c_CoM,2));    
else
    % run 1 single trajectory
    tic;
    [data] = sim_gle_nd(params,1);
    toc;
end

%% Post simulation tasks (check thermalization, reduce data)

% Calc Kinetic energy (in K)
for i=1:length(data.prtcl)
    
    [data.prtcl(i).KE, data.prtcl(i).Rot_KE] = ...
        calc_kinetic_energy(data.prtcl(i),params.prtcl(i),params.k_B,params.z_enabled,params.theta_enabled);
    
    disp(['mass =' num2str(params.prtcl(i).mass) ' ... T  = ' ...
        num2str(data.prtcl(i).KE) ' , ' num2str(data.prtcl(i).Rot_KE)])
end

% reduce data
if reduceData > 0
    data.prtcl = rmfield(data.prtcl,'r_supercell');
    data.prtcl = rmfield(data.prtcl,'p');
    if reduceData > 1
        for i=1:length(data.prtcl)
            data.prtcl(i).r = data.prtcl(i).r(:,1,:);
        end
    end
    if reduceData > 2
        data.prtcl = rmfield(data.prtcl,'r');
    end
end

%% Post computation tasks: Save, Plot, local function

if isSave
    if ~exist(data_file,'file')
        save(data_file,'params','data');    
    else
        save(data_file,'params','data','-append');
    end
    
    if isISF
        if sum(ISF2save == 1)
            save(data_file,'isf_inc_CoM','isf_inc','max_isf_inc_CoM','max_isf_inc','-append');
        end
        
        if sum(ISF2save == 2)
            save(data_file,'isf_c_CoM','isf_c','max_isf_c_CoM','max_isf_c','-append');
        end
        
        if sum(ISF2save == 3)
            save(data_file,'isf_1_CoM','isf_1','max_isf_1_CoM','max_isf_1','-append');
        end
    end
end

if toPlot && isISF
    % Plot the data:
    figure; for i=1:size(dK,2), plot(params.t_isf',real(isf_c(i,:,1))); hold on; end
    xlabel('t / ps'); ylabel('Normalised ISF'); title('Coherent ISF, 1st azimuth')
    figure; for i=1:size(dK,2), plot(params.t_isf',real(isf_inc(i,:,1))); hold on; end
    xlabel('t / ps'); ylabel('Normalised ISF'); title('Incoherent ISF, 1st azimuth')
end


