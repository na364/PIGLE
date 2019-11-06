% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

global data_file
global pigle_path

% CONFIG JOB is a script to create a collection of folders and scripts
% which then used to run PIGLE over different sets of parameters. In
% otherwords, to probe a multi-variable space.
% Define the jobname, parameters and values in supp_surface_params.m
% Then run this script (config_job), to create a folder <job_name>, and
% under it, folders <job_name>_<i> where <i> runs over all the permutations in the variable space.
% Each <job_name>_<i> folder will have a few files. For example,
% supp_surface_params.m, which holds the specific values of the
% permutation. The supp_surface_params.m file is loaded by surface_params.m
% during the execution of run_PIGLE.m
% <job_name>/<job_name>.mat is created to hold information about the job.
% Lastly, 'run_job' file is created to contain the execution lines for all
% the permutations. Following an execution of config_job.m, execute run_job.
%
% getting started:
% * Edit config_job_params (under the user interface folder, <PIGLE HOME>/UI/ )
% * run config_job - matlab will create <job_name> folder under the PIGLE
%                    root folder.
% * run 'run_job', either from MATLAB (if on a local computer), or
%               externaly (if distributed computing is enabled and set at
%               config_job_params).

run('../UI/config_job_params.m')

tmp = mfilename('fullpath');
[pigle_shell_path,~,~] = fileparts(tmp);

cd(pigle_shell_path); cd('..')
prep_environment
cd(pigle_shell_path);

create_params_file('pigle_shell_params.m',{'job_name','proj_name'},{job_name,proj_name})

% We are going to save everything under folder 'job_name'. That includes
% subfolders for every run, config_file named job_name.mat, which contains
% all the info, and (possibly) more.

if ~exist([pigle_path '/' job_name],'dir'), mkdir([pigle_path '/' job_name]); end
job_path = [pigle_path '/' job_name];
cd(job_path);

% prepare permutations for core-jobs
permutations = prepFuncs.allPerm(sweepVal{:});

Njobs = size(permutations,2);

% If <job_name> folder exists:
% Make sure it addresses the same parameter space (or throw an error).
% Find the serial number for the first 'new' (current) sub_job
% (<job_name>_i), assuming the aim is to add sub jobs (rather replacing sub
% jobs)
if exist([job_name '.mat'],'file')
    load([job_name '.mat'])
    for i=1:length(sweepParams)
        res(i)=strcmp(job_strct(1).sweepParams{i},sweepParams{i});
    end
    if sum(~res) > 0 || i < length(job_strct(1).sweepParams)
        error('A job with the same name exist already, but sweepParams are different')
    end
    start_i = length(job_strct)+1;
else
    start_i = 1;
end
    
% For each of the sub-jobs, create a sub-folder, and files with the
% new/overriding values and parameters
for i=start_i:(start_i-1)+Njobs

    job_strct(i).d_name = [job_name '_' num2str(i)];
    mkdir(job_strct(i).d_name);
    mkdir([job_strct(i).d_name '/tmp_cluster_data']);
    
    job_strct(i).perm = permutations(:,i-start_i+1);
    job_strct(i).sweepParams = sweepParams;
    
    create_params_file([job_strct(i).d_name '/supp_surface_params.m'],job_strct(i).sweepParams,job_strct(i).perm);
    create_md_serialnum_file(job_strct(i).d_name, proj_name);
    job_strct(i).data_file = data_file;
end

% save the job structure
save(job_name,'job_strct')

% create run_job file
if exist('run_job.m','file'), delete run_job.m; end
[ret, name] = system('hostname');
name = strtrim(name);

cd(pigle_shell_path)
fid  = fopen('run_job.m','w');
for i=start_i:(start_i-1)+Njobs
    str1 = ['cd ' job_path '/' job_strct(i).d_name];
    if distributed_computing
        str2 = dist_comp_cmd;
    else
        str2 = ['sub_job_path=pwd; run(''../../sweepParams/pigle_run_single_task.m'')'];
    end
    %str3 = ['sleep ' num2str(5*i)];
    fprintf(fid,'%s\n%s\n',str1,str2);
end
if distributed_computing, fprintf(fid,'%s\n','wait'); end
fclose(fid);

function create_surface_params_file(dirname,sweepParams,perm)
f_name = [dirname '/supp_surface_params.m'];
fid  = fopen(f_name,'w');
for i=1:length(sweepParams)
    fprintf(fid,'%s\n',[sweepParams{i} ' = [' num2str(perm{i}) '];']);
end
fclose(fid);
end

function create_md_serialnum_file(dirname, proj_name)
global data_file

prepFuncs.get_data_filename(proj_name);
f_name = [dirname '/md_data_file.mat'];
data_file_path = data_file;
save(f_name,'data_file_path')
end

function create_params_file(f_name,args,vals)
fid  = fopen(f_name,'w');

if ~iscell(vals)
    vals = {vals};
end

for i=1:length(args)
    if isnumeric(vals{i})
        fprintf(fid,'%s\n',[args{i} ' = [' num2str(vals{i}) '];']);
    elseif ischar(vals{i})
        fprintf(fid,'%s\n',[args{i} ' = ''' vals{i} ''';']);
    else
        error('arg value is not recognized')
    end
end
fclose(fid);
end
