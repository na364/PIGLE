% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

% PREP_ENVIRONMENT
% Add to the path required folders including pigle_path and sufolders, set
% project name, etc.

global pigle_path
global proj_name
global data_path_script_name
global isGraphicOn


%% Define user customaized functions or scripts for the simulation environment

% data_path_script_name is called (if not set to empty string), to provide 
data_path_script_name = '';%'make_data_path.m'; % set to empty if irrelevant

isGraphicOn = usejava('desktop');

%% Define path of PIGLE and Project name

pigle_path = pwd;
if ~exist('run_pigle.m','file')
   if isGraphicOn
      pigle_path = uigetdir('','Select the correct path for pigle_path');
   else
      error('Graphics off, pigle_path is not set')
   end
end

addpath(pigle_path)
addpath([pigle_path '/UI'])
addpath([pigle_path '/aux_files'])
addpath([pigle_path '/generatePES'])
addpath([pigle_path '/generatePES/PES_library'])
addpath([pigle_path '/pigle_sim'])

if ~exist('proj_name','var'), proj_name=''; end

% proj_name must be non-empty for any operation in PIGLE
if isempty(proj_name), proj_name = 'smplProj'; end

% if in graphic mode, ask whether to change proj_name
if isGraphicOn && str2num(timeinput(10,0,['proj_name = "' proj_name '" . Re-define project name?']))
    tmp = inputdlg('Type the name of the project','title1',1,{'smplProj'});
    proj_name = tmp{1};
end

disp(['The project name was set to: #########' num2str(proj_name) '#########'])

