% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

classdef prepFuncs
    %PREPFUNCS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        
        function unitcell = make_unitcell(x,y,varargin)
            %% parse varargin
            prsdArgs = inputParser;   % Create instance of inputParser class.
            prsdArgs.addParameter('numOfPrmtvCells', [1 1], @isnumeric);
            prsdArgs.addParameter('z', [], @isnumeric);
            prsdArgs.addParameter('theta', [], @isnumeric);
            prsdArgs.parse(varargin{:});

            numOfPrmtvCells = prsdArgs.Results.numOfPrmtvCells;
            z              = prsdArgs.Results.z;
            theta              = prsdArgs.Results.theta;
            
            if size(numOfPrmtvCells) ~= [1 2], error('numOfPrmtvCells is not of size [1 2]'); end

            unitcell.numOfPrmtvCells = numOfPrmtvCells; % How many primitive cells exist in the XY potential
            
            nx = x(1); xdim = x(2); x0 = x(3);
            dx = xdim/nx; x = x0:dx:xdim-dx;
            unitcell.x = x;
            celldim = xdim;
            
            ny = y(1); ydim = y(2); y0 = y(3);
            dy = ydim/ny; y = y0:dy:ydim-dy;
            unitcell.y = y;
            celldim = [celldim ydim];
            
            
            if ~isempty(z) && z(4) == 1
                nz = z(1); zdim = z(2); z0 = z(3);
                dz = zdim/nz; z = z0:dz:zdim-dz;
                unitcell.z = z;
                celldim = [celldim zdim];
            end

            if ~isempty(theta) && theta(4) == 1
                ntheta = theta(1); thetadim = theta(2); theta0 = theta(3);
                dtheta = thetadim/ntheta; theta = theta0:dtheta:thetadim-dtheta;
                unitcell.theta = theta;
                celldim = [celldim thetadim];
            end
            
            unitcell.celldim = celldim;
        end
        function params = r_init(params)
            % R_INIT 
            % The function calculates random initial positions for each
            % particle, with the exception that the potision in 'z' dimension, if
            % enabled, is set to be at the min of the potential.
            
            % Prepare an r_init of evenly spreaded the particles
            Nprtcl = sum([params.prtcl.Nprtcl]);
            celldim = params.supercell.celldim;                
            
            rx = linspace(0,celldim(1)*0.9,ceil(sqrt(Nprtcl*1.5)));
            ry = linspace(0,celldim(2)*0.9,ceil(sqrt(Nprtcl*1.5)));
            
            videoFlag = 1; % to have the particle in thw mid of the frame
            if Nprtcl == 1 && videoFlag == 1
                rx = celldim(1)/2; ry = celldim(2)/2;
            end
            
            [RX,RY] = meshgrid(rx,ry);
            rx = RX(:);
            ry = RY(:);
            randInd = randperm(length(rx));
            
            r_init(1,:)=rx(randInd);
            r_init(2,:)=ry(randInd);
                        
            % Now actually distribute the particles. If even distribution
            % is requested, take from r_init (and asign prtcl.r_init).
            for i=1:length(params.prtcl)
                % dimensions of r
                dim(1) = params.model_dim;
                dim(2) = params.prtcl(i).Nprtcl;
                randMat = rand(dim(1),dim(2));
                
                % set default r_inti values for all dimensions.
                params.prtcl(i).r_init =  randMat .* celldim';
                
                % Equally spaced distribution (overied prev assignment, if required)
                if 1
                    disp('equally distributed particles')
                    params.prtcl(i).r_init(1:2,:) = r_init(:,1:dim(2));
                    r_init(:,1:dim(2)) = [];
                end
                
                % set z_init to the bottom of the adsorption well (in 'z')
                % If 'z' is enabled, its expected to be the 3rd dimension
                if params.z_enabled > 0
                    [z_minV,z_minV_indx]=min(params.prtcl(i).pes.PotMatrix(1,1,:,1));
                    params.prtcl(i).r_init(3,:) = repmat(params.unitcell.z(z_minV_indx),1,dim(2));
                end
            end
            
        end
        
        function params = p_init(params)
            % P_INIT
            % Calculate the initial momentum for each particle, complying
            % with thermal distribution. zero_p_init=1 will force zero
            % initial potential
            
            mltply = 1;
            if params.zero_p_init == 1, mltply = 0; end
            
            for i=1:length(params.prtcl)
                % dimensions of p
                dim(1) = params.model_dim;
                dim(2) = params.prtcl(i).momenta_dimension;
                dim(3) = params.prtcl(i).Nprtcl;

                sign_mat = randi([0,1],dim(1),dim(2),dim(3)) * 2 - 1;
                rand_mat = rand(size(sign_mat))*mltply;
                
                m = params.prtcl(i).mass;
                m = repmat(m,2+double(params.z_enabled>0),1);
                if params.theta_enabled > 0
                    m = [m ; params.prtcl(i).angular_mass];
                end
                m = repmat(m,1,dim(2),dim(3));
                params.prtcl(i).p_init = sign_mat .* sqrt(-1.0 * params.k_B .* m .* params.T .* log(1.0-rand_mat));
            end

        end
        
        function get_data_filename(proj_name)
            % GET_DATA_FILENAME works out where to save the simulation
            % The function will save the path to a global variable.
            % 
            
            global data_file
            global data_path
            global verbal
            global data_path_script_name
            global sub_job_path
            
            if isempty(sub_job_path), sub_job_path = '.'; end
            
            % First priority, check if data_path is defined in md_data_file (expected for
            % sweepParams multi-run job)        
            if exist([sub_job_path '/md_data_file.mat'],'file')
                disp('loading md_data_file.mat into data_file (global var)')
                load([sub_job_path '/md_data_file.mat'],'data_file_path')
                data_file = data_file_path;
                [data_path,~,~] = fileparts(data_file);
                
            % Otherwise, if data_path_script_name is defined (see prep_environment), workout the data_path
            % using user specific script - the script is expected to define
            % data_file and data_path
            elseif exist('data_path_script_name','var') && ~isempty(data_path_script_name) && exist(data_path_script_name,'file')
                if verbal, disp('data_path_script_name exist'); end
                run(data_path_script_name);
                
            % If nothing else, assume a scratch folder (or similar), and
            % save in monotonically increasing file identifiers for md_local_<serial number>.mat
            else
                if exist('~/rds/hpc-work','dir')
                    data_path = '~/rds/hpc-work';
                elseif exist('~/scratch','dir')
                    data_path = '~/scratch';
                elseif exist('/tmp','dir')
                    data_path = '/tmp/';
                else
                    data_path = '';
                end
                
                % convert symbolic links to the actual path
                if ~ispc
                    [~,b]=system(['readlink -f ' data_path]);
                    data_path = b(1:end-1);
                end
                serial = prepFuncs.get_serial_num(data_path);
                
                serstr = sprintf('%.6d',serial);
                filename = ['md_local_' serstr '.mat'];
                data_file = [data_path '/' filename];
                if ispc && isempty(data_path), data_file = filename; end
                filetime.creation=fix(clock);
                save(data_file,'filetime');                
            end
        end

        function serial = get_serial_num(data_path)            
            tmp=dir(data_path);
            for i=1:length(tmp)
                str=tmp(i).name;
                if isempty(regexp(str,'md_local_\w*.mat')), continue; end
                C = strsplit(str,'md_local_'); c = strsplit(C{2},'.');
                d(i)=str2num(c{1});
            end
            if ~exist('d') d=0; elseif isempty(d), d=0; end
            serial = max(d)+1;
        end
        
        function paramIndMat = allPerm(varargin)
            % ALLPERM spans space of variables
            % Recurcively span all the permutations of values in cell-array
            % arguments. For example, if we have 3 parameters with 3
            % variable each, var1={val1.1,val1.2,val1.3},
            % var2={val2.1,val2.2,val2.3}, etc, paramIndMat =
            % allPerm(var1,var2,var3) will generate paramIndMat were
            % {:,j} is a combination of var1-3. So size(paramIndMat)=3x27
            % Each cell-array of values can be of strings, numbers, but not
            % both (as its meaningless). For example, varargin could be
            % {{1,2},{'aa','bb'}} but not {{1,'bb'},{'aa',2}}

            if nargin == 1
                paramIndMat = varargin{1};
                return
            end

            % After all calls are done, paramIndMat holds (in every 'call')
            % the previous permutated parameters.
            paramIndMat = prepFuncs.allPerm(varargin{1:end-1});
            L = size(paramIndMat,2);
            
            % concatenate repetitions of the previous permutated
            % parameters, as preparation to add the 'new' parameters
            paramIndMat = repmat(paramIndMat,1,length(varargin{nargin}));
            
            % repeat the new parameters in the 'amount' needed for all the
            % previous permutations
            tmp_new_params = repmat(varargin{nargin},1,L);
            if iscell(tmp_new_params) && isnumeric(tmp_new_params{1})
                tmp_new_params = cell2mat(tmp_new_params);
            end
            
            tmp = sort(tmp_new_params);
            if isnumeric(tmp), tmp = num2cell(tmp); end
            if isnumeric(paramIndMat), paramIndMat = num2cell(paramIndMat); end
            paramIndMat(nargin,:) = tmp;
        end
        
        function PotMatrix = adjust_PES(params,PotMatrix)
            % ADJUST_PES reduce the PES to the enabled dimensions
            % The function assums that if the dimensions of the PES are
            % higher than the enabled dimensions of the simulation, the PES
            % needs to be reduced. The PES is reduced at the min PES value of each inactive
            % dimention. For example, if the PES is V=V(x,y,z), but 'z'
            % is not enabled, the PES will be rediced to V=V(x,y,z1), where
            % z1 is the z-value in which the PES is minimal (so called bottom of the
            % well, calculated for [x,y]=[1,1]).
            
            if length(size(PotMatrix)) == 4
                if ~params.z_enabled
                    [zvl,indx] = min(PotMatrix(1,1,:,1));
                    PotMatrix = PotMatrix(:,:,indx,:);
                end
                if ~params.theta_enabled
                    [zvl,indx] = min(PotMatrix(1,1,1,:));
                    PotMatrix = PotMatrix(:,:,:,indx);
                end
            elseif length(size(PotMatrix)) == 3 && ~params.theta_enabled && ~params.z_enabled
                [zvl,indx] = min(PotMatrix(1,1,:));
                PotMatrix = PotMatrix(:,:,indx);
            end
            PotMatrix = squeeze(PotMatrix);
        %     PotMatrix(:,end+1,:,:)=PotMatrix(:,1,:,:);
        %     PotMatrix(end+1,:,:,:)=PotMatrix(1,:,:,:);        
        end
        
    end
    
end

