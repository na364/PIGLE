% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function conf = prepare_configuration(r_conf,varargin)
%PREPARE_CONFIGURATION of a molecule/particle
% Create a description of subgroups (or atoms) composing a
% molecule/particle.
% 
% Inputs: conf_params struct with fields
%           caseNum - Selects the type of configuration.
%                     Case 1 - top-symmetric benzene like molecule.
%
%
% Output: conf struct with the following fields
%           r_conf(i,:)      - the coordinates (x,y) or (x,y,z) of the i'th scattering center
%           form_factor(i,:) - F(dK) for the i'th scattering center
%

%% parse varargin
prsdArgs = inputParser;   % Create instance of inputParser class.
prsdArgs.addParameter('dK', [], @isnumeric);
prsdArgs.addParameter('theta_tot', [], @isnumeric);
prsdArgs.addParameter('beam_ki', [], @isnumeric);
prsdArgs.addParameter('form_factor_conf', [], @isstruct);
prsdArgs.addParameter('CoM_form_factor_conf', [], @isstruct);
prsdArgs.parse(varargin{:});

dK = prsdArgs.Results.dK;
theta_tot = prsdArgs.Results.theta_tot;
beam_ki = prsdArgs.Results.beam_ki;
form_factor_conf = prsdArgs.Results.form_factor_conf;
CoM_form_factor_conf = prsdArgs.Results.CoM_form_factor_conf;

Natoms = r_conf.Natoms;

%% Create potisions of scattering centers
switch r_conf.caseNum
    
    % Create a top symmetric molecule.
    % r0 is the radius, and Natoms is the number of atoms. These arguments
    % can be provided directly - prepare_configuration(1,r0,Natoms) or as
    % fields of a structure - prepare_configuration(1,strct).
    case 1        
        r0 = r_conf.r0;
        rotAngle=linspace(0,2*pi,Natoms+1);
        rotAngle(end)=[];
        r_conf = zeros(Natoms,3);
        for i=1:Natoms
            r_conf(i,1:2) = hlp_f.Rmat(rotAngle(i))*r0*[1 0]';    
        end
    otherwise
        warning('Configuration was NOT set-up')
end

conf.r_conf = r_conf;

%% Create form factor for each scattering center

if ~isempty(form_factor_conf)
    for j=1:Natoms
        conf_params_scat_centers = form_factor_conf.scatCntr(j);
        for i=1:size(dK,3) % for each azimuth dK might be different
            dK_norm = vecnorm(dK(:,:,i),2,2);
            [FF(:,j,i),FF_envelop(:,j,i),~] = calc_form_factor(conf_params_scat_centers,'dK',dK_norm,'beam_ki',beam_ki,'theta_tot',theta_tot);
        end
    end
    form_factor.FF = FF; form_factor.FF_envelop = FF_envelop;
end

% Create average FF.
if ~isempty(CoM_form_factor_conf)
    for i=1:size(dK,3) % for each azimuth dK might be different
        dK_norm = vecnorm(dK(:,:,i),2,2);
        [FF_CoM(:,i),FF_envelop_CoM(:,i),~] = calc_form_factor(CoM_form_factor_conf,'dK',dK_norm,'beam_ki',beam_ki,'theta_tot',theta_tot);
    end
elseif ~isempty(form_factor_conf)
    % TODO test
    FF_CoM = squeeze(mean(FF,2));
    FF_envelop_CoM = squeeze(mean(FF_envelop,2));
else
    
end
    
if exist('FF_CoM','var') % implies that structure form_factor exists
    form_factor.FF_CoM = FF_CoM; form_factor.FF_envelop_CoM = FF_envelop_CoM;
    conf.form_factor = form_factor;
end

end




