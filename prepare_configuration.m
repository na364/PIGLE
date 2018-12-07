% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function r_conf = prepare_configuration(caseNum, varargin)
%PREPARE_CONFIGURATION of a molecule/particle
% Create a description of subgroups (or atoms) composing a
% molecule/particle.
% 
% Inputs:
%           caseNum - Selects the type of configuration.
%                     Case 1 - top-symmetric benzene like molecule.
%
%           varargin - all the arguments which are case specific.
%                      Case 1 - r0 (distance from center of mass)
%                               Natoms (number of atoms).
%
% Output:
%           r_conf(i,:) - the coordinates (x,y) or (x,y,z) of the point 'i'
%

if isstruct(caseNum) && isfield(caseNum,'caseNum')
    varargin = caseNum;
    caseNum = varargin.caseNum;
    varargin.caseNum = [];
end

%% Create configuration
switch caseNum
    
    % Create a top symmetric molecule.
    % r0 is the radius, and Natoms is the number of atoms. These arguments
    % can be provided directly - prepare_configuration(1,r0,Natoms) or as
    % fields of a structure - prepare_configuration(1,strct).
    case 1
        if isstruct(varargin) && isfield(varargin,'r0') && isfield(varargin,'Natoms')
            r0 = varargin.r0; Natoms = varargin.Natoms;
        elseif length(varargin) == 2
            r0=varargin{1}; Natoms = varargin{2};
        else
            error('Failure in setting configuration (case #1)')
        end

        rotAngle=linspace(0,2*pi,Natoms+1);
        rotAngle(end)=[];
        r_conf = zeros(Natoms,3);
        for i=1:Natoms
            r_conf(i,1:2) = hlp_f.Rmat(rotAngle(i))*r0*[1 0]';    
        end
    otherwise
        warning('Configuration was NOT set-up')
end
        
end