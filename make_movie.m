% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.


function F = make_movie(params,data,varargin)
% MAKE_MOVIE from PIGLE simulation
% A function to create a movie structure from 'data' structure which is the output of run_pigle.
% It makes a use of the 'getframe' matlab function. The movie can be then
% saved.
%
% Inputs:
%    params   - pre- PIGLE-simulation configuration struct.
%    data     - result struct from PIGLE simulation
%    varargin - optional parameters, given as <param name>, <param value>
%                loop_user_limit, 500     - number of frames.
%                dt_movie, 0.1            - movie time step in ps.
%                isBallMode, 1            - should each trajectory be represented as a ball, or a point
%                rSphere, 1               - the radius of each ball, if BallMode is turned on.
%                enforce_configuration, 0 - if a configuration wasn't applied at (post) simulation
%                                           time, replace the center of mass and the azimuthal information with a
%                                           configuration of molecule.
%                k_view                   - How many subplots (and views) are displayed. Value should be of the form
%                                           {view1,view2,..}, were view1 is either "view number' as known in MATLAB,
%                                           or a vector [angle1 angle2] which defines the view
%                isf                      - the name of variable to load from the data file
%                dK_indx                  - indx of dK to plot as the trajectories evolve
%                show_trace               - trace of each particle
%
% Output
%    F - movie struct.
%

%% parse varargin
prsdArgs = inputParser;   % Create instance of inputParser class.
prsdArgs.addParameter('loop_user_limit', 500, @isnumeric);
prsdArgs.addParameter('dt_movie', 0.1, @isnumeric);
prsdArgs.addParameter('isBallMode', 1, @isnumeric);
prsdArgs.addParameter('rSphere', 0.5, @isnumeric);
prsdArgs.addParameter('enforce_configuration', 0, @isnumeric);
prsdArgs.addParameter('k_view', {[30 65]}, @iscell);
%prsdArgs.addParameter('k_view', {2,[-80 10]}, @iscell);
prsdArgs.addParameter('isf', [], @isnumeric);
prsdArgs.addParameter('dK_indx', [], @isnumeric);
prsdArgs.addParameter('show_trace', 0, @isnumeric);
prsdArgs.addParameter('pesDigitization', 10, @isnumeric);

prsdArgs.parse(varargin{:});

loop_user_limit = prsdArgs.Results.loop_user_limit;
dt_movie = prsdArgs.Results.dt_movie;
isBallMode = prsdArgs.Results.isBallMode;
rSphere = prsdArgs.Results.rSphere;
enforce_configuration = prsdArgs.Results.enforce_configuration;
k_view = prsdArgs.Results.k_view;
isf = prsdArgs.Results.isf;
dK_indx = prsdArgs.Results.dK_indx;
show_trace = prsdArgs.Results.show_trace;
pesDigitization = prsdArgs.Results.pesDigitization;

%% Calculate number of frames (=loops), the step size between frames, and reduce r_supercell to speed up
if ~isfield(data.prtcl(1),'r_supercell')
    for i=1:length(data.prtcl), data.prtcl(i).r_supercell = hlp_f.calc_r_supercell(params,i,data.prtcl(i).r); end
end
steps = ceil(dt_movie/params.isf_sample_time); max_loops = floor(length(params.t_isf)/steps);
loops = min(max_loops, loop_user_limit);
params.t(loops*steps+1:end) = [];
params.t_isf(loops*steps+1:end) = [];
if ~isempty(isf), isf(:,loops*steps+1:end) = []; end
for i=1:length(data.prtcl)
    data.prtcl(i).r_supercell(:,:,loops*steps+1:end) = [];
end

%% 
clear F h1
h=figure;

%% Create configuration
if ~isfield(data.prtcl,'conf') || ...
       ~isfield(data.prtcl(1).conf,'r_conf') || enforce_configuration
    disp('Enforcing configuration')
    caseNum =1; r0=1; Natoms = 8;
    for i=1:length(data.prtcl)
        data.prtcl(i).r_conf = prepare_configuration(caseNum,r0,Natoms);
    end
end

%% Create background image - PES for the supercell at min (well) in 'z' - by tiling the unitcell PES
% Reminder: The PES should be constructed such that each row stands for values at specific y, and varaying x
% Also, this part is badly written - needs rewriting

colors = 'mbcrgk'; colors = repmat(colors,1,20);
lineStyle = {'-','--','-.'}; lineStyle = repmat(lineStyle,1,20);

dim = params.supercell.celldim./params.unitcell.celldim;
dim = round(dim);

% If the PES has 'z' dimension with more than two 'layers' of XY, take the
% XY PES for the minimum of the PES in 'z'
if sum(params.z_enabled) > 1 && length(params.unitcell.z)>2
    [z_minV,z_minVindx]=min(params.prtcl(1).pes.PotMatrix(1,1,:,1));
else
    z_minVindx=1;
end
pes = squeeze(params.prtcl(1).pes.PotMatrix(:,:,z_minVindx,1));
pes = pes-max(max(pes));

% reduce the PES
i_ = 0;
for i=1:pesDigitization:size(pes,1)
    i_ = i_ + 1; j_ = 0;
    for j=1:pesDigitization:size(pes,2)
        j_ = j_ + 1;
        pes_(i_,j_)=pes(i,j);
    end
end

pes1 = repmat(pes_,dim(2),dim(1));

x2 = linspace(0,params.supercell.celldim(1),size(pes1,2)); % columns of pes1 are 'x'
y2 = linspace(0,params.supercell.celldim(2),size(pes1,1)); % rows of pes1 are 'y'
pes3=pes1/max(max(abs(pes)))*2;
if isempty(find(~isnan(pes3),1)), pes3 = zeros(size(pes3)); end

[x3,y3]=meshgrid(x2,y2);
fPES=scatteredInterpolant(x3(:),y3(:),pes3(:));

for k=1:length(k_view)
    subp_h(k) = subplot(length(k_view),2,k,'FontSize', 24);
    hold on; surf(x2,y2,pes3);
    shading interp; caxis([min(min(pes3)) max(max(pes3))])
end

%% Create background image - ISF
if ~isempty(isf)
    if isempty(dK_indx), dK_indx = 1:size(isf,1); end
    subp_h(length(k_view)+1) = subplot(length(k_view),2,length(k_view)+1,'FontSize', 24);
    hold on
    lgnd = {};
    for i=1:length(dK_indx)
        h2(i) = plot(params.t(1),isf(i,1),lineStyle{i},'LineWidth',4);
        lgnd = {lgnd{:},['\Delta K = ' num2str(norm(params.dK(dK_indx(i))))]};
    end
    axis([params.t_isf(1) params.t_isf(end) -0.1 1])
    title('ISF');
    xlabel('Spin Echo Time [ps]')
    ylabel('Polarization')
    legend(lgnd)
end

%% Create initial plot

% Create a sphere
[x, y, z] = sphere(10);
x=x*rSphere;
y=y*rSphere;
z=z*rSphere;

% Plot the particles (either as points of spheres)
for i=1:length(data.prtcl)
    Nprtcl = length(data.prtcl(i).r_supercell(1,:,1));
    r_conf1 = hlp_f.calc_new_r(data.prtcl(i).r_supercell(:,:,1),data.prtcl(i).conf.r_conf,params.z_enabled,params.theta_enabled);
    if ~params.z_enabled, r_conf1 = [r_conf1;zeros(1,size(r_conf1,2))]; end
    r_conf1(3,:) = r_conf1(3,:) + 3;
    if ~params.z_enabled, r_conf1(3,:) = r_conf1(3,:) + fPES(r_conf1(1,:),r_conf1(2,:)); end
    
    for k=1:length(k_view)
        subplot(length(k_view),2,k);
        if isBallMode
            for j=1:size(r_conf1,2)
                hold on; h1(i,j,k) = surf(r_conf1(1,j)+x,r_conf1(2,j)+y,r_conf1(3,j)+z);
            end
        end
        
        if show_trace
            hold on; h3(i) = plot3(r_conf1(1,:),r_conf1(2,:),r_conf1(3,:),[colors(i) '.'],'MarkerSize',10);
        end
    end
end

% Manipulating the figures, in a specific order -shading, color, light, etc. (not
% sure why, but it only works in that order!)
for k=1:length(k_view)
    subplot(length(k_view),2,k)
    shading interp;
end

for i=1:length(data.prtcl)
    for k=1:length(k_view)
        tmp = isgraphics(h1(i,:,k)); % if its not just a placeholder due to different Nprtcl in each data.partcl
        set(h1(i,tmp,k), 'FaceColor',colors(i));
    end
end

for k=1:length(k_view)
    subplot(length(k_view),2,k)
    ch(k) = camlight; lighting phong
    ch(k).Position=[150 90 120];
    view(k_view{k})
    xlabel('x / $\rm{\AA}$', 'interpreter', 'LaTex');
    ylabel('y / $\rm{\AA}$', 'interpreter', 'LaTex');
    zlabel({'Potential Energy','(Arbitrary Units)'});
    title('Supercell position')
    axis equal;
    axis([0 params.supercell.celldim(1) 0 params.supercell.celldim(2) min(min(pes3)) max(max(pes3))+5+(params.z_enabled)*10])
    
    axis manual
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    set(gcf, 'Position', get(0, 'Screensize'));
end

% Adjust size of subplots
for k=1:length(k_view)
    Position(k,:) = get(subp_h(k), 'Position');
end
% 
% set(subp_h(1), 'Position',Position(1,:).*[-0.75 0.75 1.5 1.5]);
% if length(k_view) > 1
%     set(subp_h(2), 'Position',Position(2,:).*[1.5 1 0.75 0.75]);
% end

%% Record movie
F(loops) = struct('cdata',[],'colormap',[]);
for l = 1:loops
    
  %  delete(h1)
%     surf(x2,y2,pes2); view(2)
    title(['t = ' num2str(round(params.t(l*steps))) ' [ps] out of ' num2str(round(params.t(loops*steps))) ' [ps]'])
%     shading interp
    
    for i=1:length(data.prtcl)
        Nprtcl = length(data.prtcl(i).r_supercell(1,:,1));
        r_conf1 = hlp_f.calc_new_r(data.prtcl(i).r_supercell(:,:,l*steps),data.prtcl(i).conf.r_conf,params.z_enabled,params.theta_enabled);
        if ~params.z_enabled, r_conf1 = [r_conf1;zeros(1,size(r_conf1,2))]; end
        r_conf1(3,:) = r_conf1(3,:) + 3;
        if ~params.z_enabled, r_conf1(3,:) = r_conf1(3,:) + fPES(r_conf1(1,:),r_conf1(2,:))*0.3; end
        
        for j=1:size(r_conf1,2)
            for k=1:length(k_view)
                if ~isBallMode
                    x=0;y=0;z=0;
                end
                h1(i,j,k).XData = r_conf1(1,j)+x;
                h1(i,j,k).YData = r_conf1(2,j)+y;
                h1(i,j,k).ZData = r_conf1(3,j)+z;
                if show_trace
                    h3(i).XData = [h3(i).XData r_conf1(1,:)];
                    h3(i).YData = [h3(i).YData r_conf1(2,:)];
                    h3(i).ZData = [h3(i).ZData r_conf1(3,:)];
                end
            end
        end
        
    end
    
    %% Create background image - ISF
    if ~isempty(isf)
        for i=1:length(dK_indx)
            h2(i).XData = params.t_isf(1:l*steps);
            h2(i).YData = isf(dK_indx(i),1:l*steps);
        end
    end
    
    drawnow    
    F(l) = getframe(h);
end

close(h)

disp('The movie can be PLAYED with the command:')
disp('figure; movie(gcf,F,1,50)')
disp('')
disp('The movie can be SAVED with the commands:')
disp('v = VideoWriter(''diff_20prtcl.mp4''); v.open; writeVideo(v,F); v.close; clear v')
