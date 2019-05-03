
function [FF,FF_envelop,dK] = calc_form_factor(form_factor_conf,varargin)
% CALC_FORM_FACTOR will return a form factor based on
% requested case number

%% parse varargin
prsdArgs = inputParser;   % Create instance of inputParser class.
prsdArgs.addParameter('dK', [], @isnumeric);
prsdArgs.addParameter('theta_tot', [], @isnumeric);
prsdArgs.addParameter('beam_ki', [], @isnumeric);
prsdArgs.parse(varargin{:});

dK = prsdArgs.Results.dK;
theta_tot = prsdArgs.Results.theta_tot;
beam_ki = prsdArgs.Results.beam_ki;

if ~isempty(dK) && ~isempty(theta_tot) && ~isempty(beam_ki)
    % assuming elastic scattering, beam_ki=beam_kf
    theta_i = get_theta_i(beam_ki, dK, theta_tot);
    theta_f = theta_tot - theta_i;
end

if ~isfield(form_factor_conf,'caseNum') || (isfield(form_factor_conf,'caseNum') && (isempty(form_factor_conf.caseNum) || isnan(form_factor_conf.caseNum)))
    FF=nan(size(dK)); FF_envelop=nan(size(dK)); return
end
    

switch form_factor_conf.caseNum
    case 0
        FF = ones(size(dK)); FF_envelop = ones(size(dK));
    case 1
        r0         = form_factor_conf.hemisphere_radius;
        [FF,FF_envelop,dK] = hemisphere_form_factor(beam_ki,r0,theta_i,theta_f);    
    otherwise
end

%plot(dK,FF,dK,FF_envelop,dK,FF.^2,dK,FF_envelop.^2); legend('FF', ...
% 'FF_{envelop}', 'FF^2', 'FF_{envelop} ^2'); xlabel('$\Delta
% K$','interpreter','latex'); set(gca, 'YScale', 'log')

end


function [FF,FF_envelop,dK] = hemisphere_form_factor(ki,r0,theta_i,theta_f)
    % HEMISPHERE_FORM_FACTOR is based on Lahee et al http://dx.doi.org/10.1063/1.452321
    % The fundtion returns the Fraunhofer ellastic form factor
    % assuming the interaction of the helium atom with the scattering
    % particle can be represented as a hard-wall hemisphere potential.
    %
    % Input:
    %       ki - incident wavevector (assumed to be also the final wavevector
    %       r0 - the apparent height of the adsorbate
    %       theta_i,theta_f - initial and final scattering angles with respect to the surface normal
    %
    % Outputs:
    %       FF       - the Fraunhofer form factor
    %       FF_envelop - the envelop of the FF (significally diverges from I near dK=0
    %       dK      - the calculated parallel momentum transfer (assuming ellastic scattering)
    %

    dK=ki*(sind(theta_f)-sind(theta_i));
    theta=theta_f-theta_i;
    z=ki*r0*sind(theta);
    FF=(cosd(theta)+1) .* besselj(1,z)./sind(theta);
    FF_envelop = (cosd(theta)+1) .* sqrt(2/pi./z)./sind(theta);
end


function theta_i = get_theta_i(ki, dK, theta_tot)
% Returns the theta_i equivalent to 'dK'

theta_i_vec=-89.9:0.025:89.0;
dKvec = ki*(sind(theta_tot-theta_i_vec)-sind(theta_i_vec));

theta_i = zeros(length(dK),1);
for i=1:length(dK)
    indx=[]; threshold=0.01;
    while isempty(indx)
        indx=find(dKvec-dK(i)<threshold & dKvec-dK(i)>-threshold);
        threshold = threshold+0.01;
    end
    theta_i(i) = theta_i_vec(indx(ceil(length(indx)/2)));
end

end
