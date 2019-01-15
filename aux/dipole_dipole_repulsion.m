classdef dipole_dipole_repulsion
    
    % The whole procedure is adopted from Cucchetti, A. and
            % Ying, S.C., 1999. Diffusion of Na atoms on a Cu (001)
            % surface. Physical Review B, 60(15), p.11110.
            
    properties (Constant)
        eps0 = 1 / (4*pi * 1e-7) / 299792458^2 ; % Vacuum permittivity in SI units, [Amper^2 * sec^4 /Kg / meter^3] = [Farad/meter] = [C/(V*m)] = [C*m/(V*m^2)]
        Cm_2_Debye = 299792458 * 1E21; % [Debye/(C*m)]
        eps0_Debye_Angstrm = dipole_dipole_repulsion.eps0 * dipole_dipole_repulsion.Cm_2_Debye * 1E-20; % [Debye / (V * Angstrm^2)]
        coulomb_const = 1/(4*pi*dipole_dipole_repulsion.eps0_Debye_Angstrm); %[(V * Angstrm^2) / Debye]
        
    end
    
    methods (Static)
        
        function [wrk_fnc] = topping_work_function(mu0,n0,theta,alpha)
            %TOPPING_MODEL Summary of this function goes here
            %   Detailed explanation goes here
            %   Topping, J., 1927. On the mutual potential energy of a plane network of doublets. Proceedings of the Royal Society of London. Series A, Containing Papers of a Mathematical and Physical Character, 114(766), pp.67-72.
            %
            %   Inputs:
            %           mu0   [Debye]      - adsorbate dipole moment in the zero coverage limit (could be a fit parameter)
            %           n0    [1/ML/Angstrm^2]       - adsorbate surface density at a coverage of 1 ML, given,
            %                                in the case of a hexagonal surface, (2/√3)*(1/a2), where a is the surface lattice constant in units of Angstrom;
            %           theta [ML]         - adsorbate coverage
            %           alpha [Angstrom^3] - adsorbate polarizability (could be a fit parameter)
            %
            %   Outputs:
            %           wrk_fnc [eV]       The coverage depended work function
            %
            
            mu_eff = dipole_dipole_repulsion.effective_dipole_moment(mu0,alpha,n0,theta); % [Debye]
            eps0_Debye_Angstrm = dipole_dipole_repulsion.eps0_Debye_Angstrm; % [Debye / (V * Angstrm^2)]
            
            wrk_fnc = - (1/eps0_Debye_Angstrm) * n0*theta .* mu_eff; % [(V * Angstrm^2) / Debye] * [1/Angstrm^2] * [Debye] = [V]
            
        end
        
        function mu_eff = effective_dipole_moment(mu0,alpha,n0,theta)
            %
            % Topping depolarization formula
            % Topping, J., 1927. On the mutual potential energy of a plane network of doublets. Proceedings of the Royal Society of London. Series A, Containing Papers of a Mathematical and Physical Character, 114(766), pp.67-72.
            %
            
            mu_eff = mu0./(1+9*alpha*(n0*theta).^1.5);            
        end
        
        function [f_const, f] = f_Kohn_Lau(mu,r)
            
            % Assumes a dipole-dipole repulsion. The potential is V(r) =
            % 2*mu^2/r^3 . The factor '2' comes from image charges (Kohn and Lau
            % Kohn, W. and Lau, K.H., 1976. Adatom dipole moments on metals
            % and their interactions. Solid State Communications, 18(5),
            % pp.553-555)
            %
            % Inputs: 
            %       mu [Debye]   - coverage specific EFFECTIVE dipole moment
            %       r  [Angstrm] - Seperation distance between pair of particles
            % 
            % Outputs:
            %       f_const [meV*Angstrm^3] - force constant
            %       f       [meV/Angstrm]   - force
            %
            
            if nargin == 1 % assume only mu was given
                r=1;
            end
            
            f_const_debye = 6*mu^2; % [Debye^2]
            
            coulomb_const = dipole_dipole_repulsion.coulomb_const; %[(V * Angstrm^2) / Debye]
            
            Debye_2_Cm = 1/dipole_dipole_repulsion.Cm_2_Debye; % [(C*m)/Debye]
            Debye_2_eV_Angstrm = Debye_2_Cm / 1.602176634E-19 * 1e10; % [(e*Angstrm)/Debye]
            
            f_const = coulomb_const * Debye_2_eV_Angstrm * f_const_debye; % [(V * Angstrm^2) / Debye] * [(e*Angstrm) / Debye] * [Debye^2] = [eV * Angstrm^3]
            f_const = 1000 * f_const; % [meV * Angstrm^3]
            
            f = f_const/ r.^4; %  [meV * Angstrm^3] / Angstrm^4 = [meV / Angstrm]
            
        end
        
    end
end

