function [wavenumber]=KE2k(Kinetic_Energy)
%Converts a kinetic energy in meV to a wavenumber in A^-1 for a He-3 + ion.

wavenumber=10^-10* 2.99792e8*sqrt((Kinetic_Energy*1.78266*10^-36*0.001+3.01603*1.660539*10^-27)^2-(3.01603*1.660539*10^-27)^2)/(1.05457e-34);
end