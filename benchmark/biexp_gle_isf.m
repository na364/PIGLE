% Copyright (c) 2018, Peter Townsend.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function isf=biexp_gle(DK,time,gamma1,wc1,gamma2,wc2,kb,T,m)
% Straightforward generalisation of Townsend and Ward 2018
% method to include a biexponential friction kernel, hence involves finding
% roots of a cubic denominator not quadratic..
% Intended application: benchmark PIGLE simulation results

% assumes friction kernel gamma1*exp(-wc1*t)+gamma2*exp(-wc2*t)

c3=1; % coefficients of cubic denominator
c2=wc1+wc2;
c1=wc1*wc2+gamma1*wc1+gamma2*wc2;
c0=(gamma1+gamma2)*wc1*wc2;

r=roots([c3 c2 c1 c0]);

x=r(1); y=r(2); z=r(3); % write the roots out long-hand, clearer on balance
px=(x+wc1)*(x+wc2)/((x-y)*(x-z));
py=(y+wc1)*(y+wc2)/((y-z)*(y-x));
pz=(z+wc1)*(z+wc2)/((z-x)*(z-y));

norm_exp=px*(exp(x*time)-x*time-1)/x^2+py*(exp(y*time)-y*time-1)/y^2+pz*(exp(z*time)-z*time-1)/z^2;

exponent=-kb*T*norm_exp*DK^2/m; % scale for temperature, mass etc.. and give it the right sign 

isf=exp(exponent);
