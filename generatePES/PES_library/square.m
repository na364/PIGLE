% Copyright (c) 2020, Daniel Cropper.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function [V,nx,ny,celldim, X, Y] = square(celldim, nx, ny, hollow,bridge, Tmode)

a=celldim(1);
a1=celldim(1);
a2=celldim(2);

zeta = 2*pi/a;

V=zeros(nx,ny);
X=zeros(nx,ny);
Y=zeros(nx,ny);
if Tmode==0
for ix=1:nx
    x=(ix-1)*a1/nx;
    for iy=1:ny
        y=(iy-1)*a2/ny;
        X(ix,iy)=x;
        Y(ix,iy)=y;
        V(ix,iy)=0.25*hollow*(1-cos(zeta*x))*(1-cos(zeta*y))+0.5*bridge*(1-cos(zeta*x) * cos(zeta*y));
    end
end

else
    A=zeros(nx,ny);
    B=zeros(nx,ny);
    S=zeros(nx,ny);
    for ix=1:nx
    x=(ix-1)*a1/nx;
    for iy=1:ny
        y=(iy-1)*a2/ny;
        X(ix,iy)=x;
        Y(ix,iy)=y;
        S(ix,iy)=0.25*hollow*(1-cos(zeta*x))*(1-cos(zeta*y))+0.5*bridge*(1-cos(zeta*x) * cos(zeta*y));
        A(ix,iy)=0.5*95*(x^2+y^2)*(1+0.45*(x^2+y^2));
        B(ix,iy)=exp(-1000*((x/a1)^2+(y/a2)^2)^4);
        V(ix,iy)=B(ix,iy)*A(ix,iy)+(1-B(ix,iy))*S(ix,iy);
    end
    end
end

