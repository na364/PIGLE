% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function [V,nx,ny,celldim, X, Y] = hexagonal(celldim, nx, ny, scalingF, A, p)

a=celldim(1);
a1=celldim(1);
a2=celldim(2);

zeta = 4*pi/sqrt(3)/a;
g = [0 zeta; zeta*sin(pi/3) -zeta*cos(pi/3); -zeta*sin(pi/3) -zeta*cos(pi/3)];

V=zeros(nx,ny);
X=zeros(nx,ny);
Y=zeros(nx,ny);
for ix=1:nx
    x=(ix-1)*a1/nx;
    for iy=1:ny
        y=(iy-1)*a2/ny;
        X(ix,iy)=x;
        Y(ix,iy)=y;
        r=[x,y]';
        for n=1:length(A)
            for i=1:3
                V(ix,iy)=V(ix,iy)+A(n)*cos(n*g(i,:)*r);
            end
        end
    end
end

% V = -V;
% V = V-max(max(V));
% V = V/(max(max(V))-min(min(V)));
% V = scalingF*V;

V = V+abs(min(V(:)));
V = V./max(abs(V(:)));
V = V.^p;
V = V./max(abs(V(:)));
V=V.*-1*scalingF; %negative pot
