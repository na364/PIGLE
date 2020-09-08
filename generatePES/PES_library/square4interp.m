% Copyright (c) 2020, Daniel Cropper, John Ellis.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function [V,nx,ny,celldim, X, Y] = square4interp(celldim, nx, ny, max, min, saddle, slopef)
% [V,nx,ny,celldim, X, Y] = square4interp(celldim, nx, ny, max, min, saddle, slopef)
%
% This programme is based on John% Ellis's programme hexagonal6interp.
% square4interp takes 4 points in the irreducible square of a fcc 100
% surface - the maximum (0,0), the minimum (0.5,0.5), the saddle
% point(0,0.5), and a point on the slope slopef*(0.5,0.5) - and from them
% finds the Fourier components in the first 4 shells of G vectors, and then
% from these calculates the full 2D potential. The idea is to give the max,
% min, and saddle energies and use the variable slopef to calculate the
% inbetween points. Slopef defines the point on the unit cell diagonal
% where the potential falls to the average of its min and max values. Do
% not set slopef= 1 or 0!!!! slopef is best in the range [0.4,0.6], as this
% does not cause an overshoot in the potential.
%
% copper 'a' constant - nearest neighbour distance in 100 plane = close
% packed separation
a1=celldim(1);
a2=celldim(2);
%
% set up the coordinates of the 4 reference points
xypos=zeros(4,2);
xypos(1,1)=0;
xypos(1,2)=0;
xypos(3,1)=0;
xypos(3,2)=a1/2;
xypos(2,1)=slopef*a1/2;
xypos(2,2)=slopef*a1/2;
xypos(4,1)=a1/2;
xypos(4,2)=a1/2;
% ixy=1 is max
V2D(1)=max; %128;
%ixy=4 is min
V2D(4)=min; %43;
V2D(2)=(max+min)/2;
V2D(3)=saddle;
%
% % and plot the potential
% figure(1)
% clf
% hold on
% nxy=6;
% plot3(xypos(:,1),xypos(:,2),V2D(:),'rx')
% xlabel('x')
% ylabel('y')
% title('energy of reference points')
%G1 is the distance in K space between G vectors



G1=2*pi/a1;

% now we need to set up vectors that tell you which
% G vectors are related to a particular (ix, iy) - store
% these as the x and y components of the 4 first order G
% vectors.  V1Gx means the x coordinate of the potential
% G vector for the first ring of G vectors away from the centre (Gx=0,
% Gy=0) - but later the index ishell runs 1 to 4 - 1 is the central
% (G=(0,0)) point, 2 - 4 are the first three rings of G vectors away from
% the centre: 2 has |G|=1, 3 has |G|=sqrt(2), 4 has |G|=2. All G vectors
% within a shell are symmetry related.


%Calculate the G vectors in Cartesian coordinates
V1GY(1)=1;
V1GX(1)=0;
V1GY(2)=0;
V1GX(2)=1;
V1GY(3)=-1;
V1GX(3)=0;
V1GY(4)=0;
V1GX(4)=-1;

%
%
V2GY(1)=1;
V2GX(1)=-1;
V2GY(2)=1;
V2GX(2)=1;
V2GY(3)=-1;
V2GX(3)=-1;
V2GY(4)=-1;
V2GX(4)=1;

%
V3GY(1)=2;
V3GX(1)=0;
V3GY(2)=0;
V3GX(2)=2;
V3GY(3)=0;
V3GX(3)=-2;
V3GY(4)=-2;
V3GX(4)=0;


% and then use these to create the G vectors that we will use for fitting
% the potentail - G(ishell,indexwithinshell, idirection)


% vector nGshell contains the number of G vectors within each shell


G=zeros(4,4,2);
G(1,1,1)=0;
G(1,1,2)=0;
for i=1:4
    G(2,i,1)=V1GX(i);
    G(2,i,2)=V1GY(i);
    G(3,i,1)=V2GX(i);
    G(3,i,2)=V2GY(i);
    G(4,i,1)=V3GX(i);
    G(4,i,2)=V3GY(i);
end
G=G*G1;

nGshell(1)=1;
nGshell(2)=4;
nGshell(3)=4;
nGshell(4)=4;

%
% now plot out these positions to check
% figure(41)
% clf
% hold on
% plot(G(1,1,1),G(1,1,2),'rx')
% plot(G(2,:,1),G(2,:,2),'gx')
% plot(G(3,:,1),G(3,:,2),'bx')
% plot(G(4,:,1),G(4,:,2),'ro')
%
% Now set up a matrix M(ixy,ishell) such that when it is multiplied by a
% column vector, Vcoef, which contains the (complex) coefficients of each G
% vectors (remember the coef within a shell is the same for every G vector)
%

M=zeros(4,4);
for ixy=1:4
    for ishell=1:4
        M(ixy,ishell)=0.0;
        for idum=1:nGshell(ishell)
            M(ixy,ishell)=M(ixy,ishell)+exp(1i*(xypos(ixy,1)*G(ishell,idum,1)+xypos(ixy,2)*G(ishell,idum,2)));
        
        end
    end
end
% and calculate the column of coefficients that matches the given data
% points
Vcoef=inv(M)*V2D';
%%


%nx=96;
%ny=96;
V=zeros(nx,ny);
X=zeros(nx,ny);
Y=zeros(nx,ny);
for ix=1:nx
    x=(ix-1)*a1/nx;
    for iy=1:ny
        y=(iy-1)*a2/ny;
        X(ix,iy)=x;
        Y(ix,iy)=y;
        V(ix,iy)=0.0;
        for ishell=1:4
            for idum=1:nGshell(ishell)
            V(ix,iy)=V(ix,iy)+Vcoef(ishell)*exp(1i*(x*G(ishell,idum,1)+y*G(ishell,idum,2)));
            end
        end
    end
end
% 
% figure(2)
% clf
% colorbar
% surf(X,Y,real(V))
% xlabel('x')
% ylabel('y')
% title('real(V)');
% figure(3)
% clf
% surface(X,Y,imag(V))
% xlabel('x')
% ylabel('y')
% title('imag(V)');
% 
% figure (4)
% Ampl=fft2(V);
% Ampl=fftshift(Ampl);
% clf
% surface(real(Ampl))
% figure (5)
% clf
% surface(imag(Ampl))
% %plot(Gxt,Gyt,'r.')
