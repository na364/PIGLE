% Copyright (c) 2016, John Ellis.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

function [V,nx,ny,celldim, X, Y] = hexagonal6interp(celldim, nx, ny, top,frac1, frac2,bridge,hcp,fcc)

%
% This programme takes 6 points in the irreducible triangle of a 111
% surface - top, half way between top and hcp, hcp, half way betwen top and
% fcc, bridge, fcc - and from them finds the Fourier components in the
% first 6 shells of G vectors, and then from these calculates the full 2D
% potential. The idea is to give the top, bridge, fcc, hcp energies and use
% the variable frac to calculate the inbetween points. By varying frac you
% can vary the width of the channel.
%
% copper 'a' constant - nearest neighbour distance in 111 plane
a1=celldim(1);
a2=celldim(2);
%
% set up the coordinates of the 6 reference points
xypos=zeros(6,2);
xypos(1,1)=0;
xypos(1,2)=0;
xypos(3,1)=0;
xypos(3,2)=a1*sqrt(3)/2*2/3;
xypos(2,:)=(xypos(1,:)+xypos(3,:))./2;
xypos(6,1)=a1/2;
xypos(6,2)=a1*sqrt(3)/2/3;
xypos(5,:)=(xypos(3,:)+xypos(6,:))/2;
xypos(4,:)=xypos(6,:)/2;
% ixy=1 is top
V2D(1)=top; %128;
%ixy=3 is hcp
V2D(3)=hcp; %43;
%frac=0.5;
% ixy=2 is half way between bridge and top
if 0 < frac1 && frac1 < 1
    V2D(2)=(frac1*V2D(1)+(1-frac1)*V2D(3));
else
    V2D(2)=frac1;
end
%ixy=5 is bridge
V2D(5)=bridge; %0;
%ixy=6 is fcc
V2D(6)=fcc; %28;
%ixy=4 is half way between fcc and top
if 0 < frac2 && frac2 < 1
    V2D(4)=(frac2*V2D(1)+(1-frac2)*V2D(6));
else
    V2D(4)=frac2;
end
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



G1=4*pi/sqrt(3)/a1;
% now we need to set up vectors that tell you which
% G vectors are related to a particular (ia, ib) - store
% these as the a and b components of the 6 first order G
% vectors. This is a rather poor manual way of setting them up - Zianding
% has a much nicer way - but at least this is tranparent and easy to code
% First the G vectors are given wrt axes 'a' which rises at 30 deg to
% the usual x axis and 'b' which is parallel to the usual y axis. The
% notation is a bit confused - V1Ga means the a coordinate of the potential
% G vector for the first ring of G vectors away from the centre (Gx=0,
% Gy=0) - but later the index ishell runs 1 to 6 - 1 is the central
% (G=(0,0)) point, and the second ring of 6 neighbour is split into two
% shell components (ishell=2, ishell=3) which contain G vectos related by
% symmetry within the shell. All the G Vectors in the second ring of
% neighbours (ishell=4) are symmetrt related, but the next ring are split
% into two groups (ishell=5, ishell=6)

V1Ga(1)=1;
V1Gb(1)=0;
V1Ga(2)=0;
V1Gb(2)=1;
V1Ga(3)=-1;
V1Gb(3)=1;
V1Ga(4)=-1;
V1Gb(4)=0;
V1Ga(5)=0;
V1Gb(5)=-1;
V1Ga(6)=1;
V1Gb(6)=-1;
%
%
V2Ga(1)=2;
V2Gb(1)=-1;
V2Ga(2)=1;
V2Gb(2)=1;
V2Ga(3)=-1;
V2Gb(3)=2;
V2Ga(4)=-2;
V2Gb(4)=1;
V2Ga(5)=-1;
V2Gb(5)=-1;
V2Ga(6)=1;
V2Gb(6)=-2;
%
V3Ga(1)=2;
V3Gb(1)=0;
V3Ga(2)=0;
V3Gb(2)=2;
V3Ga(3)=-2;
V3Gb(3)=2;
V3Ga(4)=-2;
V3Gb(4)=0;
V3Ga(5)=0;
V3Gb(5)=-2;
V3Ga(6)=2;
V3Gb(6)=-2;
%
% now generate the actual components in Cartesian coordinates
for idum1=1:6
    Gx1(idum1)=V1Ga(idum1)*sqrt(3)/2*G1;
    Gy1(idum1)=(V1Ga(idum1)/2+V1Gb(idum1))*G1;
    Gx2(idum1)=V2Ga(idum1)*sqrt(3)/2*G1;
    Gy2(idum1)=(V2Ga(idum1)/2+V2Gb(idum1))*G1;
    Gx3(idum1)=V3Ga(idum1)*sqrt(3)/2*G1;
    Gy3(idum1)=(V3Ga(idum1)/2+V3Gb(idum1))*G1;

end
%
% and then use these to create the G vectors that we will use for fitting
% the potentail - G(ishell,indexwithinshell, idrection)
G=zeros(6,6,2);
G(1,1,1)=0;
G(1,1,2)=0;
for idum=1:2:5
    G(2,(idum+1)/2,1)=Gx1(idum);
    G(2,(idum+1)/2,2)=Gy1(idum);
    G(3,(idum+1)/2,1)=Gx1(idum+1);
    G(3,(idum+1)/2,2)=Gy1(idum+1);
    G(4,idum,1)=Gx2(idum);
    G(4,idum,2)=Gy2(idum);
    G(4,idum+1,1)=Gx2(idum+1);
    G(4,idum+1,2)=Gy2(idum+1);
    G(5,(idum+1)/2,1)=Gx3(idum);
    G(5,(idum+1)/2,2)=Gy3(idum);
    G(6,(idum+1)/2,1)=Gx3(idum+1);
    G(6,(idum+1)/2,2)=Gy3(idum+1);
end
% vector nGshell contains the number of G vectos within each shell
nGshell(1)=1;
nGshell(2)=3;
nGshell(3)=3;
nGshell(4)=6;
nGshell(5)=3;
nGshell(6)=3;
%
% now plot out these positions to check
% figure(41)
% clf
% hold on
% plot(G(1,1,1),G(1,1,2),'rx')
% plot(G(2,1:3,1),G(2,1:3,2),'gx')
% plot(G(3,1:3,1),G(3,1:3,2),'bx')
% plot(G(4,:,1),G(4,:,2),'ro')
% plot(G(5,1:3,1),G(5,1:3,2),'go')
% plot(G(6,1:3,1),G(6,1:3,2),'bo')
%
% Now set up a matrix M(ixy,ishell) such that when it is multiplied by a
% column vector, Vcoef, which contains the (complex) coefficients of each G
% vectos (remember the coef within a shell is the same for every G vector)
%
for ixy=1:6
    for ishell=1:6
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
        for ishell=1:6
            for idum=1:nGshell(ishell)
            V(ix,iy)=V(ix,iy)+Vcoef(ishell)*exp(1i*(x*G(ishell,idum,1)+y*G(ishell,idum,2)));
            end
        end
    end
end

% figure(2)
% clf
% surface(X,Y,real(V))
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
