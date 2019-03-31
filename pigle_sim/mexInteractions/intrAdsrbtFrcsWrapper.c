
	#include "intrAdsrbtFrcs.h"
    #define N 4
    #define Nfrc 500
    #define Nperm 3

	void main () {

        double iaf[2][N]={0}, x[Nfrc], fTbltd[Nfrc][Nperm];
        double in_cutoff_r, out_cutoff_r, z_enabled, active;
        int i, j;

        double pos[2][N] = {1,2,3,4,2,1,4,3};
        for (i=0; i<Nfrc; i++)        
            x[i] = 0.01+0.1002*i;

        for (i=0; i<Nfrc; i++)
        {
            for (j=0;j<Nperm;j++)
                fTbltd[i][j] = 1;
        }

        double celldim[2] = {43.4,70.4};
        double identity[N] = {1,1,2,2};
        double f_perm[Nperm][2] = {1,1,1,2,2,2};
        double f_func[Nperm] = {1,1,1};
        in_cutoff_r = 0.1;
        out_cutoff_r = 5.42;
        z_enabled = 1;
        double freeze[N] = {0,0,0,0};
        active = 1;

        int posDims[2] = {2,N};
        int fTbltdDims[2] = {Nfrc,Nperm};       

        intrAdsrbtFrcs(iaf, pos, posDims, x, fTbltd, fTbltdDims, celldim, identity, f_perm, f_func, &in_cutoff_r, &out_cutoff_r, &z_enabled, freeze, &active);

        for (i=0;i<N;i++)
            printf( "x=%f\ty=%f\n",iaf[1][i],iaf[2][i]);
                
    }