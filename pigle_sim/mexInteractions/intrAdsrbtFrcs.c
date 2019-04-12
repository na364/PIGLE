
    #include "interp1lin.h"
    #include <math.h>
    #include <stdio.h>


void intrAdsrbtFrcs(double *iaf, const double * pos,
                const int * posDims, const double * x, const double * fTbltd, const int * fTbltdDims,
                const double * celldim, const double  * identity, const double * f_perm,
                const double * f_func, const double * in_cutoff_r, const double * out_cutoff_r,
                const double * z_enabled, const double * freeze, const double * active, const int * clisti, const int * clist)
{
    
    int               i, j, iat1, iat2, ic, perm[2], nop;
    double            sdimhx, sdimhy, dx, dy, dz, F, iaf_tmp, F_parallel, R, R2, Rinv, r, r2, rinv;

    nop = posDims[1];
    
    for (i = 0; i < posDims[0]*posDims[1]; i++) {
        iaf[i] = 0;
    }
    
    if ( active == 0 )
        return;
    
    sdimhx = celldim[0] * 0.5;
    sdimhy = celldim[1] * 0.5;
    
    for (iat1=0; iat1 < nop; iat1++) {
        
        if ( freeze[iat1] == 1 )
                continue;
        
        /*for (iat2=iat1+1; iat2 < nop; iat2++) { */
        for (ic=0; ic < clisti[iat1]; ic++) {

            iat2=clist[iat1*nop+ic];

            if ( freeze[iat2] == 1.0 )
                continue;
            
            /*
             * Find the permutation for the force function.
             * Assume that the permutations are 'minimal' and 'sorted'.
             * For ex: {[1 1],[1 2],[2 2]} for two species system (i.e, perm [2 1] is ommitted since
             * its redundent to [1 2] - making it a minimal set. Also, its a sorted set - not sorted would mean to have
             * {[1 2], [1 1], ...})
            */
            
            if (identity[iat1] < identity[iat2]) {
                perm[0] = (int) identity[iat1];
                perm[1] = (int) identity[iat2];
            } else {
                perm[0] = (int) identity[iat2];
                perm[1] = (int) identity[iat1];
            }
              
            /*
             * Now compare the 'current' premutation ('perm') to the pre-defined permutations
             * pointed by f_perm, which is a matrix with [perm1;perm2;...]
             */
            
            for (i=0; i<fTbltdDims[1]; i++) {
                if ((int) f_perm[i] == perm[0] && (int) f_perm[fTbltdDims[1]+i] == perm[1]) {
                    perm[0] = i+1; /* perm change role, it now holds the permutaion number rather than the permutation. */
                                   /* perm[0] - permutation number, perm[1] - now is meaningless */
                    break;
                }
            }
            
            /*
             * calculate dx and dy and total distances between atoms, and account
             * for wraparound of the simulation supercell
             */
            
            dx = pos[posDims[0]*iat1]  -pos[posDims[0]*iat2];
            dy = pos[posDims[0]*iat1+1]-pos[posDims[0]*iat2+1];
           
            if (posDims[0] > 2) {
                dz = pos[posDims[0]*iat1+2]-pos[posDims[0]*iat2+2];
            } else {
                dz = 0;
            }
            
            if (dx > sdimhx) {
                dx = dx-celldim[0];
            } else if (dx < -sdimhx) {
                dx=dx+celldim[0];
            }
            if (dy > sdimhy) {
                dy = dy-celldim[1];
            } else if (dy < -sdimhy) {
                dy = dy+celldim[1];
            }
            
            R2 = dx*dx+dy*dy; /* radius projection squared */
            R = sqrt(R2);
            Rinv = 1.0/R;     /* for later saving cpu time */
            if (posDims[0] > 2) {
                r2 = R2 + dz*dz;
            } else {
                r2 = R2;
            }
            r = sqrt(r2);
            rinv = 1.0/r; /* for later saving cpu time */
            
            /* if too close, keep minimal seperation */
            if (r < in_cutoff_r[0]) {
                if (r  > in_cutoff_r[0]/4) {
                    dx = dx * (in_cutoff_r[0]) * rinv;
                    dy = dy * (in_cutoff_r[0]) * rinv;
                }              
                r = in_cutoff_r[0];
                r2 = r*r;
            }

            /* Interpulate the force for r */
            interp1lin( (double *) x, (double *) fTbltd+(fTbltdDims[0]*(perm[0]-1)),fTbltdDims[0],r,&F);
                        
            F_parallel = F*R*rinv; /* F*sin_phi */
            
            /* 'x' dimention */
            iaf_tmp = F_parallel * dx*Rinv; /* F in 'x' direction = F_parall * cosine */
            iaf[posDims[0]*iat1] += iaf_tmp;
            iaf[posDims[0]*iat2] -= iaf_tmp;
                        
            /* 'y' dimention */
            iaf_tmp = F_parallel * dy*Rinv; /* F in 'y' direction = F_parall * sine */
            iaf[posDims[0]*iat1+1] += iaf_tmp;
            iaf[posDims[0]*iat2+1] -= iaf_tmp;
            
            /* 'z' dimention */
            if (posDims[0] > 2) {
                iaf_tmp = F*dz*rinv; /* F in 'z' direction = F*cos_phi */
                iaf[posDims[0]*iat1+2] += iaf_tmp;
                iaf[posDims[0]*iat2+2] -= iaf_tmp;
            }

        }
                
        /* ssPrintf("F = [%f %f]\n",iaf[posDims[0]*iat1],iaf[posDims[0]*iat1+1]); // break; // break for debuging purposes */
        /* break; // break for debuging purposes */
    }        

} /* end intrAdsrbtFrcs */
