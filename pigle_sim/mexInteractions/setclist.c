
    #include <math.h>
    #include <stdio.h>


void setclist(int * clisti, int * clist, const double * pos, const int * posDims, const double * celldim, const double * out_cutoff_r)
{
    
    int      i, j, iat1, iat2, ic, nop;
    double   R2, sdimhx, sdimhy, dx, dy;
    
    nop = posDims[1];
    
    sdimhx = celldim[0] * 0.5;
    sdimhy = celldim[1] * 0.5;
    
    /*
     * if (1) {
     *  FILE * fp;
     * fp = fopen ("file.txt", "a");        
     * fprintf(fp, "out_cutoff_r=%f\n",out_cutoff_r[0]);
     * fclose(fp);
     * }
     */
        
    for (iat1=0; iat1<nop; iat1++) {

        clisti[iat1] = 0;
        
        for (iat2=iat1+1; iat2<nop; iat2++) {
            
            /*
             * calculate dx and dy and total distances between atoms, and account
             * for wraparound of the simulation supercell
             */
            
            dx = pos[posDims[0]*iat1]  -pos[posDims[0]*iat2];
            dy = pos[posDims[0]*iat1+1]-pos[posDims[0]*iat2+1];
            
            /* Ignore 'z' dimention, assuming that the impact on the radius is not major.*/
            /*
            if (posDims[0] > 2) {
                dz = pos[posDims[0]*iat1+2]-pos[posDims[0]*iat2+2];
            } else {
                dz = 0;
            }
             */

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
            
            if (R2 < out_cutoff_r[0]*out_cutoff_r[0]) {

                clisti[iat1] = clisti[iat1] + 1;
                clist[iat1*nop+clisti[iat1]-1]=iat2;
            }
        }
    }
    
} /* end setclist */
