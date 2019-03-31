
#ifndef HEADER_INTRADSRBTFRCS
#define HEADER_INTRADSRBTFRCS

#include "interp1lin.h"
#include <math.h>
#include <stdio.h>


void intrAdsrbtFrcs(double *iaf, const double * pos,
                const int * posDims, const double * x, const double * fTbltd, const int * fTbltdDims,
                const double * celldim, const double  * identity, const double * f_perm,
                const double * f_func, const double * in_cutoff_r, const double * out_cutoff_r,
                const double * z_enabled, const double * freeze, const double * active, const int * clisti, const int * clist);


#endif
