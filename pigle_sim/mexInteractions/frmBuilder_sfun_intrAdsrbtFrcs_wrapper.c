
/*
* Include Files
*
*/
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#include "rtwtypes.h"
#endif



/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */
#include <math.h>
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 2
#define y_width 57

/*
* Create external references here.  
*
*/
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
/* extern double func(double a); */

void intrAdsrbtFrcs(double *iaf, const double * pos,
const int * posDims, const double * x, const double * fTbltd, const int * fTbltdDims,
const double * celldim, const double  * identity, const double * f_perm,
const double * f_func, const double * in_cutoff_r, const double * out_cutoff_r,
const double * z_enabled, const double * freeze, const double * active, const int * clisti, const int * clist);
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
* Output functions
*
*/
void frmBuilder_sfun_intrAdsrbtFrcs_Outputs_wrapper(const real_T *pos,
const real_T *x,
const real_T *fTbltd,
const real_T *celldim,
const real_T *identity,
const real_T *f_perm,
const real_T *f_func,
const real_T *in_cutoff_r,
const real_T *out_cutoff_r,
const real_T *z_enabled,
const real_T *freeze,
const real_T *active,
const int32_T *clisti,
const int32_T *clist,
real_T *iaf,
const real_T *nop_1, const int_T p_width0,
const real_T *size_Fint_1, const int_T p_width1,
const real_T *model_dim_1, const int_T p_width2)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
int_T posDims[2]  = {2,57};
int_T fTbltdDims[2]  = {500,3};

/*
//posDims[0] = (int_T) &model_dim;
//posDims[1] = (int_T) &nop;
//int_T             *fTbltdDims = (int_T *) size_Fint;
*/

intrAdsrbtFrcs(iaf, pos, posDims,
x, fTbltd, fTbltdDims,
celldim, identity, f_perm, f_func,
in_cutoff_r, out_cutoff_r,
z_enabled, freeze, active,clisti,clist);
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


