
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
#define y_width 1

/*
* Create external references here.  
*
*/
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
/* extern double func(double a); */

extern void setclist(int * clisti, int * clist, const double * pos, const int * posDims, const double * celldim, const double * out_cutoff_r);
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
* Output functions
*
*/
void frmBuilder_sfun_setclist_Outputs_wrapper(const real_T *pos,
const int32_T *calc_clist,
int32_T *clisti,
int32_T *clist,
const real_T *celldim, const int_T p_width0,
const real_T *out_cutoff_r, const int_T p_width1)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
int_T posDims[2]  = {2,57};

if ( calc_clist[0] == (int32_T) 0)
setclist(clisti,clist,pos,posDims,celldim,out_cutoff_r);
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


