
#include "interp1lin.h"

/*
 * Reduction of 'splint.c' from "Numerical Recipes in C" (splint, 3.3)
 * Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order),
 * and given a value of x, this routine returns a linear interpolated value y.
 */
void interp1lin(double xa[], double ya[], int n, double x, double *y)
{
    
    int klo,khi,k;
    double h,b,a;
    klo=1;
    khi=n;
    while (khi-klo > 1) {
        k=(khi+klo) >> 1;
        if (xa[k] > x) khi=k;
        else klo=k;
    } /* klo and khi now bracket the input value of x. */
    
    h=xa[khi]-xa[klo];
    if (h == 0.0) {
        fprintf(stderr,"Bad xa input to routine interp1lin"); /* The xa’s must be distinct. */
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
    }
    a=(xa[khi]-x)/h; 
    b=1-a;
    *y=a*ya[klo]+b*ya[khi];
}