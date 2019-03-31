
#ifndef HEADER_INTERP1LIN
#define HEADER_INTERP1LIN

#include <stdio.h>
#include <stdlib.h>

/*
 * Reduction of 'splint.c' from "Numerical Recipes in C" (splint, 3.3)
 * Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order),
 * and given a value of x, this routine returns a linear interpolated value y.
 */
void interp1lin(double xa[], double ya[], int n, double x, double *y);

#endif
