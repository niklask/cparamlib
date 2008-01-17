/*
 * example1.c
 *
 * Implementation of example code given in the tutorial for cparamlib.
 *
 * $Source: /home/nkarlsson/usr/cvsroot/cparamlib/examples/example1.c,v $
 * $Author: niklas $ $Date: 2008/01/17 18:44:58 $ $Revision: 1.6 $
 *
 */

#include <stdio.h>
#include <math.h>
#include "cparamlib/cparamlib.h"

int main(int argc, char* argv[])
{
    double Tp, E;    /* proton kinetic energy and gamma-ray energy */
    double s;        /* non-diffraction cross section */
    int i;
    PARAMSET params; /* struct where parameters are stored */

    printf("Example 1: total inclusive gamma-ray cross section\n");

    memset(&params, 0, sizeof(PARAMSET));

    Tp = 512000.0; /* proton kinetic energy 512 TeV */
    E = 1.0e2;     /* gamma-ray energy 100 GeV */

    gamma_param(Tp, &params);

    /* print out the values of a0,...,a8 and b0,...,b7 */
    for (i = 0; i < 9; i++) {
        printf("a%d=%f ", i, params.a[i]);
    }
    printf("\n");
    for (i = 0; i < 8; i++) {
        printf("b%d=%f ", i, params.b[i]);
    }
    printf("\n");

    s = sigma_incl_tot(ID_GAMMA, E, Tp, &params);

    /* output the results */
    printf("Tp = %.2f GeV\n", Tp);
    printf("E = %.2f GeV\n", E);
    printf("=> dsigma/dlogE = %.1f mb\n", s);

    return;
}
