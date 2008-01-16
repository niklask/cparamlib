/*
 * example1.c
 *
 * Implementation of example code given in the tutorial for cparamlib.
 *
 * $Source: /home/nkarlsson/usr/cvsroot/cparamlib/examples/example1.c,v $
 * $Author: niklas $ $Date: 2008/01/16 19:24:01 $ $Revision: 1.2 $
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

    memset(&params, 0, sizeof(PARAMSET));

    Tp = 512000.0; /* proton kinetic energy 512 TeV */
    E = 1.0e2;     /* gamma-ray energy 100 GeV */

    gamma_param(Tp, &params);

    printf("Example 1: total inclusive gamma-ray cross section\n");

    s = sigma_incl_tot(ID_GAMMA, E, Tp, &params);

    /* output the results */
    printf("Tp = %.2f GeV\n", Tp);
    printf("E = %.2f GeV\n", E);
    printf("=> sigma_incl = %10 mb\n", s);

    return;
}
