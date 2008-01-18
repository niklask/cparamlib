/*
 * example2.c
 *
 * Implementation of example code given in the tutorial for cparamlib.
 *
 * $Source: /home/nkarlsson/usr/cvsroot/cparamlib/examples/example2.c,v $
 * $Author: niklas $ $Date: 2008/01/18 23:25:49 $ $Revision: 1.1 $
 *
 */

#include <stdio.h>
#include <math.h>
#include "cparamlib/cparamlib.h"

int main(int argc, char* argv[])
{
    double Tp, E;       /* proton kinetic energy and gamma-ray energy */
    double pT;          /* transverse momentum of gamma ray */
    double* pTdist;     /* pointer to a double array for the pT distribution */
    int i;
    PARAMSET_PT params; /* struct where parameters are stored */
    FILE *file;

    printf("Example 2: pT distribution of gamma rays\n");
    printf("  Tp = %.2f GeV\n", Tp);
    printf("  E = %.2f GeV\n", E);


    /* bin the pT distribution in 250 bins from 0 GeV/c to 2.5 GeV/c,
       i.e 10 MeV/c wide bins */
    pTdist = (double*)malloc(250 * sizeof(double));    
    memset(pTdist, 0, 250 * sizeof(double));

    memset(&params, 0, sizeof(PARAMSET_PT));

    Tp = 64000.0;  /* proton kinetic energy 64 TeV */
    E = 1.0e2;     /* gamma-ray energy 100 GeV */

    gamma_pt_param(E, Tp, &params, 1);

    for (i = 0; i < 250; i++) {
        pT = 0.01*i;
        pTdist[i] = sigma_pt_tot(ID_GAMMA, pT, E, Tp, &params);
    }

    /* save to file */
    printf("  pT distribution saved in gamma_pt_dist.csv\n");

    file = fopen("gamma_pt_dist.csv", "w");
    fprintf(file, "#pT distribution of gamma rays for Tp=%.2f and Egamma=%.2f\n", Tp, E);
    for (i = 0; i < 250; i++) {
        pT = 0.01*i;
        fprintf(file, "%e %e\n", pT, pTdist[i]);
    }
    fclose(file);

    /* free allocated memory */
    free(pTdist);

    return;
}
