/* spectrum.c
 *
 * Implementation of the second example in the cparamlib tutorial. This program
 * calculates the spectrum of gamma rays due to a pencil beam of power-law
 * protons.
 *
 * $Source: /home/nkarlsson/usr/cvsroot/cparamlib/examples/spectrum.c,v $
 * $Author: niklas $ $Date: 2008/01/17 19:31:57 $ $Revision: 1.4 $
 *
 */

#include <stdio.h>
#include <math.h>
#include "cparamlib/cparamlib.h"

/* preset bins of Tp */
double Tp[43] = {512.0e3, 362.0e3, 256.0e3, 181.0e3, 128.0e3, 90.5e3, 64.0e3, 45.3e3,
                 32.0e3, 22.6e3, 16.0e3, 11.3e3, 8.0e3, 5.66e3, 4.0e3, 2.8e3, 2.0e3,
                 1.41e3, 1.0e3, 707.0, 500.0, 354.0, 250.0, 177.0, 125.0, 88.4, 62.5,
                 44.2, 31.3, 22.1, 15.6, 11.1, 7.81, 5.52, 3.91, 2.76, 1.95, 1.38,
                 0.98, 0.82, 0.69, 0.58, 0.488};

/* the following double arrays are used to store the calculated spectra */
double *spectrum;
double *spectrum_nd;
double *spectrum_diff;
double *spectrum_delta;
double *spectrum_res;

int main(int argc, char* argv[])
{
    double E;
    double Jp, pl_index;
    double s, s_nd, s_diff, s_delta, s_res;
    double width;
    int i, j;
    PARAMSET params;
    FILE *file;

    /* allocate memory for the spectrum (180 bins) */
    spectrum = (double*)malloc(180 * sizeof(double));
    spectrum_nd = (double*)malloc(180 * sizeof(double));
    spectrum_diff = (double*)malloc(180 * sizeof(double));
    spectrum_delta = (double*)malloc(180 * sizeof(double));
    spectrum_res = (double*)malloc(180 * sizeof(double));

    /* make sure they are empty */
    memset(spectrum, 0, 180 * sizeof(double));
    memset(spectrum_nd, 0, 180 * sizeof(double));
    memset(spectrum_diff, 0, 180 * sizeof(double));
    memset(spectrum_delta, 0, 180 * sizeof(double));
    memset(spectrum_res, 0, 180 * sizeof(double));

    /* make sure the parameter struct is empty */
    memset(&params, 0, sizeof(PARAMSET));

    /* pl_index is the power-law index of the proton spectrum */
    pl_index = 2.0;

    for (i = 0; i < 43; i++) {
        /* calculate parameters for this Tp */
        gamma_param(Tp[i], &params);

        /* set the width of the bin, default is 1 but increased sampling
           near pion production threshold requires reduced bin widths */
        if (i < 38)
            width = 1.0;
        else if (i == 38)
            width = 0.75;
        else if (i > 38)
            width = 0.5;

        /* calculate the proton spectrum (power law) at given Tp [TeV] */
        Jp = pow(Tp[i]*1.0e-3, -(pl_index - 1));
           
        /* calculate the inclusive cross section in each bin of Egamma */
        for (j = 0; j < 180; j++) {
            /* the gamma-ray energy is taken in 180 bins, from 1 MeV to 1000 TeV */
            E = pow(10.0, j*0.05 - 3.0);

            /* calculate individual contributions */
            s_nd = sigma_incl_nd(ID_GAMMA, E, Tp[i], &params);
            s_diff = sigma_incl_diff(ID_GAMMA, E, Tp[i], &params);
            s_delta = sigma_incl_delta(ID_GAMMA, E, Tp[i], &params);
            s_res = sigma_incl_res(ID_GAMMA, E, Tp[i], &params);
            
            /* and add them together and add to rest of the spectrum */
            s = s_nd + s_diff + s_delta + s_res;

            /* store in spectrum for index 2 */
            spectrum[j] += s*Jp*width;
            spectrum_nd[j] += s_nd*Jp*width;
            spectrum_diff[j] += s_diff*Jp*width;
            spectrum_delta[j] += s_delta*Jp*width;
            spectrum_res[j] += s_res*Jp*width;
        }
    }
    
    /* save to a file */
    file = fopen("gamma_spectrum.csv", "w");
    fprintf(file, "#spectrum due to power-law proton index %.2f\n", pl_index);
    for (i = 0; i < 180; i++) {
        s = log10(spectrum[i] + 1.0e-12) + (i*0.05 - 3.0);
        s_nd = log10(spectrum_nd[i] + 1.0e-12) + (i*0.05 - 3.0);
        s_diff = log10(spectrum_diff[i] + 1.0e-12) + (i*0.05 - 3.0);
        s_delta = log10(spectrum_delta[i] + 1.0e-12) + (i*0.05 - 3.0);
        s_res = log10(spectrum_res[i] + 1.0e-12) + (i*0.05 - 3.0);
        fprintf(file, "%e %e %e %e %e %e\n", i*0.05 - 3.0, s, s_nd, s_diff, s_delta, s_res);
    }
    fclose(file);

    /* free allocated memory */
    free(spectrum);
    free(spectrum_nd);
    free(spectrum_diff);
    free(spectrum_delta);
    free(spectrum_res);

    return;
}
