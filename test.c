/*
		test.c

		Test app for the parameterization model lib

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/test.c,v $
		$Author: niklas $ $Date: 2006/05/31 23:08:00 $ $Revision: 1.15 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/* Example 1: calculate parameter a0,...,a8 used in the non-diff inclusive 
			gamma-ray cross-section formula */
void example1(void) {
				double Tp, E;
				double f;
				double a[9];
				int i;

				Tp = 512000.0; /* proton kinetic energy 512 TeV */
				E = 1.0e2;     /* gamma-ray energy 100 GeV */

				printf("Example 1: non-diff\n");
				gamma_param_nd(Tp, a);
				for (i = 0; i < 9; i++)
								printf("a[%d] = %10e\n", i, a[i]);
				f = sigma_nd(ID_GAMMA, E, Tp, a);
				printf("Tp = 512TeV and Egamma = 100 GeV => sigma_nd = %10e mb\n", f);

				return;
}

/* Example 2: inclusive gamma-ray cross-section */
void example2(void) {
				double Tp, E;
				double f;
				double a[9];
				int i;

				Tp = 512000.0; /* proton kinetic energy 512 TeV */
				E = 1.0e2;     /* gamma-ray energy 100 GeV */

				printf("Example 2: total inclusive gamma-ray cross-section\n");
				f = sigma(ID_GAMMA, E, Tp);
				printf("total inclusive gamma-ray cross-section\n");
				printf("Tp = 512TeV, Egamma = 100 GeV => sigma = %10e mb\n", f);

				return;
}

/* Example 3: gamma-ray spectrum due to power-law protons */
void example3(void) {
				double *spectrum;
				double a[9], b[8], c[5], d[5];
				double Tp, E;
				double s, s_nd, s_diff, s_delta, s_res;
				double ind2factor;

				FILE *file;

				int i, j;

				printf("Example 3: gamma-ray spectrum from power-law protons\n");

				/* allocate memory for the spectrum (180 bins) */
				spectrum = (double*)malloc(180 * sizeof(double));
				memset(spectrum, 0, 180 * sizeof(double));

				for (i = 0; i < 41; i++) {
								/* proton kinetic energy in GeV */
								Tp = 1000.0*pow(2.0, (i - 22)/2.0);
								/* calculate parameters for this Tp */
								gamma_param_nd(Tp, a);
								gamma_param_diff(Tp, b);
								gamma_param_delta(Tp, c);
								gamma_param_res(Tp, d);
								
								ind2factor = 1.0/(Tp*1.0e-3);

								/* calculate the inclusive cross section in each bin of Egamma */
								for (j = 0; j < 180; j++) {
												s_nd = s_diff = s_delta = s_res = 0.0;
												E = pow(10.0, j*0.05 - 3.0);
												/* calculate individual contributions */
												s_nd = sigma_nd(ID_GAMMA, E, Tp, a);
												if (i > 4)
																s_diff = sigma_diff(ID_GAMMA, E, Tp, b);
												if (i < 7)
																s_delta = sigma_diff(ID_GAMMA, E, Tp, c);
												if (i > 0 && i < 7)
																s_res = sigma_diff(ID_GAMMA, E, Tp, d);
												/* and add them together and add to rest of the spectrum */
												s = s_nd + s_diff + s_delta + s_res;
												spectrum[j] += s*ind2factor;
								}
				}

				/* save to a file */
				file = fopen("gammaspectrum.csv", "w");
				fprintf(file, "gamma-ray spectrum due to power-law proton index 2.0\n");
				fprintf(file, "logE\tspectrum\n");
				for (i = 0; i < 180; i++) {
								s = log10(spectrum[i] + 1.0e-12) + (i*0.05 - 3.0);
								fprintf(file, "%e %e\n", i*0.05 - 3.0, s);
				}
				fclose(file);

				/* free allocated memory */
				free(spectrum);

				return;
}

int main(void) {
				example1();
				/*	example2();
							example3();*/

				return 0;
}

