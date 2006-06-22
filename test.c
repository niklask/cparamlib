/*
		test.c

		Test app for the parameterization model lib

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/test.c,v $
		$Author: niklas $ $Date: 2006/06/22 21:24:14 $ $Revision: 1.16 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

typedef void (*PARAM_FUNC)(double, double*);
static PARAM_FUNC paramfunc_table[7][4] = {{&gamma_param_nd, &gamma_param_diff, &gamma_param_delta, &gamma_param_res},
																																											{&elec_param_nd, &elec_param_diff, &elec_param_delta, &elec_param_res},
																																											{&posi_param_nd, &posi_param_diff, &posi_param_delta, &posi_param_res},
																																											{&nue_param_nd, &nue_param_diff, &nue_param_delta, &nue_param_res},
																																											{&numu_param_nd, &numu_param_diff, &numu_param_delta, &numu_param_res},
																																											{&antinue_param_nd, &antinue_param_diff, &antinue_param_delta, &antinue_param_res},
																																											{&antinumu_param_nd, &antinumu_param_diff, &antinumu_param_delta, &antinumu_param_res}};

char *filenames[7] = {"gammaspectrum.csv", "elecspectrum.csv", "posispectrum.csv", "nuespectrum.csv", "numuspectrum.csv", "antinuespectrum.csv", "antinumuspectrum.csv"};

double Tp[43] = {512.0e3, 362.0e3, 256.0e3, 181.0e3, 128.0e3, 90.5e3, 64.0e3, 45.3e3, 32.0e3, 22.6e3, 
																	16.0e3, 11.3e3, 8.0e3, 5.66e3, 4.0e3, 2.8e3, 2.0e3, 1.41e3, 1.0e3, 707.0, 500.0, 354.0,
																	250.0, 177.0, 125.0, 88.4, 62.5, 44.2, 31.3, 22.1, 15.6, 11.1, 7.81, 5.52, 3.91, 2.76,
																	1.95, 1.38, 0.98, 0.82, 0.69, 0.52, 0.488};

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

/* Example 3: secondary particle spectra due to power-law protons */
void example3(void) {
				double *spectrum;
				double *spectrum_nd;
				double *spectrum_diff;
				double *spectrum_delta;
				double *spectrum_res;
				double a[9], b[8], c[5], d[5];
				double E;
				double s, s_nd, s_diff, s_delta, s_res;
				double ind2factor;
				double widthfactor;

				FILE *file;

				int i, j, k;

				printf("Example 3: secondary particle spectra from power-law protons\n");

				/* allocate memory for the spectrum (180 bins) */
				spectrum = (double*)malloc(180 * sizeof(double));
				spectrum_nd = (double*)malloc(180 * sizeof(double));
				spectrum_diff = (double*)malloc(180 * sizeof(double));
				spectrum_delta = (double*)malloc(180 * sizeof(double));
				spectrum_res = (double*)malloc(180 * sizeof(double));

				for (k = 0; k < 7; k++) {
								memset(spectrum, 0, 180 * sizeof(double));
								memset(spectrum_nd, 0, 180 * sizeof(double));
								memset(spectrum_diff, 0, 180 * sizeof(double));
								memset(spectrum_delta, 0, 180 * sizeof(double));
								memset(spectrum_res, 0, 180 * sizeof(double));

								for (i = 0; i < 43; i++) {
												/* calculate parameters for this Tp */
												paramfunc_table[k][0](Tp[i], a);
												paramfunc_table[k][1](Tp[i], b);
												paramfunc_table[k][2](Tp[i], c);
												paramfunc_table[k][3](Tp[i], d);

												ind2factor = 1.0/(Tp[i]*1.0e-3);
												if (i < 38)
																widthfactor = 1.0;
												else if (i == 38)
																widthfactor = 0.75;
												else if (i > 38)
																widthfactor = 0.5;

												/* calculate the inclusive cross section in each bin of Enue */
												for (j = 0; j < 180; j++) {
																s_nd = s_diff = s_delta = s_res = 0.0;
																E = pow(10.0, j*0.05 - 3.0);
																/* calculate individual contributions */
																s_nd = sigma_nd(k, E, Tp[i], a);

																if (i < 37)
																				s_diff = sigma_diff(k, E, Tp[i], b);
																if (i > 35)
																				s_delta = sigma_delta(k, E, Tp[i], c);
																if (i > 34 && i < 41)
																				s_res = sigma_res(k, E, Tp[i], d);
																/* and add them together and add to rest of the spectrum */
																s = s_nd + s_diff + s_delta + s_res;

																/* store in spectrum for index 2 */
																spectrum[j] += s*ind2factor*widthfactor;
																spectrum_nd[j] += s_nd*ind2factor*widthfactor;
																spectrum_diff[j] += s_diff*ind2factor*widthfactor;
																spectrum_delta[j] += s_delta*ind2factor*widthfactor;
																spectrum_res[j] += s_res*ind2factor*widthfactor;
												}
								}

								/* save to a file */
								file = fopen(filenames[k], "w");
								fprintf(file, "#spectrum due to power-law proton index 2.0\n");
								fprintf(file, "logE\ttot\tnd\tdiff\tdelta\tres\n");
								for (i = 0; i < 180; i++) {
												s = log10(spectrum[i] + 1.0e-12) + (i*0.05 - 3.0);
												s_nd = log10(spectrum_nd[i] + 1.0e-12) + (i*0.05 - 3.0);
												s_diff = log10(spectrum_diff[i] + 1.0e-12) + (i*0.05 - 3.0);
												s_delta = log10(spectrum_delta[i] + 1.0e-12) + (i*0.05 - 3.0);
												s_res = log10(spectrum_res[i] + 1.0e-12) + (i*0.05 - 3.0);
												fprintf(file, "%e %e %e %e %e %e\n", i*0.05 - 3.0, s, s_nd, s_diff, s_delta, s_res);
								}
								fclose(file);
								
				}

				/* free allocated memory */
				free(spectrum);
				free(spectrum_nd);
				free(spectrum_diff);
				free(spectrum_delta);
				free(spectrum_res);

				return;
}

int main(void) {
				example1();
				example2();
				example3();

				return 0;
}

