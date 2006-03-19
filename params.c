/*
		params.c

		Test app for the parameterization model lib
		Calculates parameters as functions of Tp

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/params.c,v $
		$Author: niklas $ $Date: 2006/03/19 05:52:13 $ $Revision: 1.10 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

double Tp_list[43] = {512., 362., 256., 181., 128., 90.5, 64., 45.3, 32.,
																						22.6, 16., 11.3, 8., 5.66, 4., 2.8, 2., 1.41, 1.,
																						0.707, 0.5, 0.354, 0.25, 0.177, 0.125, 88.4e-3,
																						62.5e-3, 44.2e-3, 31.3e-3, 22.1e-3, 15.6e-3,
																						11.1e-3, 7.81e-3, 5.52e-3, 3.91e-3, 2.76e-3,
																						1.95e-3, 1.38e-3, 0.98e-3, 0.82e-3, 0.69e-3,
																						0.58e-3, 0.488e-3};

char *filenames[7][4] = {{"GammaNonDiff.dat", "GammaDiffDiss.dat", "GammaDelta1232.dat", "GammaDelta1600.dat"},
																									{"ElecNonDiff.dat", "ElecDiffDiss.dat", "ElecDelta1232.dat", "ElecDelta1600.dat"},
																									{"PosiNonDiff.dat", "PosiDiffDiss.dat", "PosiDelta1232.dat", "PosiDelta1600.dat"},
																									{"NuENonDiff.dat", "NuEDiffDiss.dat", "NuEDelta1232.dat", "NuEDelta1600.dat"},
																									{"NuMuNonDiff.dat", "NuMuDiffDiss.dat", "NuMuDelta1232.dat", "NuMuDelta1600.dat"},
																									{"AntiNuENonDiff.dat", "AntiNuEDiffDiss.dat", "AntiNuEDelta1232.dat", "AntiNuEDelta1600.dat"},
																									{"AntiNuMuNonDiff.dat", "AntiNuMuDiffDiss.dat", "AntiNuMuDelta1232.dat", "AntiNuMuDelta1600.dat"}};

void nondiff(int p) {
				double Tp;
				double f;
				double a[9];
				int i, j;
				FILE *fp;

				fp = fopen(filenames[p][0], "w");
				for (i = 0; i < 43; i++) {
								Tp = Tp_list[i]*1000.0;
								switch (p) {
												case 0: gamma_param_nd(Tp, a);
												case 1: elec_param_nd(Tp, a);
												case 2: posi_param_nd(Tp, a);
												case 3: nue_param_nd(Tp, a);
												case 4: numu_param_nd(Tp, a);
												case 5: antinue_param_nd(Tp, a);
												case 6: antinumu_param_nd(Tp, a);
								}
								fprintf(fp, "%f ", Tp_list[i]);
								for (j = 0; j < 9; j++) {
												fprintf(fp, "%10e", a[j]);
												if (j < 8)
																fprintf(fp, " ");
												else
																fprintf(fp, "\n");
								}
				}
				fclose(fp);
}

void diffdiss(int p) {
				double Tp;
				double f;
				double b[8];
				int i, j;
				FILE *fp;

				fp = fopen(filenames[p][1], "w");
				for (i = 0; i < 37; i++) {
								Tp = Tp_list[i]*1000.0;
								switch (p) {
												case 0: gamma_param_diff(Tp, b);
												case 1: elec_param_diff(Tp, b);
												case 2: posi_param_diff(Tp, b);
												case 3: nue_param_diff(Tp, b);
												case 4: numu_param_diff(Tp, b);
												case 5: antinue_param_diff(Tp, b);
												case 6: antinumu_param_diff(Tp, b);
								}
								fprintf(fp, "%f ", Tp_list[i]);
								for (j = 0; j < 8; j++) {
												fprintf(fp, "%10e", b[j]);
												if (j < 7)
																fprintf(fp, " ");
												else
																fprintf(fp, "\n");
								}
				}
				fclose(fp);
}

void delta1232(int p) {
				double Tp;
				double f;
				double c[5];
				int i, j;
				FILE *fp;

				fp = fopen(filenames[p][2], "w");
				for (i = 36; i < 43; i++) {
								Tp = Tp_list[i]*1000.0;
								switch (p) {
												case 0: gamma_param_delta(Tp, c);
												case 1: elec_param_delta(Tp, c);
												case 2: posi_param_delta(Tp, c);
												case 3: nue_param_delta(Tp, c);
												case 4: numu_param_delta(Tp, c);
												case 5: antinue_param_delta(Tp, c);
												case 6: antinumu_param_delta(Tp, c);
								}
								fprintf(fp, "%f ", Tp_list[i]);
								for (j = 0; j < 5; j++) {
												fprintf(fp, "%10e", c[j]);
												if (j < 4)
																fprintf(fp, " ");
												else
																fprintf(fp, "\n");
								}
				}
				fclose(fp);
}

void delta1600(int p) {
				double Tp;
				double f;
				double d[5];
				int i, j;
				FILE *fp;

				fp = fopen(filenames[p][3], "w");
				for (i = 35; i < 41; i++) {
								Tp = Tp_list[i]*1000.0;
								switch (p) {
												case 0: gamma_param_res(Tp, d);
												case 1: elec_param_res(Tp, d);
												case 2: posi_param_res(Tp, d);
												case 3: nue_param_res(Tp, d);
												case 4: numu_param_res(Tp, d);
												case 5: antinue_param_res(Tp, d);
												case 6: antinumu_param_res(Tp, d);
								}
								fprintf(fp, "%f ", Tp_list[i]);
								for (j = 0; j < 5; j++) {
												fprintf(fp, "%10e", d[j]);
												if (j < 4)
																fprintf(fp, " ");
												else
																fprintf(fp, "\n");
								}
				}
				fclose(fp);
}

int main(void) {
				int i;

				for (i = 0; i <= 6; i++) {
								nondiff(i);
								diffdiss(i);
								delta1232(i);
								delta1600(i);
				}

				return 0;
}

