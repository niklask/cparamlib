/*
		params.c

		Test app for the parameterization model lib
		Calculates parameters as functions of Tp

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/params.c,v $
		$Author: niklas $ $Date: 2006/01/26 22:21:16 $ $Revision: 1.1 $
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

void nondiff(void) {
				double Tp, Pp;
				double f;
				double a[9];
				int i, j;
				FILE *fp;

				fp = fopen("GammaDiffDiss.dat", "w");
				printf("results written to file!\n");
				for (i = 0; i < 43; i++) {
								Tp = Tp_list[i]*1000.0;
								Pp = sqrt(Tp*Tp + 2.0*Tp*m_p);
								gamma_param_nd(Pp, a);
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

void diffdiss(void) {
				double Tp, Pp;
				double f;
				double b[8];
				int i, j;
				FILE *fp;

				fp = fopen("GammaDiffDiss.dat", "w");
				printf("results written to file!\n");
				for (i = 0; i < 37; i++) {
								Tp = Tp_list[i]*1000.0;
								Pp = sqrt(Tp*Tp + 2.0*Tp*m_p);
								gamma_param_diff(Pp, b);
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

int main(void) {
				nondiff();
				diffdiss();

				return 0;
}

