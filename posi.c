/*
	 posi.c

		Parameter calculation for positrons

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/posi.c,v $
		$Author: niklas $ $Date: 2006/03/19 05:52:13 $ $Revision: 1.4 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/*
		Calculate parameters for positrons from non-diff
*/
void posi_param_nd(double Tp, double* a) {
    double y;
				int i;

    y = log10(Tp*0.001);

				if (Tp > 0.4) {
								a[0] = -0.79606*(y + 3.3) + 7.7496*pow(y + 3.3, 2) - 3.9326*pow(y + 3.3, 3) + 0.80202*pow(y + 3.3, 4) - 0.054994*pow(y + 3.3, 5);
								a[1] = 6.7943e-6 - 3.5345e-6*y + 6.0927e-7*y*y + 2.0219e-7*y*y*y + 5.1005e-8*y*y*y*y - 4.2622e-8*y*y*y*y*y;
								a[2] = 44.827 + 81.378*log10(0.027733*(y + 3.5)) - 1.3886e4/pow(y + 8.4417, 2);
								a[3] = 0.52010 + 0.59336*y + 0.012032*y*y - 0.0064242*y*y*y;
								a[4] = 2.1361*(y + 3.32) + 1.8514*pow(y + 3.32, 2) - 0.47572*pow(y + 3.32, 3) + 0.0032043*pow(y + 3.32, 4) + 0.0082955*pow(y + 3.32, 5);
								a[5] = 1.0845e-6 + 1.4336e-6*log10(0.0077255*(y + 4.3)) + 0.00013018/pow(y + 4.8188, 2) + 9.3601e-8*y;
								a[6] = -267.74 + 14.175*log10(0.35391*(y + 3.4)) + 64.669*y - 7.7036*y*y;
								a[7] = 138.26 - 529.84*log10(0.12467*(y + 3.9)) - 1.9869e4/pow(y + 7.6884, 2) + 1.0675*y*y;
								a[8] = -0.14707 + 0.40135*y + 0.0039899*y*y - 0.0016602*y*y*y;
				} else {
								for (i = 0; i < 9; i++)
												a[i] = 0.0;
				}
}

/*
		Calculate parameters for positrons from diff. dissoc.
*/
void posi_param_diff(double Tp, double* b) {
    double y;
				int i;

    y = log10(Tp*0.001);

				if (Tp > 1.38) {
								if (Tp > 11.0) {
												b[0] = 29.192*tanh(-0.37879*(y + 2.2)) - 3.2196*pow(y + 0.67500, 2) + 3.6687e-3*pow(y + 9.0824, 4);
												b[1] = -142.97 + 147.86*exp(-0.37194*pow((y + 1.8781)/(1.0 + 3.8389*(y + 1.8781)), 2));
												b[2] = -14.487 - 4.2223*tanh(-13.546*(y + 2.2)) + 0.00016988*pow(y + 234.65, 2);
												b[3] = -0.0036974 - 0.41976*exp(-6.1527*pow((y + 1.8194)/(1.0 + 0.99946*(y + 1.8194)), 2));
								} else {
												b[0] = 0.0;
												b[1] = 0.0;
												b[2] = 0.0;
												b[3] = 0.0;
								}
								b[4] = 1.8108 + 0.18545*pow(y + 2.95, 2) - 2.0049e-3*pow(y + 2.95, 4) + 0.85084*exp(-14.987*pow(y + 2.29 - 0.18967*(y + 2.29), 2));
								b[5] = 2.0404 - 0.51548*tanh(2.2758*(y + 1.9)) - 0.035009*(y - 6.6555);
								b[6] = 1.5258 + 1.0132*y - 0.064388*y*y - 0.0040209*y*y*y + 0.0082772*y*y*y*y;
								b[7] = 3.0551 + 3.5240*tanh(-0.36739*(y + 2.1)) - 0.13382*pow(y - 2.7718, 2);
				} else {
								for (i = 0; i < 8; i++)
												b[i] = 0.0;
				}
}

/*
		Calculate parameters for positrons from delta(1232)
*/
void posi_param_delta(double Tp, double* c) {
    double y;
				int i;

    y = log10(Tp*0.001);

    if (Tp <= 0.4 || Tp >= 2.0)
								for (i = 0; i < 5; i++)
												c[i] = 0.0;
				else {
								c[0] = 2.9841*exp(-67.857*pow(((y + 3.1272)/(1.0 + 0.22831*(y + 3.1272))), 2)) - (6.5855 + 9.6984/y - 0.41256*y*y);
								c[1] = 6.8276 + 5.2236*y + 1.4630*y*y;
								c[2] = -6.0291 - 6.4581*tanh(5.0830*(y + 2.1)) + 0.46352*y;
								c[3] = 0.59300 + 0.36093*y;
								c[4] = 0.77368 + 0.44776*y + 0.056409*y*y;
				}
}

/*
		Calculate parameters for positrons from res(1600)
*/
void posi_param_res(double Tp, double* d) {
    double y;
				int i;

    y = log10(Tp*0.001);

    if (Tp <= 0.6 || Tp >= 2.9) {
								for (i = 0; i < 5; i++)
												d[i] = 0.0;
				} else {
								d[0] = 1.9186*exp(-56.544*pow(((y + 2.9485)/(1.0 + 1.2892*(y + 2.9485))), 2)) - (0.23720 - 0.041315*y*y);
								d[1] = -4.9866 - 3.1435*y;
								d[2] = -7.0550 - 7.2165*tanh(31.033*(y + 2.1)) + 0.38541*y;
								d[3] = -2.8915 - 2.1495*y - 0.45006*y*y;
								d[4] = -1.2970 - 0.13947*y + 0.41197*y*y + 0.10641*y*y*y;
				}
}
