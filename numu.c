/*
	 numu.c

		Parameter calculation for electron neutrinos

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/numu.c,v $
		$Author: niklas $ $Date: 2006/05/31 23:08:00 $ $Revision: 1.6 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/*
		Calculate parameters for muon neutrino from non-diff
*/
void numu_param_nd(double Tp, double* a) {
    double y;
				int i;

    y = log10(Tp*0.001);

				if (Tp > 0.4) {
								a[0] = -0.63611*(y + 3.3) + 9.9015*pow(y + 3.3, 2) - 4.5897*pow(y + 3.3, 3) + 0.91778*pow(y + 3.3, 4) - 0.060724*pow(y + 3.3, 5);
								a[1] = 6.8700e-6 - 2.8245e-6*y + 7.6032e-7*y*y - 3.2953e-7*y*y*y + 7.4292e-8*y*y*y*y;
								a[2] = -240.46 + 58.405*y - 9.8556*y*y + 3.1401*y*y*y - 0.88932*y*y*y*y;
								a[3] = 0.49935 + 0.60919*y + 0.0024963*y*y - 0.0099910*y*y*y;
								a[4] = 2.5094*(y + 3.32) + 4.1350*pow(y + 3.32, 2) - 0.89534*pow(y + 3.32, 3) - 2.7577e-3*pow(y + 3.32, 4) + 0.014511*pow(y + 3.32, 5);
								a[5] = 8.2046e-7 + 1.4085e-6*log10(0.016793*(y + 4.3)) + 0.00013340/pow(y + 4.7136, 2);
								a[6] = -267.55 - 0.21018*log10(0.35217*(y + 3.4)) + 69.586*y - 9.9930*y*y;
								a[7] = 2741.8 + 222.01*log10(9.7401e-13*(y + 3.9)) - 4772.5/(y + 19.773) - 6.1001*y*y;
								a[8] = -0.11857 + 0.39072*y - 0.037813*y*y + 0.0022265*y*y*y + 0.0046931*y*y*y*y;
				} else {
								for (i = 0; i < 9; i++)
												a[i] = 0.0;
				}
}

/*
		Calculate parameters for muon neutrino from diff. dissoc.
*/
void numu_param_diff(double Tp, double* b) {
    double y;
				int i;

    y = log10(Tp*0.001);

				if (Tp > 1.38) {
								if (Tp > 11.0) {
												b[0] = 64.682*tanh(-0.34313*(y + 2.2)) - 5.5955*pow(y + 0.44754, 2) + 0.0050117*pow(y + 9.9165, 4);
												b[1] = -7.6016 + 3042.7*exp(-1.0134e4*pow((y + 2.3066)/(1.0 + 41.612*(y + 2.3066)), 2));
												b[2] = -1.4978 - 0.58163*tanh(-0.36488*(y + 1.9)) + 0.031825*(y + 2.8097) + 0.022796*pow(y - 1.8861, 2);
												b[3] = -0.0061483 - 65.799*exp(-4.8239*pow((y + 3.8835)/(1.0 + 0.53343*(y + 3.8835)), 2));
								} else {
												b[0] = 0.0;
												b[1] = 0.0;
												b[2] = 0.0;
												b[3] = 0.0;
								}
								b[4] = 2.8009 + 0.35351*pow(y + 2.95, 2) - 0.0039779*pow(y + 2.95, 4) + 1.3012*exp(-10.592*pow(y + 2.28 - 0.19149*pow(y + 2.28, 2), 2));
								b[5] = 1.8016 - 0.69847*tanh(2.8627*(y + 1.9)) - 0.015722*(y - 45.410);
								b[6] = 1.4617 + 1.0167*y - 0.078617*y*y - 0.0038336*y*y*y + 0.010141*y*y*y*y;
								b[7] = 3.5599 + 4.0041*tanh(-0.41889*(y + 2.1)) - 0.18182*pow(y - 2.4209, 2);
				} else {
								for (i = 0; i < 8; i++)
												b[i] = 0.0;
				}
}

/*
		Calculate parameters for muon neutrino from delta(1232)
*/
void numu_param_delta(double Tp, double* c) {
    double y;
				int i;

    y = log10(Tp*0.001);

    if (Tp <= 0.4 || Tp >= 2.0) {
								for (i = 0; i < 5; i++)
												c[i] = 0.0;
				} else {
								c[0] = 3.6052*exp(-60.914*pow((y + 3.1278)/(1.0 - 0.19497*(y + 3.1278)), 2)) - (0.92514 - 2.1315/y - 0.23548*y*y);
								c[1] = 95.310 + 70.497*y + 13.636*y*y;
								c[2] = -6.2158 - 6.2939*tanh(21.592*(y + 2.1)) + 0.37440*y;
								c[3] = 2.7485 + 1.1692*y;
								c[4] = -2.7568 - 1.8461*y - 0.31376*y*y;
				}
}

/*
		Calculate parameters for muon neutrino from res(1600)
*/
void numu_param_res(double Tp, double* d) {
    double y;
				int i;

    y = log10(Tp*0.001);

    if (Tp <= 0.6 || Tp >= 2.9) {
								for (i = 0; i < 5; i++)
												d[i] = 0.0;
				} else {
								d[0] = 2.5489*exp(-58.488*pow((y + 2.9509)/(1.0 + 1.3154*(y + 2.9509)), 2)) - (0.83039 + 0.34412*y);
								d[1] = 88.173 + 65.148*y + 12.585*y*y;
								d[2] = -7.0962 - 7.1690*tanh(30.890*(y + 2.1)) + 0.38032*y;
								d[3] = -4.1440 - 3.2717*y - 0.70537*y*y;
								d[4] = 2.2624 + 1.1806*y + 0.0043450*y*y - 0.043020*y*y*y;
				}
}
