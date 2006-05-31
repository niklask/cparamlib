/*
	 antinue.c

		Parameter calculation for electron anti neutrinos

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/antinue.c,v $
		$Author: niklas $ $Date: 2006/05/31 23:08:00 $ $Revision: 1.7 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/*
		Calculate parameters for electron anti neutrino from non-diff
*/
void antinue_param_nd(double Tp, double* a) {
    double y;
				int i;

    y = log10(Tp*0.001);

				if (Tp > 0.4) {
								a[0] = 0.0013113 + 0.36538*(y + 3.31) + 1.5178*pow(y + 3.31, 2) - 0.20668*pow(y + 3.31, 3) + 0.024255*pow(y + 3.31, 4);
								a[1] = -4.7833e-6 + 4.5837e-5*exp(-0.42980*(y + 3.4)) + 6.1559e-6/(y + 4.1731) + 1.1928e-6*y;
								a[2] = -245.22 + 73.223*y - 19.652*y*y + 0.83138*y*y*y + 0.71564*y*y*y*y;
								a[3] = 0.45232 + 0.52934*y + 0.010078*y*y - 0.0017092*y*y*y;
								a[4] = -0.0025734 + 0.38424*(y + 3.32) + 1.5517*pow(y + 3.32, 2) + 0.17336*pow(y + 3.32, 3) - 0.17160*pow(y + 3.32, 4) + 0.021059*pow(y + 3.32, 5);
								a[5] = 4.7673e-5 + 5.4936e-5*log10(0.0067905*(y + 4.3)) + 0.00020740/(y + 4.9772);
								a[6] = -270.30 - 114.47*log10(0.34352*(y + 3.4)) + 80.085*y - 7.9240*y*y;
								a[7] = 3271.9 - 2.9161e5/(y + 87.847) - 6.2330*y*y;
								a[8] = -0.17787 + 0.36771*y - 0.025397*y*y + 0.0019238*y*y*y + 0.0032725*y*y*y*y;
				} else {
								for (i = 0; i < 9; i++)
												a[i] = 0.0;
				}
}

/*
		Calculate parameters for electron anti neutrino from diff. dissoc.
*/
void antinue_param_diff(double Tp, double* b) {
    double y;
				int i;

    y = log10(Tp*0.001);

				if (Tp > 1.38) {
								if (Tp > 11.0) {
												b[0] = 41.307*tanh(-0.37411*(y + 2.2)) - 4.1223*pow(y + 0.55505, 2) + 0.0042652*pow(y + 9.2685, 4);
												b[1] = -132.50 + 142.12*exp(-8.0289*pow((y + 1.9196)/(1.0 + 11.530*(y + 1.9196)), 2));
												b[2] = -17.223 + 0.011285*tanh(69.746*(y + 1.9)) - 0.048233*y + 0.00025881*pow(y + 250.77, 2);
												b[3] = 8.1991 - 9.6437*exp(-45.261*pow((y + 1.9292)/(1.0 + 16.682*(y + 1.9292)), 2));
								} else {
												b[0] = 0.0;
												b[1] = 0.0;
												b[2] = 0.0;
												b[3] = 0.0;
								}
								b[4] = 0.55919 + 0.36647*(y + 2.95) + 0.056194*pow(y + 2.95, 2) + 0.49957*exp(-5.5317*pow(y + 2.2 + 0.43867*pow(y + 2.2, 2), 2));
								b[5] = 1.2544 - 0.52362*tanh(2.7638*(y + 1.9)) - 0.055837*(y - 17.638);
								b[6] = 1.4788 + 1.0278*y - 0.092852*y*y - 0.0062734*y*y*y + 0.011920*y*y*y*y;
								b[7] = 5.1651 + 5.7398*tanh(-0.37356*(y + 2.1)) - 0.22234*pow(y - 2.7889, 2);
				} else {
								for (i = 0; i < 8; i++)
												b[i] = 0.0;
				}
}

/*
		Calculate parameters for electron anti neutrino from delta(1232)
*/
void antinue_param_delta(double Tp, double* c) {
				c[0] = c[1] = c[2] = c[3] = c[4] = 0.0;
}

/*
		Calculate parameters for electron anti neutrino from res(1600)
*/
void antinue_param_res(double Tp, double* d) {
    double y;
				int i;

    y = log10(Tp*0.001);

    if (Tp <= 0.6 || Tp >= 2.9) {
								for (i = 0; i < 5; i++)
												d[i] = 0.0;
				} else {
								d[0] = 0.36459*exp(-58.210*pow((y + 2.9537)/(1.0 + 1.4320*(y + 2.9537)), 2)) - (0.11283 + 0.046244*y);
								d[1] = -9.5066 - 5.4655*y - 0.31769*y*y;
								d[2] = -7.1831 - 7.1551*tanh(30.354*(y + 2.1)) + 0.33757*y;
								d[3] = 2.7938 + 1.6992*y + 0.20161*y*y;
								d[4] = 0.61878 + 0.62371*y + 0.18913*y*y + 0.019118*y*y*y;
				}
}
