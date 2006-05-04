/*
	 elec.c

		Parameter calculation for electrons

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/elec.c,v $
		$Author: niklas $ $Date: 2006/05/04 20:44:33 $ $Revision: 1.6 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/*
		Calculate parameters for electrons from non-diff
*/
void elec_param_nd(double Tp, double* a) {
    double y;
				int i;

    y = log10(Tp*0.001);

				if (Tp > 0.4) {
								a[0] = -0.018639*(y + 3.3) + 2.4315*pow(y + 3.3, 2) - 0.57719*pow(y + 3.3, 3) + 0.063435*pow(y + 3.3, 4);
								a[1] = 7.1827e-6 - 3.5067e-6*y + 1.3264e-6*y*y - 3.3481e-7*y*y*y + 2.3551e-8*y*y*y*y + 3.4297e-9*y*y*y*y*y;
								a[2] = 563.91 - 362.18*log10(2.7187*(y + 3.4)) - 2.8924e4/pow(y + 7.9031, 2);
								a[3] = 0.52684 + 0.57717*y + 0.0045336*y*y - 0.0089066*y*y*y;
								a[4] = 0.36108*(y + 3.32) + 1.6963*pow(y + 3.32, 2) - 0.074456*pow(y + 3.32, 3) - 0.071455*pow(y + 3.32, 4) + 0.010473*pow(y + 3.32, 5);
								a[5] = 9.7387e-5 + 7.8573e-5*log10(0.0036055*(y + 4.3)) + 0.00024660/(y + 4.9390) - 3.8097e-7*y*y;
								a[6] = -273.00 - 106.22*log10(0.34100*(y + 3.4)) + 89.037*y - 12.546*y*y;
								a[7] = 432.53 - 883.99*log10(0.19737*(y + 3.9)) - 4.1938e4/pow(y + 8.5518, 2);
								a[8] = -0.12756 + 0.43478*y - 0.0027797*y*y -0.0083074*y*y*y;
				} else {
								for (i = 0; i < 9; i++)
												a[i] = 0.0;
				}
}

/*
		Calculate parameters for electrons from diff. dissoc.
*/
void elec_param_diff(double Tp, double* b) {
    double y;
				int i;

    y = log10(Tp*0.001);

				if (Tp > 1.38) {
								if (Tp > 5.52) {
												b[0] = 0.20463*tanh(-6.2370*(y + 2.2)) - 0.16362*pow(y + 1.6878, 2) + 3.5183e-4*pow(y + 9.6400, 4);
												b[1] = 1.6537 + 3.8530*exp(-3.2027*pow((y + 2.0154)/(1.0 + 0.62779*(y + 2.0154)), 2));
												b[2] = -10.722 + 0.082672*tanh(1.8879*(y + 2.1)) + 0.00014895*pow(y + 256.63, 2);
												b[3] = -0.023752 - 0.51734*exp(-3.3087*pow((y + 1.9877)/(1.0 + 0.40300*(y + 1.988)), 2));
								} else {
												b[0] = 0.0;
												b[1] = 0.0;
												b[2] = 0.0;
												b[3] = 0.0;
								}
								b[4] = 0.94921 + 0.12280*pow(y + 2.9, 2) - 7.1585e-4*pow(y + 2.9, 4) + 0.52130*log10(y + 2.9);
								b[5] = -4.2295 - 1.0025*tanh(9.0733*(y + 1.9)) - 0.11452*(y - 62.382);
								b[6] = 1.4862 + 0.99544*y - 0.042763*y*y - 0.0040065*y*y*y + 0.0057987*y*y*y*y;
								b[7] = 6.2629 + 6.9517*tanh(-0.36480*(y + 2.1)) - 0.26033*pow(y - 2.8542, 2);
				} else {
								for (i = 0; i < 8; i++)
												b[i] = 0.0;
				}
}

/*
		Calculate parameters for electrons from delta(1232)
*/
void elec_param_delta(double Tp, double* c) {
				c[0] = c[1] = c[2] = c[3] = c[4] = 0.0;
}

/*
		Calculate parameters for electrons from res(1600)
*/
void elec_param_res(double Tp, double* d) {
    double y;
				int i;

    y = log10(Tp*0.001);

    if (Tp <= 0.6 || Tp >= 2.9) {
								for (i = 0; i < 5; i++)
												d[i] = 0.0;
				} else {
								d[0] = 0.37790*exp(-56.826*pow((y + 2.9537)/(1.0 + 1.5221*(y + 2.9537)), 2)) - 0.059458 + 0.0096583*y*y;
								d[1] = -5.5135 - 3.3988*y;
								d[2] = -7.1209 - 7.1850*tanh(30.801*(y + 2.1)) + 0.35108*y;
								d[3] = -6.7841 - 4.8385*y - 0.91523*y*y;
								d[4] = -134.03 - 139.63*y - 48.316*y*y - 5.5526*y*y*y;
				}
}
