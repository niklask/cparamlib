/*
	 elec.c

		Parameter calculation for electrons

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/elec.c,v $
		$Author: niklas $ $Date: 2006/01/28 04:17:36 $ $Revision: 1.2 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/*
		Calculate parameters for electron from non-diff
*/
void elec_param_nd(double Pp, double* a) {
    double Etot, Tp, y;
				int i;

				Etot	= sqrt(Pp*Pp + m_p*m_p);
    Tp = Etot - m_p;
    y = log10(Tp*0.001);

				/* These are the original parameter values
							note that a0 and a4 have been multiplied with 1/(10^0.05 - 1) for unit conversion
							-0.00227427738719 0.296685099602 -0.0704276934266 0.00774029362947
							7.18268483979e-06 -3.50668619831e-06 1.32643742745e-06 -3.34813563541e-07 2.35509229896e-08 3.42968631273e-09
							563.909851074 -362.181091309 2.71871185303 -28924.9960938 7.90309762955
							0.526838362217 0.577174305916 0.00453363452107 -0.00890664290637
							0.0440583191812 0.206974059343 -0.00908503122628 -0.00871878303587 0.0012778629316
							9.73868664005e-05 7.85731754149e-05 0.00360554642975 0.000246596493525 4.9389667511 -3.80967208002e-07
							-272.979888916 -106.217224121 0.341000974178 89.0373687744 -12.5462293625
							432.526733398 -883.994812012 0.197373971343 -41937.7460938 8.55176448822
							-0.127555832267 0.434775620699 -0.00277966680005 -0.0083074150607
				*/
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
		Calculate parameters for electron from diff. dissoc.
*/
void elec_param_diff(double Pp, double* b) {
    double Etot, Tp, y;
				int i;

				Etot	= sqrt(Pp*Pp + m_p*m_p);
    Tp = Etot - m_p;
    y = log10(Tp*0.001);

				/* These are the original parameter values
							note that b0 and b4 have been multiplied with 1/(10^0.05 - 1) for unit conversion
							0.0249691866338 -6.23704004288 -0.0199651997536 1.68776345253 4.29302090197e-05 9.64004135132
							1.65368855 3.85295248032 3.2026617527 -2.01541686058 0.627792596817
							-10.7221393585 0.0826715230942 1.88791573048 0.00014895023196 255.631759644
							-0.0237520020455 0.517341017723 3.30872654915 -1.98768746853 0.40300488472
							0.115821219981 0.0149836251512 -8.73467724887e-05 0.06360848248
							-4.22947454453 -1.00249087811 9.07331371307 -0.114516332746 -62.3815078735
							1.4861805439 0.995444893837 -0.0427628308535 -0.00400647893548 0.00579872075468
							6.26293468475 6.95169401169 -0.363795995712 -0.260329335928 -2.85419464111
				*/

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
		Calculate parameters for electron from delta(1232)
*/
void elec_param_delta(double Pp, double* c) {
				c[0] = c[1] = c[2] = c[3] = c[4] = 0.0;
}

/*
		Calculate parameters for gamma-ray from res(1600)
*/
void elec_param_res(double Pp, double* d) {
    double Etot, Tp, y;
				int i;

				Etot	= sqrt(Pp*Pp + m_p*m_p);
    Tp = Etot - m_p;
    y = log10(Tp*0.001);

				/* These are the original parameter values
							note that d0 have been multiplied with 1/(10^0.05 - 1) for unit conversion
							0.046110983938 56.8264007568 -2.95373129845 1.52206134796 0.00725625408813 -0.00117849593516
							-5.51345300674 -3.39880275726
							-7.12091302872 -7.18501472473 30.8012924194 0.351079583168
							-6.78410625458 -4.83850049973 -0.915232002735
							-134.033050537 -139.628967285 -48.3161201477 -5.55262756348
				*/

    if (Tp <= 0.6 || Tp >= 2.9) {
								for (i = 0; i < 5; i++)
												d[i] = 0.0;
				} else {
								d[0] = 0.37790*exp(-56.826*pow((y + 2.9537)/(1.0 + 1.5221*(y + 2.9537)), 2)) - 0.059458 + 0.0096583*y*y;
								d[1] = -5.5135 - 3.3988*y;
								d[2] = -7.1209 - 7.1850*tanh(30.801*(y + 2.1)) + 0.35108*y;
								d[3] = -6.7841 - 4.8485*y - 0.91523*y*y;
								d[4] = -134.03 - 139.63*y - 48.316*y*y - 5.5526*y*y*y;
				}
}
