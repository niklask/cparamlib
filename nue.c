/*
	 nue.c

		Parameter calculation for electron neutrinos

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/nue.c,v $
		$Author: niklas $ $Date: 2006/01/31 19:40:51 $ $Revision: 1.1 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/*
		Calculate parameters for electron neutrino from non-diff
*/
void nue_param_nd(double Pp, double* a) {
    double Etot, Tp, y;
				int i;

				Etot	= sqrt(Pp*Pp + m_p*m_p);
    Tp = Etot - m_p;
    y = log10(Tp*0.001);

				/* These are the original parameter values
							note that a0 and a4 have been multiplied with 1/(10^0.05 - 1) for unit conversion
				*/
				if (Tp > 0.4) {
								a[0] =
								a[1] =
								a[2] =
								a[3] =
								a[4] =
								a[5] =
								a[6] =
								a[7] =
								a[8] =
				} else {
								for (i = 0; i < 9; i++)
												a[i] = 0.0;
				}
}

/*
		Calculate parameters for electron neutrino from diff. dissoc.
*/
void nue_param_diff(double Pp, double* b) {
    double Etot, Tp, y;
				int i;

				Etot	= sqrt(Pp*Pp + m_p*m_p);
    Tp = Etot - m_p;
    y = log10(Tp*0.001);

				/* These are the original parameter values
							note that b0 and b4 have been multiplied with 1/(10^0.05 - 1) for unit conversion
				*/

				if (Tp > 1.38) {
								if (Tp > 11.0) {
												b[0] =
												b[1] =
												b[2] =
												b[3] =
								} else {
												b[0] = 0.0;
												b[1] = 0.0;
												b[2] = 0.0;
												b[3] = 0.0;
								}
								b[4] =
								b[5] =
								b[6] =
								b[7] =
				} else {
								for (i = 0; i < 8; i++)
												b[i] = 0.0;
				}
}

/*
		Calculate parameters for electron neutrino from delta(1232)
*/
void nue_param_delta(double Pp, double* c) {
    double Etot, Tp, y;
				int i;

				Etot	= sqrt(Pp*Pp + m_p*m_p);
    Tp = Etot - m_p;
    y = log10(Tp*0.001);

				/* These are the original parameter values
							note that d0 have been multiplied with 1/(10^0.05 - 1) for unit conversion
				*/

    if (Tp <= 0.4 || Tp >= 2.0) {
								for (i = 0; i < 5; i++)
												c[i] = 0.0;
				} else {
								c[0] =
								c[1] =
								c[2] =
								c[3] =
								c[4] =
				}
}

/*
		Calculate parameters for electron neutrino from res(1600)
*/
void nue_param_res(double Pp, double* d) {
    double Etot, Tp, y;
				int i;

				Etot	= sqrt(Pp*Pp + m_p*m_p);
    Tp = Etot - m_p;
    y = log10(Tp*0.001);

				/* These are the original parameter values
							note that d0 have been multiplied with 1/(10^0.05 - 1) for unit conversion
				*/

    if (Tp <= 0.6 || Tp >= 2.9) {
								for (i = 0; i < 5; i++)
												d[i] = 0.0;
				} else {
								d[0] =
								d[1] =
								d[2] =
								d[3] =
								d[4] =
				}
}
