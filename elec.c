/*
	 elec.c

		Parameter calculation for electrons

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/elec.c,v $
		$Author: niklas $ $Date: 2006/01/18 20:16:59 $ $Revision: 1.1 $
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
}

/*
		Calculate parameters for electron from delta(1232)
*/
void elec_param_delta(double Pp, double* c) {
				c[0] = c[1] = c[2] = c[3] = c[4] = 0.0
}

/*
		Calculate parameters for gamma-ray from res(1600)
*/
void gamma_param_res(double Pp, double* d) {
    double Etot, Tp, y;
				int i;

				Etot	= sqrt(Pp*Pp + m_p*m_p);
    Tp = Etot - m_p;
    y = log10(Tp*0.001);
}
