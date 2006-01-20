/*
		flux.c

		Main part of cparamlib; methods for flux calculations

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/flux.c,v $
		$Author: niklas $ $Date: 2006/01/20 21:48:08 $ $Revision: 1.4 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/* Table 2 of Kamae et al. */
double L_MAX[7] = {0.96, 0.96, 0.94, 0.98, 0.98, 0.94, 0.98};
double W_NDL[7] = {15.0, 20.0, 15.0, 15.0, 15.0, 20.0, 15.0};
double W_NDH[7] = {44.0, 45.0, 47.0, 42.0, 40.0, 45.0, 40.0};

typedef void (*PARAM_FUNC)(double, double*);
static PARAM_FUNC paramfunc_table[4] = {&gamma_param_nd, &gamma_param_diff, &gamma_param_delta, &gamma_param_res};

/*
		Calculate flux from non-diff process
	*/
double flux_nd(int particle, double E, double Tp, double* a) {
				double Wl, Wh, Lmin, Lmax;
				double x, y;
				double flux;
				double factor, r_factor;

				/* init some variables, given in table 2 */
				Lmin = -2.6;
				Lmax = L_MAX[particle]*log10(Tp);
				Wl = W_NDL[particle];
				Wh = W_NDH[particle];

				/* calculate the log of E and Tp */
				x = log10(E);
				y = log10(Tp);

				/* calculate the flux due to non-diffractive process for given gamma-ray energy */
				flux = a[0]*exp(-a[1]*pow(x - a[3] + a[2]*pow(x - a[3], 2), 2)) + 
								   a[4]*exp(-a[5]*pow(x - a[8] + a[6]*pow(x - a[8], 2) + a[7]*pow(x - a[8], 3), 2));
				/* factor is the kinematic limit function as in paper */
    factor = (1.0/(1.0 + exp(Wl*(Lmin - x))))*(1.0/(1.0 + exp(Wh*(x - Lmax))));			
				flux = flux*factor;
				
				if (flux < 0.0)
								flux = 0.0;
				
				/* renormalization */
				if (Tp <= 1.95) {
								r_factor = 3.05*exp(-107.0*pow((y + 3.25)/(1.0 + 8.08*(y + 3.25)), 2));
								flux = flux*r_factor;
				}

				return flux;
}

/*
		Calculate from diff. dissoc. process
	*/
double flux_diff(int particle, double E, double Tp, double* b) {
				double Wdiff, Lmax;
				double x, y;
				double flux;
				double factor;

				/* init some variables */
				Lmax = log10(Tp);
				Wdiff = 75.0;

				/* calculate the log of E and Tp */
				x = log10(E);
				y = log10(Tp);

				/* calculate the flux due to diffractive process for given gamma-ray energy */
    flux = b[0]*exp(-b[1]*pow((x - b[2])/(1.0 + b[3]*(x - b[2])), 2)) +
								   b[4]*exp(-b[5]*pow((x - b[6])/(1.0 + b[7]*(x - b[6])), 2));
				/* factor is the kinematic limit function as in paper */
				factor = 1.0/(1.0 + exp(Wdiff*(x - Lmax)));
				flux = flux*factor;

				if (flux < 0.0)
								flux = 0.0;

				return flux;
}

/*
		Calculate gamma-ray flux from either of the two resonance processes
	*/
double flux_res(int particle, double E, double Tp, double* c) {
				double Wdiff, Lmax;
				double x, y;
				double flux;
				double factor;

				/* init some variables */
				Lmax = log10(Tp);
				Wdiff = 75.0;

				/* calculate the log of E and Tp */
				x = log10(E);
				y = log10(Tp);

				/* calculate the flux due to resonance process for given gamma-ray energy */
    flux = c[0]*exp(-c[1]*pow((x - c[2])/(1.0 + c[3]*(x - c[2]) + c[4]*pow(x - c[2], 2)), 2));
				/* factor is the kinematic limit function as in paper */
				factor = 1.0/(1.0 + exp(Wdiff*(x - Lmax)));
				flux = flux*factor;

				if (flux < 0.0)
								flux = 0.0;

				return flux;
}

/*
		Calculate flux from all processes
*/
double flux(int particle, double E, double Pp) {
    double Etot, Tp;
				double f_tot;
				double f_nd, f_diff, f_d1232, f_r1600;

				double a[9];
				double b[8];
				double c[5];
				double d[5];

				f_tot = f_nd = f_diff = f_d1232 = f_r1600 = 1.0;

				Etot	= sqrt(Pp*Pp + m_p*m_p);
    Tp = Etot - m_p;

				/* calculate parameters */
				paramfunc_table[0](Pp, a);
				paramfunc_table[1](Pp, b);
				paramfunc_table[2](Pp, c);
				paramfunc_table[3](Pp, d);

				/* calculate fluxes */
				f_nd = flux_nd(particle, E, Tp, a);
				f_diff = flux_diff(particle, E, Tp, b);
				f_d1232 = flux_res(particle, E, Tp, c);
				f_r1600 = flux_res(particle, E, Tp, d);

				f_tot = f_nd + f_diff + f_d1232 + f_r1600;

				return f_tot;
}
