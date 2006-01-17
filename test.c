/*
		test.c

		Test app for the parameterization model lib

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/test.c,v $
		$Author: niklas $ $Date: 2006/01/17 19:57:51 $ $Revision: 1.4 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

int main(void) {
				double Tp, Pp, E;
				double f;
				double a[9];
				int i;

				Tp = 512000.0;
				Pp = 512000.938;
				E = 1.0e2;

				/* Example 1: calculate non-diff params and flux */
				printf("Example 1: non-diff");
				gamma_param_nd(Pp, a);
				for (i = 0; i < 8; i++)
								printf("a[%d] = %10e\n", i, a[i]);
				f = flux_nd(ID_GAMMA, E, Tp, a);
				printf("flux_nd(Tp=512TeV) = %10e\n", f);

				/* Example 2: total flux */
				printf("Example 2: total flux");
				f = flux(ID_GAMMA, E, Pp);
				printf("total flux(Tp=512TeV) = %10e\n", f);

				return 0;
}

