/*
		test.c

		Test app for the parameterization model lib

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/test.c,v $
		$Author: niklas $ $Date: 2006/03/16 22:57:58 $ $Revision: 1.9 $
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

				/* Example 1: calculate parameter a0,...,a8 used in the non-diff inclusive 
							gamma-ray cross-section formula */
				printf("Example 1: non-diff\n");
				gamma_param_nd(Pp, a);
				for (i = 0; i < 8; i++)
								printf("a[%d] = %10e\n", i, a[i]);
				f = sigma_nd(ID_GAMMA, E, Tp, a);
				printf("sigma_nd(Tp=512TeV) = %10e\n", f);

				/* Example 2: inclusive gamma-ray cross-section */
				printf("\n");
				printf("Example 2: total inclusive gamma-ray cross-section\n");
				f = sigma(ID_GAMMA, E, Pp);
				printf("total inclusive gamma-ray cross-section(Tp=512TeV) = %10e\n", f);

				return 0;
}

