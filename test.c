#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

int main(void) {
				double Tp, Pp, E;
				double f;
				double a[8];
				int i;

				Tp = 512000.0;
				Pp = 512000.938;
				E = 1.0e2;

				gamma_param_diff(Pp, a);

				for (i = 0; i < 8; i++)
								printf("a[%d] = %10e\n", i, a[i]);

				f = flux_diff(ID_GAMMA, E, Tp, a);

				printf("flux_diff(Tp=512TeV) = %10e\n", f);

				return 0;
}

