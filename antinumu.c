/*
	 antinumu.c

		Parameter calculation for muon anti neutrinos

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/antinumu.c,v $
		$Author: niklas $ $Date: 2006/02/02 22:27:37 $ $Revision: 1.2 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/*
		Calculate parameters for muon anti neutrino from non-diff
*/
void antinumu_param_nd(double Pp, double* a) {
    double Etot, Tp, y;
				int i;

				Etot	= sqrt(Pp*Pp + m_p*m_p);
    Tp = Etot - m_p;
    y = log10(Tp*0.001);

				/* These are the original parameter values
							note that a0 and a4 have been multiplied with 1/(10^0.05 - 1) for unit conversion
							-0.185991585255 1.23324215412 -0.526211738586 0.0977133139968 -0.00594520336017
							-2.62967587332e-05 9.38576922636e-05 0.323838979006 7.78205503593e-06 4.05602025986 7.61485125622e-06 -8.40908967348e-07
							-243.62197876 59.3739776611 -5.73555803299 1.98147809505 -1.04779481888
							0.508068501949 0.602210879326 0.00341200479306 -0.0111390100792
							0.323135674 0.544017136097 -0.155495405197 0.0142262140289 0.000371871079551
							9.11014865324e-07 1.38803295613e-06 0.0169981326908 0.000125834150822 4.77067899704
							-272.112670898 53.4774284363 0.355306923389 56.0406951904 -6.08761930466
							6431.83984375 893.920410156 5.70126656996e-09 2103.59521484 5.6739859581 -6.1124920845
							-0.111196398735 0.381441265345 -0.0401276387274 0.00474844034761 0.00547071173787
				*/
				if (Tp > 0.4) {
								a[0] = -1.5243*(y + 3.3) + 10.107*pow(y + 3.3, 2) - 4.3126*pow(y + 3.3, 3) + 0.80081*pow(y + 3.3, 4) - 0.048724*pow(y + 3.3, 5);
								a[1] = -2.6297e-5 + 9.3858e-5*exp(-0.32384*(y + 3.4)) + 7.7821e-6/(y + 4.0560) + 7.6149e-6*y - 8.4091e-7*y*y;
								a[2] = -243.62 + 59.374*y - 5.7356*y*y + 1.9815*y*y*y - 1.0478*y*y*y*y;
								a[3] = 0.50807 + 0.60221*y + 0.0034120*y*y - 0.011139*y*y*y;
								a[4] = 2.6483*(y + 3.32) + 4.4585*pow(y + 3.32, 2) - 1.2744*pow(y + 3.32, 3) + 0.11659*pow(y + 3.32, 4) + 0.0030477*pow(y + 3.32, 5);
								a[5] = 9.1101e-7 + 1.3880e-6*log10(0.016998*(y + 4.3)) + 0.00012583/pow(y + 4.7707, 2);
								a[6] = -272.11 + 53.477*log10(0.35531*(y + 3.4)) + 56.041*y - 6.0876*y*y;
								a[7] = 6431.8 + 893.92*log10(5.7013e-9*(y + 3.9)) + 2103.6/(y + 5.6740) - 6.1125*y*y;
								a[8] = -0.11120 + 0.38144*y - 0.040128*y*y + 0.0047484*y*y*y + 0.0054707*y*y*y*y;
				} else {
								for (i = 0; i < 9; i++)
												a[i] = 0.0;
				}
}

/*
		Calculate parameters for muon anti neutrino from diff. dissoc.
*/
void antinumu_param_diff(double Pp, double* b) {
    double Etot, Tp, y;
				int i;

				Etot	= sqrt(Pp*Pp + m_p*m_p);
    Tp = Etot - m_p;
    y = log10(Tp*0.001);

				/* These are the original parameter values
							note that b0 and b4 have been multiplied with 1/(10^0.05 - 1) for unit conversion
							8.59378910065 -0.35816475749 -0.815034985542 0.522732317448 0.000801161513664 9.52664852142
							-8.11448192596 7685.99902344 44046.3203125 -2.21900486946 81.1047286987
							-1.3094547987 0.071270160377 0.00754632754251 0.0677590221167 5.34327173233 -0.00442051235586 -1.86827158928
							0.0821493342519 2190.05322266 533.749938965 -2.83626127243 7.09758710861
							0.33603516221 0.0413145460188 -0.000418203737354 0.142508044839 10.4077215195 -0.189221933484 
							2.18170928955 -0.595835030079 2.7053797245 -0.010909171775 -14.8804807663
							1.45912253857 1.02747344971 -0.0749489739537 -0.0060396399349 0.00975683052093
							3.76090764999 4.28426742554 -0.371478229761 -0.164793148637 -2.76527428627
				*/

				if (Tp > 1.38) {
								if (Tp > 11.0) {
												b[0] = 70.430*tanh(-0.35816*(y + 2.2)) - 6.6796*pow(y + 0.52273, 2) + 0.0065659*pow(y + 9.5266, 4);
												b[1] = -8.1145 + 7686.0*exp(-44046*pow((y + 2.2190)/(1.0 + 81.105*(y + 2.2190)), 2));
												b[2] = -1.3095 + 0.071270*tanh(0.0075463*(y + 1.9)) + 0.067759*(y + 5.3433) - 0.0044205*pow(y - 1.8683, 2);
												b[3] = 0.08215 - 2190.1*exp(-533.75*pow((y + 2.8363)/(1.0 + 7.0976*(y + 2.8363)), 2));
								} else {
												b[0] = 0.0;
												b[1] = 0.0;
												b[2] = 0.0;
												b[3] = 0.0;
								}
								b[4] = 2.7540 + 0.33859*pow(y + 2.95, 2) - 0.0034274*pow(y + 2.95, 4) + 1.1679*exp(-10.408*pow(y + 2.28 - 0.18922*pow(y + 2.2, 2), 2));
								b[5] = 2.1817 - 0.59584*tanh(2.7054*(y + 1.9)) - 0.010909*(y - 14.880);
								b[6] = 1.4591 + 1.0275*y - 0.074949*y*y - 0.0060396*y*y*y + 0.0097568*y*y*y*y;
								b[7] = 3.7609 + 4.2843*tanh(-0.37148*(y + 2.1)) - 0.16479*pow(y - 2.7653, 2);
				} else {
								for (i = 0; i < 8; i++)
												b[i] = 0.0;
				}
}

/*
		Calculate parameters for muon anti neutrino from delta(1232)
*/
void antinumu_param_delta(double Pp, double* c) {
				c[0] = c[1] = c[2] = c[3] = c[4] = 0.0;
}

/*
		Calculate parameters for muon anti neutrino from res(1600)
*/
void antinumu_param_res(double Pp, double* d) {
    double Etot, Tp, y;
				int i;

				Etot	= sqrt(Pp*Pp + m_p*m_p);
    Tp = Etot - m_p;
    y = log10(Tp*0.001);

				/* These are the original parameter values
							note that d0 have been multiplied with 1/(10^0.05 - 1) for unit conversion
							0.0444865971804 58.2102165222 -2.95368266106 1.43199265003 0.0137672061101 0.00564260641113
							-9.50664615631 -5.46553182602 -0.317691653967
							-7.18312692642 -7.15508651733 30.3540401459 0.337569147348
							2.79380726814 1.6991943121 0.201611116529
							0.618783593178 0.623712420464 0.189134046435 0.0191184338182
				*/

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
