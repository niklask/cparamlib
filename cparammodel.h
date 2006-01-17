#ifndef _GAMMA_H_
#define _GAMMA_H_


/* particle ids */
#define ID_GAMMA 0
#define ID_ELECTRON 1
#define ID_POSITRON 2
#define ID_NUE 3
#define ID_NUMU 4
#define ID_ANTINUE 5
#define ID_ANTINUMU 6

/* particle masses */
#define m_p 0.938

void gamma_param_nd(double Pp, double* a);
void gamma_param_diff(double Pp, double* b);
void gamma_param_delta(double Pp, double* c);
void gamma_param_res(double Pp, double* d);

double flux_nd(int particle, double E, double Tp, double* a);
double flux_diff(int particle, double E, double Tp, double* b);
double flux_res(int particle, double E, double Tp, double* c);

double flux(int particle, double E, double Pp);

#endif
