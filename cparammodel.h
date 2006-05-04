/*
		cparammodel.h

		Header file for cparamlib

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/cparammodel.h,v $
		$Author: niklas $ $Date: 2006/05/04 20:45:09 $ $Revision: 1.7 $
*/

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

void gamma_param_nd(double Tp, double* a);
void gamma_param_diff(double Tp, double* b);
void gamma_param_delta(double Tp, double* c);
void gamma_param_res(double Tp, double* d);

void elec_param_nd(double Tp, double* a);
void elec_param_diff(double Tp, double* b);
void elec_param_delta(double Tp, double* c);
void elec_param_res(double Tp, double* d);

void posi_param_nd(double Tp, double* a);
void posi_param_diff(double Tp, double* b);
void posi_param_delta(double Tp, double* c);
void posi_param_res(double Tp, double* d);

void nue_param_nd(double Tp, double* a);
void nue_param_diff(double Tp, double* b);
void nue_param_delta(double Tp, double* c);
void nue_param_res(double Tp, double* d);

void antinue_param_nd(double Tp, double* a);
void antinue_param_diff(double Tp, double* b);
void antinue_param_delta(double Tp, double* c);
void antinue_param_res(double Tp, double* d);

void numu_param_nd(double Tp, double* a);
void numu_param_diff(double Tp, double* b);
void numu_param_delta(double Tp, double* c);
void numu_param_res(double Tp, double* d);

void antinumu_param_nd(double Tp, double* a);
void antinumu_param_diff(double Tp, double* b);
void antinumu_param_delta(double Tp, double* c);
void antinumu_param_res(double Tp, double* d);

double sigma_nd(int particle, double E, double Tp, double* a);
double sigma_diff(int particle, double E, double Tp, double* b);
double sigma_delta(int particle, double E, double Tp, double* c);
double sigma_res(int particle, double E, double Tp, double* d);

double sigma(int particle, double E, double Tp);

#endif
