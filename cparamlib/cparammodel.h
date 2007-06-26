/*

 * cparammodel.h
 *
 * Header file for cparamlib
 *
 * $Source: /home/nkarlsson/usr/cvsroot/cparamlib/cparamlib/Attic/cparammodel.h,v $
 * $Author: niklas $ $Date: 2007/06/26 17:17:53 $ $Revision: 1.6 $
 *
 */

#ifndef CPARAMLIB_H
#define CPARAMLIB_H

/* 
 * Definition of particle ids
 */
#define ID_GAMMA 0
#define ID_ELECTRON 1
#define ID_POSITRON 2
#define ID_NUE 3
#define ID_NUMU 4
#define ID_ANTINUE 5
#define ID_ANTINUMU 6

/*
 * struct Definitions
 */

/* 
 * PARAMSET
 *
 * Contains four double arrays used for calculating inclusive cross sections 
 * (dsigma/dlogE) according to Kamae et al. (2006). Each double array
 * corresponds to a parameter set in the paper.
 */
typedef struct {
    double a[9];
    double b[8];
    double c[5];
    double d[5];
} PARAMSET;

/* 
 * PARAMSET_PT
 *
 * Contains doubles and double arrays used for calculating pT distributions
 * (differential cross section d^2sigma/dlogEdpT) according to Karlsson & 
 * Kamae (2007).
 */
typedef struct {
    double a0, a1;
    double a1i[4];
    double b0, b1, b2;
    double b1i[4], b2i[4];
    double c0, c1, c2;
    double c1i[4], c2i[4];

    PARAMSET params;
} PARAMSET_PT;

/* 
 * API Definitions
 *  
 * Listed below are all the functions implemented in this library
 */

/*
 * Gamma-rays
 */
void gamma_param_nd(double Tp, PARAMSET* params);
void gamma_param_diff(double Tp, PARAMSET* params);
void gamma_param_delta(double Tp, PARAMSET* params);
void gamma_param_res(double Tp, PARAMSET* params);
void gamma_param(double Tp, PARAMSET* params);

void gamma_pt_param_nr(double E, double Tp, PARAMSET_PT* pt_params, int flag);
void gamma_pt_param_delta(double E, double Tp, PARAMSET_PT* pt_params, int flag);
void gamma_pt_param_res(double E, double Tp, PARAMSET_PT* pt_params, int flag);
void gamma_pt_param(double E, double Tp, PARAMSET_PT* pt_params, int flag);

/*
 * Electrons
 */
void elec_param_nd(double Tp, PARAMSET* params);
void elec_param_diff(double Tp, PARAMSET* params);
void elec_param_delta(double Tp, PARAMSET* params);
void elec_param_res(double Tp, PARAMSET* params);
void elec_param(double Tp, PARAMSET* params);

/*
 * Positrons
 */
void posi_param_nd(double Tp, PARAMSET* params);
void posi_param_diff(double Tp, PARAMSET* params);
void posi_param_delta(double Tp, PARAMSET* params);
void posi_param_res(double Tp, PARAMSET* params);
void posi_param(double Tp, PARAMSET* params);

/*
 * Electron neutrinos
 */
void nue_param_nd(double Tp, PARAMSET* params);
void nue_param_diff(double Tp, PARAMSET* params);
void nue_param_delta(double Tp, PARAMSET* params);
void nue_param_res(double Tp, PARAMSET* params);
void nue_param(double Tp, PARAMSET* params);

/*
 * Electron anti-neutrinos
 */
void antinue_param_nd(double Tp, PARAMSET* params);
void antinue_param_diff(double Tp, PARAMSET* params);
void antinue_param_delta(double Tp, PARAMSET* params);
void antinue_param_res(double Tp, PARAMSET* params);
void antinue_param(double Tp, PARAMSET* params);

/*
 * Muon neutrinos
 */
void numu_param_nd(double Tp, PARAMSET* params);
void numu_param_diff(double Tp, PARAMSET* params);
void numu_param_delta(double Tp, PARAMSET* params);
void numu_param_res(double Tp, PARAMSET* params);
void numu_param(double Tp, PARAMSET* params);

/*
 * Muon anti-neutrinos
 */
void antinumu_param_nd(double Tp, PARAMSET* params);
void antinumu_param_diff(double Tp, PARAMSET* params);
void antinumu_param_delta(double Tp, PARAMSET* params);
void antinumu_param_res(double Tp, PARAMSET* params);
void antinumu_param(double Tp, PARAMSET* params);

/*
 * Calculation of the inclusive cross section, dsigma/dlogE, for any given secondary
 * particle and any given energy of that secondary particle
 *
 * Proton kinetic energy Tp is passed as well as a pointer to a struct holding the
 * parameters (at that Tp)
 */
double sigma_incl_nd(int particle, double E, double Tp, PARAMSET* params);
double sigma_incl_diff(int particle, double E, double Tp, PARAMSET* params);
double sigma_incl_delta(int particle, double E, double Tp, PARAMSET* params);
double sigma_incl_res(int particle, double E, double Tp, PARAMSET* params);
double sigma_incl_tot(int particle, double E, double Tp, PARAMSET* params);

/*
 * Calculation of the differential cross section dsigma/dlogEdpT, for any gamma-ray energy
 *
 * Proton kinetic energy Tp is passed as well as a pointer to a struct keeping all the 
 * parameters (at given Tp and log(E))
 */
double sigma_pt_nr(int particle, double pT, double E, double Tp, PARAMSET_PT* pt_params);
double sigma_pt_delta(int particle, double pT, double E, double Tp, PARAMSET_PT* pt_params);
double sigma_pt_res(int particle, double pT, double E, double Tp, PARAMSET_PT* pt_params);
double sigma_pt_tot(int particle, double pT, double E, double Tp, PARAMSET_PT* pt_params);

/*
 * Calculation of the p-p cross section for a given proton momentum
 */
double sigma_pp_nd(double Pp);
double sigma_pp_diff(double Pp);
double sigma_pp_delta(double Pp);
double sigma_pp_res(double Pp);

#endif
