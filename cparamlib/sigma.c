/*
 * sigma.c
 *
 * Main part of cparamlib; methods for calculations of differential
 * cross sections. Functions for inclusive cross sections as well 
 * as kinematic cutoff functions are given in Kamae et al. (2006)
 * and functions for pT distributions are given in Karlsson and
 * Kamae (2007).
 *
 * 10/11/2006: On request we have added functions for calculating
 * sigma_pp for the for components as given in the paper. We also
 * changed sigma to sigma_incl.
 *
 * 05/01/2007: Functionality to calculate pT distributions added
 *
 * 05/22/2007: Changed all functions to take a pointer to a struct
 * containing parameters rather than taking pointers to double
 * arrays.
 *
 * $Source: /home/nkarlsson/usr/cvsroot/cparamlib/cparamlib/sigma.c,v $
 * $Author: niklas $ $Date: 2007/06/26 17:17:53 $ $Revision: 1.7 $
 *
 */

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/* 
 * Table 2 of Kamae et al. 
 */
double L_MAX[7] = {0.96, 0.96, 0.94, 0.98, 0.98, 0.94, 0.98};
double W_NDL[7] = {15.0, 20.0, 15.0, 15.0, 15.0, 20.0, 15.0};
double W_NDH[7] = {44.0, 45.0, 47.0, 42.0, 40.0, 45.0, 40.0};

/*
 * Calculate inclusive cross section from the non-diff process
 */
double sigma_incl_nd(int particle, double E, double Tp, PARAMSET* params)
{
    double Wl, Wh, Lmin, Lmax;
    double x, y;
    double xa3, xa8;
    double pow1, pow2;
    double sigma;
    double cutoff, r_factor;

    /* check whether params is a null pointer or not */
    if (params == NULL)
        return 0;

    /* calculate log(E) and log(Tp) */
    x = log10(E);
    y = log10(Tp*0.001);

    /* init some variables, given in table 2 */
    Lmin = -2.6;
    Lmax = L_MAX[particle]*(y + 3.0);
    Wl = W_NDL[particle];
    Wh = W_NDH[particle];

    /* calculate the flux due to non-diffractive process for given gamma-ray energy */
    /* 06/12/06: replaced code involving pow() */
    xa3 = x - params->a[3];
    pow1 = xa3*(1 + params->a[2]*xa3);
    xa8 = x - params->a[8];
    pow2 = xa8*(1 + xa8*(params->a[6] + params->a[7]*xa8));
    sigma = params->a[0]*exp(-params->a[1]*pow1*pow1) + params->a[4]*exp(-params->a[5]*pow2*pow2);

    /* cutoff is the kinematic limit function as in the paper */
    cutoff = (1/(1 + exp(Wl*(Lmin - x))))*(1/(1 + exp(Wh*(x - Lmax))));   
    sigma = sigma*cutoff;
    
    if (sigma < 0)
        sigma = 0;
    
    /* renormalization
       this is different for each particle, thus we must use if statements
     */
    r_factor = 1;
    switch (particle) {
        /* gamma */
        case ID_GAMMA:
            if (Tp <= 1.95) {
                pow1 = (y + 3.25)/(1 + 8.08*(y + 3.25));
                r_factor = 3.05*exp(-107.0*pow1*pow1);
            } else
                r_factor = 1.01;
            break;
        /* electron */
        case ID_ELECTRON:
            if (Tp <= 15.6) {
                pow1 = (y + 3.26)/(1 + 9.21*(y + 3.26));
                r_factor = 3.63*exp(-106*pow1*pow1) + y*(-0.182 - 0.175*y);
            } else
                r_factor = 1.01;
            break;
        /* positron */
        case ID_POSITRON:
            if (Tp <= 5.52) {
                pow1 = (y + 3.25)/(1 + 10.4*(y + 3.25));
                r_factor = 2.22*exp(-98.9*pow1*pow1);
            }
            break;
        /* electron neutrino */
        case ID_NUE:
            if (Tp <= 7.81) {
                pow1 = (y + 3.26)/(1 + 6.56*(y + 3.26));
                r_factor = 0.329*exp(-247*pow1*pow1) + y*(-0.959 - 0.229*y);
            }
            break;
        /* muon neutrino */
        case ID_NUMU:
            if (Tp <= 15.6) {
                pow1 = (y + 3.25)/(1 + 8.38*(y + 3.25));
                r_factor = 2.23*exp(-93.4*pow1*pow1) + y*(-0.376 - 0.121*y);
            }
            break;
        /* electron anti-neutrino */
        case ID_ANTINUE:
            if (Tp <= 15.6) {
                pow1 = (y + 3.27)/(1 + 6.59*(y + 3.27));
                r_factor = 2.67*exp(-45.7*pow1*pow1) + y*(-0.301 - 0.208*y);
            }
            break;
        /* muon anti-neutrino */
        case ID_ANTINUMU:
            if (Tp <= 15.6) {
                pow1 = (y + 3.25)/(1 + 8.34*(y + 3.25));
                r_factor = 2.56*exp(-107*pow1*pow1) + y*(-0.385 - 0.125*y);
            }
            break;
    }
    sigma = sigma*r_factor;

    return sigma;
}

/*
 * Calculate inclusive cross section from the diffraction dissociation process
 */
double sigma_incl_diff(int particle, double E, double Tp, PARAMSET* params)
{
    double Wdiff, Lmax;
    double x, y;
    double pow1, pow2;
    double sigma;
    double cutoff;

    /* check whether params is a null pointer or not */
    if (params == NULL)
        return 0;

    /* calculate log(E) and log(Tp) */
    x = log10(E);
    y = log10(Tp*0.001);

    /* init some variables */
    Lmax = y + 3.0;
    Wdiff = 75.0;

    /* calculate the sigma due to diffractive process for given gamma-ray energy */
    /* 06/12/06: replaced code involving pow() */
    pow1 = (x - params->b[2])/(1 + params->b[3]*(x - params->b[2]));
    pow2 = (x - params->b[6])/(1 + params->b[7]*(x - params->b[6]));
    sigma = params->b[0]*exp(-params->b[1]*pow1*pow1) + params->b[4]*exp(-params->b[5]*pow2*pow2);

    /* cutoff is the kinematic limit function as in the paper */
    cutoff = 1/(1 + exp(Wdiff*(x - Lmax)));
    sigma = sigma*cutoff;

    if (sigma < 0)
        sigma = 0;

    return sigma;
}

/*
 * Calculate inclusive cross section from the Delta(1232) resonance
 */
double sigma_incl_delta(int particle, double E, double Tp, PARAMSET* params)
{
    double Wdiff, Lmax;
    double x, y;
    double xc2;
    double pow;
    double sigma;
    double cutoff;

    /* calculate log(E) and log(Tp) */
    x = log10(E);
    y = log10(Tp*0.001);

    /* init some variables */
    Lmax = y + 3.0;
    Wdiff = 75.0;

    /* calculate the sigma due to resonance process for given gamma-ray energy */
    /* 06/12/06: replaced code involving pow() */
    xc2 = x - params->c[2];
    pow = xc2/(1 + xc2*(params->c[3] + params->c[4]*xc2));
    sigma = params->c[0]*exp(-params->c[1]*pow*pow);

    /* cutoff is the kinematic limit function as in the paper */
    cutoff = 1/(1 + exp(Wdiff*(x - Lmax)));
    sigma = sigma*cutoff;

    if (sigma < 0)
        sigma = 0;

    return sigma;
}

/*
 * Calculate inclusive cross section from the res(1600) resonance
 */
double sigma_incl_res(int particle, double E, double Tp, PARAMSET* params)
{
    double Wdiff, Lmax;
    double x, y;
    double xd2;
    double pow;
    double sigma;
    double cutoff;

    /* check whether params is a null pointer or not */
    if (params == NULL)
        return 0;

    /* calculate log(E) and log(Tp) */
    x = log10(E);
    y = log10(Tp*0.001);

    /* init some variables */
    Lmax = y + 3.0;
    Wdiff = 75.0;

    /* calculate the sigma due to resonance process for given gamma-ray energy */
    /* 06/12/06: replaced code involving pow() */
    xd2 = x - params->d[2];
    pow = xd2/(1 + xd2*(params->d[3] + params->d[4]*xd2));
    sigma = params->d[0]*exp(-params->d[1]*pow*pow);

    /* cutoff is the kinematic limit function as in the paper */
    cutoff = 1/(1 + exp(Wdiff*(x - Lmax)));
    sigma = sigma*cutoff;

    if (sigma < 0)
        sigma = 0;

    return sigma;
}

/*
 * Calculate total inclusive cross section (sum of all processes)
 */
double sigma_incl_tot(int particle, double E, double Tp, PARAMSET* params)
{
    double f_tot;
    double f_nd, f_diff, f_delta, f_res;

    /* check whether params is a null pointer or not */
    if (params == NULL)
        return 0;

    f_tot = f_nd = f_diff = f_delta = f_res = 0;

    /* calculate sigmas */
    f_nd = sigma_incl_nd(particle, E, Tp, params);
    f_diff = sigma_incl_diff(particle, E, Tp, params);
    f_delta = sigma_incl_delta(particle, E, Tp, params);
    f_res = sigma_incl_res(particle, E, Tp, params);

    f_tot = f_nd + f_diff + f_delta + f_res;

    return f_tot;
}

/*
 * Calculate the differential cross section dsigma/dlogEdpT from non-res process
 * (i.e. non-diffractive + diffractive)
 */
double sigma_pt_nr(int particle, double pT, double E, double Tp, PARAMSET_PT* pt_params)
{
    double x;
    double x1, Lp, W;
    double sigma;
    double cutoff;

    /* check whether params is a null pointer or not */
    if (pt_params == NULL)
        return 0;

    sigma = 0;
    W = 75.0;
    Lp = 0;

    x = log10(E);

    if (particle == ID_GAMMA) {
        sigma = pt_params->a0*pT*exp(-pT/pt_params->a1);
        
        /* cutoff is the kinematic limit function as in the paper */
        if (x > 0.5) {
            Lp = 2.5;
        } else if (x > -1.0) {
            x1 = x + 1.0;
            Lp = -0.793 + exp(0.271*x1 + 0.363*x1*x1);
        } else {
            Lp = 0.09755 + 0.670*exp(1.81*x);
        }
        cutoff = 1/(1 + exp(W*(pT - Lp)));
        sigma = sigma*cutoff;

        if (sigma < 0)
            sigma = 0;
    } 
    
    return sigma;
}

/*
 * Calculate the differential cross section dsigma/dlogEdpT from the Delta(1232) resonance
 */
double sigma_pt_delta(int particle, double pT, double E, double Tp, PARAMSET_PT* pt_params)
{
    double x;
    double pow;
    double x1, Lp, W;
    double sigma;
    double cutoff;

    /* check whether params is a null pointer or not */
    if (pt_params == NULL)
        return 0;

    sigma = 0;
    W = 75.0;
    Lp = 0;

    x = log10(E);

    if (particle == ID_GAMMA) {
        pow = pT - pt_params->b1;
        sigma = pt_params->b0*pT*exp(-pow*pow/pt_params->b2);

        /* cutoff is the kinematic limit function as in the paper */
        if (x > 0.5) {
            Lp = 2.5;
        } else if (x > -1.0) {
            x1 = x + 1.0;
            Lp = -0.793 + exp(0.271*x1 + 0.363*x1*x1);
        } else {
            Lp = 0.09755 + 0.670*exp(1.81*x);
        }
        cutoff = 1/(1 + exp(W*(pT - Lp)));

        sigma = sigma*cutoff;

        if (sigma < 0)
            sigma = 0;
    }

    return sigma;
}

/*
 * Calculate the differential cross section dsigma/dlogEdpT from the res(1600) resonance
 */
double sigma_pt_res(int particle, double pT, double E, double Tp, PARAMSET_PT* pt_params)
{
    double x;
    double pow;
    double x1, Lp, W;
    double sigma;
    double cutoff;

    /* check whether params is a null pointer or not */
    if (pt_params == NULL)
        return 0;

    sigma = 0;
    W = 75.0;
    Lp = 0;

    x = log10(E);

    if (particle == ID_GAMMA) {
        pow = pT - pt_params->c1;
        sigma = pt_params->c0*pT*exp(-pow*pow/pt_params->c2);

        /* cutoff is the kinematic limit function as in the paper */
        if (x > 0.5) {
            Lp = 2.5;
        } else if (x > -1.0) {
            x1 = x + 1.0;
            Lp = -0.793 + exp(0.271*x1 + 0.363*x1*x1);
        } else {
            Lp = 0.09755 + 0.670*exp(1.81*x);
        }
        cutoff = 1/(1 + exp(W*(pT - Lp)));
        sigma = sigma*cutoff;

        if (sigma < 0)
            sigma = 0;
    }

    return sigma;
}

/*
 * Calculate the differential cross section dsigma/dlogEdp from all processes
 */
double sigma_pt_tot(int particle, double pT, double E, double Tp, PARAMSET_PT* pt_params)
{
    double f_tot;
    double f_nr, f_delta, f_res;

    /* check whether params is a null pointer or not */
    if (pt_params == NULL)
        return 0;

    f_tot = f_nr = f_delta = f_res = 0;

    /* calculate sigmas */
    f_nr = sigma_pt_nr(particle, pT, E, Tp, pt_params);
    f_delta = sigma_pt_delta(particle, pT, E, Tp, pt_params);
    f_res = sigma_pt_res(particle, pT, E, Tp, pt_params);
    
    f_tot = f_nr + f_delta + f_res;

    return f_tot;
}

/*
 * Calculate non-diffractive inelastic p-p cross section
 * as given by equation 1 in Kamae et al. (2006)
 */
double sigma_pp_nd(double Pp)
{
    double a[8] = {0.1176, 0.3829, 23.10, 6.454, -5.764, -23.63, 94.75, 0.02667};
    double b[2] = {11.34, 23.72};
    double c[3] = {28.5, -6.133, 1.464};
    double x, sigma;

    x = log10(Pp);
    sigma = 0;

    if ((Pp >= 1) && (Pp < 1.3)) {
        sigma = 0.57*pow(x/a[0], 1.2)*(a[2] + x*x*(a[3] + a[4]*x) + a[5]*exp(-a[6]*(x + a[7])*(x + a[7])));
    } else if ((Pp >= 1.3) && (Pp < 2.4)) {
        sigma = (b[0]*abs(a[1] - x) + b[1]*abs(a[0] - x))/(a[1] - a[0]);
    } else if ((Pp >= 2.4) && (Pp <= 10.0)) {
        sigma = a[2] + x*x*(a[3] + a[4]*x) + a[5]*exp(-a[6]*(x + a[7])*(x + a[7]));
    } else if (Pp > 10.0) {
        sigma = c[0] + c[1]*x + c[2]*x*x;
    }

    return sigma;
}

/*
 * Calculate diffractive inelastic p-p cross section
 * as given by equation 2 in Kamae et al. (2006)
 */
double sigma_pp_diff(double Pp)
{
    double d[7] = {0.3522, 0.1530, 1.498, 2.0, 30.0, 3.155, 1.042};
    double e[2] = {5.922, 1.632};

    double x, sigma;

    x = log10(Pp);
    sigma = 0;

    if ((Pp >= 2.25) && (Pp < 3.2)) {
        sigma = sqrt((x-d[0])/d[1])*(d[2] + d[3]*log10(d[4]*(x - 0.25)) + x*x*(d[5] - d[6]*x));
    } else if ((Pp >= 3.2) && (Pp <= 100.0)) {
        sigma = d[2] + d[3]*log10(d[4]*(x - 0.25)) + x*x*(d[5] - d[6]*x);
    } else if (Pp > 100.0) {
        sigma = e[0] + e[1]*x;
    }

    return sigma;
}

/*
 * Calculate delta(1232) inelastic p-p cross section
 * as given by equation 3 in Kamae et al. (2006)
 */
double sigma_pp_delta(double Pp)
{
    double f[5] = {0.0834, 9.5, -5.5, 1.68, 3134.0};

    double Ep, sigma;

    /* proton mass 0.938 GeV/c^2 gives mass squared  0.879844 (GeV/c^2)^2 */
    Ep = sqrt(Pp*Pp + 0.879844);
    sigma = 0;

    if ((Ep >= 1.4) && (Ep < 1.6)) {
        sigma = f[0]*pow(Ep, 10);
    } else if ((Ep >= 1.6) && (Ep < 1.8)) {
        sigma = f[1]*exp(-f[2]*(Ep - f[3])*(Ep - f[3]));
    } else if ((Ep >= 1.8) && (Ep <= 10.0)) {
        sigma = f[4]*pow(Ep, -10);
    }

    return sigma;
}

/*
 * Calculate res(1600) inelastic p-p cross section
 * as given by equation 4 in Kamae et al. (2006)
 */
double sigma_pp_res(double Pp)
{
    double g[5] = {0.0004257, 4.5, -7.0, 2.1, 503.5};

    double Ep, sigma;

    /* proton mass 0.938 GeV/c^2 gives mass squared  0.879844 (GeV/c^2)^2 */
    Ep = sqrt(Pp*Pp + 0.879844);
    sigma = 0;

    if ((Ep >= 1.6) && (Ep < 1.9)) {
        sigma = g[0]*pow(Ep, 14);
    } else if ((Ep >= 1.9) && (Ep < 2.3)) {
        sigma = g[1]*exp(-g[2]*(Ep - g[3])*(Ep - g[3]));
    } else if ((Ep >= 2.3) && (Ep <= 20.0)) {
        sigma = g[4]*pow(Ep, -6);
    }

    return sigma;
}
