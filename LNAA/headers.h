
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_sf_log.h>      
#include <unistd.h>
#include <time.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>


typedef struct struct_data {
  double *y, *x,*yy;
  int L, M,N,*NoORF,*NoSUM,*NoTIME,SHIFTlmn,maxy,maxNoTIME;
  char latentfilename[50];
} struct_data;

typedef struct struct_MH {
  double hP,hnu,hK,hr,accept_P,accept_nu,accept_K,accept_r;
} struct_MH;

typedef struct struct_para {
  double
    *K_lm,            *tau_K_l,
    *r_lm,            *tau_r_l,

    *K_o_l,            sigma_K_o,
    *r_o_l,            sigma_r_o,
    *nu_l,             sigma_nu,

    K_p,
    r_p,
    nu_p,
    P,
    tau_K_p, sigma_tau_K,
    tau_r_p,sigma_tau_r,
    *zeta_lm,
    sigma_zeta,zeta_p,

    *upsilon_lm,            *tau_upsilon_l,
    *upsilon_o_l,            sigma_upsilon_o,
    upsilon_p,
    tau_upsilon_p, sigma_tau_upsilon;


} struct_para;

typedef struct struct_priors {
  double
    eta_K_o,                psi_K_o,
    eta_r_o,                psi_r_o,
    eta_nu,                 psi_nu,

    K_mu,                   eta_K_p,
    r_mu,                   eta_r_p,
    nu_mu,                  eta_nu_p,
    P_mu,                   eta_P,
    df,
    tau_K_mu,tau_r_mu,
    eta_tau_K,eta_tau_r,
    eta_tau_K_p,eta_tau_r_p,
    psi_tau_K,psi_tau_r,

    eta_zeta,                 psi_zeta,
    zeta_mu,                  eta_zeta_p,

    eta_upsilon_o, psi_upsilon_o,
    upsilon_mu, eta_upsilon_p,
    tau_upsilon_mu, eta_tau_upsilon,
    eta_tau_upsilon_p, psi_tau_upsilon;

} struct_priors;


