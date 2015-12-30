
#include "./sersic.h"
#include "../component.h"
#include "../recipe.h"
#include "../utils.h"

#include "stdlib.h"
#include "string.h"
#include "math.h"

// useful parameters
#define IRS (par[ 0]) // surface density at the effective radius
#define XC  (par[ 1]) // center: x
#define YC  (par[ 2]) // center: y
#define PHI (par[ 3]) // pitch angle
#define RS  (par[ 4]) // effective radius
#define SN  (par[ 5]) // sersic n
#define Q   (par[ 6]) // axis ratio b/a
#define C   (par[ 7]) // shape parameter C
#define AMC (par[ 8])
#define ASC (par[ 9])
#define ZMC (par[10])
#define ZSC (par[11])
#define AMK (par[12])
#define ASK (par[13])
#define ZMK (par[14])
#define ZSK (par[15])

#define Sq(X) ((X) * (X))

const char * sersic_par_names[] =
  {
    "I_e", "X_c", "Y_c", "phi", "R_s", "n", "q", "c",
    "Mu_Age_c", "Std_Age_c", "Mu_M_c", "Std_M_c",
    "Mu_Age_k", "Std_Age_k", "Mu_M_k", "Std_M_k"
  };

// make a sersic component
component *
cps_sersic_2d(double * par, double * par_lim, int * is_const, const char * name)
{
  // par: I_Rs, Xc, Yc, Phi, Rs, n, q, c, Amc, Asc, Zmc, Zsc, Amk, Ask, Zmk, Zsk
  int I_par, N_fp = 0;

  // make new component
  component * cp_t = make_base_component(NPAR_SERSIC2D);

  // parameters into component
  for(I_par = 0; I_par < NPAR_SERSIC2D; ++ I_par)
    *(cp_t -> par + I_par) = *(par + I_par),
    *(cp_t -> is_const + I_par) = *(is_const + I_par);

  // range of parameters
  for(I_par = 0; I_par < 2 * NPAR_SERSIC2D; ++ I_par)
    *(cp_t -> par_lim + I_par) = *(par_lim + I_par);

  // number of free params
  for(I_par = 0; I_par < NPAR_SERSIC2D; ++ I_par)
    if(! *(is_const + I_par)) ++ N_fp;
  cp_t -> N_fp = N_fp;

  // give names to component and parameters
  strcpy(cp_t -> name, name),
  cp_t -> par_name = (char **) sersic_par_names;

  // set function headers
  cp_t -> sigma = & _cps_sersic_2d_sigma;
  cp_t -> recipe = & _cps_sersic_2d_recipe;

  return cp_t;
}

// surface density of a 2d sersic disk
double
_cps_sersic_2d_sigma(double x, double y, double const * par)
{
  // rotate back
  double r_t = sqrt(Sq(x) + Sq(y));
  double cos_rt = x / r_t,  sin_rt = y / r_t,
         cos_pt = cos(PHI), sin_pt = sin(PHI);
  if(!isfinite(cos_rt)) cos_rt = 1., sin_rt = 0.;
  double x_t = r_t * (cos_rt * cos_pt + sin_rt * sin_pt), // cos(rt - pt)
         y_t = r_t * (sin_rt * cos_pt - cos_rt * sin_pt); // sin(rt - pt)

  // generalized radius
  double r = pow( pow(fabs(x_t - XC),     C)
                + pow(fabs(y_t - YC) / Q, C), 1. / C);

  // calculate dn
  double n = SN, n_sq = Sq(SN);
  double dn = 3. * n - (1. / 3.) + (8. / 1215.) / n
            + (184. / 229635.) / n_sq + (1048. / 31000725.) / (n * n_sq)
            - (17557576. / 1242974068875.) / Sq(n_sq);

  // surface density
  return IRS * exp(-dn * (pow(r, 1. / n) - 1.));
}

int
_cps_sersic_2d_sigma_arr(int N_pt, double * x, double * y,
    double const * par, double * sigma)
{
  int I_pt;

  // rotate back
  double * r_t = TALLOC(double, N_pt),
         * cos_rt = TALLOC(double, N_pt),
         * sin_rt = TALLOC(double, N_pt),
         * x_t = TALLOC(double, N_pt),
         * y_t = TALLOC(double, N_pt);
  double cos_pt = cos(PHI), sin_pt = sin(PHI);

  FOREACH(I_pt, N_pt)
    r_t[I_pt] = sqrt(Sq(x[I_pt]) + Sq(y[I_pt]));

  FOREACH(I_pt, N_pt)
    {
      if(r_t[I_pt] != 0.)
        cos_rt[I_pt] = x[I_pt] / r_t[I_pt],
        sin_rt[I_pt] = y[I_pt] / r_t[I_pt];
      else cos_rt[I_pt] = 1., sin_rt[I_pt] = 0.;
    }

  FOREACH(I_pt, N_pt)
    {
      x_t[I_pt] = r_t[I_pt] * (cos_rt[I_pt] * cos_pt + sin_rt[I_pt] * sin_pt);
      y_t[I_pt] = r_t[I_pt] * (sin_rt[I_pt] * cos_pt - cos_rt[I_pt] * sin_pt);
    }

  // generalized radius
  double * r = TALLOC(double, N_pt);

  FOREACH(I_pt, N_pt)
    r[I_pt] = pow( pow(fabs(x_t[I_pt] - XC),     C)
                 + pow(fabs(y_t[I_pt] - YC) / Q, C), 1. / C);

  // surface density
  double n = SN, n_sq = Sq(SN);
  double dn = 3. * n - (1. / 3.) + (8. / 1215.) / n
            + (184. / 229635.) / n_sq + (1048. / 31000725.) / (n * n_sq)
            - (17557576. / 1242974068875.) / Sq(n_sq);

  FOREACH(I_pt, N_pt)
    sigma[I_pt] = IRS * exp(-dn * (pow(r[I_pt], 1. / n) - 1.));

  // free objects
  free(r_t), free(cos_rt), free(sin_rt), free(x_t), free(y_t), free(r);
}

// recipe of a 2d sersic component
int
_cps_sersic_2d_recipe(double x, double y,
    recipe * rcp_t, double const * par)
{
  // rotate back
  double r_t = sqrt(Sq(x) + Sq(y));
  double cos_rt = x / r_t,  sin_rt = y / r_t,
         cos_pt = cos(PHI), sin_pt = sin(PHI);
  double x_t = r_t * (cos_rt * cos_pt + sin_rt * sin_pt), // cos(rt - pt)
         y_t = r_t * (sin_rt * cos_pt - cos_rt * sin_pt); // sin(rt - pt)

  // generalized radius
  double r = pow( pow(fabs(x_t - XC),     C)
                + pow(fabs(y_t - YC) / Q, C), 1. / C);

  // calculate dn
  double n = SN, n_sq = Sq(SN);
  double dn = 3. * n - (1. / 3.) + (8. / 1215.) / n
            + (184. / 229635.) / n_sq + (1048. / 31000725.) / (n * n_sq)
            - (17557576. / 1242974068875.) / Sq(n_sq);

  // surface density
  double sigma = IRS * exp(-dn * (pow(r, 1. / n) - 1.));

  // grid points for recipe
  int I_age, I_Z, N_age = rcp_t -> N_age, N_Z = rcp_t -> N_Z;
  double * age_ax = rcp_t -> age_ax, * Z_ax = rcp_t -> Z_ax;

  // mean and std. of age and metallicity
  double c_age = AMC + AMK * r, s_age_sq = Sq(ASC + ASK * r),
         c_Z   = ZMC + ZMK * r, s_Z_sq   = Sq(ZSC + ZSK * r);

  // make recipe!
  double * rc = rcp_t -> rcp;
  for(I_age = 0; I_age < N_age; ++ I_age)
    for(I_Z = 0; I_Z < N_Z; ++ I_Z)
      *(rc + I_age * N_Z + I_Z) =
            exp(-Sq(age_ax[I_age] - c_age) / (2. * s_age_sq))
          * exp(-Sq(Z_ax[I_Z] - c_Z) / (2. * s_Z_sq));

  // normalize to sigma
  int N_size = N_age * N_Z, I_px;
  double sum_t = 0;
  for(I_px = 0; I_px < N_size; ++ I_px) sum_t += rc[I_px];

  double K_ft = sigma / sum_t;
  for(I_px = 0; I_px < N_size; ++ I_px) rc[I_px] *= K_ft;

  return 0;
}

// useful parameters
#undef IRS
#undef XC
#undef YC
#undef PHI
#undef RS
#undef SN
#undef Q
#undef C
#undef AMC
#undef ASC
#undef ZMC
#undef ZSC
#undef AMK
#undef ASK
#undef ZMK
#undef ZSK

#undef Sq
