
#include "../sersic.h"
#include "../component.h"
#include "../recipe.h"

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

  // misc
  strcpy(cp_t -> name, name);
}

// surface density of a 2d sersic disk
double
_cps_sersic_2d_sigma(double x, double y, double const * par)
{
  // rotate back
  double r_t = sqrt(Sq(x) + Sq(y);
  double cos_rt = x / r_t,  sin_rt = y / r_t,
         cos_pt = cos(PHI), sin_pt = sin(PHI);
  double x_t = r_t * (cos_rt * cos_pt + sin_rt * sin_pt), // cos(rt - pt)
         y_t = r_t * (sin_rt * cos_pt - cos_rt * sin_pt); // sin(rt - pt)

  // generalized radius
  double r = pow( pow(abs(x_t - XC),     C)
                + pow(abs(y_t - YC) / Q, C), 1. / C);

  // calculate dn
  double n = SN, n_sq = Sq(SN);
  double dn = 3. * n - (1. / 3.) + (8. / 1215.) / n
            + (184. / 229635.) / n_sq + (1048. / 31000725.) / (n * n_sq)
            - (17557576. / 1242974068875.) / Sq(n_sq);

  // surface density
  return IRS * exp(-dn * (pow(r, 1. / n) - 1.));
}

// recipe of a 2d sersic component
double
_cps_expdisk_2d_recipe(double x, double y,
    const double * rcp, double const * par)
{
  // rotate back
  double r_t = sqrt(Sq(x) + Sq(y);
  double cos_rt = x / r_t,  sin_rt = y / r_t,
         cos_pt = cos(PHI), sin_pt = sin(PHI);
  double x_t = r_t * (cos_rt * cos_pt + sin_rt * sin_pt), // cos(rt - pt)
         y_t = r_t * (sin_rt * cos_pt - cos_rt * sin_pt); // sin(rt - pt)

  // generalized radius
  double r = pow( pow(abs(x_t - XC),     C)
                + pow(abs(y_t - YC) / Q, C), 1. / C);

  // calculate dn
  double n = SN, n_sq = Sq(SN);
  double dn = 3. * n - (1. / 3.) + (8. / 1215.) / n
            + (184. / 229635.) / n_sq + (1048. / 31000725.) / (n * n_sq)
            - (17557576. / 1242974068875.) / Sq(n_sq);

  // surface density
  double sigma = IRS * exp(-dn * (pow(r, 1. / n) - 1.));

  // grid points for recipe
  int I_age, I_Z, N_age = rcp_t -> N_age, N_Z = rcp_t -> Z;
  double * age_ax = rcp_t -> age_ax, * Z_ax = rcp_t -> Z_ax;

  // mean and std. of age and metallicity
  double c_age = AMC + AMK * r, s_age_sq = Sq(ASC + ASK * r),
         c_Z   = ZMC + ZMK * r, s_Z_sq   = Sq(ZSC + ZSK * r);

  // make recipe!
  double * rc = rcp_t -> rcp;
  for(I_age = 0; I_age < N_age; ++ I_age)
    for(I_Z = 0; I_Z < N_Z; ++ I_Z)
      *(rc + I_age * N_Z + I_Z) =
            exp(-Sq(age_ax[I_age]) / (2. * s_age_sq))
          * exp(-Sq(Z_ax[I_Z]) / (2. * s_Z_sq));

  // normalize to sigma
  int N_size = N_age * N_Z; I_px;
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
