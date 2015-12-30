
#include "expdisk.h"
#include "../component.h"
#include "../recipe.h"
#include "../utils.h"

#include "stdlib.h"
#include "string.h"
#include "math.h"

#include "stdio.h"

// useful definitions
#define IC  (par[ 0])
#define XC  (par[ 1])
#define YC  (par[ 2])
#define PHI (par[ 3])
#define A   (par[ 4])
#define Q   (par[ 5])
#define C   (par[ 6])
#define AMC (par[ 7])
#define ASC (par[ 8])
#define ZMC (par[ 9])
#define ZSC (par[10])
#define AMK (par[11])
#define ASK (par[12])
#define ZMK (par[13])
#define ZSK (par[14])

#define Sq(X) ((X) * (X))

const char * expdisk_par_names[] =
  {
    "I_c", "X_c", "Y_c", "phi", "a", "q", "c",
    "Mu_Age_c", "Std_Age_c", "Mu_M_c", "Std_M_c",
    "Mu_Age_k", "Std_Age_k", "Mu_M_k", "Std_M_k"
  };

component *
cps_expdisk_2d(double * par, double * par_lim, int * is_const,
    const char * name)
{
  // par: S_c, x_c, y_c, phi, a, Q, c, Amc, Asc, Zmc, Zsc, Amk, Ask, Zmk, Zsk
  int I_par, N_fp = 0;

  // make new component
  component * cp_t = make_base_component(NPAR_EXPDISK2D);

  // parameters into component
  for(I_par = 0; I_par < NPAR_EXPDISK2D; ++ I_par)
    *(cp_t -> par + I_par) = *(par + I_par),
    *(cp_t -> is_const + I_par) = *(is_const + I_par);

  // range of parameters
  for(I_par = 0; I_par < 2 * NPAR_EXPDISK2D; ++ I_par)
    *(cp_t -> par_lim + I_par) = *(par_lim + I_par);

  // number of free params
  for(I_par = 0; I_par < NPAR_EXPDISK2D; ++ I_par)
    if(! *(is_const + I_par)) ++ N_fp;
  cp_t -> N_fp = N_fp;

  // set component nane and param name
  strcpy(cp_t -> name, name);
  cp_t -> par_name = (char **) expdisk_par_names;

  // set function headers
  cp_t -> sigma = & _cps_expdisk_2d_sigma;
  cp_t -> recipe = & _cps_expdisk_2d_recipe;

  return cp_t;
}

// surface density of a 2d exponential disk
double
_cps_expdisk_2d_sigma(double x, double y, double const * par)
{
  // rotate back
  double r_t = sqrt(Sq(x) + Sq(y));
  double cos_rt = x / r_t,  sin_rt = y / r_t,
         cos_pt = cos(PHI), sin_pt = sin(PHI);
  if (!isfinite(cos_rt)) cos_rt = 1., sin_rt = 0.;
  double x_t = r_t * (cos_rt * cos_pt + sin_rt * sin_pt), // cos(rt - pt)
         y_t = r_t * (sin_rt * cos_pt - cos_rt * sin_pt); // sin(rt - pt)

  // generalized radius
  double r = pow( pow(fabs(x_t - XC),     C)
                + pow(fabs(y_t - YC) / Q, C), 1. / C);

  // surface density
  return exp(-r / A) * IC;
}

int
_cps_expdisk_2d_sigma_arr(int N_pt, double * x, double * y,
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

  // print_head(r_t, 9);

  FOREACH(I_pt, N_pt)
    {
      if(r_t[I_pt] != 0.)
        cos_rt[I_pt] = x[I_pt] / r_t[I_pt],
        sin_rt[I_pt] = y[I_pt] / r_t[I_pt];
      else cos_rt[I_pt] = 1., sin_rt[I_pt] = 0.;
    }

  /*
  print_head(cos_rt, 9);
  print_head(sin_rt, 9);
  */

  FOREACH(I_pt, N_pt)
    {
      x_t[I_pt] = r_t[I_pt] * (cos_rt[I_pt] * cos_pt + sin_rt[I_pt] * sin_pt);
      y_t[I_pt] = r_t[I_pt] * (sin_rt[I_pt] * cos_pt - cos_rt[I_pt] * sin_pt);
    }

  /*
  print_head(x_t, 9);
  print_head(y_t, 9);
  */

  // generalized radius
  double * r = TALLOC(double, N_pt);

  FOREACH(I_pt, N_pt)
    r[I_pt] = pow( pow(fabs(x_t[I_pt] - XC),     C)
                 + pow(fabs(y_t[I_pt] - YC) / Q, C), 1. / C);

  //print_head(r, 9);
  //printf("Xc: %f, Yc: %f, C: %f, Q: %f\n", XC, YC, C, Q);

  //save_array(r, 257 * 257, "xt.tmp");
  //draw_image("xt.tmp", 257, 257);

  // surface density
  FOREACH(I_pt, N_pt) sigma[I_pt] = exp(-r[I_pt] / A) * IC;

  //print_head(sigma, 9);

  // free objects
  free(r_t), free(cos_rt), free(sin_rt), free(x_t), free(y_t), free(r);
}

// recipe of a 2d exponential disk
int
_cps_expdisk_2d_recipe(double x, double y,
    recipe * rcp_t, double const * par)
{
  // rotate back
  double r_t = sqrt(Sq(x) + Sq(y));
  double cos_rt = x / r_t,  sin_rt = y / r_t,
         cos_pt = cos(PHI), sin_pt = sin(PHI);
  if (!isfinite(cos_rt)) cos_rt = 1., sin_rt = 0.;
  double x_t = r_t * (cos_rt * cos_pt + sin_rt * sin_pt), // cos(rt - pt)
         y_t = r_t * (sin_rt * cos_pt - cos_rt * sin_pt); // sin(rt - pt)

  // generalized radius and surface density
  double r = pow( pow(fabs(x_t - XC),     C)
                + pow(fabs(y_t - YC) / Q, C), 1. / C);
  double sigma = exp(-r / A) * IC;

  // grid points for recipe
  int I_age, I_Z, N_age = rcp_t -> N_age, N_Z = rcp_t -> N_Z;
  double * age_ax = rcp_t -> age_ax, * Z_ax = rcp_t -> Z_ax;

  // mean and std. of age and metallicity
  double c_age = AMC + AMK * r, s_age_sq = Sq(ASC + ASK * r),
         c_Z   = ZMC + ZMK * r, s_Z_sq   = Sq(ZSC + ZSK * r);

  /* //debug
  printf("c_age = %f, s_age_sq = %f, c_Z = %f, s_Z_sq = %f\n",
         c_age, s_age_sq, c_Z, s_Z_sq);
  */

  print_head(age_ax, N_age),
  print_head(Z_ax, N_Z);

  // make recipe!
  double * rc = rcp_t -> rcp;
  for(I_age = 0; I_age < N_age; ++ I_age)
    for(I_Z = 0; I_Z < N_Z; ++ I_Z)
      rc[I_age * N_Z + I_Z] =
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

// undefine them
#undef IC
#undef XC
#undef YC
#undef PHI
#undef A
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
// end of expdisk.c
