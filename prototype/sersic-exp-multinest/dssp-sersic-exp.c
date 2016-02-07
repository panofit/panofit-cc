
#include "stddef.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "float.h"
#include "string.h"
#include "time.h"

#include "multinest.h"

typedef double realn;

#define TALLOC(TYPE, NUM) ((TYPE *) malloc(sizeof(TYPE) * (NUM)))
#define FOREACH(IDX, NUM) for(IDX = 0; IDX < (NUM); ++ IDX)
#define Sq(X) ((X) * (X))

struct _sersic_exp_wspac
{
  // param range & lower limit
  double par_h[15], par_a[15];

  // observed data.
  int Nx, Ny, Nwl;
  realn * X_ax, * Y_ax, * wl_ax;
  realn * Val, * Err; int * Msk;

  // SSP model
  int N_age; realn * age_ax, * ssp;

  // PSF settings
  int use_psf;
  double psf_radius;

  // pre-computed
  int N_sp, N_vol;
  realn * X_axt, * Y_axt;

  // variables for model generation
  realn * Xt, * Yt, * Rt, * cos_t, * sin_t;
  realn * Xs, * Ys, * Ms, * Ks, * Xe, * Ye, * Me, * Ke;
  realn * Fs, * Fe;

  // mock dataws
  realn * Val_fit, * chisq;

  // save name
  char dump_fname[128];
};

int print_array(realn * arr, int N)
{
  int I;
  FOREACH(I, N) printf("% 3u, % f\n", I, arr[I]);
  return 0;
}

int save_array(realn * arr, int N, const char * filename)
{
  FILE * fp = fopen(filename, "wb");
  if (fp == NULL) return -1;
  fwrite (arr, sizeof(realn), N, fp);
  fclose (fp);
  return 0;
}

int _axbsch(int N, realn * X, realn v)
{
  // indices for the left/right/middle point
  int i = 0, j = N - 1;
  int k = (i + j) / 2;

  // out-of-range check
  if(v <  X[i]) return 0;
  if(v >= X[j]) return N;

  // start binary search
  while(j - i > 1)
  {
    if(v <  X[k]) j = k;
    else i = k;
    k = (i + j) / 2;
  }

  // return the right point idx
  return j;
}


// a quick and dirty 2d Gaussian filter
int
gaussian_filter(realn * img, int Nx, int Ny, realn sigma, int N_steps)
{
  realn   delta, bscl, pscl, * ptr;
  int     N_sp = Nx * Ny, It, Ix, Iy, Ist;

  if(sigma <= 0 || N_steps < 0) return -1;

  double lambda = (sigma * sigma) / (2.0 * N_steps);
  double d_delta = (1. + 2. * lambda - sqrt(1. + 4. * lambda)) / (2. * lambda);

  delta = (realn) d_delta;
  bscl  = (realn) (1. / (1. - d_delta));
  pscl  = (realn) (pow(d_delta / lambda, 2 * N_steps));

  FOREACH(Iy, Ny) FOREACH(Ist, N_steps)
    {
      ptr = img + Nx * Iy;
      ptr[0] *= bscl; for(Ix = 1; Ix < Nx; Ix++) ptr[Ix] += delta * ptr[Ix - 1];
      ptr[Ix = Nx - 1] *= bscl; for(; Ix > 0; Ix--) ptr[Ix - 1] += delta * ptr[Ix];
    }

  FOREACH(Ix, Nx) FOREACH(Ist, N_steps)
    {
      ptr = img + Ix;
      ptr[0] *= bscl; for(It = Nx; It < N_sp; It += Nx) ptr[It] += delta * ptr[It - Nx];
      ptr[It = N_sp - Nx] *= bscl; for(; It > 0; It -= Nx) ptr[It - Nx] += delta * ptr[It];
    }

  FOREACH(It, N_sp) img[It] *= pscl;

  return 0;
}


/*
  Parameters follow different symbols than the python version.

  Sersic component:
    Ns, Rs, Ks, Ps, Qs, Cs, As: Sersic idx, effective radius, surface brightness at
        the effective radius, position angle, axis ratio b/a, isophot shape C, log SSP age

  Exponential component:
    Ke, Re, Pe, Qe, Ce, Ae: surface brightness at the center, scale length, position angle,
        axis ratio b/a, isophot shape, log SSP age

  void LogLike(double *Cube, int *ndim, int *npars, double *lnew, void *wspac)
*/
void _sersic_exp_lnp(double * p, int * N_dim, int * N_pars, double * L_new, void * wspac)
{
  int It, Iw, Ix, Iy;

  // unpack workspace
  struct _sersic_exp_wspac * ws = (struct _sersic_exp_wspac *) wspac;

  // X, Y and lambda axes.
  int Nx = ws -> Nx, Ny = ws -> Ny, Nwl = ws -> Nwl;
  realn * X_ax = ws -> X_ax, * Y_ax = ws -> Y_ax, * wl_ax = ws -> wl_ax;

  // datacube: flux, error, mask bit
  realn * Val = ws -> Val, * Err = ws -> Err;
  int * Msk = ws -> Msk;

  // used to convert unit cube to param space (assuming flat prior)
  double * par_a = ws -> par_a, * par_h = ws -> par_h;

  // fitting result
  realn * Val_fit = ws -> Val_fit, * chisq = ws -> chisq;

  // SSP model
  int N_age = ws -> N_age;
  realn * age_ax = ws -> age_ax, * ssp = ws -> ssp;

  // PSF model
  int use_psf = ws -> use_psf;
  double psf_radius = ws -> psf_radius;

  // pre-computed
  int N_sp = ws -> N_sp, N_vol = ws -> N_vol;
  realn * X_axt = ws -> X_axt, * Y_axt = ws -> Y_axt;

  // variables for model generation
  realn * Xt = ws -> Xt, * Yt = ws -> Yt, * Rt = ws -> Rt,
      * cos_t = ws -> cos_t, * sin_t = ws -> sin_t;
  realn * Xs = ws -> Xs, * Ys = ws -> Ys, * Ms = ws -> Ms, * Ks = ws -> Ks,
      * Xe = ws -> Xe, * Ye = ws -> Ye, * Me = ws -> Me, * Ke = ws -> Ke;
  realn * Fs = ws -> Fs, * Fe = ws -> Fe;

  // convert from unit cube to physical parameters, linear scaling
  FOREACH(It, 15) p[It] = (ws -> par_a[It]) + p[It] * (ws -> par_h[It]);

  // unpack parameters
  realn Xc = p[ 0], Yc = p[ 1], Ns = p[ 2], Rs = p[ 3], Is = p[ 4], Ps = p[ 5],
        Qs = p[ 6], Cs = p[ 7], As = p[ 8], Ie = p[ 9], Re = p[10], Pe = p[11],
        Qe = p[12], Ce = p[13], Ae = p[14];

  // calculate shifted coordinates
  FOREACH(It, N_sp) Xt[It] = X_axt[It] - Xc; //v
  FOREACH(It, N_sp) Yt[It] = Y_axt[It] - Yc; //v
  FOREACH(It, N_sp) Rt[It] = sqrt(Xt[It] * Xt[It] + Yt[It] * Yt[It]);
  FOREACH(It, N_sp) cos_t[It] = Xt[It] / Rt[It]; //v
  FOREACH(It, N_sp) sin_t[It] = Yt[It] / Rt[It]; //v

  // calculate the sersic component
  realn Ns_sq = Ns * Ns;
  realn d_n = -(3. * Ns - 0.3333333333333333 + 0.006584362139917695 / Ns
            + 0.0008012715831645873 / Ns_sq + 3.3805660996637985e-05 / (Ns_sq * Ns)
            - 1.4125456386947105e-05 / (Ns_sq * Ns_sq)); // FIXME: Ref?
  realn cos_s = cos(Ps), sin_s = sin(Ps), Is_t = pow(2.512, -Is), inv_Cs = 1. / Cs, inv_Ns = 1. / Ns;
  FOREACH(It, N_sp) Xs[It] = Rt[It] * (cos_t[It] * cos_s + sin_t[It] * sin_s); //v
  FOREACH(It, N_sp) Ys[It] = Rt[It] * (sin_t[It] * cos_s - cos_t[It] * sin_s); //v
  FOREACH(It, N_sp) Ms[It] = pow(pow(fabs(Xs[It]), Cs) + pow(fabs(Ys[It]) / Qs, Cs), inv_Cs);
  FOREACH(It, N_sp) Ks[It] = Is_t * exp(d_n * (pow(Ms[It] / Rs, inv_Ns) - 1.));

  // calculate the exponential component
  realn cos_e = cos(Pe), sin_e = sin(Pe), Ie_t = pow(2.512, -Ie), inv_Ce = 1. / Ce;
  FOREACH(It, N_sp) Xe[It] = Rt[It] * (cos_t[It] * cos_e + sin_t[It] * sin_e); //v
  FOREACH(It, N_sp) Ye[It] = Rt[It] * (sin_t[It] * cos_e - cos_t[It] * sin_e); //v
  FOREACH(It, N_sp) Me[It] = pow(pow(fabs(Xe[It]), Ce) + pow(fabs(Ye[It]) / Qe, Ce), inv_Ce);
  FOREACH(It, N_sp) Ke[It] = Ie_t * exp(-Me[It] / Re);

  int id_a, id_b; realn u_t, v_t; ptrdiff_t os_a, os_b, os_t;

  // interpolate two SSPs
  id_b = _axbsch(N_age, age_ax, As);
  if(id_b % N_age == 0) u_t = 1., v_t = 0., id_a = (id_b == N_age)?(N_age - 1):(0), id_b = id_a;
  else id_a = id_b - 1, u_t = (age_ax[id_b] - As) / (age_ax[id_b] - age_ax[id_a]), v_t = 1. - u_t;
  os_a = id_a * Nwl, os_b = id_b * Nwl;
  FOREACH(It, Nwl) Fs[It] = ssp[os_a + It] * u_t + ssp[os_b + It] * v_t; //v

  id_b = _axbsch(N_age, age_ax, Ae);
  if(id_b % N_age == 0) u_t = 1., v_t = 0., id_a = (id_b == N_age)?(N_age - 1):(0), id_b = id_a;
  else id_a = id_b - 1, u_t = (age_ax[id_b] - Ae) / (age_ax[id_b] - age_ax[id_a]), v_t = 1. - u_t;
  os_a = id_a * Nwl, os_b = id_b * Nwl;
  FOREACH(It, Nwl) Fe[It] = ssp[os_a + It] * u_t + ssp[os_b + It] * v_t; //v

  // find the final ws
  FOREACH(Iw, Nwl)
    FOREACH(It, N_sp) //v
      Val_fit[Iw * N_sp + It] = Fe[Iw] * Ke[It] + Fs[Iw] * Ks[It];

  // apply gaussian filter
  if(use_psf)
    FOREACH(Iw, Nwl)
      gaussian_filter(Val_fit + Iw * N_sp, Nx, Ny, psf_radius, 3);

  // calculate a chi-square
  FOREACH(It, N_vol) chisq[It]  = Val_fit[It] - Val[It]; //v
  FOREACH(It, N_vol) chisq[It] /= Err[It]; //v
  FOREACH(It, N_vol) chisq[It] *= chisq[It]; //v

  // calculate chi-sq sum
  realn chisq_sum = 0.;
  FOREACH(It, N_vol) chisq[It] *= Msk[It]; //v
  FOREACH(It, N_vol) chisq_sum += chisq[It];

  printf("% e  ", -chisq_sum * 0.5);
  FOREACH(It, 15) printf("% .3e  ", p[It]);
  printf("\n");

  // return ln_prob
  * L_new = -0.5 * chisq_sum;
}

void _sersic_exp_dump(int * N_samples, int * N_live, int * N_par, double ** phys_live, double ** posterior, double ** param_constr, double * max_ln_L, double * ln_Z, double * IS_ln_Z, double * ln_Z_err, void * wspac)
{
  // unpack workspace
  struct _sersic_exp_wspac * ws = (struct _sersic_exp_wspac *) wspac;

  // open and write posterior distributions
  char save_filename[128];
  sprintf(save_filename, "%s-samples", ws -> dump_fname);
  FILE * fpa = fopen(save_filename, "wb");
  if(fpa == NULL) {perror("Unable to dump file!"); return;}
  fwrite(* posterior, sizeof(double), (* N_samples) * (* N_par + 2), fpa); // written in fortran order!
  fclose(fpa);

  // open and write best-fitting models
  sprintf(save_filename, "%s-params", ws -> dump_fname);
  FILE * fpb = fopen(save_filename, "wb");
  if(fpb == NULL) {perror("Unable to dump file!"); return;}
  fwrite(* param_constr, sizeof(double), 4 * (* N_par), fpb); // written in fortran order!
  fclose(fpb);

  // printf("dumper: posterior sample and parameters saved.\n");

  return;
}

struct _sersic_exp_wspac * read_sersic_exp_wspac(const char * filename)
{
  // try to open the file
  FILE * fp = fopen(filename, "rb");
  if(fp == NULL) return NULL; // not accessable, return NULL

  // allocate space and read data
  struct _sersic_exp_wspac * ws = TALLOC(struct _sersic_exp_wspac, 1);

  // read array dimensions
  fread (& (ws -> Nx),  sizeof(int), 1, fp),
  fread (& (ws -> Ny),  sizeof(int), 1, fp),
  fread (& (ws -> Nwl), sizeof(int), 1, fp);

  // read axes
  ws -> X_ax  = TALLOC(realn, ws -> Nx);
  fread (ws -> X_ax, sizeof(realn), ws -> Nx, fp);
  ws -> Y_ax  = TALLOC(realn, ws -> Ny);
  fread (ws -> Y_ax, sizeof(realn), ws -> Ny, fp);
  ws -> wl_ax = TALLOC(realn, ws -> Nwl);
  fread (ws -> wl_ax, sizeof(realn), ws -> Nwl, fp);

  ws -> N_sp = (ws -> Nx) * (ws -> Ny);
  ws -> N_vol = (ws -> N_sp) * (ws -> Nwl);

  // read dataws to fit
  ws -> Val = TALLOC(realn, ws -> N_vol);
  fread (ws -> Val, sizeof(realn), ws -> N_vol, fp);
  ws -> Err = TALLOC(realn, ws -> N_vol);
  fread (ws -> Err, sizeof(realn), ws -> N_vol, fp);
  ws -> Msk = TALLOC(int, ws -> N_vol);
  fread (ws -> Msk, sizeof(int), ws -> N_vol, fp);

  // read the SSP model
  fread (& (ws -> N_age), sizeof(int), 1, fp);
  ws -> age_ax = TALLOC(realn, ws -> N_age);
  fread (ws -> age_ax, sizeof(realn), ws -> N_age, fp);
  ws -> ssp = TALLOC(realn, (ws -> N_age) * (ws -> Nwl));
  fread (ws -> ssp, sizeof(realn), (ws -> N_age) * (ws -> Nwl), fp);

  // read par range
  fread (ws -> par_a, sizeof(realn), 15, fp),
  fread (ws -> par_h, sizeof(realn), 15, fp);

  // read psf settings
  fread (& (ws -> use_psf), sizeof(int), 1, fp),
  fread (& (ws -> psf_radius), sizeof(double), 1, fp);

  // allocate mem
  ws -> Fs      = TALLOC(realn, ws -> Nwl  ),
  ws -> Fe      = TALLOC(realn, ws -> Nwl  ),
  ws -> Xs      = TALLOC(realn, ws -> N_sp ),
  ws -> Ys      = TALLOC(realn, ws -> N_sp ),
  ws -> Ms      = TALLOC(realn, ws -> N_sp ),
  ws -> Ks      = TALLOC(realn, ws -> N_sp ),
  ws -> Xe      = TALLOC(realn, ws -> N_sp ),
  ws -> Ye      = TALLOC(realn, ws -> N_sp ),
  ws -> Me      = TALLOC(realn, ws -> N_sp ),
  ws -> Ke      = TALLOC(realn, ws -> N_sp ),
  ws -> Xt      = TALLOC(realn, ws -> N_sp ),
  ws -> Yt      = TALLOC(realn, ws -> N_sp ),
  ws -> Rt      = TALLOC(realn, ws -> N_sp ),
  ws -> X_axt   = TALLOC(realn, ws -> N_sp ),
  ws -> Y_axt   = TALLOC(realn, ws -> N_sp ),
  ws -> cos_t   = TALLOC(realn, ws -> N_sp ),
  ws -> sin_t   = TALLOC(realn, ws -> N_sp ),
  ws -> chisq   = TALLOC(realn, ws -> N_vol),
  ws -> Val_fit = TALLOC(realn, ws -> N_vol);

  // fill X_axt and Y_axt
  int It;
  FOREACH(It, ws -> N_sp)
    ws -> X_axt[It] = ws -> X_ax[It % (ws -> Nx)],
    ws -> Y_axt[It] = ws -> Y_ax[(It / (ws -> Nx)) % (ws -> Ny)];

  fclose(fp);

  // write dunp filename
  sprintf(ws -> dump_fname, "%s-sav", filename);

  return ws;
}

int free_wspac(struct _sersic_exp_wspac * ws)
{
  if(ws == NULL) return 0;

  free (ws -> X_ax), free (ws -> Y_ax), free (ws -> wl_ax),
  free (ws -> Val), free (ws -> Err), free (ws -> Msk),
  free (ws -> age_ax), free (ws -> ssp),
  free (ws -> X_axt), free (ws -> Y_axt),
  free (ws -> Xt), free (ws -> Yt), free (ws -> Rt),
  free (ws -> cos_t), free (ws -> sin_t),
  free (ws -> Xs), free (ws -> Ys), free (ws -> Ms), free (ws -> Ks),
  free (ws -> Xe), free (ws -> Ye), free (ws -> Me), free (ws -> Ke),
  free (ws -> Fs), free (ws -> Fe),
  free (ws -> Val_fit), free (ws -> chisq),

  free (ws); ws = NULL;
  return 0;
}

int main(int argc, char ** argv)
{
  // check argc
  if (argc != 4) {perror("Run with: ./dssp-sersic-exp [wspac_file] [0] [0]"); return -1;}

  // import the workspace
  struct _sersic_exp_wspac * wspac = read_sersic_exp_wspac(argv[1]);
  if (wspac == NULL) {perror("Can't access the workspace file"); return -1;}

  int i;

  // setup and run the fitting module
  int     imp_sampl       =  1;
  int     is_multi_modal  =  1;
  int     is_const_eff    =  1;
  int     N_live_pts      =  1500;
  double  sampl_eff       =  0.8;
  double  tol             =  0.99; // important!
  int     N_dims          =  15;
  int     N_pars          =  15;
  int     N_clpar         =  15;
  int     update_intv     =  100;
  double  Z_tol           = -1E90;
  int     N_max_modes     =  100;
  int *   is_periodic     =  TALLOC(int, N_dims);
  char    root_path[100]  = "testrun-masswt-";
  int     rnd_seed        = -1;
  int     is_verbose      =  1;
  int     if_resume       =  1;
  int     if_write_file   =  1;
  int     if_init_MPI     =  0;
  double  ln_zero         = -DBL_MAX;
  int     N_iters_max     =  200000;

  // phi_s and phi_e are periodic
  FOREACH(i, N_dims) is_periodic[i] = 0;
  is_periodic[5] = 1, is_periodic[11] = 1;

  // make file name
  // sprintf(root_path, "%s-run-", argv[1]);
  // printf("Files write to: %s\n", root_path);


  run(imp_sampl, is_multi_modal, is_const_eff, N_live_pts, tol, sampl_eff,
      N_dims, N_pars, N_clpar, N_max_modes, update_intv, Z_tol,
      root_path, rnd_seed, is_periodic, is_verbose, if_resume, if_write_file, if_init_MPI,
      ln_zero, N_iters_max, _sersic_exp_lnp, _sersic_exp_dump, (void *) wspac);


  // a set of "good" parameters, for debug

  /*
  double p_opt[15] = {-1.383e+00, -2.430e+00, 5.176e+00, 2.411e+01, -4.643e+00, 4.948e-01,
                      6.910e-01, 1.788e+00, 6.014e-01, -8.481e+00, 1.079e+01, 4.124e-02,
                      7.278e-01, 1.893e+00, 2.807e-01 };

  // back convert to unit cube
  FOREACH(i, 15) p_opt[i] = (p_opt[i] - (wspac -> par_a[i])) / (wspac -> par_h[i]);

  // evaluate likelihood
  double L_tmp;
  _sersic_exp_lnp(p_opt, & N_dims, & N_pars, & L_tmp, (void *) wspac);
  printf("Getting ln_L:%f\n", L_tmp);

  // dump workspace for checking
  save_array(wspac -> ssp, 165 * 184, "ssp.tmp");
  save_array(wspac -> X_axt, wspac -> N_sp, "X_axt.tmp");
  save_array(wspac -> Y_axt, wspac -> N_sp, "Y_axt.tmp");
  save_array(wspac -> Xt, wspac -> N_sp, "Xt.tmp");
  save_array(wspac -> Yt, wspac -> N_sp, "Yt.tmp");
  save_array(wspac -> cos_t, wspac -> N_sp, "cos_t.tmp");
  save_array(wspac -> sin_t, wspac -> N_sp, "sin_t.tmp");
  save_array(wspac -> Ks, wspac -> N_sp, "Ks.tmp");
  save_array(wspac -> Ke, wspac -> N_sp, "Ke.tmp");
  save_array(wspac -> Xs, wspac -> N_sp, "Xs.tmp");
  ssave_array(wspac -> Ys, wspac -> N_sp, "Ys.tmp");
  save_array(wspac -> Xe, wspac -> N_sp, "Xe.tmp");
  save_array(wspac -> Ye, wspac -> N_sp, "Ye.tmp");
  save_array(wspac -> Ms, wspac -> N_sp, "Ms.tmp");
  save_array(wspac -> Me, wspac -> N_sp, "Me.tmp");
  save_array(wspac -> Val_fit, wspac -> N_vol, "val_fit.tmp");

  print_array(wspac -> X_ax, wspac -> Nx);
  print_array(wspac -> Y_ax, wspac -> Ny);
  */

  // clean
  free (is_periodic);
  free_wspac (wspac);
}
