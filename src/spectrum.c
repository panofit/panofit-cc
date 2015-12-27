
#include "utils.h"
#include "spectrum.h"

#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

#include "stdio.h"
#include "stdlib.h"
#include "string.h"

spectrum *
make_empty_spectrum(spec_lib * lib_t)
{
  // allocate empty spectrum
  spectrum * sp_t = TALLOC(spectrum, 1);

  // set properties
  sp_t -> N_spx = lib_t -> N_spx,
  sp_t -> wl = lib_t -> wl,
  sp_t -> is_shared_wl = 1;

  // allocate space for flux and err
  sp_t -> flux = TALLOC(double, lib_t -> N_spx),
  sp_t -> err  = TALLOC(double, lib_t -> N_spx),
  sp_t -> mask = TALLOC(int, lib_t -> N_spx);

  // fill with zeros
  int I_spx;
  //for(I_spx = 0; I_spx < lib_t -> N_spx; ++ I_spx)
  FOREACH(I_spx, (lib_t -> N_spx))
    sp_t -> flux[I_spx] = 0.,
    sp_t -> err[I_spx]  = 0.,
    sp_t -> mask[I_spx] = 1 ;

  // return the pointer
  return sp_t;
}

int
free_spectrum(spectrum * sp_t)
{
  free(sp_t -> flux),
  free(sp_t -> err);

  if(! (sp_t -> is_shared_wl))
    free(sp_t -> wl);

  sp_t = NULL;

  return 0;
}

spec_lib *
load_spec_lib_raw(const char * spec_file)
{
  // open file,
  FILE * spl_fp = fopen(spec_file, "rb");
  if(spl_fp == NULL) return NULL;

  // values to read
  char name[128]; // name of the spec lib
  int N_age, N_Z, N_spx; // grid pts in age, metallicity, N of spectral px.

  // read header
  fread(name, sizeof(char), 128, spl_fp),
  fread(& N_Z, sizeof(int), 1, spl_fp),
  fread(& N_age, sizeof(int), 1, spl_fp),
  fread(& N_spx, sizeof(int), 1, spl_fp);

  // allocate spec lib object
  spec_lib * lib_t = TALLOC(spec_lib, 1);
  lib_t -> age_ax  = TALLOC(double, N_age),
  lib_t -> Z_ax    = TALLOC(double, N_Z),
  lib_t -> wl      = TALLOC(double, N_spx);
  lib_t -> data    = TALLOC(double, N_age * N_Z * N_spx);

  // read arrays
  fread(lib_t -> age_ax, sizeof(double), N_age, spl_fp),
  fread(lib_t -> Z_ax,   sizeof(double), N_Z, spl_fp),
  fread(lib_t -> wl,     sizeof(double), N_spx, spl_fp),
  fread(lib_t -> data,   sizeof(double), N_age * N_Z * N_spx, spl_fp);

  // write properties
  lib_t -> N_age = N_age,
  lib_t -> N_Z = N_Z,
  lib_t -> N_spx = N_spx;
  strcpy(lib_t -> name, name);

  // return object
  return lib_t;
}

int free_spec_lib(spec_lib * lib_t)
{
  free(lib_t -> wl),
  free(lib_t -> Z_ax),
  free(lib_t -> age_ax),
  free(lib_t -> data);
  lib_t =  NULL;

  return 0;
}

spectrum *
make_empty_spectrum_as(spectrum * sp_i)
{
  int I_wl;

  // allocate empty
  spectrum * sp_t = TALLOC(spectrum, 1);
  sp_t -> flux = TALLOC(double, sp_i -> N_spx),
  sp_t -> err = TALLOC(double, sp_i -> N_spx),
  sp_t -> mask = TALLOC(int, sp_i -> N_spx);

  // copy parameters
  sp_t -> N_spx = sp_i -> N_spx;

  // wavelength
  sp_t -> shared_wl = sp_i -> shared_wl;
  if(! (sp_i -> shared_wl))
    {
      sp_t -> wl = TALLOC(double, sp_i -> N_spx);
      FOREACH(I_wl, sp_i -> N_spx)
        sp_t -> wl[I_wl] = sp_i -> wl[I_wl];
    }

  // copy mask
  FOREACH(I_wl, sp_i -> N_spx)
    sp_t -> mask[I_wl] = sp_i -> mask[I_wl];

  return sp_t;
}

int
resample_spectrum(spectrum * sp_t, spectrum sp_i)
{
  // create intepolator
  gsl_interp_accel * ac_t = gsl_interp_accel_alloc();
  gsl_spline * spl_t = gsl_spline_alloc(gsl_interp_cspline, sp_i -> N_spx);
  gsl_spline_init(spl_t, sp_i -> wl, sp_i -> flux, sp_i -> N_spx);

  // do interpolation
  int I_spx;
  FOREACH(I_spx, sp_t -> N_spx)
    sp_t -> flux[I_spx] = gsl_spline_eval(spl_t, sp_t -> wl[I_spx], ac_t);

  // fill mask area
  // TODO

  // free interpolator
  gsl_spline_free(spl_t);
  gsl_interp_accel_free(ac_t);
}

spectrum *
specfilter_gaussian_evenspacing(spectrum * sp_i, double lw, int N_spx_c)
{
  // make empty spectrum
  spectrum * sp_t = make_empty_spectrum_as(sp_i);

  // do gaussian filter
  // TODO, 1d horizontal gaussian smoothing.

  return sp_t;
}
