
#include "utils.h"
#include "spectrum.h"

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
  sp_t -> err = TALLOC(double, lib_t -> N_spx);

  // fill with zeros
  int I_spx;
  for(I_spx = 0; I_spx < lib_t -> N_spx; ++ I_spx)
    sp_t -> flux[I_spx] = 0., sp_t -> err[I_spx] = 0.;
    
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

  return 0;
}
