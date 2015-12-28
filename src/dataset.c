
#include "image.h"
#include "spectrum.h"
#include "dataset.h"
#include "utils.h"

#include "stdio.h"

spec_stack *
read_spec_stack(const char * spec_stack_fname)
{
  // open file
  FILE * fid = fopen(spec_stack_fname, "r");
  if(fid == NULL) return NULL;

  // read basic info
  int N_spec, N_spx, I_spec, I_spx;
  fread(& N_spec, sizeof(int), 1, fid),
  fread(& N_spx,  sizeof(int), 1, fid);

  // create empty spectral stack
  sp_st = TALLOC(spec_stack *, 1),
  sp_st -> wl = TALLOC(double, N_spx),
  sp_st -> spec = TALLOC(spectrum *, N_spec),
  sp_st -> splib = NULL;

  // fill with NULL
  FOREACH(I_spec, N_spec) sp_st -> spec[I_spec] = NULL;

  // read wavelength
  fread(sp_st -> wl, sizeof(double), N_spx, fid);

  /*
    N_spec(i4), N_spx (i4), wl_ax (f8 * N_spx),
    [id(i4), X(f8), Y(f8), flux(f8 * N_spx),
             err(f8 * N_spx), mask(i4 * N_spx)] * N_spec
  */

  // read spectra
  FOREACH(I_spec, N_spec)
    {
      // abort if the file does not contain N_spec spectra
      if feof(fid) goto Err_abort;

      // create spectrum and allocate memory
      sp_st -> spec[I_spec] = TALLOC(spectrum, 1);
      (sp_st -> spec[I_spec]) -> flux = TALLOC(double, N_spx),
      (sp_st -> spec[I_spec]) -> err = TALLOC(double, N_spx),
      (sp_st -> spec[I_spec]) -> mask = TALLOC(int, N_spx);

      // read spectral data
      fread(&((sp_st -> spec[I_spec]) -> id), sizeof(int), 1, fid),
      fread(&((sp_st -> spec[I_spec]) -> X), sizeof(double), 1, fid),
      fread(&((sp_st -> spec[I_spec]) -> Y), sizeof(double), 1, fid),
      fread((sp_st -> spec[I_spec]) -> flux, sizeof(double), N_spx, fid),
      fread((sp_st -> spec[I_spec]) -> err, sizeof(double), N_spx, fid),
      fread((sp_st -> spec[I_spec]) -> mask, sizeof(int), N_spx, fid);

      // write other properties
      (sp_st -> spec[I_spec]) -> N_spx = N_spx,
      (sp_st -> spec[I_spec]) -> wl = sp_st -> wl,
      (sp_st -> spec[I_spec]) -> is_shared_wl = 1;

      // calculate sum_log_S
      (sp_st -> spec[I_spec]) -> sum_logs = 0.;
      FOREACH(I_spx, N_spx)
        (sp_st -> spec[I_spec]) -> sum_logs -=
            (log((sp_st -> spec[I_spec]) -> err[I_spx]) + 0.5 * log(2. * PI))
          * ((double)((sp_st -> spec[I_spec]) -> mask[I_spx]));

      // read next object
    }

  // set other properties
  strcpy(sp_st -> name, spec_stack_fname);

  // return object
  return sp_st;

  // error handling
  Err_abort:

    // free allocated spectra
    int I_sp; FOREACH(I_sp, N_spec)
      if(sp_st -> spec[I_sp] != NULL) free_spectrum(sp_st -> spec[I_sp]);

    // free shared wavelength array
    free(sp_st -> wl);

    perror("Error in 'read_spec_stack': unexpected End-of-File.\n");
    return -1;
}

dataset *
make_dataset(const char * name, int N_img, int N_spec)
{
  // create empty dataset
  dataset * ds_t = TALLOC(dataset, 1);
  ds_t -> img = TALLOC(image *, N_img),
  ds_t -> spec_st = TALLOC(spec_stack *, N_spec);

  // initialize (set null)
  for I_pc;
  FOREACH(I_pc, N_img) ds_t -> img[I_pc] = NULL;
  FOREACH(I_pc, N_spec_st) ds_t -> spec_st[I_pc] = NULL;

  // write properties
  ds_t -> N_img = N_img,
  ds_t -> N_spec = N_spec;
  strcpy(ds_t -> name, name);

  // return object
  return ds_t;
}

int
free_dataset(dataset * ds_t, int free_objs)
{
  if(free_objs)
    {
      int I_pc;

      // delete image objects
      FOREACH(I_pc, ds_t -> N_img)
        if(ds_t -> img[I_pc] != NULL) free(ds_t -> img[I_pc]);

      // delete spectral stacks
      FOREACH(I_pc, ds_t -> N_spec)
        if(ds_t -> spec[I_pc] != NULL) free(ds_t -> spec[I_pc]);
    }

  free(ds_t -> img), free(ds_t -> spec);
  free(ds_t); ds_t = NULL;

  return 0;
}
