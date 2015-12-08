
#include "synthesis.h"
#include "recipe.h"

#include "assert.h"

spectrum *
synthesis(recipe * rcp_t, spec_lib * lib_t, double x, double y);
{
  // make a new spectrum
  spectrum * sp_t = make_empty_spectrum(lib_t);

  // call noalloc
  synthesis_noalloc(rcp_t, lib_t, x, y, sp_t);

  // return spectrum
  return sp_t;
}

int
synthesis_noalloc(recipe * rcp_t, spec_lib * lib_t,
    double x, double y, spectrum * sp_t)
{
  // sanity check
  assert(lib_t -> N_spx == sp_t -> N_spx),
  assert(rcp_t -> N_age == lib_t -> N_age),
  assert(rcp_t -> N_Z == lib_t -> N_Z);

  // useful variables
  int I_spx, I_age, I_Z, N_spx = lib_t -> N_spx,
      N_age = lib_t -> N_age, N_Z = lib_t -> N_Z;

  // library, recipe and target
  double * spl = lib_t -> data, * rcp = rcp_t -> rcp;
      * flux = sp_t -> flux, * err = sp_t -> err;

  // fill flux and err with zeros
  for(I_spx = 0; I_spx < N_spx; ++ I_spx) flux[I_spx] = 0.;
  for(I_spx = 0; I_spx < N_spx; ++ I_spx) err[I_spx] = 0.;

  // loop over all age and Z /* TODO performance test ! */
  double * spl_i, wt_i;
  for(I_age = 0; I_age < N_age; ++ I_age)
    for(I_Z = 0; I_Z < N_Z; ++ I_Z)
      {
        wt_i = rcp[I_Z + I_age * N_Z],
        spl_i = spl + N_spx * (I_age * N_Z + I_Z);
        for(I_spx = 0; I_spx < N_spx; ++ I_spx)
          flux[I_spx] += spl_i[I_spx] * wt_i;
      }

  // done,
  return 0;
}
