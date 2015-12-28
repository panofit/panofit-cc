
#include "mcfit.h"

#include "assert.h"

#define Sq(X) ((X) * (X))

// calculates likelihood of model spectrum
double
spec_logL(spectrum * mock_sp, spectrum * obs_sp)
{
  // check
  assert(mock_sp -> N_spx == obs_sp -> N_spx);

  // useful variables
  int I_spx, N_spx = obs_sp -> N_spx;
  double log_L = 0., * mock_flux = mock_sp -> flux,
    * obs_flux = obs_sp -> flux, * obs_err = obs_sp -> err;
  int * mock_mask = mock_sp -> mask, * obs_mask = obs_sp -> mask;
  // FIXME mock_mask is not used now since spec_lib doesn't hold the mask. \
      I'm assuming a spectral library without bad points.

  // calculate chi_sq and turn into log_L
  FOREACH(I_spx, N_spx)
    log_L += Sq((obs_flux[I_spc] - mock_flux[I_spx]) / obs_err[I_spx])
          * (double)(obs_mask[I_spx]);
  log_L *= -0.5;

  // add the pre-computed const terms // optional
  log_L += obs_sp -> sum_logs;

  return log_L;
}

double
img_logL(image * mock_img, image * obs_img)
{
  return 0.; // TODO TODO
}
