
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

// generate a trial solution for MCMC fitting. \
   ps_d for destination, ps_s for previous solution, T for temperature.
int
sample_param_set(param_set * ps_d, param_set * ps_s, double T)
{
  // for the same model?
  if(ps_d -> N_cps == ps_s -> N_cps) return -1;
  if(ps_d -> shared_with == ps_s -> shared_with) return -4;

  // temporatory parameters
  double par_t; int is_valid = 0, I_cp, I_par

  // iterate over components
  FOREACH(I_cp, N_cps)
    // generate params for this component until a valid set was found
    do
      {
        // iterate over parameters
        FOREACH(I_par, ps_s -> N_par[I_cp])
          {
            // skip if this is a const param
            if(ps_s -> is_const[I_par]) continue;

            // generate a new value
            do par_t = fabs((ps_s -> is_const[I_par * 2 + 1])
                           -(ps_s -> is_const[I_par * 2])) * rand_gauss() * T
                     + (ps_s -> par[I_cp]);

            // make sure the value is valid
            while((par_t > (ps_s -> is_const[I_par * 2])) &&
                  (par_t < (ps_s -> is_const[I_par * 2 + 1])));
          }
      }
    while(); // TODO add boundary check function header in param set
}
