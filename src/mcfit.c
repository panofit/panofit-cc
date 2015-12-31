
#include "utils.h"
#include "mcfit.h"
#include "spectrum.h"
#include "image.h"
#include "model.h"

#include "assert.h"
#include "math.h"
#include "stdio.h"

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
    log_L += Sq((obs_flux[I_spx] - mock_flux[I_spx]) / obs_err[I_spx])
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
  if(ps_d -> N_cps != ps_s -> N_cps) return -1;
  if(ps_d -> shared_with != ps_s -> shared_with) return -4;

  // temporatory parameters
  double par_t; int is_valid = 0, I_cp, I_par, N_cps = ps_s -> N_cps;

  // iterate over components
  FOREACH(I_cp, N_cps)
    // generate params for this component until a valid set was found
    do
      {
        // iterate over parameters
        FOREACH(I_par, ps_s -> N_par[I_cp])
          {
            // skip if this is a const param
            if((ps_s -> is_const[I_cp])[I_par]) continue;

            // generate a new value
            do par_t = fabs(((ps_s -> par_lim[I_cp])[I_par * 2 + 1])
                     -((ps_s -> par_lim[I_cp])[I_par * 2])) * rand_gauss() * T
                     + (ps_s -> par[I_cp])[I_par];

            // make sure the value is valid
            while((par_t > (ps_s -> is_const[I_cp])[I_par * 2]) &&
                  (par_t < (ps_s -> is_const[I_cp])[I_par * 2 + 1]));

            // write par_t to param set
            (ps_d -> par[I_cp])[I_par] = par_t;
          }
      }
    while((ps_d -> check[I_cp])(ps_d -> par[I_cp]));
}

/*
  Run Monte-carlo fitting
    [IN]
      model * m_t: configuration of the model
      dataset * ds_t: dataset (spectra, images) to fit,
      param_set * ps_t: param set associated with the model
      int N_steps: number of steps for the chain
    [OUT]
      param_set * ps_best: best-fitting parameter set,
      FILE * chain_rec: the MCMC sampling set
      FILE * log_rec: used to record something.
    [RET]
      1: err abort
      0: success.
*/
int
metropolis_sampling(model * m_t, dataset * ds_t, param_set * ps_t, int N_steps,
    double T, param_set * ps_best, FILE * chain_rec, FILE * log_rec)
{
  // sanity check

  // useful variables
  int I_step, I_sp_st, N_sp_st = ds_t -> N_spec_st;

  // get prepared
  param_set * ps_current = make_param_set_for(m_t),
            * ps_trial = make_param_set_for(m_t);

  FOREACH(I_step, N_steps)
    {
      // generate new solution from existing solution
      sample_param_set(ps_new, ps_t, T);

      // test likelihood: spectra, iterate over spec stacks
      FOREACH(I_sp_st, N_sp_st)
    }
    // generate new solution
    // test likelihood
    // accept or reject?
      // accept: write log,
      //         if best solution?
    // next step
}
