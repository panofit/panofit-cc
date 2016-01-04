
#include "utils.h"
#include "mcfit.h"
#include "spectrum.h"
#include "image.h"
#include "model.h"
#include "recipe.h"
#include "synthesis.h"
#include "mcfit.h"
#include "bayesian.h"

#include "assert.h"
#include "math.h"
#include "stdio.h"

#define Sq(X) ((X) * (X))

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
  int I_step, I_sp, I_sp_st, N_sp_st = ds_t -> N_spec_st, N_spx;
  spec_st * stack_t; spectrum * spec_t; recipe * rcp_t;
  int if_accept;

  // parameter sets
  param_set * ps_current = make_param_set_for(m_t),
            * ps_trial = make_param_set_for(m_t),
            * ps_best = make_param_set_for(m_t);
  param_set * ps_param = get_assoc_param_set(m_t);

  // corresponding likelihood
  double ln_L_current = -DBL_MAX, ln_L_trial, ln_L_best;

  // initialize: current = ps_param
  copy_param_set(ps_current, ps_param);

  // MCMC loop
  FOREACH(I_step, N_steps)
    {
      // new step, clear temp variables.
      ln_L_trial = 0;

      // generate new solution from existing solution
      sample_param_set(ps_trial, ps_current, T);

      // write new param set to the model
      copy_param_set(ps_param, ps_trial);

      // test likelihood: spectra, iterate over spec stacks
      FOREACH(I_sp_st, N_sp_st)
        {
          // which stack?
          stack_t = ds_t -> spec_st[I_sp_st];
          N_spec = stack_t -> N_spec, N_spx = stack_t -> N_spx;

          // make mock spectrum object
          spectrum * mock_sp_t =
              make_empty_spectrum_as((stack_t -> spec[I_sp])[0]);

          // iterate over individual spectra
          FOREACH(I_sp, N_spec)
            {
              // which spectrum?
              spec_t = stack_t -> spec[I_sp];

              // generate recipe for it.
              sample_recipe_noalloc(m_t, spec_t -> X, spec_t -> Y, rcp_t);

              // generate mock spectrum at this point
              synthesis_noalloc(rcp_t, stack_t -> splib, mock_sp_t);

              // calculate likelihood
              ln_L_trial += spec_ln_L(mock_sp_t, spec_t);
            }

          // TODO free mock spectrum object
          free_spectrum(mock_sp_t);
        }

      // iterate over images;
      // FOREACH(...)

      // accept or not?
      if(ln_L_trial > ln_L_current) if_accept = 1;
      else if_accept = rand_uniform() < exp(ln_L_trial - ln_L_current);

      if(if_accept)
        {
          // DEBUG
          print_param_set(param_set * param_trial, 1);

          // TODO: write chain

          // TODO: write log

          // copy the new solution into current
          copy_param_set(param_current, param_trial);
          ln_L_current = ln_L_trial;
        }
      else
        {
          // DEBUG
          print_param_set(param_set * param_trial, 1);

          // TODO: write log

          // rewind I_step
          -- I_step;
        }

    }
}
