
#include "utils.h"
#include "model.h"
#include "component.h"

#include "stdlib.h"
#include "stdarg.h"
#include "string.h"

#include "assert.h"

model *
make_model(const char * name, int N_cps, ...)
{
  // alloc model object
  model * m_t = TALLOC(model, 1);
  m_t -> cps = TALLOC(component *, N_cps);

  // assign components
  int I_cp;
  va_list cps_list;
  va_start(cps_list, N_cps);
  for(I_cp = 0; I_cp < N_cps; ++ I_cp)
    *(m_t -> cps + I_cp) = va_arg(cps_list, component *);
  va_end(cps_list);

  // do other things
  strcpy(m_t -> name, name);
  m_t -> N_cps = N_cps;

  // return model
  return m_t;
}

int
free_model(model * m_t, int free_cps)
{
  if(free_cps)
    {
      int I_cp, N_cps = m_t -> N_cps;
      for(I_cp = 0; I_cp < N_cps; ++ I_cp)
        ((m_t -> cps[I_cp]) -> kill)(m_t -> cps[I_cp]);
    }

  // free componen and return
  free(m_t -> cps), free(m_t);

  return 0;
}

// stripe a parameter set from the model,\
   parameter arrays are shared with the model.
param_set *
get_assoc_param_set(model * m_t)
{
  // make param_set object
  param_set * ps_t = TALLOC(param_set, 1);
  ps_t -> N_par = TALLOC(int, m_t -> N_cps),
  ps_t -> par = TALLOC(double *, m_t -> N_cps),
  ps_t -> par_lim = TALLOC(double *, m_t -> N_cps),
  ps_t -> is_const = TALLOC(int *, m_t -> N_cps);

  // write number of parameters and set pointers
  int I_cp, N_cps = m_t -> N_cps;
  FOREACH(I_cp, N_cps)
    ps_t -> N_par[I_cp] = m_t -> cps[I_cp] -> N_par,
    ps_t -> par[I_cp] = m_t -> cps[I_cp] -> par,
    ps_t -> par_lim[I_cp] = m_t -> cps[I_cp] -> par_lim,
    ps_t -> is_const[I_cp] = m_t -> cps[I_cp] -> is_const;

  // set other properties
  ps_t -> N_cps = N_cps;
  ps_t -> is_shared = 1;

  // return
  return ps_t;
}

// make a deep copy of the parameter set \
   parameter arrays are indepent with the model
param_set *
make_param_set_for(model * m_t)
{
  // allocate parameter set
  param_set * ps_t = TALLOC(param_set, 1);
  ps_t -> N_par = TALLOC(int, m_t -> N_cps),
  ps_t -> par = TALLOC(double *, m_t -> N_cps),
  ps_t -> par_lim = TALLOC(double *, m_t -> N_cps),
  ps_t -> is_const = TALLOC(int *, m_t -> N_cps);

  // FIXME
}

/*
  TODO: in this way?
    param_set * make_param_set_for (model *),
    param_set * get_assoc_param_set (model *),
    int read_param_from (param_set *, param_set *),
    int apply_param_set (param_set *, param_set *)
*/

int
free_param_set(param_set * ps_t)
{
  free(ps_t -> N_par), free(ps_t -> par), free(ps_t -> par_lim),
  free(ps_t -> is_const); free(ps_t);
  return 0;
}
