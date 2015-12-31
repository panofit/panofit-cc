
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
  ps_t -> is_const = TALLOC(int *, m_t -> N_cps),
  ps_t -> par_name = TALLOC(char **, m_t -> N_cps),
  ps_t -> name = TALLOC(char *, m_t -> N_cps),
  ps_t -> check = malloc(sizeof(int (**)(double const *)) * (m_t -> N_cps));

  // write number of parameters and set pointers
  int I_cp, N_cps = m_t -> N_cps;
  FOREACH(I_cp, N_cps)
    ps_t -> N_par[I_cp] = m_t -> cps[I_cp] -> N_par,
    ps_t -> par[I_cp] = m_t -> cps[I_cp] -> par,
    ps_t -> par_lim[I_cp] = m_t -> cps[I_cp] -> par_lim,
    ps_t -> is_const[I_cp] = m_t -> cps[I_cp] -> is_const,
    ps_t -> par_name[I_cp] = (m_t -> cps[I_cp]) -> par_name,
    ps_t -> name[I_cp] = (m_t -> cps[I_cp]) -> name,
    ps_t -> check[I_cp] = (m_t -> cps[I_cp]) -> check;

  // set other properties
  ps_t -> N_cps = N_cps,
  ps_t -> is_shared = 1,
  ps_t -> shared_with = m_t;

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
  ps_t -> is_const = TALLOC(int *, m_t -> N_cps),
  ps_t -> par_name = TALLOC(char **, m_t -> N_cps),
  ps_t -> name = TALLOC(char *, m_t -> N_cps),
  ps_t -> check = malloc(sizeof(int (**)(double const *)) * (m_t -> N_cps));

  // allocate memory and copy parameters
  int I_cp, N_cps = m_t -> N_cps;
  FOREACH(I_cp, N_cps)
    ps_t -> N_par[I_cp] = m_t -> cps[I_cp] -> N_par,
    ps_t -> par[I_cp] = TALLOC(double, m_t -> cps[I_cp] -> N_par),
    ps_t -> par_lim[I_cp] = TALLOC(double, 2 * (m_t -> cps[I_cp] -> N_par)),
    ps_t -> is_const[I_cp] = TALLOC(int, m_t -> cps[I_cp] -> N_par);

  FOREACH(I_cp, N_cps)
    memcpy(ps_t -> par[I_cp], m_t -> cps[I_cp] -> par,
        (m_t -> cps[I_cp] -> N_par) * sizeof(double)),
    memcpy(ps_t -> par_lim[I_cp], m_t -> cps[I_cp] -> par_lim,
        (m_t -> cps[I_cp] -> N_par) * sizeof(double) * 2),
    memcpy(ps_t -> is_const[I_cp], m_t -> cps[I_cp] -> is_const,
        (m_t -> cps[I_cp] -> N_par) * sizeof(int));

  // copy name of param
  FOREACH(I_cp, N_cps)
    ps_t -> name[I_cp] = (m_t -> cps[I_cp]) -> name,
    ps_t -> par_name[I_cp] = (m_t -> cps[I_cp]) -> par_name,
    ps_t -> check[I_cp] = (m_t -> cps[I_cp]) -> check;

  // set other properties
  ps_t -> N_cps = N_cps,
  ps_t -> is_shared = 0,
  ps_t -> shared_with = m_t;

  // return
  return ps_t;
}

int
copy_param_set(param_set * ps_d, param_set * ps_i)
{
  // check if destination is independent
  if(ps_d -> is_shared) return 1;
  if(ps_d -> N_cps != ps_i -> N_cps) return -1;
  // if(ps_d -> shared_with != ps_i -> shared_with) return -4;

  // make deep copy
  int I_cp, N_cps = ps_i -> N_cps;
  FOREACH(I_cp, N_cps)
    {
      // if number of parameter match
      if(ps_i -> N_par[I_cp] != ps_d -> N_par[I_cp]) return -2;

      memcpy(ps_d -> par[I_cp], ps_d -> par[I_cp],
          (ps_i -> N_par[I_cp]) * sizeof(double)),
      memcpy(ps_d -> par_lim[I_cp], ps_d -> par_lim[I_cp],
          (ps_i -> N_par[I_cp]) * sizeof(double) * 2),
      memcpy(ps_d -> is_const[I_cp], ps_d -> is_const[I_cp],
          (ps_i -> N_par[I_cp]) * sizeof(int));
    }

  // set other paramters
  FOREACH(I_cp, N_cps)
    ps_d -> name[I_cp] = ps_i -> name[I_cp],
    ps_d -> par_name[I_cp] = ps_i -> par_name[I_cp],
    ps_d -> check[I_cp] = ps_i -> check[I_cp];

  return 0;
}


int
free_param_set(param_set * ps_t)
{
  // deep-copy or associated
  if(!(ps_t -> is_shared))
    {
      int I_cp, N_cps = ps_t -> N_cps;
      FOREACH(I_cp, N_cps)
        free(ps_t -> par[I_cp]), free(ps_t -> par_lim[I_cp]),
        free(ps_t -> is_const[I_cp]);
    }

  free(ps_t -> N_par), free(ps_t -> par), free(ps_t -> par_lim),
  free(ps_t -> is_const), free(ps_t -> par_name), free(ps_t -> check);
  free(ps_t);
  return 0;
}
