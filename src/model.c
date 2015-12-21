
#include "utils.h"
#include "model.h"

#include "stdarg.h"
#include "string.h"

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
        (*((m_t -> cps + I_cp) -> kill))(m_t -> cps + I_cp);
    }

  // free componen and return
  free(m_t -> cps), free(m_t);
  
  return 0;
}

param_set *
get_param_set(model * m_t)
{
  // make param_set object
  param_set * ps_t = TALLOC(param_set, 1);
  ps_t -> N_par = TALLOC(int, m_t -> N_cps);
  ps_t -> par = TALLOC(double *, m_t -> N_cps);
  ps_t -> par_lim = TALLOC(double *, m_t -> N_cps);

  // write number of parameters and set pointers
  int I_cp, N_cps = m_t -> N_cps;
  for(I_cp = 0; I_cp < N_cps; ++ I_cp)
    *(ps_t -> N_par + I_cp) = (*(m_t -> cps + I_cp)) -> N_par,
    *(ps_t -> par + I_cp) = (*(m_t -> cps + I_cp)) -> par,
    *(ps_t -> par_lim + I_cp) = (*(m_t -> cps + I_cp)) -> par_lim;

  // set other properties
  ps_t -> N_cps = N_cps;

  // return
  return ps_t;
}
