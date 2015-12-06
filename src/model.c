
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
free_model(model *, int free_cps)
{
  return 0;
}
