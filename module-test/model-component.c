
#include "../src/model.h"
#include "../src/component.h"

#include "stdio.h"

int
main(int argc, char ** argv)
{
  // make and delete single component // passed
  if(0)
    {
      component * cp_t = make_base_component(10);
      printf("N_par: %u, N_fp: %u\n", cp_t -> N_par, cp_t -> N_fp);
      free_component(cp_t);
    }

  // make and delete model // passed
  if(1)
    {
      component * cp_1 = make_base_component(8);
      component * cp_2 = make_base_component(12);

      model * m_t = make_model("Simple model", 2, cp_1, cp_2);

      free_model(m_t, 0);
      free_component(cp_1);
      free_component(cp_2);
    }
}
