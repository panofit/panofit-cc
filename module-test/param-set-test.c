
#include "../src/cps/sersic.h"
#include "../src/cps/expdisk.h"
#include "../src/component.h"
#include "../src/model.h"
#include "../src/utils.h"
#include "../src/spectrum.h"
#include "../src/recipe.h"

#include "stdio.h"
#include "stdlib.h"

int
main(int argc, char ** argv)
{
  // allocate, copy and delete parameter sets // passed
  if(0)
    {
      double disk_par[] = {1., 0., 0., 0., 1., 0.5, 2.,
                           8., 0.2, 0., 0.2, 0., 0., 0., 0.};
      double disk_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 0.5, 0.5, 2.,
                           2., 8., 8., 0.2, 0.2, 0., 0., 0.2, 0.2, 0., 0.,
                           0., 0., 0., 0., 0., 0.};
      int disk_var[]    = {1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1};

      component * cp_disk =
        cps_expdisk_2d(disk_par, disk_lim, disk_var, "disk");

      double bulge_par[] = {1., 0., 0., 0., 1., 2., 0.5, 2.,
                            8., 0.2, 0., 0.2, 0., 0., 0., 0.};
      double bulge_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 2., 2.,
                            0.5, 0.5, 2., 2., 8., 8., 0.2, 0.2, 0., 0.,
                            0.2, 0.2, 0., 0., 0., 0., 0., 0., 0., 0.};
      int bulge_var[]    = {1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1};

      component * cp_bulge =
        cps_sersic_2d(bulge_par, bulge_lim, bulge_var, "bulge");

      // make a model
      model * m_t = make_model("Simple model", 2, cp_bulge, cp_disk);

      // get shared param set
      param_set * ps_t = get_assoc_param_set(m_t);

      // try print parameters
      print_param_set(ps_t, 0);

      // copy parameters
      param_set * ps_w = make_param_set_for(m_t);

      // try print parameters again
      print_param_set(ps_w, 1);

      // free this deep-copy
      free_param_set(ps_w);

      free_param_set(ps_t);
      free_model(m_t, 1);
    }

  // test monte-carlo sampling
  if(1)
    {
      double disk_par[] = {1., 0., 0., 0., 1., 0.5, 2.,
                           8., 0.2, 0., 0.2, 0., 0., 0., 0.};
      double disk_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 0.5, 0.5, 2.,
                           2., 8., 8., 0.2, 0.2, 0., 0., 0.2, 0.2, 0., 0.,
                           0., 0., 0., 0., 0., 0.};
      int disk_var[]    = {1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1};

      component * cp_disk =
        cps_expdisk_2d(disk_par, disk_lim, disk_var, "disk");

      double bulge_par[] = {1., 0., 0., 0., 1., 2., 0.5, 2.,
                            8., 0.2, 0., 0.2, 0., 0., 0., 0.};
      double bulge_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 2., 2.,
                            0.5, 0.5, 2., 2., 8., 8., 0.2, 0.2, 0., 0.,
                            0.2, 0.2, 0., 0., 0., 0., 0., 0., 0., 0.};
      int bulge_var[]    = {1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1};

      component * cp_bulge =
        cps_sersic_2d(bulge_par, bulge_lim, bulge_var, "bulge");

      // make a model
      model * m_t = make_model("Simple model", 2, cp_bulge, cp_disk);

      // get shared param set
      param_set * ps_t = get_assoc_param_set(m_t);
      param_set * ps_w = make_param_set_for(m_t);

      // generate new solutions

      free_param_set(ps_w), free_param_set(ps_t);
      free_model(m_t, 1);
    }
}
