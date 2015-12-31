
#include "../src/cps/sersic.h"
#include "../src/cps/expdisk.h"
#include "../src/component.h"
#include "../src/model.h"
#include "../src/mcfit.h"
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

  // test rand_gauss
  if(0)
    {
      int I, N = 20;
      FOREACH(I, N) printf("%f\n", rand_gauss());
    }

  // test monte-carlo sampling
  if(1)
    {
      int I_par;

      double disk_par[] = {1., 0., 0., 0., 1., 0.5, 2.,
                           8., 0.2, 0., 0.2, 0., 0., 0., 0.};
      double disk_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 0.5, 0.5, 2.,
                           2., 8., 8., 0.2, 0.2, 0., 0., 0.2, 0.2, 0., 0.,
                           0., 0., 0., 0., 0., 0.};
      int disk_var[]    = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

      FOREACH(I_par, NPAR_EXPDISK2D)
        disk_lim[I_par * 2] = disk_par[I_par] - 0.5,
        disk_lim[I_par * 2 + 1] = disk_par[I_par] + 0.5;

      component * cp_disk =
        cps_expdisk_2d(disk_par, disk_lim, disk_var, "disk");

      double bulge_par[] = {1., 0., 0., 0., 1., 2., 0.5, 2.,
                            8., 0.2, 0., 0.2, 0., 0., 0., 0.};
      double bulge_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 2., 2.,
                            0.5, 0.5, 2., 2., 8., 8., 0.2, 0.2, 0., 0.,
                            0.2, 0.2, 0., 0., 0., 0., 0., 0., 0., 0.};
      int bulge_var[]    = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

      FOREACH(I_par, NPAR_SERSIC2D)
        bulge_lim[I_par * 2] = bulge_par[I_par] - 0.5,
        bulge_lim[I_par * 2 + 1] = bulge_par[I_par] + 0.5;

      component * cp_bulge =
        cps_sersic_2d(bulge_par, bulge_lim, bulge_var, "bulge");

      // make a model
      model * m_t = make_model("Simple model", 2, cp_bulge, cp_disk);

      // get shared param set
      param_set * ps_t = get_assoc_param_set(m_t);

      // make deep copy
      param_set * ps_w = make_param_set_for(m_t);

      // generate new solutions
      sample_param_set(ps_w, ps_t, 0.01);

      // print new parameter set
      print_param_set(ps_t, 0);
      print_param_set(ps_w, 0);

      free_param_set(ps_w), free_param_set(ps_t);
      free_model(m_t, 1);
    }
}
