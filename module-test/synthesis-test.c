
#include "../src/utils.h"
#include "../src/component.h"
#include "../src/model.h"
#include "../src/spectrum.h"
#include "../src/recipe.h"
#include "../src/synthesis.h"

#include "../src/cps/expdisk.h"

int
main(int argc, char ** argv)
{
  // make a synthetic spectrum // passed.
  if(1)
  {
    // read spec lib
    spec_lib * lib_t = load_spec_lib_raw("mock-speclib.dat");

    // make proper model
    double cp_par[] = {1., 0., 0., 0., 1., 0.5, 2.,
                       10., 1., 0., 0.2, 0., 0., 0., 0.};
    double cp_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 0.5, 0.5, 2.,
                       2., 8., 8., 0.2, 0.2, 0., 0., 0.2, 0.2, 0., 0.,
                       0., 0., 0., 0., 0., 0.};
    int    cp_var[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    component * cp_t = cps_expdisk_2d(cp_par, cp_lim, cp_var, "disk");
    model * m_t = make_model("model", 1, cp_t);

    // generate a recipe
    recipe * rc_t = make_empty_recipe(lib_t);
    sample_recipe_noalloc(m_t, 1., 1., rc_t);

    // do synthesis
    spectrum * sp_t = synthesis(rc_t, lib_t);

    // print them.
    print_arr(sp_t -> flux, sp_t -> N_spx);
  }
  return 0;
}
