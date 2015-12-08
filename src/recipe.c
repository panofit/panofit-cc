
#include "model.h"
#include "recipe.h"

recipe *
make_empty_recipe(spec_lib * lib_t)
{
  // allocate recipe
  recipe * rcp_t = TALLOC(recipe, 1);
  rcp_t -> rcp = TALLOC(double, N_age_pts * N_Z_pts);

  // set to zero
  int I_pt, N_pts = N_age_pts * N_Z_pts;
  for(I_pt = 0; I_pt < N_pts; ++ I_pt)
    *(rcp_t -> rcp + I_pt) = 0.;

  // set properties
  rcp_t -> x = 0., rcp -> y = 0., rcp -> z = 0.;
  rcp_t -> N_age = lib_t -> N_age,
  rcp_t -> N_Z = lib_t -> N_Z;

  // return object
  return rcp_t;
}

int
free_recipe(recipe * rcp_t)
{
  free(rcp_t -> rcp);
  free(rcp_t);

  return 0;
}

recipe *
sample_recipe(model * m_t, double x, double y)
{
  // make recipe object
  recipe * rcp_t = TALLOC(recipe, 1);

  // calculate using noalloc function
  sample_recipe_noalloc(m_t, x, y, rcp_t);

  // return object
  return rcp_t;
}

int
sample_recipe_noalloc(model * m_t, double x, double y, recipe * rcp_t)
{
  // loop over components

  // for each component, estimate recipe,

  // find sum,
}
