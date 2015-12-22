
#include "spectrum.h"
#include "model.h"
#include "recipe.h"
#include "utils.h"

#include "stdlib.h"

recipe *
make_empty_recipe(spec_lib * lib_t)
{
  int I_pt, N_pts = (lib_t -> N_age) * (lib_t -> N_Z);

  // allocate recipe
  recipe * rcp_t = TALLOC(recipe, 1);
  rcp_t -> rcp   = TALLOC(double, N_pts);

  // set to zero
  FOREACH(I_pt, N_pts) rcp_t -> rcp[I_pt] = 0.;

  // set properties
  rcp_t -> x = 0., rcp_t -> y = 0., rcp_t -> z = 0.;
  rcp_t -> N_age = lib_t -> N_age,
  rcp_t -> N_Z = lib_t -> N_Z;

  rcp_t -> age_ax = lib_t -> age_ax,
  rcp_t -> Z_ax = lib_t -> Z_ax;

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
