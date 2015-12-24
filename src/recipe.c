
#include "spectrum.h"
#include "model.h"
#include "recipe.h"
#include "utils.h"
#include "component.h"

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

recipe *
make_empty_recipe_as(recipe * rc_t)
{
  int I_pt, N_pts = (rc_t -> N_age) * (rc_t -> N_Z);

  // allocate recipe
  recipe * rcp_t = TALLOC(recipe, 1);
  rcp_t -> rcp   = TALLOC(double, N_pts);

  // set to zero
  FOREACH(I_pt, N_pts) rcp_t -> rcp[I_pt] = 0.;

  // set properties
  rcp_t -> x = 0., rcp_t -> y = 0., rcp_t -> z = 0.;
  rcp_t -> N_age = rc_t -> N_age,
  rcp_t -> N_Z   = rc_t -> N_Z;

  rcp_t -> age_ax = rc_t -> age_ax,
  rcp_t -> Z_ax = rc_t -> Z_ax;

  return rcp_t;
}

// co-add recipes. rc_s = rc_s + rc_i
int
recipe_coadd(recipe * rc_s, recipe * rc_i)
{
  // check!
  int I_pt, N_Z = rc_s -> N_Z, N_age = rc_s -> N_age;
  if((N_Z != rc_i -> N_Z) || (N_age != rc_i -> N_age)) return -1;

  // co-add
  int N_size = N_Z * N_age;
  FOREACH(I_pt, N_size) rc_s -> rcp[I_pt] += rc_i -> rcp[I_pt];

  // return
  return 0;
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
  int I_cp, N_cps = m_t -> N_cps;

  // empty recipe (for co-adding)
  recipe * rcp_i = make_empty_recipe_as(rcp_t);

  // for each component, estimate recipe and co-add them
  for(I_cp = 0; I_cp < N_cps; ++ I_cp)
    ((m_t -> cps[I_cp]) -> recipe)(x, y, rcp_i, (m_t -> cps[I_cp]) -> par),
    recipe_coadd(rcp_t, rcp_i);

  return 0;
}
