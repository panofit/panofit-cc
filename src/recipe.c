
#include "recipe.h"

extern int N_age_pts, N_Z_pts;

recipe
make_empty_recipe()
{
  // allocate recipe
  recipe * rcp_t = TALLOC(recipe, 1);
  rcp_t -> rcp = TALLOC(double, N_age_pts * N_Z_pts);

  // set to zero
  int I_pt, N_pts = N_age_pts * N_Z_pts;
  for(I_pt = 0; I_pt < N_pts; ++ I_pt)
    *(rcp_t -> rcp + I_pt) = 0.;

  rcp_t -> x = 0., rcp -> y = 0., rcp -> z = 0.;
}

int
free_recipe(recipe * rcp_t)
{
  free(rcp_t -> rcp);
  free(rcp_t);

  return 0;
}
