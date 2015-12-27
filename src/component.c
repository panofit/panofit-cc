
#include "utils.h"
#include "component.h"

#include "stdlib.h"
#include "float.h"

// make an empty component
component *
make_base_component(int N_par)
{
  int I_par;

  // allocate space
  component * cp_t = TALLOC(component, 1);

  // parameters
  cp_t -> N_par = N_par;
  cp_t -> N_fp = N_par;
  cp_t -> par = TALLOC(double, N_par);
  cp_t -> par_lim = TALLOC(double, N_par * 2);
  cp_t -> is_const = TALLOC(int, N_par);

  // all params are free by default
  for(I_par = 0; I_par < N_par; ++ I_par)
    (cp_t -> is_const)[I_par] = 0;

  // params are zero by default.
  for(I_par = 0; I_par < N_par; ++ I_par)
    (cp_t -> par)[I_par] = 0.;

  // ranges are unlimited by default
  for(I_par = 0; I_par < N_par; ++ I_par)
    (cp_t -> par_lim)[I_par * 2] = -DBL_MAX,
    (cp_t -> par_lim)[I_par * 2 + 1] = DBL_MAX;

  // bounds check
  cp_t -> check = & cps_bounds_ok;

  // init and kill
  cp_t -> init = & cps_uni_init;
  cp_t -> kill = & cps_uni_kill;

  // sigma and recipe
  cp_t -> sigma = & cps_uni_sigma;
  cp_t -> recipe = & cps_uni_recipe;
}

int free_component(component * cp_t)
{
  (cp_t -> kill)(cp_t); free(cp_t); cp_t = NULL;
  return 0;
}

// universal init and kill func.
int cps_uni_init(component * cp_t) {return 0;}

int cps_uni_kill(component * cp_t)
{
  free(cp_t -> par),
  free(cp_t -> par_lim),
  free(cp_t -> is_const);

  return 0;
}

// dummy function;
double cps_uni_sigma(double x, double y, double const * par) {return -1.;}

int cps_uni_recipe(double x, double y,
    recipe * rcp, double const * par) {return 0;}

// dummy function for bounds check. always pass.
int cps_bounds_ok(double const * par) {return 0; /* banana! */}

// compiled, Dec 23, 2015
