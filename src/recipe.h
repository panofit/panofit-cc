
#ifndef _RECIPE_H
#define _RECIPE_H

typedef struct _spec_lib spec_lib;
typedef struct _model model;

typedef struct _recipe
{
  // coordinate (maybe unnecesary)
  double x, y, z;

  // grid points in age and Z
  int N_age, N_Z;
  double * age_ax, * Z_ax;

  // data
  double * rcp;

} recipe;

recipe * make_empty_recipe(spec_lib *);
int free_recipe();

recipe * sample_recipe(model *, double, double);
int sample_recipe_noalloc(model *, double, double, recipe *);

#endif // _RECIPE_H
