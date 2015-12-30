
#ifndef _COMPONENT_H
#define _COMPONENT_H

typedef struct _component component;
typedef struct _recipe recipe;

struct _component
{
  // info
  char name[128];
  int type;

  // parameters
  int N_par, N_fp; // N of params, N of free params
  double * par, * par_lim;
  int * is_const;
  char ** par_name;

  // workspace (if needed)
  void * ws;

  // calculate surface density and stellar recipe
  double (* sigma)(double, double, double const *);
  int (* recipe)(double, double, recipe *, double const *);

  // boundary check
  int (* check)(double const *);

  // initialize and finalize
  int (* init)(component *);
  int (* kill)(component *);
};

component * make_base_component(int);
int free_component(component *);

int cps_uni_init(component *);
int cps_uni_kill(component *);
double cps_uni_sigma(double, double, double const *);
int cps_uni_recipe(double, double, recipe *, double const *);
int cps_bounds_ok(double const *);

#endif // _COMPONENT_H
