
#ifndef _MODEL_H
#define _MODEL_H

typedef struct _component component;

typedef struct _model
{
  // name of the model
  char name[128];

  // number of components
  int N_cps;

  // list of components
  component ** cps;

} model;

typedef struct _param_set
{
  // number of components
  int N_cps;

  // name of components
  char ** name;

  // number of parameters, for each component
  int * N_par;

  // pointers to parameters
  double ** par;
  double ** par_lim;
  int    ** is_const;
  char  *** par_name;

  // parameter checkers
  int (** check)(double const *);

  // are parameter arrays shared with the model?
  int is_shared;

  // pointer to the associated model
  model * shared_with;

} param_set;

model * make_model(const char * name, int N_cps, ...);
int free_model(model *, int);
param_set * get_assoc_param_set(model *);
param_set * make_param_set_for(model *);
int copy_param_set(param_set *, param_set *);


#endif // _MODEL_H
