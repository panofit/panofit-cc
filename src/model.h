
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

  // number of parameters, for each component
  int * N_par;

  // pointers to parameters
  double ** par;
  double ** par_lim;
  int ** is_const;

} param_set;

model * make_model(const char * name, int N_cps, ...);
int free_model(model *, int);
param_set * get_param_set(model * m_t);

#endif // _MODEL_H
