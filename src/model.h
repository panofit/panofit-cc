
#ifndef _MODEL_H
#define _MODEL_H

typedef struct _component component;

typedef struct _model
{
  char name[128];

  int N_cps;
  component ** cps;
} model;

model * make_model(int N_cps, ...);


#endif // _MODEL_H
