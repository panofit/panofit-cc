
#ifndef _SYNTHESIS_H
#define _SYNTHESIS_H

typedef struct _spectrum spectrum;
typedef struct _spec_lib spec_lib;
typedef struct _model    model;
typedef struct _recipe   recipe;

spectrum * synthesis(recipe * rcp_t, spec_lib * lib_t, double x, double y);
int synthesis_noalloc(recipe * rcp_t, spec_lib * lib_t,
    double x, double y, spectrum * sp_t);

#endif // _SYNTHESIS_H
