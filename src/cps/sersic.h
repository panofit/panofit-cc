
#ifndef SERSIC_H
#define SERSIC_H

#define NPAR_SERSIC2D 16

typedef struct _component component;
typedef struct _recipe recipe;

component * cps_sersic_2d(double *, double *, int *, const char *);
double _cps_sersic_2d_sigma(double, double, double const *);
int _cps_sersic_2d_recipe(double, double, recipe *, double const *);

#endif
