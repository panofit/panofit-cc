
#ifndef EXPDISK_H

#define NPAR_EXPDISK2D 15

typedef struct _component component;
typedef struct _recipe recipe;

component * cps_expdisk_2d(double *, double *, int *, const char *);

double _cps_expdisk_2d_sigma(double, double, double const *);
int _cps_expdisk_2d_sigma_arr(int, double *, double *,
    double const *, double *);

int _cps_expdisk_2d_recipe(double, double, recipe *, double const *);

#define EXPDISK_H
#endif
