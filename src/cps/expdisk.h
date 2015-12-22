
#ifndef EXPDISK_H

#define NPAR_EXPDISK2D 15

component * cps_expdisk_2d(double *, double *, int *, const char *);
double _cps_expdisk_2d_sigma(double, double, double const *);
double _cps_expdisk_2d_recipe(double, double, const double *, double const *);

#define EXPDISK_H
#endif
