
#ifndef SERSIC_H
#define SERSIC_H

#define NPAR_SERSIC2D 16

component * cps_sersic_2d(double *, double *, int *, const char *);
double _cps_sersic_2d_sigma(double, double, double const *);
double _cps_sersic_2d_recipe(double, double, const double *, double const *);

#endif
