
#ifndef _MCFIT_H

typedef struct _spectrum spectrum;
typedef struct _image image;
typedef struct _param_set param_set;

// likelihood functions
double spec_logL(spectrum * mock_sp, spectrum * obs_sp);
double img_logL(image * mock_img, image * obs_img);

// monte-carlo sampling
int sample_param_set(param_set * ps_d, param_set * ps_s, double T);

#define _MCFIT_H
#endif
