
#ifndef _BAYESIAN_H
#define _BAYESIAN_H

typedef struct _spectrum spectrum;
typedef struct _image image;

double spec_ln_L(spectrum * mock_sp, spectrum * obs_sp);
double img_ln_L(image * mock_img, image * obs_img);

#endif
