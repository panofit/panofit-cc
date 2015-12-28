
#ifndef _DATASET_H
#define _DATASET_H

typedef struct _spec_lib spec_lib;
typedef struct _spectrum spectrum;

// a batch of spectra.
typedef struct _spec_stack
{
  // a fancy name
  char name[256];

  // N of spectra, N of spectral pixels
  int N_spec, N_spx;

  // wavelength array
  double * wl;

  // array of spectra
  spectrum ** spec;

  // corresponding spec_lib
  spec_lib * splib;

} spec_stack;

typedef struct _dataset
{
  // basic information
  char name[256];

  // number of images and spectral stacks
  int N_img, N_spec_st;

  // images
  image ** img;

  // spectra
  spec_stack ** spec_st;

  // other stuff
} dataset;

#endif
