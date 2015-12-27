
#ifndef _DATASET_H
#define _DATASET_H

typedef struct _dataset
{
  // basic information
  char name[256];
  int N_img, N_spec;

  // images
  image ** img;

  // spectra
  spectrum ** spec;

  // other stuff
} dataset;

#endif
