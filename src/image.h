
#ifndef _IMAGE_H
#define _IMAGE_H

typedef struct _image
{
  // basic info of this image
  char band[16]; // name of the band (not used as identifier)
  int idx_band;  // index of the band in the photometric library

  // axis
  int X_pts, Y_pts; // grid points in X and Y
  double * X_ax, * Y_ax; // axes, includes register info.

  // image data
  double * data, * err;
  int * mask;

  // other info?
} image;

#endif
