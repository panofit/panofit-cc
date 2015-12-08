
#ifndef _SPECTRUM_H
#define _SPECTRUM_H

// spectral library object
typedef struct _spec_lib
{
  // name
  char name[128];

  // N of log_age bins, N of metallicity bins, N of spectral pixels
  int N_age, N_Z, N_spx;

  // grid axes of log_age and Z
  double * age_ax, * Z_ax;

  // wavelength / wavenumber
  double * wl;

  // data hold there
  double * data;

} spec_lib;

// spectrum object
typedef struct _spectrum spectrum;

struct _spectrum
{
  // number of spectral pixels
  int N_spx;

  // array to hold flux and err
  double * flux, * err;

  // wavelength / wavenumber
  double * wl;

  // sharing wl array with other spectra?
  int is_shared_wl; // avoid multi-free
};

spectrum * make_empty_spectrum(spec_lib *);
int free_spectrum(spectrum *);

spec_lib * load_spec_lib_raw(const char *);
int free_spec_lib(spec_lib *);

#endif
