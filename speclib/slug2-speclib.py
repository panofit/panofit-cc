#!/usr/bin/python

'''
  Make spectral libs for CALIFA V500/V1200 datacubes.
  YJ Qin, Dec 2015 @ Shanghai
'''

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import interp1d

# parameters of the original spectral lib
N_spx_raw = 448
N_age_div, N_M_div = 41, 32
age_min, age_max = 1.e6, 1.e10

# age axis in the original library
age_alpha = np.log(age_max / age_min) / (N_age_div - 1.)
age_ax = age_min * np.exp(age_alpha * np.arange(N_age_div))
age_ax = np.log10(age_ax)

# metallicity axis
M_min, M_max = 0.05, 2.0
M_ax = np.linspace(np.log10(M_min), np.log10(M_max), N_M_div)

# directory of the original spec library
spec_dir = '/home/qinyj/workspace/panofit/slug2/my_output/integrated_results/'

# make spectral library for a specific dataset (V500/V1200)
def _make_spec_lib_for(sample_cube_fname, output_speclib_fname):

  spec_lib = np.zeros((N_age_div, N_M_div, N_spx_raw))

  spec_files = [('SSP_10GYR_005Z_spec-cut.npz', np.log10(0.05)),
                ('SSP_10GYR_02Z_spec-cut.npz',  np.log10(0.2)),
                ('SSP_10GYR_04Z_spec-cut.npz',  np.log10(0.4)),
                ('SSP_10GYR_1Z_spec-cut.npz',   np.log10(1.)),
                ('SSP_10GYR_2Z_spec-cut.npz',   np.log10(2.))]

  spec_lib_raw = np.zeros((N_spx_raw, len(spec_files), N_age_div))
  M_raw = np.zeros(len(spec_files))

  for I_ft, frec_t in enumerate(spec_files):

    # print 'Processing file', I_ft, ':', frec_t[0]

    ssp_t = np.load(spec_dir + frec_t[0])
    wl_t = ssp_t['wl']; spec_t = ssp_t['spec']

    M_raw[I_ft] = frec_t[1]
    spec_lib_raw[:, I_ft, :] = np.nanmean(spec_t, axis = -1)

  # print 'Getting array:', spec_lib_raw.shape

  # create interpolator for metallicity
  M_intp = interp1d(M_raw, spec_lib_raw, kind = 'quadratic', axis = 1)

  # interpolate metallicity
  for Im, Mt in enumerate(M_ax):

    # interpolate a metallicity.
    spec_t = M_intp(Mt) # < in order (N_spx_raw, age)
    spec_t = np.swapaxes(spec_t, 0, 1)
    spec_lib[:, Im, :] = spec_t

  # open the example datacube to make new wavelength axis
  hlst = fits.open(sample_cube_fname)
  N_spx_cube = hlst[0].header['NAXIS3']
  wl_ax = hlst[0].header['CRVAL3'] + hlst[0].header['CDELT3'] * np.arange(N_spx_cube)
  hlst.close()

  # the third cube, for the resampled library
  flux_intp = interp1d(wl_t, spec_lib, kind = 'linear', axis = 2, assume_sorted = True)
  spec_lib_intp = flux_intp(wl_ax)

  # DEBUG
  np.save(output_speclib_fname + '.dbg', spec_lib_intp)

  # write the output file
  with open(output_speclib_fname, "wb") as fid:

    # write name
    name_str = np.empty(shape = 128, dtype = 'c')
    name_str[:] = "Mock spectral library\0\0\0"
    name_str.tofile(fid)

    # write header
    header = np.array([N_M_div, N_age_div, N_spx_cube]).astype('i4')
    header.tofile(fid)

    # write age/metallicity/wavelength axes
    age_ax.astype('f8').tofile(fid)
    M_ax.astype('f8').tofile(fid)
    wl_ax.astype('f8').tofile(fid)

    # write the spec lib
    spec_lib_intp.tofile(fid)

  return 0

if __name__ == "__main__":

  _make_spec_lib_for("./califa_sample/NGC0001.V1200.rscube.fits", "./NGC0001/speclib-V1200.dat")
  _make_spec_lib_for("./califa_sample/NGC0001.V500.rscube.fits", "./NGC0001/speclib-V500.dat")

# lo6zYxwD1vuKfBDezcfj
