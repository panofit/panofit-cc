#!/usr/bin/python

'''
   Generates pre-computed spectra.
'''

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import interp1d

import cfg

def _make_spec_lib():

  spec_dir = '/home/qinyj/workspace/panofit/slug2/my_output/integrated_results/'
  spec_lib = np.zeros((cfg.age_div, cfg.M_div, 448, 2))

  spec_files = [('SSP_10GYR_005Z_spec-cut.npz', np.log10(0.05)),
                ('SSP_10GYR_02Z_spec-cut.npz',  np.log10(0.2)),
                ('SSP_10GYR_04Z_spec-cut.npz',  np.log10(0.4)),
                ('SSP_10GYR_1Z_spec-cut.npz',   np.log10(1.)),
                ('SSP_10GYR_2Z_spec-cut.npz',   np.log10(2.))]

  spec_lib_raw = np.zeros((448, len(spec_files), cfg.age_div, 2))
  M_raw = np.zeros(len(spec_files))

  for I_ft, frec_t in enumerate(spec_files):

    print 'Processing file', I_ft, ':', frec_t[0]

    ssp_t = np.load(spec_dir + frec_t[0])
    wl_t = ssp_t['wl']; spec_t = ssp_t['spec']

    M_raw[I_ft] = frec_t[1]
    spec_lib_raw[:, I_ft, :, 0] = np.nanmean(spec_t, axis = -1)
    spec_lib_raw[:, I_ft, :, 1] = np.nanstd(spec_t, axis = -1)

  print 'Getting array:', spec_lib_raw.shape

  # create interpolator.
  M_intp = interp1d(M_raw, spec_lib_raw, kind = 'quadratic', axis = 1)

  for Im, Mt in enumerate(cfg.M_ax_log):

    # interpolate a metallicity.
    spec_t = M_intp(Mt) # < in order (448, age, 2)
    spec_t = np.swapaxes(spec_t, 0, 1)

    spec_lib[:, Im, :, :] = spec_t

  print 'Getting spec_lib:', spec_lib.shape

  with open("mock-speclib.dat", "wb") as fid:

    name_str = np.empty(shape = 128, dtype = 'c')
    name_str[:] = "Mock spectral library\0\0\0"
    name_str.tofile(fid)

    header = np.array([cfg.M_div, cfg.age_div, 448]).astype('i4')
    header.tofile(fid)

    cfg.age_ax_log.tofile(fid)
    cfg.M_ax_log.tofile(fid)
    wl_t.astype('f8').tofile(fid)

    (spec_lib[:, :, :, 0]).tofile(fid)

if __name__ == '__main__': _make_spec_lib()
