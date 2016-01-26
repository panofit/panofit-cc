#!/usr/bin/python

'''
  Functions for spectra resampling
'''

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d

# move the spectra to rest frame
def restframe(wl, flux, err, mask, wl_new, z, mask_nan = False):

  # rest-frame wavelength
  wl_restf = wl / (1. + z)

  # flux interpolator
  flux_atd = InterpolatedUnivariateSpline(wl, flux).antiderivative() # FIXME bounds?
  flux_spl = InterpolatedUnivariateSpline(wl_restf, flux_atd(wl)).derivative()

  # error interpolator (not the right way...)
  err_intp = InterpolatedUnivariateSpline(wl_restf[mask == 1], err[mask == 1], k = 1, ext = 3)
  # mask_spl = InterpolatedUnivariateSpline(wl_restf, mask.astype('f8'), k = 1)

  # interpolate
  flux_new, err_new = flux_spl(wl_new), err_intp(wl_new)
  mask_new = np.ones(flux_new.size, dtype = 'i4')

  # mask bad points
  badpx_restwl = wl_restf[mask == 0]
  wlt_a, wlt_b = 2. * wl_new[0] - wl_new[1], 2. * wl_new[-1] - wl_new[-2]
  badpx_restwl = badpx_restwl[(badpx_restwl > wlt_a) & (badpx_restwl < wlt_b)]
  id_badpts = np.unique(np.searchsorted(0.5 * (wl_new[1:] + wl_new[:-1]), badpx_restwl))

  # mask bad pixels as NaN?
  if mask_nan:
    flux_new[id_badpts], err_new[id_badpts] = np.nan, np.nan
    return flux_new, err_new

  # if not, fill with large values.
  else:
    err[id_badpts], mask_new[id_badpts] = 1.e10, 0
    return flux_new, err_new, mask_new

def rebin(wl, flux, err, mask, wl_new, mask_nan = False):

  # construct anti-deriv
  flux_atd = InterpolatedUnivariateSpline(wl, flux).antiderivative()
  flux_spl = InterpolatedUnivariateSpline(wl_new, flux_atd(wl_new)).derivative()

  err_intp = InterpolatedUnivariateSpline(wl[mask == 1], err[mask == 1], k = 1, ext = 3)
  # mask_spl = InterpolatedUnivariateSpline(wl, mask.astype('f8'), k = 1)

  # evaluate
  flux_new, err_new = flux_spl(wl_new), err_intp(wl_new)
  mask_new = np.ones(flux_new.size, dtype = 'i4')

  # mask bad points
  badpx_wl = wl[mask == 0]
  wlt_a, wlt_b = 2. * wl_new[0] - wl_new[1], 2. * wl_new[-1] - wl_new[-2]
  badpx_wl = badpx_wl[(badpx_wl > wlt_a) & (badpx_wl < wlt_b)]
  id_badpts = np.unique(np.searchsorted(0.5 * (wl_new[1:] + wl_new[:-1]), badpx_wl))

  # mask bad pixels as NaN?
  if mask_nan:
    flux_new[mask_new != 1], err_new[mask_new != 1] = np.nan, np.nan
    return flux_new, err_new

  # if not, fill with large values.
  else:
    err[mask_new < 0.5], mask_new[mask_new < 0.5] = 1.e10, 0
    return flux_new, err_new, mask_new.astype('i4')

# perform tests
if __name__ == "__main__":

  import fsps
  import matplotlib.pyplot as plt

  if 1: # test case #1

    # make a spectrum for testing
    wl_ida, wl_idb = 369, 557
    sp = fsps.StellarPopulation(sfh = 0, sf_start = 0.)
    flux_sp, mass_sp, L_sp = sp.ztinterp(0., 10., peraa = True)
    flux_sp, wl_sp = flux_sp[wl_ida: wl_idb + 1], sp.wavelengths[wl_ida: wl_idb + 1]

    # rebin and draw
    N_wlpts = 100
    wl_new = np.linspace(wl_sp[0], wl_sp[-1], N_wlpts)
    flux_new, err_new, mask_new = rebin(wl_sp, flux_sp, \
        flux_sp * 0.1, np.ones(wl_sp.shape), wl_new)

    # draw
    plt.plot(wl_sp, np.log10(flux_sp), c = 'r', label = "Original")
    plt.plot(wl_new, np.log10(flux_new), c = 'b', label = "Rebinned")
    plt.xlabel("wavelength [A]"), plt.ylabel("log Flux"), plt.legend(loc = 'best')
    plt.show()
