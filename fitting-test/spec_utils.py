#!/usr/bin/python

'''
  Functions for spectra resampling
'''

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d

def rebin(wl, flux, err, mask, wl_new, mask_nan = False):

  # construct anti-deriv
  #'''
  flux_atd = InterpolatedUnivariateSpline(wl, flux).antiderivative()
  flux_spl = InterpolatedUnivariateSpline(wl_new, flux_atd(wl_new)).derivative()
  #'''
  # flux_spl = InterpolatedUnivariateSpline(wl, flux).antiderivative().derivative() # -> fig2

  err_intp = InterpolatedUnivariateSpline(wl, err)
  mask_spl = interp1d(wl, mask, kind = 'nearest', copy = False,
      bounds_error = False, fill_value = 0., assume_sorted = True)

  # evaluate
  flux_new, err_new, mask_new = flux_spl(wl_new), err_intp(wl_new), mask_spl(wl_new)

  # correct bad pixels
  if mask_nan:
    flux_new[mask_new == 0] = np.nan
    err_new[mask_new == 0]  = np.nan
    return flux_new, err_new

  else: return flux_new, err_new, mask_new

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
        np.zeros_like(wl_sp), np.ones(wl_sp.shape), wl_new)

    # draw
    plt.plot(wl_sp, np.log10(flux_sp), c = 'r', label = "Original")
    plt.plot(wl_new, np.log10(flux_new), c = 'b', label = "Rebinned")
    plt.xlabel("wavelength [A]"), plt.ylabel("log Flux"), plt.legend(loc = 'best')
    plt.show()
