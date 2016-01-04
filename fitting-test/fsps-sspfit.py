#!/usr/bin/python

'''
  Match spectra of NGC 1 / NGC 171 with python-fsps
  (as a proof-of-concept, to evaluate the effect of redshift and LOSVD)
  YJ Qin, Jan 2016 @ Shanghai

  Two parameters to optimize: ssp_age, ssp_Z
'''

# specify the name of the data cube
datacube_fname = "./califa_sample/NGC0001.V500.rscube.fits"

import numpy as np
from scipy.optimize import fmin_tnc # bounded solver for multidimensional minimization
from scipy.interpolate import interp1d # resampling the spectrum (not the optimal solution!)
import matplotlib.pyplot as plt

from astropy.io import fits
import fsps

import sys
import itertools as itt

# re-bin a spectrum
def rebin(wl, flux, err, out_wl):

  # FIXME: turn simple spline interpolation to realistic rebinning.

  # create interpolator
  spec_in = np.vstack((flux, err)).T
  sp_intp = interp1d(wl, spec_in, kind = 'linear', axis = 0, copy = False,
      bounds_error = False, fill_value = 0., assume_sorted = True)

  # interplate flux and err at new wavelengths
  spec_out = sp_intp(out_wl)
  flux_out, err_out = spec_out.T

  # correct bad pixels
  err_out[err_out == 0.] = 1.e9
  err_out[~np.isfinite(err_out)] = 1.e9

  return flux_out, err_out

# object function to optimize
def _min_objfc(par, flux_obs, err_obs, id_a, id_b, sp):

  # unpack parameters
  ssp_age, ssp_Z, flux_scale = par

  # interpolate spectrum, cut and rescale
  flux_model, mass, lbol = sp.ztinterp(ssp_Z, np.power(10., ssp_age), peraa = True)
  flux_model = flux_model[id_a: id_b + 1]
  flux_model = flux_model / np.nanmean(flux_model) * np.power(10., flux_scale)

  # print ssp_age, ssp_Z, flux_scale, id_a, id_b
  # plt.plot(np.log10(flux_obs)), plt.plot(np.log10(flux_model)), plt.show()

  # calculate chi_sq and return
  chi_sq = ((flux_obs - flux_model) / err_obs) ** 2
  chi_sq[~np.isfinite(chi_sq)] = 0.

  # print chi_sq.sum()

  return chi_sq.sum()

# match a single spectrum
def match_spectrum(wl, flux, err):

  # make stellar population
  sp = fsps.StellarPopulation(sfh = 0, sf_start = 0.)
  # FIXME: make sp object somewhere else.

  # determine the range of lambda
  id_a, id_b = np.searchsorted(sp.wavelengths, [wl[0], wl[-1]], side = 'left')
  wl_new = (sp.wavelengths)[id_a: id_b + 1]

  # rebin the spectrum
  flux_new, err_new = rebin(wl, flux, err, wl_new)

  # call minimization solver
  par_opt, N_eval, ret = fmin_tnc(_min_objfc, (0., 0., 1.), args = (flux_new, err_new, id_a, id_b, sp),
      approx_grad = True, bounds = [(-3, 1.1), (-1.98, 0.2), (-10., 10.)])

  # return best solution
  return par_opt

# main program
if __name__ == "__main__":

  # read spectrum from CALIFA cubes
  hlst = fits.open(datacube_fname)

  # read basic information
  N_ra, N_dec = hlst[0].header['NAXIS1'], hlst[0].header['NAXIS2']
  N_spx = hlst[0].header['NAXIS3']

  # make wavelength axis, in A
  wl_ax = hlst[0].header['CRVAL3'] + hlst[0].header['CDELT3'] * np.arange(N_spx)

  # RA and Dec axes, in asec, origin at the galaxy center
  ra_ax  = -hlst[0].header['CRPIX1'] + hlst[0].header['CDELT1'] * np.arange(N_ra)
  dec_ax = -hlst[0].header['CRPIX2'] + hlst[0].header['CDELT2'] * np.arange(N_dec)

  # data, err and mask
  flux, err, mask = hlst[0].data, hlst[1].data, hlst[2].data

  # fitting output
  ssp_age, ssp_Z = np.zeros([N_ra, N_dec]), np.zeros([N_ra, N_dec])

  # number of valid pixels fitted.
  N_valid = 0

  # for each spaxel
  for I_ra, I_dec in itt.product(range(N_ra), range(N_dec)):

    # if not valid, fill with nan
    if np.sum(flux[:, I_dec, I_ra]) == 0.:
      ssp_age[I_ra, I_dec], ssp_Z[I_ra, I_dec] = float("NaN"), float("NaN")

    # do fitting
    else:
      flux_pt = flux[:, I_dec, I_ra].astype('f8')
      err_pt = err[:, I_dec, I_ra].astype('f8')
      opt = match_spectrum(wl_ax, flux_pt, err_pt)
      ssp_age[I_ra, I_dec], ssp_Z[I_ra, I_dec] = opt[0], opt[1]

      N_valid += 1

      print "N_valid:", N_valid, " I_ra:", I_ra, " I_dec:", I_dec, " log(Age):", opt[0], " log(Z/Zs):", opt[1]

  # DEBUG
  # plt.ylim(None, 5.), plt.plot(wl_ax, flux_ct), plt.plot(wl_ax, err_ct), plt.show()

  # save results
  np.save("ssp_age.dat", ssp_age)
  np.save("ssp_Z.dat", ssp_Z)
