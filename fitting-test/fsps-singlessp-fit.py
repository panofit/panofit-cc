#!/usr/bin/python

'''
  Fit CALIFA V500 spectral cube of NGC 1 with single component (proof of concept)
  YJ Qin, Jan 2016 @ Shanghai

  Eight parameters to optimize: sersic_n, sersic_re, sersic_Ie, sersic_phi,
                                sersic_q, sersic_c, ssp_age, ssp_Z
  with a bounded minimization solver.

  First run [tnc solver]:
    n (!)            Re               Ie               phi
    1.00000000e-02   1.06920563e+01  -1.24847128e+01   4.10142121e-01
    7.47263830e-01   2.55921267e+00  -8.22526725e-01  -4.88497142e-02
    q                c                age              Z

  2nd run [tnc solver]: (with fixed log_Z = 0.)
    n                Re               Ie               phi
    0.94896813       6.65533725      -11.81804128      0.03088737
    1.               2.00342504       0.1750717        0.
    q (!)            c                age              Z
'''

# specify the name of the data cube
datacube_path = "./califa_sample/NGC0001.V500.rscube.fits"

import numpy as np
from scipy.optimize import fmin_tnc # bounded solver for multidimensional minimization
from scipy.interpolate import interp1d, UnivariateSpline # resampling the spectrum (not the optimal solution!)
import matplotlib.pyplot as plt

from astropy.io import fits
import fsps

import sys
import itertools as itt
import gc

# rebin the spectrum, TODO: with a better method.
def rebin(wl, flux, err, wl_new):

  # create interpolator
  spec_in = np.vstack((flux, err)).T
  sp_intp = interp1d(wl, spec_in, kind = 'linear', axis = 0, copy = False,
      bounds_error = False, fill_value = 0., assume_sorted = True)

  # interplate flux and err at new wavelengths
  spec_out = sp_intp(wl_new)
  flux_out, err_out = spec_out.T

  # correct bad pixels
  err_out[err_out == 0.] = 1.e10
  err_out[~np.isfinite(err_out)] = 1.e10

  return flux_out, err_out

# object function for fitting
def _min_objfc(par, spec_list, sp, wl_cut_idx):

  # unpack parameters to fit
  sersic_n, sersic_re, sersic_Ie, sersic_phi, \
      sersic_q, sersic_c, ssp_age, ssp_Z = par
  id_a, id_b = wl_cut_idx

  # DEBUG
  # ssp_Z = 0.

  # interpolate spectrum
  flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(ssp_Z, np.power(10., ssp_age), peraa = True)
  flux_ssp = flux_ssp[id_a: id_b + 1] * 1.e15

  print "Target function called with param:", par

  # calculate d_n
  ek1, ek2, ek3, ek4, ek5, ek6 = 3., -1. / 3., 8. / 1215., \
      184. / 229635., 1048. / 31000725., -17557576. / 1242974068875.
  d_n = ek1 * sersic_n + ek2 + ek3 / sersic_n + ek4 / (sersic_n ** 2) \
      + ek5 / (sersic_n ** 3) + ek6 / (sersic_n ** 4)

  # array of single-spectrum chi_sq
  chi_sq = np.zeros([len(spec_list)], dtype = 'f8')

  # loop over spectra
  for id_spec, spec_t in enumerate(spec_list):

    # unpack spectrum
    X_i, Y_i, flux_i, err_i = spec_t

    # calculate local "surface brightness"
    r_i = np.sqrt(X_i ** 2 + Y_i ** 2)
    cos_ri, sin_ri = X_i / r_i, Y_i / r_i
    if np.any(~np.isfinite(cos_ri)): cos_ri, sin_ri = 1., 0.
    cos_pt, sin_pt = np.cos(sersic_phi), np.sin(sersic_phi)

    X_w = r_i * (cos_ri * cos_pt + sin_ri * sin_pt)
    Y_w = r_i * (sin_ri * cos_pt - cos_ri * sin_pt)
    M_w = np.power( np.power(np.abs(X_w), sersic_c)
                  + np.power(np.abs(Y_w) / sersic_q, sersic_c), 1. / sersic_c)

    I_i = np.power(10., sersic_Ie) * np.exp(-d_n * (np.power(M_w / sersic_re, 1. / sersic_n) - 1.))

    # find chi_sq of this spectra
    chi_sq_i = ((flux_i - flux_ssp * I_i) / err_i) ** 2
    chi_sq_i[~np.isfinite(chi_sq_i)] = 0.

    # write chi_sq of this spectrum
    chi_sq[id_spec] = chi_sq_i.sum()

  # return a total chi_sq
  chi_sq_sum = np.sum(chi_sq)

  print "  Target function exit with", chi_sq_sum, "\n"

  return chi_sq_sum

def fit_datacube(filename):

  # read spectrum from CALIFA cubes
  hlst = fits.open(filename)

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

  # valid pixels
  valid_spaxels = []

  # for each spaxel, if valid, add into list
  for I_ra, I_dec in itt.product(range(N_ra), range(N_dec)):
    if np.sum(flux[:, I_dec, I_ra]) != 0.: valid_spaxels.append((I_ra, I_dec))

  # how many valid spaxels?
  print len(valid_spaxels), "spectra are valid in this cube."

  # make stellar pop
  sp = fsps.StellarPopulation(sfh = 0, sf_start = 0.)

  # cut the valid range of wavelength
  id_a, id_b = np.searchsorted(sp.wavelengths, [wl_ax[0], wl_ax[-1]], side = 'left')
  wl_new = (sp.wavelengths)[id_a: id_b + 1]

  # interpolate the cube.
  flux_intp = interp1d(wl_ax, flux, kind = 'linear', axis = 0, copy = False,
      bounds_error = False, fill_value = 0., assume_sorted = True)
  err_intp = interp1d(wl_ax, err, kind = 'nearest', axis = 0, copy = False,
      bounds_error = False, fill_value = 0., assume_sorted = True)
  flux_new = flux_intp(wl_new)
  err_new = err_intp(wl_new)

  # close fits file.
  hlst.close()

  # force delete objects.
  del flux_intp, err_intp, flux, err, mask
  gc.collect()

  # make the list of spectra to fit.
  spec_list = []
  for I_ra, I_dec in valid_spaxels:
    spec_list.append((ra_ax[I_ra], dec_ax[I_dec], \
        flux_new[:, I_dec, I_ra], err_new[:, I_dec, I_ra]))

  # call fitting routine.
  par_opt, N_eval, ret = fmin_tnc(_min_objfc, (2., 1., 1., 0., 0.5, 2., 0., 0.),
      args = (spec_list, sp, (id_a, id_b)), approx_grad = True,
      bounds = [(1.e-2, 8.), (1.e-2, 4.e1), (-16., 16.), (-np.pi, np.pi),
                (1.e-1, 1.), (1.e-1, 4.), (-3, 1.1), (-1.98, 0.2)])

  print par_opt

if __name__ == "__main__": fit_datacube(datacube_path)
