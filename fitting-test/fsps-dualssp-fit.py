#!/usr/bin/python

'''
  Fit CALIFA V500 spectral cube of NGC 1 with two components (proof of concept)
  YJ Qin, Jan 2016 @ Shanghai

  Parameters to optimize: sersic_n, sersic_re, sersic_Ie,  sersic_phi,
                          sersic_q, sersic_c,  sersic_age,
                          exp_Ic,   exp_h,     exp_phi,    exp_q,
                          exp_c,    exp_age,
  with a bounded minimization solver.

  First run [tnc solver]: (looks reasonable)
    1.20189845   5.99849347 -12.0130215   -0.31569187   0.51323838
    1.99219028   1.09990431 -11.29316191  12.61844848  -0.14629537
    0.57732247   2.31568646   0.15574498

  Second run [tnc solver]:
    0.86428252   2.7846034  -15.36964374  -2.22694922   0.80212827
    1.06718218   1.01979302 -11.08551644   9.79371104  -2.06361907
    0.98261746   2.42270186   0.32113472
'''

# specify the name of the data cube
datacube_path = "./califa_sample/NGC0001.V500.rscube.fits"

# specify the solver
min_solver = "tnc"

# import everything
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin_tnc, differential_evolution # bounded multivar minimization solver
from scipy.interpolate import interp1d
import emcee

import fsps # python-fsps
from astropy.io import fits

import sys, gc
import itertools as itt

# object function for fitting
def _min_objfc(par, spec_list, sp, wl_cut_idx):

  # unpack parameters to fit
  sersic_n, sersic_re, sersic_Ie, sersic_phi, sersic_q, sersic_c, sersic_age, \
      exp_Ic, exp_h, exp_phi, exp_q, exp_c, exp_age = par
  id_a, id_b = wl_cut_idx

  print "Target function called with param:", par

  # calculate d_n for sersic component
  ek1, ek2, ek3, ek4, ek5, ek6 = 3., -1. / 3., 8. / 1215., \
      184. / 229635., 1048. / 31000725., -17557576. / 1242974068875.
  d_n = ek1 * sersic_n + ek2 + ek3 / sersic_n + ek4 / (sersic_n ** 2) \
      + ek5 / (sersic_n ** 3) + ek6 / (sersic_n ** 4)

  # array of single-spectrum chi_sq
  chi_sq = np.zeros([len(spec_list)], dtype = 'f8')

  # interpolate ssp
  flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(0., np.power(10., sersic_age), peraa = True)
  flux_s = flux_ssp[id_a: id_b + 1] * 1.e15

  flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(0., np.power(10., exp_age), peraa = True)
  flux_e = flux_ssp[id_a: id_b + 1] * 1.e15

  # loop over spectra
  for id_spec, spec_t in enumerate(spec_list):

    # print "\r", id_spec, "...",
    # sys.stdout.flush()

    # unpack spectrum
    X_i, Y_i, flux_i, err_i = spec_t

    # pos on the sky
    r_i = np.sqrt(X_i ** 2 + Y_i ** 2)
    cos_ri, sin_ri = X_i / r_i, Y_i / r_i
    if np.any(~np.isfinite(cos_ri)): cos_ri, sin_ri = 1., 0.

    # calculate "surface brightness" of sersic component
    cos_ps, sin_ps = np.cos(sersic_phi), np.sin(sersic_phi)
    X_s = r_i * (cos_ri * cos_ps + sin_ri * sin_ps)
    Y_s = r_i * (sin_ri * cos_ps - cos_ri * sin_ps)
    M_s = np.power( np.power(np.abs(X_s), sersic_c)
                  + np.power(np.abs(Y_s) / sersic_q, sersic_c), 1. / sersic_c)
    I_s = np.power(10., sersic_Ie) * np.exp(-d_n * (np.power(M_s / sersic_re, 1. / sersic_n) - 1.))

    # calculate "surface brightness" of exponential component
    cos_pe, sin_pe = np.cos(exp_phi), np.sin(exp_phi)
    X_e = r_i * (cos_ri * cos_pe + sin_ri * sin_pe)
    Y_e = r_i * (sin_ri * cos_pe - cos_ri * sin_pe)
    M_e = np.power( np.power(np.abs(X_e), exp_c)
                  + np.power(np.abs(Y_e) / exp_q, exp_c), 1. / exp_c)
    I_e = np.power(10., exp_Ic) * np.exp(-np.abs(M_e / exp_h))

    # model spectrum
    flux_model = flux_s * I_s + flux_e * I_e

    # find chi_sq of this spectra
    chi_sq_i = ((flux_i - flux_model) / err_i) ** 2
    chi_sq_i[~np.isfinite(chi_sq_i)] = 0.

    # write chi_sq of this spectrum
    chi_sq[id_spec] = chi_sq_i.sum()

  # return a total chi_sq
  chi_sq_sum = np.sum(chi_sq)

  print "  Target function gets chi_sq:", chi_sq_sum, "\n"

  return chi_sq_sum

# log_L function, for MCMC sampling.
def _log_prob(par, spec_list, sp, wl_cut_idx, bounds):

  # check bounds: within bounds, return log_L, else return -inf
  if np.all(np.logical_and(np.asarray(bounds)[:, 0] < np.asarray(par), \
                           np.asarray(par) < np.asarray(bounds)[:, 1])):
    return -0.5 * _min_objfc(par, spec_list, sp, wl_cut_idx)
  else: return -np.inf

# *almost* the same as single-ssp case.
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

  '''
    sersic_n, sersic_re, sersic_Ie,  sersic_phi,
    sersic_q, sersic_c,  sersic_age,
    exp_Ic,   exp_h,     exp_phi,    exp_q,
    exp_c,    exp_age,
  '''

  # call fitting routine.

  init_guess = (2., 1., 1., 0., 0.5, 2., 1., 1., 1., 0., 0.5, 2., 0.)
  min_bounds = [(1.e-2, 8.), (1.e-2, 4.e1), (-16., 16.), (-np.pi, np.pi),
                (1.e-1, 1.), (1.e-1, 4.),   (-3, 1.1),
                (-16., 16.), (1.e-2, 4.e1), (-np.pi, np.pi), (1.e-1, 1.),
                (1.e-1, 4.), (-3, 1.1)]

  if min_solver == "tnc":
    par_opt, N_eval, ret = fmin_tnc(_min_objfc, init_guess,
        args = (spec_list, sp, (id_a, id_b)), approx_grad = True,
        bounds = min_bounds)
    print par_opt

  elif min_solver == "diffevo":
    res = differential_evolution(_min_objfc, bounds = min_bounds,
        args = (spec_list, sp, (id_a, id_b)))
    print res.x

  elif min_solver == "mcmc":

    raise RuntimeError("MCMC solver doesn't work well.")

    N_dim, N_walkers = 13, 32
    MC_sampler = emcee.EnsembleSampler(N_walkers, N_dim, _log_prob,
        args = (spec_list, sp, (id_a, id_b), min_bounds))
    MC_sampler.run_mcmc(init_guess, 1)

    samples = MC_sampler.chain
    np.save("mc_sampling.dat", samples)

if __name__ == "__main__": fit_datacube(datacube_path)
