#!/usr/bin/python

'''
  Fit CALIFA V500 spectral cube of NGC 1 with single component (proof of concept)
  YJ Qin, Jan 2016 @ Shanghai

  Eight parameters to optimize: sersic_n, sersic_re, sersic_Ie, sersic_phi,
                                sersic_q, sersic_c, ssp_age_c, ssp_age_k
  with a bounded minimization solver.

'''

# FIXME it doesn't work now...

# specify the name of the data cube
datacube_path = "./califa_sample/NGC0001.V500.rscube.fits"

mcsolver = "mc" # or pt

from spec_utils import *

import numpy as np
from scipy.optimize import fmin_tnc # bounded solver for multidimensional minimization
from scipy.interpolate import interp1d, interp2d, griddata, RectBivariateSpline

from astropy.io import fits
import fsps, emcee

import sys, time, gc
import itertools as itt

'''
# NOTE depricated. use spec_utils.rebin!
def rebin(wl, flux, err, mask, wl_new):

  raise RuntimeError("'rebin' depricated, use spec_utils.rebin instead.")

  # do interpolation
  flux_new = griddata(wl, flux, wl_new, method = 'linear')
  err_new  = griddata(wl, err,  wl_new, method = 'linear')
  mask_new = griddata(wl, mask, wl_new, method = 'nearest')

  # correct bad pixels
  flux_new[mask_new == 0] = np.nan
  err_new[mask_new == 0]  = np.nan

  return flux_new, err_new
'''

# object function for fitting
def _min_objfc(par, spec_list, spec_lib, ax_age, ax_Z):

  # unpack parameters to fit
  sersic_n, sersic_re, sersic_Ie, sersic_phi, \
      sersic_q, sersic_c, ssp_age_c, ssp_age_k = par

  # interpolate spectrum
  # NOTE: replaced with my own bilinear interpolator.
  # flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(ssp_Z, np.power(10., ssp_age), peraa = True)
  # flux_ssp = flux_ssp[id_a: id_b + 1] * 1.e15

  # print "Target function called with param:", par

  # calculate d_n
  ek1, ek2, ek3, ek4, ek5, ek6 = 3., -1. / 3., 8. / 1215., \
      184. / 229635., 1048. / 31000725., -17557576. / 1242974068875.
  d_n = ek1 * sersic_n + ek2 + ek3 / sersic_n + ek4 / (sersic_n ** 2) \
      + ek5 / (sersic_n ** 3) + ek6 / (sersic_n ** 4)

  # array of single-spectrum chi_sq
  chi_sq = np.zeros([len(spec_list)], dtype = 'f8')
  del_pc = np.zeros([len(spec_list)], dtype = 'f8')
  del_pd = np.zeros([len(spec_list)], dtype = 'f8')

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

    # inerpolate spectrum
    age_t = ssp_age_c + ssp_age_k * (M_w / sersic_re)
    id_a, id_b, w = None, np.searchsorted(ax_age, age_t), 1.

    if id_b == 0: id_a, id_b, w = 0, 0, 1.
    elif id_b == ax_age.size: id_a, id_b, w = ax_age.size - 1, ax_age.size - 1, 1.
    else: id_a = id_b - 1; w = (age_t - ax_age[id_a]) / (ax_age[id_b] - ax_age[id_a])

    flux_m = I_i * (spec_lib[id_b] * w + spec_lib[id_a] * (1. - w))

    # find chi_sq of this spectra
    chi_sq_i = ((flux_i - flux_m) / err_i) ** 2
    chi_sq_i[~np.isfinite(chi_sq_i)] = 0.

    # write chi_sq of this spectrum
    chi_sq[id_spec] = chi_sq_i.sum()

    # calculate average precentage diff
    del_pc[id_spec] = np.nanmean(np.abs((flux_i - flux_m) / flux_i))
    del_pd[id_spec] = np.nanmean(flux_i / flux_m)

  # return a total chi_sq
  chi_sq_sum = np.sum(chi_sq)

  # print "  Target function exit with", chi_sq_sum, "\n"

  print "% e "%chi_sq_sum,
  print "% e "%(np.nanmean(del_pc)),
  print "% e "%(np.nanmean(del_pd)),
  for i in par: print "% e "%i,
  print " "

  return chi_sq_sum

# log likelihood
def _log_prob(par, spec_list, spec_lib, ax_age, ax_Z, bounds):

  # check bounds: within bounds, return log_L, else return -inf
  if np.all(np.logical_and(np.asarray(bounds)[:, 0] < np.asarray(par), \
                           np.asarray(par) < np.asarray(bounds)[:, 1])):
    return -0.5 * _min_objfc(par, spec_list, spec_lib, ax_age, ax_Z)
  else: return -np.inf

def fit_datacube(filename, out_fname):

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

  # DBG
  '''
  flux[mask == 0] = np.nan
  print np.nanmean(flux)
  sys.exit()
  # '''

  # valid pixels
  valid_spaxels = []

  # for each spaxel, if valid, add into list
  for I_ra, I_dec in itt.product(range(N_ra), range(N_dec)):
    if np.sum(flux[:, I_dec, I_ra]) != 0.: valid_spaxels.append((I_ra, I_dec))

  # how many valid spaxels?
  print len(valid_spaxels), "spectra are valid i-n this cube."

  # make stellar pop
  sp = fsps.StellarPopulation(sfh = 0, sf_start = 0.)

  # cut the valid range of wavelength
  id_a, id_b = np.searchsorted(sp.wavelengths, [wl_ax[0], wl_ax[-1]], side = 'left')
  wl_new = (sp.wavelengths)[id_a: id_b + 1]

  # prepare spec lib: spec_lib, ax_age, ax_Z
  ax_age = np.linspace(-3., 1.1, 41)
  ax_Z = np.zeros((1)) # reserved for future use

  # interpolate
  spec_lib = np.zeros((ax_age.size, wl_new.size))
  for i_age, v_age in enumerate(ax_age):
    flux_t, m_t, l_t = sp.ztinterp(0., np.power(10., v_age), peraa = True)
    spec_lib[i_age, :] = flux_t[id_a: id_b + 1]

  # interpolate the cube.
  # NOTE rebin the observed spectra to match wl axis of speclib.
  # do this later when collecting valid spaxels
  '''
  flux_intp = interp1d(wl_ax, flux, kind = 'linear', axis = 0, copy = False,
      bounds_error = False, fill_value = 0., assume_sorted = True)
  err_intp = interp1d(wl_ax, err, kind = 'nearest', axis = 0, copy = False,
      bounds_error = False, fill_value = 0., assume_sorted = True)
  flux_new = flux_intp(wl_new)
  err_new = err_intp(wl_new)
  '''

  # close fits file.
  hlst.close()

  # force delete objects.
  # del flux_intp, err_intp, flux, err, mask
  # gc.collect()

  # make the list of spectra to fit.
  spec_list = []
  for I_ra, I_dec in valid_spaxels:

    # rebin the spectrum
    flux_i, err_i = rebin(wl_ax, flux[:, I_dec, I_ra], \
        err[:, I_dec, I_ra], mask[:, I_dec, I_ra], wl_new, mask_nan = True)

    # put into the list.
    spec_list.append((ra_ax[I_ra], dec_ax[I_dec], flux_i, err_i))

  # call fitting routine. NOTE TNC solver not used.
  '''
  par_opt, N_eval, ret = fmin_tnc(_min_objfc, (2., 1., 1., 0., 0.5, 2., 0., 0.),
      args = (spec_list, sp, (id_a, id_b)), approx_grad = True,
      bounds = [(1.e-2, 8.), (1.e-2, 4.e1), (-16., 16.), (-np.pi, np.pi),
                (1.e-1, 1.), (1.e-1, 4.), (-3, 1.1), (-1.98, 0.2)])
  print par_opt
  '''

  # call MCMC
  init_guess = (1., 1.e0, 1., 0.0, 0.75, 2., 1.0, -0.02)
  min_bounds = [(1.e-2, 1.e1), (1.e-2, 8.e1), (-5., 5.), (-np.pi, np.pi),
                (1. / 5., 1.), (1.e-1, 4.  ), (-3, 1.1), (-1.4, 0.)]

  # DEBUG test call
  #print _log_prob(init_guess, spec_list, spec_lib, ax_age, ax_Z, min_bounds)
  #sys.exit()

  # with "normal" monte-carlo
  if mcsolver == "mc":

    N_dim, N_walkers = len(init_guess), len(init_guess) * 2
    MC_sampler = emcee.EnsembleSampler(N_walkers, N_dim, _log_prob,
        args = (spec_list, spec_lib, ax_age, ax_Z, min_bounds))

    init_sol = np.outer(np.ones(N_walkers), np.array(init_guess)) \
             + 2.5e-3 * np.random.randn(N_walkers * N_dim).reshape((N_walkers, N_dim))

    t0 = time.clock()
    MC_sampler.run_mcmc(init_sol, 2500)
    t1 = time.clock()
    print "takes", t1 - t0, "sec to run."

    samples = MC_sampler.chain
    np.save(out_fname, samples)

  # with parallel tempering
  elif mcsolver == "pt":

    N_temps = 4
    N_dim, N_walkers = len(init_guess), len(init_guess) * 2
    # 8 * 8 * 16 = 512

    _log_pri = lambda par, spec_list, spec_lib, ax_age, ax_Z, bounds: 0.

    MC_sampler = emcee.PTSampler(N_temps, N_walkers, N_dim, _log_prob, _log_pri,
        threads = 1, betas = None, a = 2., Tmax = None,
        loglargs = (spec_list, spec_lib, ax_age, ax_Z, min_bounds),
        logpargs = (spec_list, spec_lib, ax_age, ax_Z, min_bounds))

    init_sol = np.multiply.outer(np.ones(N_temps),\
        np.outer(np.ones(N_walkers), np.array(init_guess)))
    init_sol += 2.5e-1 * np.random.randn(init_sol.size).reshape(init_sol.shape)

    t0 = time.clock()
    MC_sampler.run_mcmc(init_sol, 1000)
    t1 = time.clock()

    print "takes", t1 - t0, "sec to run."

    samples = MC_sampler.chain
    np.save(out_fname, samples)

if __name__ == "__main__": fit_datacube(datacube_path, "fsps-singlessp-agegradient-test.dat")
