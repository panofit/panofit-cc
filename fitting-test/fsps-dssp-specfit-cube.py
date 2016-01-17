#!/usr/bin/python

'''
  Fit CALIFA V500 spectral cube of NGC 1 with two components (proof of concept)
  YJ Qin, Jan 2016 @ Shanghai

  Parameters to optimize: sersic_n,   sersic_re,  sersic_Ie,  sersic_phi, sersic_q
                          sersic_c,   sersic_age, exp_Ic,     exp_h,      exp_phi,
                          exp_q,      exp_c,      exp_age,
  with a Monte-carlo solver.

  # Note: the mock model here is generated as cubes.
'''

# specify the name of the data cube
datacube_path = "./califa_sample/NGC0001.V500.rscube.fits"
source_redshift = 0.0146275
PSF_size = 2.

# specify the solver
min_solver = "mcmc"

# import everything
import numpy as np
import matplotlib.pyplot as plt
import emcee

import fsps # python-fsps
from astropy.io import fits

from scipy.ndimage.filters import gaussian_filter

import sys, gc, time
import itertools as itt

from spec_utils import *

def generate_mock_datacube(par, sp, ra_ax, dec_ax, id_a, id_b):

  # unpack param
  sersic_n, sersic_re, sersic_Ie, sersic_phi, sersic_q, sersic_c, sersic_age, \
      exp_Ic, exp_h, exp_phi, exp_q, exp_c, exp_age = par

  # calculate structural information
  ek1, ek2, ek3, ek4, ek5, ek6 = 3., -1. / 3., 8. / 1215., \
      184. / 229635., 1048. / 31000725., -17557576. / 1242974068875.
  d_n = ek1 * sersic_n + ek2 + ek3 / sersic_n + ek4 / (sersic_n ** 2) \
      + ek5 / (sersic_n ** 3) + ek6 / (sersic_n ** 4)

  # interpolate ssp
  flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(0., np.power(10., sersic_age), peraa = True)
  flux_s = flux_ssp[id_a: id_b + 1]
  flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(0., np.power(10., exp_age), peraa = True)
  flux_e = flux_ssp[id_a: id_b + 1]

  # mock cube, follow CALIFA style.
  cube_t = np.zeros([id_b - id_a + 1, dec_ax.size, ra_ax.size])

  # loop over spectra,
  for I_ra, v_ra in enumerate(ra_ax):
    for I_dec, v_dec in enumerate(dec_ax):

      # calculate structural param
      r_i = np.sqrt(v_ra ** 2 + v_dec ** 2)
      cos_ri, sin_ri = v_ra / r_i, v_dec / r_i
      if np.any(~np.isfinite(cos_ri)): cos_ri, sin_ri = 1., 0.

      # calculate "surface brightness" of sersic component
      cos_ps, sin_ps = np.cos(sersic_phi), np.sin(sersic_phi)
      X_s = r_i * (cos_ri * cos_ps + sin_ri * sin_ps)
      Y_s = r_i * (sin_ri * cos_ps - cos_ri * sin_ps)
      M_s = np.power( np.power(np.abs(X_s), sersic_c)
                    + np.power(np.abs(Y_s) / sersic_q, sersic_c), 1. / sersic_c)
      I_s = np.power(2.512, -sersic_Ie) * np.exp(-d_n * (np.power(M_s / sersic_re, 1. / sersic_n) - 1.))

      # calculate "surface brightness" of exponential component
      cos_pe, sin_pe = np.cos(exp_phi), np.sin(exp_phi)
      X_e = r_i * (cos_ri * cos_pe + sin_ri * sin_pe)
      Y_e = r_i * (sin_ri * cos_pe - cos_ri * sin_pe)
      M_e = np.power( np.power(np.abs(X_e), exp_c)
                    + np.power(np.abs(Y_e) / exp_q, exp_c), 1. / exp_c)
      I_e = np.power(2.512, -exp_Ic) * np.exp(-np.abs(M_e / exp_h))

      # model spectrum
      cube_t[:, I_dec, I_ra] = flux_s * I_s + flux_e * I_e

  return cube_t

# object function for fitting
def _min_objfc(par, cube, err, sp, ra_ax, dec_ax, wl_ida, wl_idb):

  t0 = time.clock()

  # make the cube.
  cube_m = generate_mock_datacube(par, sp, ra_ax, dec_ax, wl_ida, wl_idb)

  t1 = time.clock()

  # Gaussian blurring
  cube_new = np.zeros(cube_m.shape)
  for I_wl in xrange(cube_m.shape[0]):
    cube_new[I_wl, :, :] = gaussian_filter(cube_m[I_wl, :, :], PSF_size)

  t2 = time.clock()

  # calculate chi_sq
  chi_sq = np.nanmean(((cube_m - cube) / err) ** 2)

  # print chi_sq and param
  print "% E"%chi_sq,
  for vpar in par: print "% E"%vpar,
  print ' '

  return chi_sq

# log_L function, for MCMC sampling.
def _log_prob(par, cube, err, sp, ra_ax, dec_ax, wl_ida, wl_idb, bounds):

  # check bounds: within bounds, return log_L, else return -inf
  if np.all(np.logical_and(np.asarray(bounds)[:, 0] < np.asarray(par), \
                           np.asarray(par) < np.asarray(bounds)[:, 1])):
    return -0.5 * _min_objfc(par, cube, err, sp, ra_ax, dec_ax, wl_ida, wl_idb)
  else: return -np.inf

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

  # data, err and mask, clear bad pixels
  flux, err, mask = hlst[0].data, hlst[1].data, hlst[2].data

  # DEBUG
  # sys.exit()

  # pre-processing of the data cube.
  mask[err > 1.e5] = 0                          # mask zero-value pixels
  mask[np.log10(np.abs(flux / err)) > 2.4] = 0  # mask abnormal SNR

  # make stellar pop
  sp = fsps.StellarPopulation(sfh = 0, sf_start = 0.)

  # cut the valid range of wavelength
  id_a, id_b = np.searchsorted(sp.wavelengths,
      [wl_ax[0] / (1. + source_redshift), wl_ax[-1] / (1. + source_redshift)])
  wl_new = (sp.wavelengths)[id_a: id_b + 1]

  # new cubes for fitting
  flux_new = np.zeros((wl_new.size, dec_ax.size, ra_ax.size))
  err_new  = np.zeros((wl_new.size, dec_ax.size, ra_ax.size))
  mask_new = np.zeros((wl_new.size, dec_ax.size, ra_ax.size), dtype = 'i4')

  # de-redshift and rebin
  for i_ra in np.arange(N_ra):
    for i_dec in np.arange(N_dec):
      flux_i, err_i, mask_i = restframe(wl_ax, flux[:, i_dec, i_ra],
          err[:, i_dec, i_ra], mask[:, i_dec, i_ra], wl_new,
          source_redshift, mask_nan = False)
      flux_new[:, i_dec, i_ra] = flux_i
      err_new[:, i_dec, i_ra] = err_i
      mask_new[:, i_dec, i_ra] = mask_i

  # close fits file.
  hlst.close()

  # force delete objects.
  del flux, err, mask
  gc.collect()

  # pre-processing of spectral cubes
  flux_new[mask_new < 0.5] = np.nan
  err_new[mask_new < 0.5]  = np.nan

  # DEBUG PLOT
  ''' # works!
  plt.imshow(flux_new[:, 36, :], interpolation = 'nearest', aspect = 'auto')
  plt.show()
  sys.exit()
  #'''

  '''
    Parameters:
      sersic_n, sersic_re, sersic_Ie,  sersic_phi,
      sersic_q, sersic_c,  sersic_age,
      exp_Ic,   exp_h,     exp_phi,    exp_q,
      exp_c,    exp_age,
  '''

  # call fitting routine.

  # init condition and bounds

  init_guess = (2.0, 12., -8.0, 0.00, 1.00, 2.00, 1.0, -8., 15., 0.85, 0.75, 2.0, 0.3) # for N0001

  # for N7716
  '''
  init_guess = (2.0,  # sersic n
                12.,  # sersic re
               -9.0,  # sersic Ie (log-2.512-based)
                0.00, # sersic phi
                1.00, # sersic q
                2.00, # sersic c
                1.0,  # sersic age
               -9.0,  # exp Ic
                12.,  # exp h
                0.85, # exp phi
                0.75, # exp q
                2.0,  # exp c
                0.3)  # exp_age
  '''

  min_bounds = [(1.e-1, 8.),                # sersic n
                (1.e0, 3.5e1),              # sersuc re
                (-16., 0.),                 # sersic Ie
                (-np.pi / 2., np.pi / 2.),  # sersic phi
                (2.e-1, 1.),                # sersic q
                (1., 3.),                   # sersic c
                (-3, 1.1),                  # sersic age
                (-16., 0.),                 # exp Ic
                (1.e-1, 3.5e1),             # exp h
                (-np.pi / 2., np.pi / 2.),  # exp phi
                (2.e-1, 1.),                # exp q
                (1., 3.),                   # exp c
                (-3, 1.1)]                  # exp age

  if min_solver == "tnc":
    par_opt, N_eval, ret = fmin_tnc(_min_objfc, init_guess,
        args = (flux_new, err_new, sp, ra_ax, dec_ax, id_a, id_b),
        approx_grad = True, bounds = min_bounds, messages = 0)
    print par_opt

  elif min_solver == "diffevo":
    res = differential_evolution(_min_objfc, bounds = min_bounds,
        args = (flux_new, err_new, sp, ra_ax, dec_ax, id_a, id_b))
    print res.x

  elif min_solver == "mcmc":

    N_dim, N_walkers = len(init_guess), len(init_guess) * 2

    # read init condition
    '''
    init_guess = np.outer(np.ones(N_walkers), np.array(init_guess))
    init_guess = init_guess + np.random.randn(init_guess.size).reshape(init_guess.shape) * 2.5e-1
    #'''

    # DEBUG
    #'''
    chain_pt = np.load("./fsps-dssp-specfit-cube-NGC0001V500-PSF2/NGC0001.V500.rscube.fits.chain-02.npy")
    init_guess = chain_pt[:, -1, :]
    #'''

    MC_sampler = emcee.EnsembleSampler(N_walkers, N_dim, _log_prob,
        args = (flux_new, err_new, sp, ra_ax, dec_ax, id_a, id_b, min_bounds))
    t0 = time.clock()
    MC_sampler.run_mcmc(init_guess, 1500)
    t1 = time.clock()

    print "Takes", t1 - t0, "sec to run."

    samples = MC_sampler.chain
    np.save(datacube_path.split('/')[-1] + '.chain', samples)

  elif min_solver == "pt":

    N_temps = 4
    N_dim, N_walkers = len(init_guess), len(init_guess) * 2
    # 8 * 8 * 16 = 512

    _log_pri = lambda par, flux_new, err_new, sp, ra_ax, dec_ax, id_a, id_b, min_bounds: 0.

    MC_sampler = emcee.PTSampler(N_temps, N_walkers, N_dim, _log_prob, _log_pri,
        threads = 1, betas = None, a = 2., Tmax = None,
        loglargs = (flux_new, err_new, sp, ra_ax, dec_ax, id_a, id_b, min_bounds),
        logpargs = (flux_new, err_new, sp, ra_ax, dec_ax, id_a, id_b, min_bounds))

    init_sol = np.load("./NGC2410.V500.rscube.fits.ptchain.npy")[:, :, -1, :]

    '''
    init_sol = np.multiply.outer(np.ones(N_temps),\
        np.outer(np.ones(N_walkers), np.array(init_guess)))
    init_sol += 2.5e-1 * np.random.randn(init_sol.size).reshape(init_sol.shape)
    '''

    t0 = time.clock()
    MC_sampler.run_mcmc(init_sol, 1500)
    t1 = time.clock()

    print "takes", t1 - t0, "sec to run."

    samples = MC_sampler.chain
    np.save(datacube_path.split('/')[-1] + '.ptchain', samples)

if __name__ == "__main__": fit_datacube(datacube_path)
