#!/usr/bin/python

'''
  Generate a mock dual-ssp disk + bulge model and fit it.
  YJ Qin, Jan 2016 @ Shanghai

  Parameters to optimize: sersic_n,   sersic_re,  sersic_Ie,  sersic_phi, sersic_q
                          sersic_c,   sersic_age, exp_Ic,     exp_h,      exp_phi,
                          exp_q,      exp_c,      exp_age,
  with a bounded minimization solver.
'''

'''
  TNC solver, No noise:
    [ 2.74259793  3.03021036  1.4389679   0.10002821  0.59519267  1.83170749
      0.95330637  0.83443975  9.66147203  0.32600399  0.5365104   1.75574811
      0.32194648]
    Model: (2.75, 4., 1., 0.3, 0.75, 2.25, 1.0, 0.8, 10., 0.3, 0.5, 2.0, 0.3)
'''

'''
  TNC solver, No noise
    2.232633e+00 7.266384e+00 4.836127e-01 2.821076e-01 6.958582e-01 1.922107e+00
    1.021880e+00 2.006268e+00 7.515663e+00 8.528878e-01 4.988663e-01 1.971169e+00
    3.068533e-01
    Model: (2.75, 7.5, 0.3, 0.3, 0.75, 2.25, 1.0, 2., 7.5, 0.85, 0.5, 2.0, 0.3)
'''

# the "correct answer"
par_opt = (2.75, 7.5, 0.3, 0.3, 0.75, 2.25, 1.0, 2., 7.5, 0.85, 0.5, 2.0, 0.3)

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

# important stuff
N_ra, N_dec = 32, 32
ra_ax = np.linspace(-30., 30., N_ra)
dec_ax = np.linspace(-30., 30., N_dec)

# wavelength cut
wl_ida, wl_idb = 369, 557

# SNR (global, constant)
SNR = 20.

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
  cube_t = np.zeros([wl_idb - wl_ida + 1, dec_ax.size, ra_ax.size])

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
      I_s = np.power(10., sersic_Ie) * np.exp(-d_n * (np.power(M_s / sersic_re, 1. / sersic_n) - 1.))

      # calculate "surface brightness" of exponential component
      cos_pe, sin_pe = np.cos(exp_phi), np.sin(exp_phi)
      X_e = r_i * (cos_ri * cos_pe + sin_ri * sin_pe)
      Y_e = r_i * (sin_ri * cos_pe - cos_ri * sin_pe)
      M_e = np.power( np.power(np.abs(X_e), exp_c)
                    + np.power(np.abs(Y_e) / exp_q, exp_c), 1. / exp_c)
      I_e = np.power(10., exp_Ic) * np.exp(-np.abs(M_e / exp_h))

      # model spectrum
      cube_t[:, I_dec, I_ra] = flux_s * I_s + flux_e * I_e
      #cube_t[:, I_dec, I_ra] = np.log10(I_s + I_e)

  return cube_t

# object function to minimize
def _min_objfc(par, cube, sp, ra_ax, dec_ax, wl_ida, wl_idb):

  # unpack parameters to fit
  sersic_n, sersic_re, sersic_Ie, sersic_phi, sersic_q, sersic_c, sersic_age, \
      exp_Ic, exp_h, exp_phi, exp_q, exp_c, exp_age = par

  # calculate d_n for sersic component
  ek1, ek2, ek3, ek4, ek5, ek6 = 3., -1. / 3., 8. / 1215., \
      184. / 229635., 1048. / 31000725., -17557576. / 1242974068875.
  d_n = ek1 * sersic_n + ek2 + ek3 / sersic_n + ek4 / (sersic_n ** 2) \
      + ek5 / (sersic_n ** 3) + ek6 / (sersic_n ** 4)

  # array of single-spectrum chi_sq
  chi_sq = np.zeros((ra_ax.size, dec_ax.size), dtype = 'f8')

  # interpolate ssp
  flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(0., np.power(10., sersic_age), peraa = True)
  flux_s = flux_ssp[wl_ida: wl_idb + 1]

  flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(0., np.power(10., exp_age), peraa = True)
  flux_e = flux_ssp[wl_ida: wl_idb + 1]

  # loop over spectra
  for I_ra, v_ra in enumerate(ra_ax):
    for I_dec, v_dec in enumerate(dec_ax):

      # pos on the sky
      r_i = np.sqrt(v_ra ** 2 + v_dec ** 2)
      cos_ri, sin_ri = v_ra / r_i, v_dec / r_i
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
      flux_obs = cube[:, I_dec, I_ra]
      err_obs = 1#flux_obs / SNR

      # find chi_sq of this spectra
      chi_sq_i = ((flux_obs - flux_model) / err_obs) ** 2
      chi_sq_i[~np.isfinite(chi_sq_i)] = 0. # mask out nan and inf

      # write chi_sq of this spectrum
      chi_sq[I_ra, I_dec] = chi_sq_i.sum()

  # return a total chi_sq
  chi_sq_sum = np.sum(chi_sq)

  # print "  Target function gets chi_sq:", chi_sq_sum, "\n"
  print "%e"%chi_sq_sum,

  #print "Target function called with param:", par
  for ipar in par: print "%e"%ipar,
  print " "

  return chi_sq_sum

# log_L function, for MCMC sampling.
def _log_prob(par, cube, sp, ra_ax, dec_ax, wl_ida, wl_idb, bounds):

  # check bounds: within bounds, return log_L, else return -inf
  if np.all(np.logical_and(np.asarray(bounds)[:, 0] < np.asarray(par), \
                           np.asarray(par) < np.asarray(bounds)[:, 1])):
    return -0.5 * _min_objfc(par, cube, sp, ra_ax, dec_ax, wl_ida, wl_idb)
  else: return -np.inf

def fit_mock_datacube(cube, sp, ra_ax, dec_ax, wl_ida, wl_idb):

  # par_opt = (2.75, 4., 1., 0.3, 0.75, 2.25, 1.0, 0.8, 10., 0.3, 0.5, 2.0, 0.3)

  #init_guess = (2., 1., 1., 0., 0.5, 2., 1., 1., 1., 0., 0.5, 2., 0.)
  #init_guess = (2.5, 3.75, 1.25, 0.3, 0.7, 2., 0.975, 0.75, 8.5, 0.3125, 0.5125, 1.875, 0.3125)
  init_guess = (2.75, 7.5, 0.3, 0.3, 0.75, 2.25, 1.0, 2., 7.5, 0.85, 0.5, 2.0, 0.3)
  min_bounds = [(1.e-2, 8.), (1.e-2, 8.e1), (-16., 16.), (-np.pi, np.pi),
                (1.e-1, 1.), (1.e-1, 4.),   (-3, 1.1),
                (-16., 16.), (1.e-2, 3.e2), (-np.pi, np.pi), (1.e-1, 1.),
                (1.e-1, 4.), (-3, 1.1)]

  # add some noise?
  cube = cube + np.random.randn(cube.size).reshape(cube.shape) * np.mean(cube) * 1.e-2

  #solver = "TNC"
  solver = "MCMC"

  if solver == 'TNC':

    par_fit, N_eval, ret = fmin_tnc(_min_objfc, init_guess,
        args = (cube, sp, ra_ax, dec_ax, wl_ida, wl_idb), approx_grad = True,
        bounds = min_bounds, messages = 0)
    return par_fit

  if solver == 'MCMC':

    import emcee, time

    N_dim, N_walkers = 13, 32
    MC_sampler = emcee.EnsembleSampler(N_walkers, N_dim, _log_prob,
        args = (cube, sp, ra_ax, dec_ax, wl_ida, wl_idb, min_bounds))

    init_sol = np.outer(np.ones(N_walkers), np.array(init_guess)) \
             + 2.5e-3 * np.random.randn(N_walkers * N_dim).reshape((N_walkers, N_dim))

    t0 = time.clock()
    MC_sampler.run_mcmc(init_sol, 4500)
    t1 = time.clock()
    print "takes", t1 - t0, "sec to run."

    samples = MC_sampler.chain
    np.save("mc_sampling_more_noisy-4500.dat", samples)

if __name__ == "__main__":

  # make stellar pop
  sp = fsps.StellarPopulation(sfh = 0, sf_start = 0.)
  cube = generate_mock_datacube(par_opt, sp, ra_ax, dec_ax, wl_ida, wl_idb)

  print np.mean(cube)

  # fit it.
  #par_fit = fit_mock_datacube(cube, sp, ra_ax, dec_ax, wl_ida, wl_idb)

  #print par_opt
  #print par_fit

  # plt.imshow(np.rot90((cube[ 5, :, :])), interpolation = 'nearest'), plt.colorbar(), plt.show()
  # plt.imshow(np.rot90((cube[-5, :, :])), interpolation = 'nearest'), plt.colorbar(), plt.show()
