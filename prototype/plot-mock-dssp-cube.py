#!/usr/bin/python

'''
  Make figures of mock dual-ssp cubes.
  YJ Qin, Jan. 2016 @ Shanghai
'''
# import everything
import numpy as np
import matplotlib.pyplot as plt
import emcee

import fsps # python-fsps
from astropy.io import fits

import sys, gc
import itertools as itt

# important stuff
N_ra, N_dec = 78, 72
ra_ax  = np.linspace(-38., 39., 78)
dec_ax = np.linspace(-36., 35., 72)
wl_ida, wl_idb = 369, 557

# the function is identical to that in fsps-mock-dssp-specfit.py
# may not match the latest version. please check.
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

if __name__ == "__main__":

  '''
  Parameters to optimize: sersic_n,   sersic_re,  sersic_Ie,  sersic_phi, sersic_q
                          sersic_c,   sersic_age, exp_Ic,     exp_h,      exp_phi,
                          exp_q,      exp_c,      exp_age,
  '''

  par_opt = (2.0, # sersic n
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
  wl_shift = 0
  wl_ida += wl_shift; wl_idb += wl_shift
  '''

  # make a mock cube
  sp = fsps.StellarPopulation(sfh = 0, sf_start = 0., zred = 0.0)
  cube = generate_mock_datacube(par_opt, sp, ra_ax, dec_ax, wl_ida, wl_idb)
  wl_ax = sp.wavelengths[wl_ida: wl_idb + 1]

  # some Gaussian PSF?
  #'''
  from scipy.ndimage.filters import gaussian_filter
  cube_new = np.zeros(cube.shape)
  for I_wl in xrange(cube.shape[0]):
    cube_new[I_wl, :, :] = gaussian_filter(cube[I_wl, :, :], 2.)
    # cube_new[I_wl, :, :] = cube[I_wl, :, :]
  #'''

  # plot it.
  cube_sl = cube_new[:, 36, :]

  fig = plt.figure(figsize = (18., 4.))
  plt_rg = [wl_ax[0], wl_ax[-1], ra_ax[0], ra_ax[-1]]

  ax1 = fig.add_subplot(1, 1, 1)
  pl1 = ax1.imshow(np.rot90(np.log10(cube_sl)), interpolation = 'nearest',
      extent = plt_rg, aspect = 'auto', vmin = -3., vmax = 0.5)
  plt.colorbar(pl1, ax = ax1)

  plt.show()
