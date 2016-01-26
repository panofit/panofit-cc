#!/usr/bin/python

'''
  Models
'''

import numpy as np
import itertools as itt

def fsps_dssp_exp_sersic_cf(par, sp, ra_ax, dec_ax, id_a, id_b):

  # unpack param
  x_c, y_c, sersic_n, sersic_re, sersic_Ie, sersic_phi, sersic_q, sersic_c, sersic_age, \
      exp_Ic, exp_h, exp_phi, exp_q, exp_c, exp_age = par

  # calculate structural information
  ek1, ek2, ek3, ek4, ek5, ek6 = 3., -1. / 3., 8. / 1215., \
      184. / 229635., 1048. / 31000725., -17557576. / 1242974068875.
  d_n = ek1 * sersic_n + ek2 + ek3 / sersic_n + ek4 / (sersic_n ** 2) \
      + ek5 / (sersic_n ** 3) + ek6 / (sersic_n ** 4)

  # interpolate ssp
  flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(0., np.power(10., sersic_age), peraa = True)
  flux_s = flux_ssp[id_a: id_b]
  flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(0., np.power(10., exp_age), peraa = True)
  flux_e = flux_ssp[id_a: id_b]

  # mock cube, follow CALIFA style.
  ra_pts, dec_pts = ra_ax.size, dec_ax.size
  I_exp, I_sersic = np.zeros((ra_pts, dec_pts)), np.zeros((ra_pts, dec_pts))
  cube_t = np.zeros([id_b - id_a, dec_pts, ra_pts])

  # calculate pos
  dec_t, ra_t = np.meshgrid(dec_ax - y_c, ra_ax - x_c)
  r_t = np.sqrt(ra_t ** 2 + dec_t ** 2)
  cos_t, sin_t = ra_t / r_t, dec_t / r_t

  # mask central singularity
  cos_t[~np.isfinite(cos_t)] = 1.
  sin_t[~np.isfinite(sin_t)] = 0.

  # calculate t
  cos_s, sin_s = np.cos(sersic_phi), np.sin(sersic_phi)
  X_s = r_t * (cos_t * cos_s + sin_t * sin_s)
  Y_s = r_t * (sin_t * cos_s - cos_t * sin_s)
  M_s = np.power( np.power(np.abs(X_s), sersic_c)
                + np.power(np.abs(Y_s) / sersic_q, sersic_c), 1. / sersic_c)
  I_s = np.power(2.512, -sersic_Ie) * np.exp(-d_n * (np.power(M_s / sersic_re, 1. / sersic_n) - 1.))

  # calculate the exponential component
  cos_e, sin_e = np.cos(exp_phi), np.sin(exp_phi)
  X_e = r_t * (cos_t * cos_e + sin_t * sin_e)
  Y_e = r_t * (sin_t * cos_e - cos_t * sin_e)
  M_e = np.power( np.power(np.abs(X_e), exp_c)
                + np.power(np.abs(Y_e) / exp_q, exp_c), 1. / exp_c)
  I_e = np.power(2.512, -exp_Ic) * np.exp(-np.abs(M_e / exp_h))

  # iterate over the cube and co-add
  for I_ra, I_dec in itt.product(range(ra_pts), range(dec_pts)):
    cube_t[:, I_dec, I_ra] = flux_s * I_s[I_ra, I_dec] + flux_e * I_e[I_ra, I_dec]
  # TODO: maybe using ravel + outer can boost speed?

  return cube_t

def fsps_ssp_sersic_agrad_cf(par, sp, ra_ax, dec_ax, id_a, id_b):

  # unpack parameters to fit
  x_c, y_c, sersic_n, sersic_re, sersic_Ie, sersic_phi, \
      sersic_q, sersic_c, ssp_age_c, ssp_age_k = par

  # calculate d_n
  ek1, ek2, ek3, ek4, ek5, ek6 = 3., -1. / 3., 8. / 1215., \
      184. / 229635., 1048. / 31000725., -17557576. / 1242974068875.
  d_n = ek1 * sersic_n + ek2 + ek3 / sersic_n + ek4 / (sersic_n ** 2) \
      + ek5 / (sersic_n ** 3) + ek6 / (sersic_n ** 4)

  # construct arrays
  ra_pts, dec_pts = ra_ax.size, dec_ax.size
  I_sersic = np.zeros((ra_pts, dec_pts))
  cube_t = np.zeros([id_b - id_a + 1, dec_pts, ra_pts])

  # calculate pos
  dec_t, ra_t = np.meshgrid(dec_ax - y_c, ra_ax - x_c)
  r_t = np.sqrt(ra_t ** 2 + dec_t ** 2)
  cos_t, sin_t = ra_t / r_t, dec_t / r_t

  # mask central singularity
  cos_t[~np.isfinite(cos_t)] = 1.
  sin_t[~np.isfinite(sin_t)] = 0.

  # calculate surface brightness and local age
  cos_s, sin_s = np.cos(sersic_phi), np.sin(sersic_phi)
  X_s = r_t * (cos_t * cos_s + sin_t * sin_s)
  Y_s = r_t * (sin_t * cos_s - cos_t * sin_s)
  M_s = np.power( np.power(np.abs(X_s), sersic_c)
                + np.power(np.abs(Y_s) / sersic_q, sersic_c), 1. / sersic_c)
  I_s = np.power(2.512, -sersic_Ie) * np.exp(-d_n * (np.power(M_s / sersic_re, 1. / sersic_n) - 1.))
  A_s = ssp_age_k * M_s / sersic_re + ssp_age_c

  # interpolate the age gradient
  for I_ra, I_dec in itt.product(range(ra_pts), range(dec_pts)):
    flux_i, M_i, L_i = sp.ztinterp(0., np.power(10., A_s[I_ra, I_dec]), peraa = True)
    cube_t[:, I_dec, I_ra] = flux_i[id_a: id_b] * I_s[I_ra, I_dec]

  return cube_t
