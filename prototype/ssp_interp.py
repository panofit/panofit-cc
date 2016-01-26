#!/usr/bin/python

import numpy as np

def ssp_interp(splib, age_ax, Z_ax, age, Z):

  '''
    Bilinear interpolation of SSP age/metallicity
  '''

  # number of points
  N_age, N_Z = age_ax.size, Z_ax.size

  # determine the weights
  age_ida, age_idb, age_w = None, np.searchsorted(age_ax, age), 1.
  if age_idb == 0: age_ida, age_idb, age_w = 0, 0, 1.
  elif age_idb == N_age: age_ida, age_idb, age_w = N_age, N_age, 1.
  else: age_ida, age_w = age_idb - 1, (age_ax[wl_idb] - age) / (age_ax[wl_idb] - age_ax[wl_idb - 1])

  Z_ida, Z_idb, Z_w = None, np.searchsorted(Z_ax, Z), 1.
  if Z_idb == 0: Z_ida, Z_idb, Z_w = 0, 0, 1.
  elif Z_idb == N_Z: Z_ida, Z_idb, Z_w = N_Z, N_Z, 1.
  else: Z_ida, Z_w = Z_idb - 1, (Z_ax[wl_idb] - Z) / (Z_ax[wl_idb] - Z_ax[wl_idb - 1])

  return splib[age_ida, Z_ida] * age_w * Z_w + splib[age_ida, Z_idb] * age_w * (1. - Z_w) \
       + splib[age_idb, Z_ida] * (1. - age_w) * Z_w + splib[age_idb, Z_idb] * (1. - age_w) * (1. - Z_w)

def ssp_age_interp(splib, age_ax, age):

  '''
    Linear interpolation of SSP age
  '''

  # determine the weights
  N_age = age_ax.size
  age_ida, age_idb, age_w = None, np.searchsorted(age_ax, age), 1.
  if age_idb == 0: age_ida, age_idb, age_w = 0, 0, 1.
  elif age_idb == N_age: age_ida, age_idb, age_w = N_age, N_age, 1.
  else: age_ida, age_w = age_idb - 1, (age_ax[wl_idb] - age) / (age_ax[wl_idb] - age_ax[wl_idb - 1])

  return splib[age_ida] * w_age + splib[age_idb] * (1. - w_age)
