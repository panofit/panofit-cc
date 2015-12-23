#!/usr/bin/python

# this file provides parameters to generate mock spec library.
# most options are meaningless.

import numpy as np

'''
  Re-sample the model along z-axis.
'''
z_div, z_max = 3001, 125.
z_ax = np.linspace(-z_max, z_max, z_div)
dz = 2. * z_max / (z_div - 1.)

'''
  Re-sample the age distribution.
'''
age_div = 41      # 50 points in log interval
age_min, age_max = 1.e6, 1.e10    # 1 Myr to 10 Gyr
age_alpha = np.log(age_max / age_min) / (age_div - 1.)
age_ax = age_min * np.exp(age_alpha * np.arange(age_div))
age_ax_log = np.log10(age_ax)

'''
  Re-sample the metallicity distribution.
'''
M_div = 32
M_min, M_max = 0.05, 2.0
M_ax_log = np.linspace(np.log10(M_min), np.log10(M_max), M_div)
M_ax = np.power(10., M_ax_log)

'''
  Spec and phot library
'''
N_phot_bands = 5
phot_bands = ['SDSS u', 'SDSS g', 'SDSS r', 'SDSS i', 'SDSS z']

N_spec_px = 448
spec_wl = None
# It works.
