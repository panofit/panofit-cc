#!/usr/bin/python

'''
  Visualize the result of fsps-dualssp-fit.py
  YJ Qin, Jan 2016 @ Shanghai

  Parameters: sersic_n, sersic_re, sersic_Ie,  sersic_phi,
              sersic_q, sersic_c,  sersic_age,
              exp_Ic,   exp_h,     exp_phi,    exp_q,
              exp_c,    exp_age,
'''

# specify the name of the data cube
datacube_path = "./califa_sample/NGC0001.V500.rscube.fits"

# best-fitting parameter
par_opt = (1.20189845, 5.99849347, -12.0130215, -0.31569187, 0.51323838,
           1.99219028, 1.09990431, -11.29316191, 12.61844848, -0.14629537,
           0.57732247, 2.31568646, 0.15574498)

# import everything
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin_tnc # bounded multivar minimization solver
from scipy.interpolate import interp1d

import fsps # python-fsps
from astropy.io import fits

import sys, gc
import itertools as itt

# TODO
