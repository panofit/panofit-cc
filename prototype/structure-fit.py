#!/usr/bin/python

'''
  Generate a bulge-disk model and fit it. (without color/stellar pop info)
  YJ Qin, Jan 2016 @ Shanghai

  Parameters to optimize: sersic_n,   sersic_re,  sersic_Ie,  sersic_phi, sersic_q,   sersic_c,
                          exp_Ic,     exp_h,      exp_phi,    exp_q,      exp_c
  with a bounded minimization solver.
'''

import numpy as np
import matplotlib.pyplot as plt
import emcee, time

from scipy.optimize import fmin_tnc

# global variables

# Trial setup 1.
# N_x, N_y, box_size = 256, 256, 3. * 60.

N_x, N_y, box_size = 553, 553, 30. # mock-dssp-fit
X_ax = np.linspace(-box_size / 2., box_size / 2., N_x)
Y_ax = np.linspace(-box_size / 2., box_size / 2., N_y)

# degree to radians
deg2rad = lambda x: np.pi * x / 180.
flux2mag = lambda arr: np.log(arr) / np.log(2.512)

def _make_image(par, ax_x, ax_y):

  # unpack parameters
  sersic_n, sersic_re, sersic_Ie, sersic_phi, sersic_q, sersic_c, \
      exp_Ic, exp_h, exp_phi, exp_q, exp_c = par

  # calculate structural information
  ek1, ek2, ek3, ek4, ek5, ek6 = 3., -1. / 3., 8. / 1215., \
      184. / 229635., 1048. / 31000725., -17557576. / 1242974068875.
  d_n = ek1 * sersic_n + ek2 + ek3 / sersic_n + ek4 / (sersic_n ** 2) \
      + ek5 / (sersic_n ** 3) + ek6 / (sersic_n ** 4)

  # make image array
  Y_t, X_t = np.meshgrid(ax_y, ax_x)
  X_pts, Y_pts = ax_x.size, ax_y.size
  Img = np.zeros_like(X_t)

  # calculate pos
  r_t = np.sqrt(X_t ** 2 + Y_t ** 2)
  cos_t, sin_t = X_t / r_t, Y_t / r_t

  # mask central singularity
  cos_t[~np.isfinite(cos_t)] = 1.
  sin_t[~np.isfinite(sin_t)] = 0.

  # calculate sersic compnent
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

  # DBG: print bulge-to-total ratio
  # print np.sum(I_s) / np.sum(I_e + I_s)

  # return image
  return I_s + I_e

def _imfit_fc(par, img, err, ax_x, ax_y):

  # generate mock image
  img_mock = _make_image(par, ax_x, ax_y)

  # calculate and return chi^2
  chi_sq = np.nanmean(((img_mock - img) / err) ** 2)

  # print something
  #'''
  print "% e |"%(chi_sq,),
  for vpar in par: print "% e "%(vpar,),
  print " "
  #'''

  return chi_sq

def _imfit_lnprob(par, img, err, ax_x, ax_y, bounds):

  # check bounds
  if np.all(np.logical_and(np.asarray(bounds)[:, 0] < np.asarray(par), \
                           np.asarray(par) < np.asarray(bounds)[:, 1])):
    return -0.5 * _imfit_fc(par, img, err, ax_x, ax_y)
  else: return -np.inf

def fit_image(img, err, ax_x, ax_y, out_fname, solver = "mcmc"):

  # set bounds of fitting
  ''' # Test case: 1
  par_init = (2., 10, 18., 0.5, 0.75, 2., 18., 35., 0.5, 0.45, 2.0)
  par_bds  = ((0., 4.), (1.e-2, 60.), (10., 24.), (0., np.pi), (1.e-3, 1.), (1., 3.),
              (10., 24.), (1.e-2, 60.), (0., np.pi), (1.e-3, 1.), (1., 3.))
  '''
  par_init = (2.75, 7.5, -2.25, 0.3, 0.75, 2.25, -4., 7.5, 0.85, 0.5, 2.0)
  par_bds  = ((0., 4.), (1.e-1, 25.), (-5., 0.), (-np.pi / 2., np.pi / 2.), (1.e-3, 1.), (1., 3.),
              (-5., 0.), (1.e-1, 25.), (-np.pi / 2., np.pi / 2.), (1.e-3, 1.), (1., 3.))

  if solver == "mcmc":

    # set mcmc
    N_dim, N_walkers = len(par_init), len(par_init) * 2
    MC_sampler = emcee.EnsembleSampler(N_walkers, N_dim, _imfit_lnprob,
        args = (img, err, ax_x, ax_y, par_bds), threads = 4)

    init_sol = np.outer(np.ones(N_walkers), np.array(par_init)) \
             + 2.5e-2 * np.random.randn(N_walkers * N_dim).reshape((N_walkers, N_dim))

    # print _imfit_lnprob(par_init, img, ax_x, ax_y, par_bds)

    t0 = time.clock()
    MC_sampler.run_mcmc(init_sol, 5500)
    t1 = time.clock()
    print "takes", t1 - t0, "sec to run."

    samples = MC_sampler.chain
    np.save(out_fname, samples)

  elif solver == "tnc":

    # set some arbitary init condition
    par_init = np.mean(np.asarray(par_bds), axis = 1)

    par_fit, N_eval, ret = fmin_tnc(_imfit_fc, par_init,
        args = (img, err, ax_x, ax_y), approx_grad = True,
        bounds = par_bds, messages = 0)
    np.save(out_fname, par_fit)
    return par_fit

  else: raise RuntimeError("Invalid solver.")

if __name__ == "__main__":

  import sys

  # model parameter (modified from NGC 936, Table 2., Kim+2014a)
  par_id = "sersic_n sersic_re sersic_Ie sersic_phi sersic_q sersic_c \
            exp_Ic exp_h exp_phi exp_q exp_c"
  # TEST CASE 1
  # par_opt = (2., 10, 18., deg2rad(30.), 0.75, 2., 18., 35., deg2rad(30.), 0.45, 2.0)
  par_opt = (2.75, 7.5, -2.25, 0.3, 0.75, 2.25, -4., 7.5, 0.85, 0.5, 2.0)

  # generate mock image
  img = _make_image(par_opt, X_ax, Y_ax)
  #err = img / 50. + img.min() # SNR Model 1: nearly constant, disk highlighted.
  err = np.sqrt(np.round(img / img.min())) * img.min()

  # plot mock image
  '''
  img_mag = flux2mag(img); img_mag -= img_mag.max()
  plt_ext = [X_ax[0] / 60., X_ax[-1] / 60., Y_ax[0] / 60., Y_ax[-1] / 60.]
  plt.contourf(np.rot90(img_mag), 64, extent = plt_ext, cmap = "gray_r", origin = "upper")
  #plt.imshow(np.rot90(img_mag), extent = plt_ext, cmap = "gray_r")
  plt.colorbar()
  plt.show()

  plt.contourf(np.rot90(np.log10(img / err)), 64, extent = plt_ext, cmap = "gray_r", origin = "upper")
  plt.colorbar()
  plt.show()

  sys.exit()
  #'''

  # fit mock image
  fit_image(img, err, X_ax, Y_ax, "structure-fit-test-dssp", solver = "tnc")

  # plot results (MCMC solver)
  '''
  chain = np.load("structure-fit-test-dssp.npy")
  #chain = chain.reshape((-1, len(par_opt)))
  chain = (chain[:, 500:, :]).reshape((-1, len(par_opt)))

  import corner
  fig = corner.corner(chain, truths = par_opt, labels = par_id.split())
  plt.savefig("structure-fit-test-553-mockdssp.npy.pdf")
  '''
