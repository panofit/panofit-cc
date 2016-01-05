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
#''' # 1st run, TNC
par_opt = (1.20189845, 5.99849347, -12.0130215, -0.31569187, 0.51323838,
           1.99219028, 1.09990431, -11.29316191, 12.61844848, -0.14629537,
           0.57732247, 2.31568646, 0.15574498)
#'''

''' # 2nd run, TNC
par_opt = (0.86428252, 2.7846034, -15.36964374, -2.22694922, 0.80212827,
           1.06718218, 1.01979302, -11.08551644, 9.79371104, -2.06361907,
           0.98261746, 2.42270186, 0.32113472)
#'''

''' # 3rd run, DIFF-EVO (large area of invalid solutions!)
par_opt = (0.01170914, 6.89123059, 10.92779896, 3.02595948, 0.68800141,
           0.12277318, -0.77716597, -10.09690454, 2.38735419, 1.95948294,
           0.31298667, 3.12340962, 1.1)
#'''

# import everything
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin_tnc # bounded multivar minimization solver
from scipy.interpolate import interp1d

import fsps # python-fsps
from astropy.io import fits

import sys, gc
import itertools as itt

# return a mock spectrum
def _mock_spec_cube(ra_ax, dec_ax, par, sp, wl_cut_idx):

  # unpack parameters
  sersic_n, sersic_re, sersic_Ie, sersic_phi, sersic_q, sersic_c, sersic_age, \
      exp_Ic, exp_h, exp_phi, exp_q, exp_c, exp_age = par
  id_a, id_b = wl_cut_idx

  # calculate structural information
  ek1, ek2, ek3, ek4, ek5, ek6 = 3., -1. / 3., 8. / 1215., \
      184. / 229635., 1048. / 31000725., -17557576. / 1242974068875.
  d_n = ek1 * sersic_n + ek2 + ek3 / sersic_n + ek4 / (sersic_n ** 2) \
      + ek5 / (sersic_n ** 3) + ek6 / (sersic_n ** 4)

  # interpolate ssp
  flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(0., np.power(10., sersic_age), peraa = True)
  flux_s = flux_ssp[id_a: id_b + 1] * 1.e15

  flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(0., np.power(10., exp_age), peraa = True)
  flux_e = flux_ssp[id_a: id_b + 1] * 1.e15

  # construct output datacube
  cube_t = np.zeros(shape = (id_b - id_a + 1, dec_ax.size, ra_ax.size), dtype = 'f8')

  # loop over spectra,
  for I_ra, v_ra in enumerate(ra_ax):

    print "\r", "at I_ra =", I_ra, "...",
    sys.stdout.flush()

    for I_dec, v_dec in enumerate(dec_ax):

      # calculate "surface brightness" for the disk and the bulge
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
      cube_t[:, I_dec, I_ra] = flux_model

  return cube_t

# main program.
if __name__ == "__main__":

  # read spectrum from CALIFA cubes
  hlst = fits.open(datacube_path)

  # read basic information
  N_ra, N_dec = hlst[0].header['NAXIS1'], hlst[0].header['NAXIS2']
  N_spx = hlst[0].header['NAXIS3']

  # make wavelength axis, in A
  wl_ax = hlst[0].header['CRVAL3'] + hlst[0].header['CDELT3'] * np.arange(N_spx)

  # RA and Dec axes, in asec, origin at the galaxy center
  ra_ax  = -hlst[0].header['CRPIX1'] + hlst[0].header['CDELT1'] * np.arange(N_ra)
  dec_ax = -hlst[0].header['CRPIX2'] + hlst[0].header['CDELT2'] * np.arange(N_dec)
  ra_c, dec_c = int(hlst[0].header['CRPIX1']), int(hlst[0].header['CRPIX2'])

  # data, err and mask
  flux, err, mask = hlst[0].data, hlst[1].data, hlst[2].data

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

  # generate the best-fitting model cube
  cube = _mock_spec_cube(ra_ax, dec_ax, par_opt, sp, (id_a, id_b))

  # = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  # make plots

  # mock long-slit
  slit_fit, slit_obs = cube[:, dec_c, :], flux_new[:, dec_c, :]
  plt_ext = [wl_new[0], wl_new[-1], ra_ax[0], ra_ax[-1]]

  fig = plt.figure()
  ax1, ax2 = fig.add_subplot(2, 1, 1), fig.add_subplot(2, 1, 2)
  ax1.imshow(np.rot90(slit_obs), interpolation = 'nearest', extent = plt_ext, aspect = 'auto')
  ax2.imshow(np.rot90(slit_fit), interpolation = 'nearest', extent = plt_ext, aspect = 'auto')
  ax2.set_xlabel("Wavelength [A]"), ax2.set_ylabel("Position [asec]")
  plt.show()

  # "images" at three wavelengths
  idx_t = np.searchsorted(wl_new, [4500, 5500, 7000])
  im450_obs, im550_obs, im700_obs = np.swapaxes(np.log10(flux_new), 1, 2)[idx_t, :, :]
  im450_fit, im550_fit, im700_fit = np.swapaxes(np.log10(cube), 1, 2)[idx_t, :, :]

  fig = plt.figure()
  plt_ext = [ra_ax[0], ra_ax[-1], dec_ax[0], dec_ax[-1]]

  ax1, ax4 = fig.add_subplot(2, 3, 1), fig.add_subplot(2, 3, 4)
  ax1.imshow(np.rot90(im450_obs), extent = plt_ext, interpolation = 'nearest')
  ax4.imshow(np.rot90(im450_fit), extent = plt_ext, interpolation = 'nearest')

  ax2, ax5 = fig.add_subplot(2, 3, 2), fig.add_subplot(2, 3, 5)
  ax2.imshow(np.rot90(im550_obs), extent = plt_ext, interpolation = 'nearest')
  ax5.imshow(np.rot90(im550_fit), extent = plt_ext, interpolation = 'nearest')

  ax3, ax6 = fig.add_subplot(2, 3, 3), fig.add_subplot(2, 3, 6)
  ax3.imshow(np.rot90(im700_obs), extent = plt_ext, interpolation = 'nearest')
  ax6.imshow(np.rot90(im700_fit), extent = plt_ext, interpolation = 'nearest')

  plt.show()
