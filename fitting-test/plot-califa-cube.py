#!/usr/bin/python

'''
  Make figures of CALIFA V500/V1200 cubes.
  YJ Qin, Jan. 2016 @ Shanghai
'''

import sys
import itertools as itt

import numpy as np
from astropy.io import fits

import matplotlib.pyplot as plt

def load_califa_cube(fname):

  # load fits cube
  hlst = fits.open(fname)

  print hlst.info()

  # see info
  N_ra, N_dec = hlst[0].header['NAXIS1'], hlst[0].header['NAXIS2']
  N_spx = hlst[0].header['NAXIS3']

  # wavelength axis, in A
  wl_ax = hlst[0].header['CRVAL3'] + hlst[0].header['CDELT3'] * np.arange(N_spx)

  # RA and Dec axes, in asec
  ra_ax = -hlst[0].header['CRPIX1'] + hlst[0].header['CDELT1'] * np.arange(N_ra)
  dec_ax = -hlst[0].header['CRPIX2'] + hlst[0].header['CDELT2'] * np.arange(N_dec)

  # data, err and mask
  flux, err, mask = hlst[0].data, hlst[1].data, hlst[2].data

  # mask bad pixels
  #flux[mask == 0] = np.nan
  err[mask == 0]  = np.nan

  # mask zeros
  #flux[err > 1.e8] = np.nan
  err[err > 1.e8]  = np.nan

  # return values
  return flux, err, ra_ax, dec_ax, wl_ax

if __name__ == "__main__":

  flux, err, ra_ax, dec_ax, wl_ax = load_califa_cube("./califa_sample/NGC0001.V500.rscube.fits")

  # cut slice
  id_c = np.searchsorted(dec_ax, 0.)
  flux_sl, err_sl = flux[:, id_c, :], err[:, id_c, :]

  # plot
  #'''
  fig = plt.figure()
  plt_rg = [wl_ax[0], wl_ax[-1], ra_ax[0], ra_ax[-1]]

  ax1 = fig.add_subplot(3, 1, 1)
  pl1 = ax1.imshow(np.rot90(flux_sl), interpolation = 'nearest', extent = plt_rg, aspect = 'auto')
  plt.colorbar(pl1, ax = ax1)

  ax2 = fig.add_subplot(3, 1, 2)
  pl2 = ax2.imshow(np.rot90(err_sl), interpolation = 'nearest', extent = plt_rg, aspect = 'auto')
  plt.colorbar(pl2, ax = ax2)

  ax3 = fig.add_subplot(3, 1, 3)
  pl3 = ax3.imshow(np.rot90(np.abs(flux_sl / err_sl)), interpolation = 'nearest',
             extent = plt_rg, aspect = 'auto', vmin = np.power(10., -1.), vmax = np.power(10., 2.2))
  plt.colorbar(pl3, ax = ax3)

  plt.show()
  #'''

  # plot the central pixel
  '''
  idc_ra  = np.searchsorted(ra_ax,  0.)
  idc_dec = np.searchsorted(dec_ax, 0.)

  flux_t, err_t = flux[:, idc_dec, idc_ra], err[:, idc_dec, idc_ra]
  plt.xlim(wl_ax[0], wl_ax[-1])
  plt.plot(wl_ax, flux_t); plt.plot(wl_ax, err_t)
  plt.show()
  '''

  '''
  for I_d in np.arange(dec_ax.size - idc_dec):

    flux_t, err_t = flux[:, idc_dec + I_d, idc_ra], err[:, idc_dec + I_d, idc_ra]
    plt.xlim(wl_ax[0], wl_ax[-1])
    plt.plot(wl_ax, flux_t); plt.plot(wl_ax, err_t)
    plt.show()
  '''

  ''' # SNR plot
  snr_t = np.abs(flux_sl / err_sl).ravel()
  snr_t[~np.isfinite(snr_t)] = 1.
  plt.hist(np.log10(snr_t), bins = 128, histtype = 'step')
  plt.show()
  #'''
