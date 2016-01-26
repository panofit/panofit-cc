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
  flux, err, mask = hlst[0].data, hlst[1].data, hlst[3].data

  # mask bad pixels
  # flux[mask == 0] = np.nan
  # err[mask == 0]  = np.nan

  # mask zeros
  # flux[err > 1.e8] = np.nan
  # err[err > 1.e8]  = np.nan

  # return values
  return flux, err, mask, ra_ax, dec_ax, wl_ax

if __name__ == "__main__":

  # flux, err, ra_ax, dec_ax, wl_ax = load_califa_cube("./califa_sample/NGC0001.V500.rscube.fits")
  flux, err, mask, ra_ax, dec_ax, wl_ax = load_califa_cube("/home/qinyj/workspace/panofit/califa_sample/NGC4185.V500.rscube.fits")

  ''' 2916
  d_ra, d_dec, d_r = -4.5, 11.5, 3.
  for i_ra, i_dec in itt.product(range(ra_ax.size), range(dec_ax.size)):
    if np.sqrt((ra_ax[i_ra] - d_ra) ** 2 + (dec_ax[i_dec] - d_dec) ** 2) < d_r:
      flux[:, i_dec, i_ra] = np.nan
  #'''

  ''' # 4185
  d_ra, d_dec, d_r = -12., -1., 3.5
  for i_ra, i_dec in itt.product(range(ra_ax.size), range(dec_ax.size)):
    if np.sqrt((ra_ax[i_ra] - d_ra) ** 2 + (dec_ax[i_dec] - d_dec) ** 2) < d_r: #pass
      flux[:, i_dec, i_ra] = np.nan
  d_ra, d_dec, d_r = 3., 30., 2.5
  for i_ra, i_dec in itt.product(range(ra_ax.size), range(dec_ax.size)):
    if np.sqrt((ra_ax[i_ra] - d_ra) ** 2 + (dec_ax[i_dec] - d_dec) ** 2) < d_r: #pass
      flux[:, i_dec, i_ra] = np.nan
  #'''

  '''
  plt.imshow(np.rot90(flux[:, 36, :]), interpolation = 'nearest', aspect = 'auto')
  plt.colorbar(); plt.show()
  plt.imshow(np.rot90(err[:, 36, :]), interpolation = 'nearest', aspect = 'auto')
  plt.colorbar(); plt.show()
  plt.imshow(np.rot90(mask[:, 36, :]), interpolation = 'nearest', aspect = 'auto')
  plt.colorbar(); plt.show()
  '''

  # print ra/dec axes
  '''
  print ra_ax
  print dec_ax
  #'''

  # cut slice
  id_c = np.searchsorted(dec_ax, 0.)
  flux_sl, err_sl = flux[:, id_c, :], err[:, id_c, :]
  print id_c

  # plot single-panel, log
  '''
  fig = plt.figure(figsize = (18., 4.))
  plt_rg = [wl_ax[0], wl_ax[-1], ra_ax[0], ra_ax[-1]]

  flux_sl = np.abs(flux_sl)
  ax1 = fig.add_subplot(1, 1, 1)
  pl1 = ax1.imshow(np.rot90(np.log10(flux_sl)), interpolation = 'nearest',
      extent = plt_rg, aspect = 'auto', vmin = -3., vmax = 0.5)
  plt.colorbar(pl1, ax = ax1)

  plt.show()
  #'''

  # plot three-panel
  '''
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

  # plot flux sum
  #'''
  fig = plt.figure(figsize = (12., 12.))
  plt_rg = [ra_ax[0], ra_ax[-1], dec_ax[0], dec_ax[-1]]

  ax1 = fig.add_subplot(1, 1, 1)
  pl1 = ax1.imshow(np.abs(np.rot90(np.swapaxes(flux.sum(axis = 0), 0, 1))), interpolation = 'nearest', extent = plt_rg, aspect = 'equal')
  plt.colorbar(pl1, ax = ax1)

  plt.show()
  #'''

  # plot three images
  '''
  fig = plt.figure(figsize = (18., 5.))

  id1, id2, id3 = np.searchsorted(wl_ax, [4500., 5500., 7000.])

  ax1 = fig.add_subplot(1, 3, 1)
  pl1 = ax1.imshow(np.rot90(np.swapaxes(flux[id1, :, :], 0, 1)), interpolation = 'nearest',
      extent = [ra_ax[0], ra_ax[-1], dec_ax[0], dec_ax[-1]], aspect = 'equal')
  plt.colorbar(pl1, ax = ax1)

  ax2 = fig.add_subplot(1, 3, 2)
  pl2 = ax2.imshow(np.rot90(np.swapaxes(flux[id2, :, :], 0, 1)), interpolation = 'nearest',
      extent = [ra_ax[0], ra_ax[-1], dec_ax[0], dec_ax[-1]], aspect = 'equal')
  plt.colorbar(pl2, ax = ax2)

  ax3 = fig.add_subplot(1, 3, 3)
  pl3 = ax3.imshow(np.rot90(np.swapaxes(flux[id3, :, :], 0, 1)), interpolation = 'nearest',
      extent = [ra_ax[0], ra_ax[-1], dec_ax[0], dec_ax[-1]], aspect = 'equal')
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
  #'''

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

  # test de-redshift // it works.
  '''
  wl_ida, wl_idb = 369, 557
  idc_ra  = np.searchsorted(ra_ax,  0.)
  idc_dec = np.searchsorted(dec_ax, 0.)
  flux_t, err_t = flux[:, idc_dec, idc_ra], err[:, idc_dec, idc_ra]

  import fsps
  from spec_utils import *
  sp = fsps.StellarPopulation(sfh = 0, sf_start = 0., zred = 0.1)
  flux_w, m_w, lbol_w = sp.ztinterp(0., 10., peraa = True)
  flux_w = flux_w[wl_ida: wl_idb + 1]
  wl_w = sp.wavelengths[wl_ida: wl_idb + 1]

  # de-redshift
  flux_dr, err_dr = restframe(wl_ax, flux_t, err_t, \
      np.ones(wl_ax.size), wl_w, 0.01463, mask_nan = True)

  # plt.xlim(wl_ax[0], wl_ax[-1])
  # plt.plot(wl_ax, flux_t)
  plt.plot(wl_w, flux_dr / np.nanmean(flux_dr))
  plt.plot(wl_w, flux_w / np.nanmean(flux_w))
  plt.show()
  #'''
