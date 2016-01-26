#!/usr/bin/python

'''
  Function to load datacubes, do rebinning and mask bad pixels
  YJ Qin, Jan 2016, Shanghai
'''

import itertools as itt
import numpy as np
from astropy.io import fits
from spec_utils import *

def read_califa_datacube(filename):

  '''
    Load CALIFA V500/V1200 datacubes
  '''

  # read spectrum from CALIFA cubes
  hlst = fits.open(filename)

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
  flux, err, mask = np.copy(hlst[0].data), np.copy(hlst[1].data), 1 - np.copy(hlst[3].data)
  # my definition, mask = 1 for a valid pixel.

  # close fits file.
  hlst.close()

  return flux, err, mask, wl_ax, ra_ax, dec_ax

def rebin_to_restframe(flux, err, mask, wl, wl_new, source_z = 0., mask_nan = True):

  '''
    Convert the datacube to restframe
  '''

  # shape of the original cube
  N_wl, N_dec, N_ra = flux.shape

  # cut the valid range of wavelength
  id_a, id_b = np.searchsorted(wl_new,
      [wl[0] / (1. + source_z), wl[-1] / (1. + source_z)])
  wl_cut = wl_new[id_a: id_b]

  # rebinned cubes
  flux_t = np.zeros((wl_cut.size, N_dec, N_ra)) + np.nan
  err_t  = np.zeros((wl_cut.size, N_dec, N_ra)) + np.nan
  mask_t = np.zeros((wl_cut.size, N_dec, N_ra), dtype = 'i4') + np.nan

  # iterate over spaxels and rebin
  for i_ra, i_dec in itt.product(np.arange(N_ra), np.arange(N_dec)):

    # if outside the hexagon, set nan.
    if np.sum(mask[:, i_dec, i_ra]) == 0: continue

    # extract from the cube
    flux_i, err_i, mask_i = \
        flux[:, i_dec, i_ra], err[:, i_dec, i_ra], mask[:, i_dec, i_ra]

    flux_i, err_i, mask_i = restframe(wl, flux_i, err_i, mask_i,
        wl_cut, source_z, mask_nan = False)

    # put into the new cube
    flux_t[:, i_dec, i_ra], err_t[:, i_dec, i_ra], \
        mask_t[:, i_dec, i_ra] = flux_i, err_i, mask_i

  # mask bad pixels as NaN if necessary
  if mask_nan:
    flux_t[mask_t == 0] = np.nan
    err_t[mask_t == 0]  = np.nan

  return flux_t, err_t, mask_t, wl_cut, id_a, id_b
