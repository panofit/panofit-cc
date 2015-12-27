#!/usr/bin/python

'''
  extract CALIFA data cubes for fitting.
  YJ Qin, Dec 2015

  Format for output file:
    N_spec(i4), N_spx (i4), wl_ax (f8 * N_spx),
    [id(i4), X(f8), Y(f8), flux(f8 * N_spx), err(f8 * N_spx), mask(i4 * N_spx)] * N_spec
'''

import sys
import itertools as itt

import numpy as np
from astropy.io import fits

import matplotlib.pyplot as plt

def extract_califa_cube(fname, output_fname):

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

  # plot to check
  # plt.imshow(np.rot90(flux[:, :, int(hlst[0].header['CRPIX1'])]), interpolation = 'nearest', aspect = 'auto'), plt.show()

  # count valid pixels
  id_valid = []
  for I_ra, I_dec in itt.product(range(N_ra), range(N_dec)):
    if np.sum(flux[:, I_dec, I_ra]) != 0.: id_valid.append((I_ra, I_dec))

  '''
  print "N. of valid pixels", len(id_valid)
  sys.exit()
  #'''

  # create output file
  with open(output_fname, "wb") as fp:

    # spectrum index
    idx = np.array(0, dtype = 'i4')

    # write header
    np.asarray(len(id_valid), dtype = 'i4').tofile(fp)
    np.asarray(N_spx, dtype = 'i4').tofile(fp)

    # write wl_ax
    wl_ax.astype('f8').tofile(fp)

    # loop over valid pixels
    for I_ra, I_dec in id_valid:

      # write id
      idx.tofile(fp)

      # write X, Y
      np.asarray(ra_ax[I_ra], dtype = 'f8').tofile(fp)
      np.asarray(dec_ax[I_dec], dtype = 'f8').tofile(fp)

      # write flux, err, mask
      flux[:, I_dec, I_ra].astype('f8').tofile(fp)
      err[:, I_dec, I_ra].astype('f8').tofile(fp)
      mask[:, I_dec, I_ra].astype('i4').tofile(fp)

if __name__ == '__main__':

  extract_califa_cube("./califa_sample/NGC0001.V1200.rscube.fits", "./NGC0001/NGC0001.V1200.dat")
  extract_califa_cube("./califa_sample/NGC0001.V500.rscube.fits",  "./NGC0001/NGC0001.V500.dat")

# ioAlWyKDXfyV1MBSN4m4
