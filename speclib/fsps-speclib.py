#!/usr/bin/python

'''
  Make spectral libs for CALIFA V500/V1200 datacubes. (with fsps)
  YJ Qin, Jan 2016 @ Shanghai

  [!] Not working yet, just some tests.
'''

import fsps
import numpy as np
import matplotlib.pyplot as plt

'''
  Test: play with python-fsps
'''
if __name__ == "__main__":

  # test
  sp = fsps.StellarPopulation(sfh = 0, sf_start = 0.)

  # plot ssp at different ages
  #'''
  plt.xlim(np.log(91.), np.log(10000.))
  plt.ylim(-25., -10.)
  plt.title("SSP Spectra with FSPS (vel space)")
  plt.xlabel("ln lambda [A]"), plt.ylabel("log flux")

  # loop over ages
  for age_i in [1.e-3, 1.e-2, 1.e-1, 1.e0, 1.e1]:

    wl, flux = sp.get_spectrum(tage = age_i)
    plt.plot(np.log(wl), np.log10(flux), label = str(age_i) + " Gyr")

  plt.legend(loc = 'best')
  plt.show()
  #'''

  # wavelength range
  '''
  plt.plot(np.log10(wl))
  plt.xlabel('Index'), plt.ylabel("log(lambda) [A]")
  plt.title("Wavelength axis")
  plt.show()
  '''

  # resolution
  '''
  plt.plot(np.log10(wl[:-1]), np.log10(np.diff(wl)))
  plt.xlabel('log(lambda) [A]'), plt.ylabel("log(Delta lambda) [A]")
  plt.title("Resolution in the wavelength axis")
  plt.show()
  '''

  # interpolate spectrum in log(Z/Z_sun), Age space.
  import time
  t0 = time.clock()
  flux, mass, lbol = sp.ztinterp(0., 12., peraa = True)
  delta_t = time.clock() - t0
  print "takes", delta_t, "sec to interpolate."
  plt.plot(np.log10(sp.wavelengths), np.log10(flux)), plt.show()
