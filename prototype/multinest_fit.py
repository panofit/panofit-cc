
import sys
import numpy as np
import pymultinest, fsps

from models     import *
from datacube   import *
from ssp_interp import *

# suppress divide-by-zero warnings
np.seterr(all = 'ignore')

# "PSF" settings, for debug
PSF_size = None
from scipy.ndimage.filters import gaussian_filter

# conceptual test of a specific model
def fsps_dssp_exp_sersic_cf_fit(filename, redshift, descr):

  # load datacube
  flux, err, mask, wl_ax, ra_ax, dec_ax = read_califa_datacube(filename)

  # mask zero-error points
  mask[np.log10(np.abs(flux / err)) > 2.4] = 0  # mask abnormal SNR
  mask[err > 1.e5] = 0 # mask zero-value pixels

  # initialize stellar population model
  sp = fsps.StellarPopulation(sfh = 0, sf_start = 0.)

  # shift to rest frame and rebin to the wl of SPS backend
  flux, err, mask, wl_new, id_a, id_b = \
    rebin_to_restframe(flux, err, mask, wl_ax, sp.wavelengths, source_z = redshift)

  # limit of parameters
  par_lim = [ (-5., 5.), (-5., 5.), # x_c, y_c
              (1.e-1, 8.), (1.e0, 3.5e1), (-16., 0.), # sersic n, Re, Ie
              (-np.pi / 2., np.pi / 2.), (2.e-1, 1.), (1., 3.), # sersic phi, q, c
              (-3, 1.1), # sersic age
              (-16., 0.), (1.e-1, 3.5e1), # exp Ic, Rh
              (-np.pi / 2., np.pi / 2.), (2.e-1, 1.), (1., 3.), # exp phi, q, c
              (-3, 1.1)] # exp age
  par_lim = np.array(par_lim)
  par_names = ["Xc", "Yc", "Sersic n", "Sersic Re", "Sersic Ie", "Sersic phi", "Sersic b/a",
               "Sersic c", "Sersic log(Age)", "Disk Ic", "Disk Rh", "Disk phi",
               "Disk b/a", "Disk c", "Disk log(Age)"]
  n_param = len(par_names)

  # extract parameters
  def _unpack_param(cube, ndim):
    return tuple(cube[i] for i in xrange(ndim))

  # the cube-mapping function. Note the "cube" here is not the datacube.
  def _prior(cube, ndim, nparams):

    # uniform mapping for all parameters, assuming independent
    for i in xrange(ndim):
      cube[i] = par_lim[i, 0] + (par_lim[i, 1] - par_lim[i, 0]) * cube[i]

  def _ln_like(cube, ndim, nparams):

    # unpack from lp_c_double to tuple of float n
    par = _unpack_param(cube, ndim)

    # make model, Note: argument "cube" is not a datacube
    cube_t = fsps_dssp_exp_sersic_cf(par, sp, ra_ax, dec_ax, id_a, id_b)

    # do "PSF" stuff
    assert PSF_size != None
    if PSF_size != 0.:
      for I_wl in xrange(cube_t.shape[0]):
        cube_t[I_wl, :, :] = gaussian_filter(cube_t[I_wl, :, :], PSF_size)

    # find chi-square and return likelihood
    _ln_like_t = -0.5 * np.nansum(((cube_t - flux) / err) ** 2)

    # DEBUG print
    print "% .6E"%(_ln_like_t,),
    for i_par in par: print "% .3E"%(i_par,),
    print " "

    return _ln_like_t

  # run sampling
  pymultinest.run(_ln_like, _prior, n_param, outputfiles_basename = filename + '-' + descr + '-', \
      resume = False, verbose = True, const_efficiency_mode = True, sampling_efficiency = 1., n_live_points = 1500)

if __name__ == "__main__":

  data_cube = "./NGC0001.V500.rscube.fits"
  source_z  = 0.015147
  PSF_size  = 2.
  fsps_dssp_exp_sersic_cf_fit(data_cube, source_z, "")
