#!/usr/bin/python

'''
  Visualize the result of fsps-dualssp-fit.py
  YJ Qin, Jan 2016 @ Shanghai

  Parameters: sersic_n, sersic_re, sersic_Ie,  sersic_phi,
              sersic_q, sersic_c,  sersic_age,
              exp_Ic,   exp_h,     exp_phi,    exp_q,
              exp_c,    exp_age,
'''

# import everything
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import fsps # python-fsps
from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import interp1d

import corner

import sys
import itertools as itt

par_names = ["Sersic n", "Bulge Re", "Bulge Ie", "Bulge phi", "Bulge b/a",
             "Bulge c", "Bulge log(Age)", "Disk Ic", "Disk R_h", "Disk phi",
             "Disk b/a", "Disk c", "Disk log(Age)"]

def _read_spec_cube(datacube_path):

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
  flux, err, mask = np.copy(hlst[0].data), np.copy(hlst[1].data), np.copy(hlst[2].data)

  # close fits file.
  hlst.close()

  return flux, err, mask, wl_ax, ra_ax, dec_ax

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
  flux_s = flux_ssp[id_a: id_b + 1]

  flux_ssp, mass_ssp, lbol_ssp = sp.ztinterp(0., np.power(10., exp_age), peraa = True)
  flux_e = flux_ssp[id_a: id_b + 1]

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
      I_s = np.power(2.512, -sersic_Ie) * np.exp(-d_n * (np.power(M_s / sersic_re, 1. / sersic_n) - 1.))

      # calculate "surface brightness" of exponential component
      cos_pe, sin_pe = np.cos(exp_phi), np.sin(exp_phi)
      X_e = r_i * (cos_ri * cos_pe + sin_ri * sin_pe)
      Y_e = r_i * (sin_ri * cos_pe - cos_ri * sin_pe)
      M_e = np.power( np.power(np.abs(X_e), exp_c)
                    + np.power(np.abs(Y_e) / exp_q, exp_c), 1. / exp_c)
      I_e = np.power(2.512, -exp_Ic) * np.exp(-np.abs(M_e / exp_h))

      # model spectrum
      flux_model = flux_s * I_s + flux_e * I_e
      cube_t[:, I_dec, I_ra] = flux_model

  return cube_t

def _density_peak(chain):

  # not the best solution. just returns the mean value.
  return np.mean(chain, axis = 0)

# main program.
if __name__ == "__main__":

    # get parameters:
  datacube_path, chain_path = sys.argv[1],sys.argv[2]
  source_redshift, PSF_size = float(sys.argv[3]), float(sys.argv[4])
  basename = (datacube_path.split('/')[-1]).split('.')
  basename = '-'.join(basename)

  # load the datacube
  flux, err, mask, wl_ax, ra_ax, dec_ax = _read_spec_cube(datacube_path)

  # get center of the galaxy
  ra_ic, dec_ic = np.argmin(np.abs(ra_ax)), np.argmin(np.abs(dec_ax))

  # make stellar pop
  sp = fsps.StellarPopulation(sfh = 0, sf_start = 0.)

  # cut the valid range of wavelength
  id_a, id_b = np.searchsorted(sp.wavelengths,
      [wl_ax[0] / (1. + source_redshift), wl_ax[-1] / (1. + source_redshift)], side = 'left')
  wl_new = (sp.wavelengths)[id_a: id_b + 1]

  # load the chain and determine the best solution
  chain = (np.load(chain_path)[:, 500:, :]).reshape((-1, 13))
  chain = chain[chain[:,  4] <= 1., :]
  chain = chain[chain[:, 10] <= 1., :]
  par_opt = _density_peak(chain)
  print par_opt

  # generate the best-fitting model cube
  cube = _mock_spec_cube(ra_ax, dec_ax, par_opt, sp, (id_a, id_b))

  # Gaussian filter
  for I_wl in xrange(cube.shape[0]):
    cube[I_wl, :, :] = gaussian_filter(cube[I_wl, :, :], PSF_size)

  # DEBUG
  # sys.exit()

  # = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  # make plots

  # mock long-slit
  #'''
  slit_fit, slit_obs = np.log10(np.abs(cube))[:, dec_ic, :], np.log10(np.abs(flux))[:, dec_ic, :]
  pltrg_obs = [wl_ax[0], wl_ax[-1], ra_ax[0], ra_ax[-1]]
  pltrg_fit = [wl_new[0], wl_new[-1], ra_ax[0], ra_ax[-1]]

  # determine the range of plotting
  imgs = np.concatenate((slit_fit.ravel(), slit_obs.ravel()))
  imgs = imgs[np.isfinite(imgs)].ravel()

  #   create normalized cumulative histogam
  hist_t, ed_t = np.histogram(imgs, bins = 512, normed = True,
      range = (imgs.min() - 0.1, imgs.max() + 0.1))
  hist_t = np.cumsum(hist_t); hist_t /= hist_t.max()
  ed_t = 0.5 * (ed_t[:-1] + ed_t[1:])

  #   find proper min and max points
  plt_vlim_itpfc = interp1d(hist_t, ed_t, bounds_error = False)
  vcmin, vcmax = plt_vlim_itpfc(5.e-3), plt_vlim_itpfc(1. - 5.e-3)

  fig = plt.figure(figsize = (15., 6.))
  ax1, ax2 = fig.add_subplot(2, 1, 1), fig.add_subplot(2, 1, 2)
  pl1 = ax1.imshow(np.rot90(slit_obs), interpolation = 'nearest',
            extent = pltrg_obs, aspect = 'auto', cmap = "jet", vmin = vcmin, vmax = vcmax)
  pl2 = ax2.imshow(np.rot90(slit_fit), interpolation = 'nearest',
            extent = pltrg_fit, aspect = 'auto', cmap = "jet", vmin = vcmin, vmax = vcmax)
  ax2.set_xlabel("Wavelength [A]"), ax2.set_ylabel("Position [asec]")
  dv1 = make_axes_locatable(ax1); ca1 = dv1.append_axes("right", size = "1%", pad = 0.05)
  dv2 = make_axes_locatable(ax2); ca2 = dv2.append_axes("right", size = "1%", pad = 0.05)
  cb1 = plt.colorbar(pl1, cax = ca1); cb2 = plt.colorbar(pl2, cax = ca2)

  plt.savefig(basename + '-slit.pdf', bbox_inches = 'tight')
  plt.clf(); plt.close()
  #plt.show()
  #'''

  # "images" at three wavelengths
  wls = np.array([4500, 5500, 7000])
  idx_t, idx_w = np.searchsorted(wl_new, wls), np.searchsorted(wl_ax, wls / (1. + source_redshift))
  im450_obs, im550_obs, im700_obs = np.swapaxes(np.log10(np.abs(flux)), 1, 2)[idx_w, :, :]
  im450_fit, im550_fit, im700_fit = np.swapaxes(np.log10(np.abs(cube)), 1, 2)[idx_t, :, :]

  im450_fit[~np.isfinite(im450_obs)] = np.nan
  im550_fit[~np.isfinite(im550_obs)] = np.nan
  im700_fit[~np.isfinite(im700_obs)] = np.nan

  # determine the range of plotting
  imgs = np.hstack((im450_fit, im450_obs, im550_fit, im550_obs, im700_fit, im700_obs))
  imgs = imgs[np.isfinite(imgs)].ravel()

  #   create normalized cumulative histogam
  hist_t, ed_t = np.histogram(imgs, bins = 512, normed = True,
      range = (imgs.min() - 0.1, imgs.max() + 0.1))
  hist_t = np.cumsum(hist_t); hist_t /= hist_t.max()
  ed_t = 0.5 * (ed_t[:-1] + ed_t[1:])

  #   find proper min and max points
  plt_vlim_itpfc = interp1d(hist_t, ed_t, bounds_error = False)
  vcmin, vcmax = plt_vlim_itpfc(5.e-3), plt_vlim_itpfc(1. - 5.e-3)

  # DEBUG
  '''
  plt.hist(imgs, bins = 128, histtype = 'step', normed = True)
  plt.axvline(x = vcmin); plt.axvline(x = vcmax)
  plt.show()
  #sys.exit()
  '''

  fig = plt.figure(figsize = (15., 7.5))
  plt_ext = [ra_ax[0], ra_ax[-1], dec_ax[0], dec_ax[-1]]

  ax1, ax4 = fig.add_subplot(2, 3, 1), fig.add_subplot(2, 3, 4)
  ax1.set_title("4500 A"); ax1.set_ylabel("Dec [asec]")
  ax4.set_ylabel("Dec [asec]"); ax4.set_xlabel("RA [asec]")
  pl1 = ax1.imshow(np.rot90(im450_obs), extent = plt_ext,
      interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet')
  pl4 = ax4.imshow(np.rot90(im450_fit), extent = plt_ext,
      interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet')
  dv1 = make_axes_locatable(ax1); ca1 = dv1.append_axes("right", size = "3%", pad = 0.05)
  dv4 = make_axes_locatable(ax4); ca4 = dv4.append_axes("right", size = "3%", pad = 0.05)
  cb1 = plt.colorbar(pl1, cax = ca1); cb4 = plt.colorbar(pl4, cax = ca4)

  ax2, ax5 = fig.add_subplot(2, 3, 2), fig.add_subplot(2, 3, 5)
  ax2.set_title("5500 A"); ax5.set_xlabel("RA [asec]")
  pl2 = ax2.imshow(np.rot90(im550_obs), extent = plt_ext,
      interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet')
  pl5 = ax5.imshow(np.rot90(im550_fit), extent = plt_ext,
      interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet')
  dv2 = make_axes_locatable(ax2); ca2 = dv2.append_axes("right", size = "3%", pad = 0.05)
  dv5 = make_axes_locatable(ax5); ca5 = dv5.append_axes("right", size = "3%", pad = 0.05)
  cb2 = plt.colorbar(pl2, cax = ca2); cb5 = plt.colorbar(pl5, cax = ca5)

  ax3, ax6 = fig.add_subplot(2, 3, 3), fig.add_subplot(2, 3, 6)
  ax3.set_title("7000 A"); ax6.set_xlabel("RA [asec]")
  pl3 = ax3.imshow(np.rot90(im700_obs), extent = plt_ext,
      interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet')
  pl6 = ax6.imshow(np.rot90(im700_fit), extent = plt_ext,
      interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet')
  dv3 = make_axes_locatable(ax3); ca3 = dv3.append_axes("right", size = "3%", pad = 0.05)
  dv6 = make_axes_locatable(ax6); ca6 = dv6.append_axes("right", size = "3%", pad = 0.05)
  cb3 = plt.colorbar(pl3, cax = ca3); cb6 = plt.colorbar(pl6, cax = ca6)

  plt.savefig(basename + '-imgs.pdf', bbox_inches = 'tight')
  #plt.show()
  plt.clf(), plt.close()

  fig = corner.corner(chain, labels = par_names, truths = par_opt)
  fig.savefig(basename + '-prob.pdf', bbox_inches = 'tight')
