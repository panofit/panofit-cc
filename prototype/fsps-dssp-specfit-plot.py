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
from scipy.optimize import fmin, minimize

from spec_utils import *

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
  cube_sersic = np.zeros(shape = (id_b - id_a + 1, dec_ax.size, ra_ax.size), dtype = 'f8')
  cube_exp = np.zeros(shape = (id_b - id_a + 1, dec_ax.size, ra_ax.size), dtype = 'f8')

  # loop over spectra,
  for I_ra, v_ra in enumerate(ra_ax):
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
      cube_exp[:, I_dec, I_ra] = flux_e * I_e
      cube_sersic[:, I_dec, I_ra] = flux_s * I_s

  return cube_sersic, cube_exp

def _density_peak_min_tgfc(par_t, ch_t, par_std):

  dist = np.zeros_like(ch_t)
  for i_dim in range(ch_t.shape[1]):
    if i_dim == 3 or i_dim == 9: # two dims with periodic
      dist[:, i_dim] = np.fmod(ch_t[:, i_dim] - par_t[i_dim], (np.pi / 2.) / par_std[i_dim])
    else: dist[:, i_dim] = ch_t[:, i_dim] - par_t[i_dim]

  return np.exp(-dist ** 2 / 0.001).mean()

def _density_peak(chain):

  return np.median(chain, axis = 0)

  # not the best solution. just returns the mean value.
  par_mean, par_std = np.mean(chain, axis = 0), np.std(chain, axis = 0)
  chain_scaled = chain / np.outer(np.ones(chain.shape[0]), par_std)
  mean_scaled = par_mean / par_std

  # find bottom using optimization
  ''' # Fixed: using periodic boundary for sersic_phi and exp_phi (Jan 18, 2016)
  fmin_tgfc = lambda par_t, ch_t: np.exp(-(ch_t - np.outer(np.ones(ch_t.shape[0]), par_t)) ** 2 / 0.001).mean()
  fmin_res = minimize(fmin_tgfc, mean_scaled, args = (chain_scaled,),)
  '''

  # find bottom using optimization
  fmin_res = minimize(_density_peak_min_tgfc, mean_scaled, args = (chain_scaled, par_std))

  return fmin_res.x * par_std

# main program.
if __name__ == "__main__":

  # get parameters:
  datacube_path, chain_path = sys.argv[1],sys.argv[2]
  source_redshift, PSF_size = float(sys.argv[3]), float(sys.argv[4])
  basename = (chain_path.split('/')[-1]).split('.')
  basename = '-'.join(basename)

  # load the datacube
  print "Loading datacube..."
  flux, err, mask, wl_ax, ra_ax, dec_ax = _read_spec_cube(datacube_path)

  # get center of the galaxy
  ra_ic, dec_ic = np.argmin(np.abs(ra_ax)), np.argmin(np.abs(dec_ax))

  # make stellar pop
  print "Initializing stellar population..."
  sp = fsps.StellarPopulation(sfh = 0, sf_start = 0.)

  # cut the valid range of wavelength
  id_a, id_b = np.searchsorted(sp.wavelengths,
      [wl_ax[0] / (1. + source_redshift), wl_ax[-1] / (1. + source_redshift)], side = 'left')
  wl_new = (sp.wavelengths)[id_a: id_b + 1]

  # redshift correction and rebinning of the original cube
  print "Rebinning..."
  flux_rf = np.zeros((wl_new.size, dec_ax.size, ra_ax.size))
  for i_ra, i_dec in itt.product(np.arange(ra_ax.size), np.arange(dec_ax.size)):
    flux_i, err_i = restframe(wl_ax, flux[:, i_dec, i_ra],
        err[:, i_dec, i_ra], mask[:, i_dec, i_ra], wl_new, source_redshift, mask_nan = True)
    flux_rf[:, i_dec, i_ra] = flux_i

  # load the chain and determine the best solution
  '''
  print "Finding optimal solution..."
  chain = np.load(chain_path)
  if chain.ndim == 3: # normal ensemble sampler
    chain = (chain[:, 500:, :]).reshape((-1, 13))
  elif chain.ndim == 4: # parallel tempering
    chain = (chain[:, :, 500:, :]).reshape((-1, 13))
  else: raise RuntimeError("Wrong MCMC chain dimension.")
  chain = chain[chain[:,  4] <= 1., :]
  chain = chain[chain[:, 10] <= 1., :]
  par_opt = _density_peak(chain)
  print par_opt
  #''' # DEBUG
  '''
  par_opt = (0.86925408, 7.89477418, -8.82317326, -0.4356502, 0.75141268,
             1.93851449, 0.55313915, -9.90119334, 12.76594508, -0.20848418,
             0.88245159, 2.15420082, 0.59015511)
  #'''

  param = np.fromfile(chain_path, dtype = 'f8').reshape((15, 4)).T
  if np.isfinite(np.sum(param[3, :])): par_opt = param[3, :]
  elif np.isfinite(np.sum(param[0, :])): par_opt = param[0, :]
  elif np.isfinite(np.sum(param[2, :])): par_opt = param[2, :]
  else: print "Bad param file."; sys.exit()

  # generate the best-fitting model cube
  print "Making best-fitting model..."
  cube_sersic, cube_exp = _mock_spec_cube(ra_ax, dec_ax, par_opt, sp, (id_a, id_b))

  # Gaussian filter
  print "Applying PSF..."
  for I_wl in xrange(cube_sersic.shape[0]):
    cube_sersic[I_wl, :, :] = gaussian_filter(cube_sersic[I_wl, :, :], PSF_size)
    cube_exp[I_wl, :, :] = gaussian_filter(cube_exp[I_wl, :, :], PSF_size)

  # total flux from two components
  cube = cube_sersic + cube_exp

  # = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  # make plots
  print "Start making plots."

  # two-panel long-slit plot
  if True:

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
              extent = pltrg_obs, aspect = 'auto', cmap = "RdYlBu_r", vmin = vcmin, vmax = vcmax)
    pl2 = ax2.imshow(np.rot90(slit_fit), interpolation = 'nearest',
              extent = pltrg_fit, aspect = 'auto', cmap = "RdYlBu_r", vmin = vcmin, vmax = vcmax)
    ax1.set_ylabel("Observed"), ax2.set_ylabel("Best-fitting"); ax2.set_xlabel('Wavelength [A]')
    dv1 = make_axes_locatable(ax1); ca1 = dv1.append_axes("right", size = "1%", pad = 0.05)
    dv2 = make_axes_locatable(ax2); ca2 = dv2.append_axes("right", size = "1%", pad = 0.05)
    cb1 = plt.colorbar(pl1, cax = ca1); cb2 = plt.colorbar(pl2, cax = ca2)
    cb1.set_label("log flux"); cb2.set_label("log flux")

    plt.savefig(basename + '-slit.pdf', bbox_inches = 'tight')
    plt.clf(); plt.close()

  # five-panel long-slit plot
  if True:

    slit_raw = np.log10(np.abs(flux))[:, dec_ic, :] # observed, unbinned original
    slit_bin = np.log10(np.abs(flux_rf))[:, dec_ic, :] # observed, binned
    slit_fit = np.log10(np.abs(cube))[:, dec_ic, :] # best-fitting
    slit_sersic = np.log10(np.abs(cube_sersic))[:, dec_ic, :] # best-fitting, sersic
    slit_exp = np.log10(np.abs(cube_exp))[:, dec_ic, :] # best-fitting, exponential
    slit_delta = flux_rf[:, dec_ic, :] - cube[:, dec_ic, :] # difference
    pltrg_obs = [wl_ax[0], wl_ax[-1], ra_ax[0], ra_ax[-1]]
    pltrg_fit = [wl_new[0], wl_new[-1], ra_ax[0], ra_ax[-1]]

    #   determine the range of original/best-fitting plots
    imgs = np.concatenate((slit_fit.ravel(), slit_raw.ravel()))
    imgs = imgs[np.isfinite(imgs)].ravel()
    hist_t, ed_t = np.histogram(imgs, bins = 512, normed = True,
        range = (imgs.min() - 0.1, imgs.max() + 0.1))
    hist_t = np.cumsum(hist_t); hist_t /= hist_t.max()
    ed_t = 0.5 * (ed_t[:-1] + ed_t[1:])
    plt_vlim_itpfc = interp1d(hist_t, ed_t, bounds_error = False)
    vcmin, vcmax = plt_vlim_itpfc(5.e-3), plt_vlim_itpfc(1. - 5.e-3)

    #   determine the color range of diff plot
    vdmin, vdmax = -np.abs(slit_delta.max()), np.abs(slit_delta.max())

    #   make figure
    fig = plt.figure(figsize = (15., 15.))

    #   plot unbinned original data and the binned
    ax1, ax2 = fig.add_subplot(6, 1, 1), fig.add_subplot(6, 1, 2)
    pl1 = ax1.imshow(np.rot90(slit_raw), interpolation = 'nearest',
              extent = pltrg_obs, aspect = 'auto', cmap = 'RdYlBu_r', vmin = vcmin, vmax = vcmax)
    pl2 = ax2.imshow(np.rot90(slit_bin), interpolation = 'nearest',
              extent = pltrg_fit, aspect = 'auto', cmap = 'RdYlBu_r', vmin = vcmin, vmax = vcmax)
    dv1 = make_axes_locatable(ax1); ca1 = dv1.append_axes("right", size = "1%", pad = 0.05)
    dv2 = make_axes_locatable(ax2); ca2 = dv2.append_axes("right", size = "1%", pad = 0.05)
    cb1 = plt.colorbar(pl1, cax = ca1); cb2 = plt.colorbar(pl2, cax = ca2)

    #   plot best fitting model, sersic and exponential components
    ax3, ax4, ax5 = fig.add_subplot(6, 1, 3), fig.add_subplot(6, 1, 4), fig.add_subplot(6, 1, 5)
    pl3 = ax3.imshow(np.rot90(slit_fit), interpolation = 'nearest',
              extent = pltrg_fit, aspect = 'auto', cmap = 'RdYlBu_r', vmin = vcmin, vmax = vcmax)
    pl4 = ax4.imshow(np.rot90(slit_sersic), interpolation = 'nearest',
              extent = pltrg_fit, aspect = 'auto', cmap = 'RdYlBu_r', vmin = vcmin, vmax = vcmax)
    pl5 = ax5.imshow(np.rot90(slit_exp), interpolation = 'nearest',
              extent = pltrg_fit, aspect = 'auto', cmap = 'RdYlBu_r', vmin = vcmin, vmax = vcmax)
    dv3 = make_axes_locatable(ax3); ca3 = dv3.append_axes("right", size = "1%", pad = 0.05)
    dv4 = make_axes_locatable(ax4); ca4 = dv4.append_axes("right", size = "1%", pad = 0.05)
    dv5 = make_axes_locatable(ax5); ca5 = dv5.append_axes("right", size = "1%", pad = 0.05)
    cb3 = plt.colorbar(pl3, cax = ca3); cb4 = plt.colorbar(pl4, cax = ca4); cb5 = plt.colorbar(pl5, cax = ca5)

    # plot difference
    ax6 = fig.add_subplot(6, 1, 6)
    pl6 = ax6.imshow(np.rot90(slit_delta), interpolation = 'nearest',
              extent = pltrg_fit, aspect = 'auto', cmap = "seismic", vmin = vdmin, vmax = vdmax)
    dv6 = make_axes_locatable(ax6); ca6 = dv6.append_axes("right", size = "1%", pad = 0.05)
    cb6 = plt.colorbar(pl6, cax = ca6)

    # axis labels
    ax1.set_ylabel("Original"), ax2.set_ylabel("Rest-frame rebinned"), ax3.set_ylabel("Best-fitting")
    ax4.set_ylabel("Bulge"), ax5.set_ylabel("Disk"), ax6.set_ylabel("Residue")
    ax6.set_xlabel("Wavelength [A]")
    cb1.set_label("log flux"), cb2.set_label("log flux"), cb3.set_label("log flux")
    cb4.set_label("log flux"), cb5.set_label("log flux"), cb6.set_label("$\Delta$ flux")

    plt.savefig(basename + '-slitm.pdf', bbox_inches = 'tight')
    plt.clf(); plt.close()

  # "images" at three wavelengths
  if True:

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

    fig = plt.figure(figsize = (15., 7.5))
    plt_ext = [ra_ax[0], ra_ax[-1], dec_ax[0], dec_ax[-1]]

    ax1, ax4 = fig.add_subplot(2, 3, 1), fig.add_subplot(2, 3, 4)
    ax1.set_title("Rest-frame 4500A"); ax1.set_ylabel("Observed")
    ax4.set_ylabel("Dec [asec]"); ax4.set_xlabel("RA [asec]")
    pl1 = ax1.imshow(np.rot90(im450_obs), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet')
    pl4 = ax4.imshow(np.rot90(im450_fit), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet')
    dv1 = make_axes_locatable(ax1); ca1 = dv1.append_axes("right", size = "3%", pad = 0.05)
    dv4 = make_axes_locatable(ax4); ca4 = dv4.append_axes("right", size = "3%", pad = 0.05)
    cb1 = plt.colorbar(pl1, cax = ca1); cb4 = plt.colorbar(pl4, cax = ca4)

    ax2, ax5 = fig.add_subplot(2, 3, 2), fig.add_subplot(2, 3, 5)
    ax2.set_title("Rest-frame 5500A"); ax5.set_xlabel("RA [asec]")
    pl2 = ax2.imshow(np.rot90(im550_obs), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet')
    pl5 = ax5.imshow(np.rot90(im550_fit), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet')
    dv2 = make_axes_locatable(ax2); ca2 = dv2.append_axes("right", size = "3%", pad = 0.05)
    dv5 = make_axes_locatable(ax5); ca5 = dv5.append_axes("right", size = "3%", pad = 0.05)
    cb2 = plt.colorbar(pl2, cax = ca2); cb5 = plt.colorbar(pl5, cax = ca5)

    ax3, ax6 = fig.add_subplot(2, 3, 3), fig.add_subplot(2, 3, 6)
    ax3.set_title("Rest-frame 7000A"); ax6.set_xlabel("Best-fitting")
    pl3 = ax3.imshow(np.rot90(im700_obs), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet')
    pl6 = ax6.imshow(np.rot90(im700_fit), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet')
    dv3 = make_axes_locatable(ax3); ca3 = dv3.append_axes("right", size = "3%", pad = 0.05)
    dv6 = make_axes_locatable(ax6); ca6 = dv6.append_axes("right", size = "3%", pad = 0.05)
    cb3 = plt.colorbar(pl3, cax = ca3); cb6 = plt.colorbar(pl6, cax = ca6)

    plt.savefig(basename + '-imgs.pdf', bbox_inches = 'tight')
    plt.clf(), plt.close()

  # images at three wagelengths, 5 panels
  if True:

    wls = np.array([4500, 5500, 7000])
    idx_t, idx_w = np.searchsorted(wl_new, wls), np.searchsorted(wl_ax, wls / (1. + source_redshift))
    im450_obs, im550_obs, im700_obs = np.swapaxes(np.log10(np.abs(flux)), 1, 2)[idx_w, :, :]
    im450_fit, im550_fit, im700_fit = np.swapaxes(np.log10(np.abs(cube)), 1, 2)[idx_t, :, :]
    im450_s, im550_s, im700_s = np.swapaxes(np.log10(cube_sersic), 1, 2)[idx_t, :, :]
    im450_e, im550_e, im700_e = np.swapaxes(np.log10(cube_exp), 1, 2)[idx_t, :, :]
    im450_dt, im550_dt, im700_dt = im450_obs - im450_fit, im550_obs - im550_fit, im700_obs - im700_fit

    # mask bad points
    bad_450, bad_550, bad_700 = ~np.isfinite(im450_obs), ~np.isfinite(im550_obs), ~np.isfinite(im700_obs)
    im450_fit[bad_450] = np.nan; im550_fit[bad_550] = np.nan; im700_fit[bad_700] = np.nan
    im450_e[bad_450] = np.nan; im550_e[bad_550] = np.nan; im700_e[bad_700] = np.nan
    im450_s[bad_450] = np.nan; im550_s[bad_550] = np.nan; im700_s[bad_700] = np.nan
    im450_dt[bad_450] = np.nan; im550_dt[bad_550] = np.nan; im700_dt[bad_700] = np.nan

    # determine range of plotting
    imgs = np.hstack((im450_fit, im450_obs, im550_fit, im550_obs, im700_fit, im700_obs))
    imgs = imgs[np.isfinite(imgs)].ravel()
    hist_t, ed_t = np.histogram(imgs, bins = 512, normed = True,
        range = (imgs.min() - 0.1, imgs.max() + 0.1))
    hist_t = np.cumsum(hist_t); hist_t /= hist_t.max()
    ed_t = 0.5 * (ed_t[:-1] + ed_t[1:])
    plt_vlim_itpfc = interp1d(hist_t, ed_t, bounds_error = False)
    vcmin, vcmax = plt_vlim_itpfc(5.e-3), plt_vlim_itpfc(1. - 5.e-3)

    # determine range of differnece
    dt_img = np.hstack((im450_dt, im550_dt, im700_dt))
    dt_img = dt_img[np.isfinite(dt_img)].ravel()
    hist_t, ed_t = np.histogram(dt_img, bins = 512, normed = True,
        range = (dt_img.min() - 0.1, dt_img.max() + 0.1))
    hist_t = np.cumsum(hist_t); hist_t /= hist_t.max()
    ed_t = 0.5 * (ed_t[:-1] + ed_t[1:])
    plt_dlim_itpfc = interp1d(hist_t, ed_t, bounds_error = False)
    vd_rg = np.abs(plt_dlim_itpfc([5.e-2, 1. - 5.e-2]))
    vdmin, vdmax = -vd_rg.min(), vd_rg.min()

    # make figure
    fig = plt.figure(figsize = (15., 22.5))
    plt_ext = [ra_ax[0], ra_ax[-1], dec_ax[0], dec_ax[-1]]

    # plot original
    ax1, ax2, ax3 = fig.add_subplot(5, 3, 1), fig.add_subplot(5, 3, 2), fig.add_subplot(5, 3, 3)
    pl1 = ax1.imshow(np.rot90(im450_obs), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet' )
    pl2 = ax2.imshow(np.rot90(im550_obs), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet' )
    pl3 = ax3.imshow(np.rot90(im700_obs), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet' )
    dv1 = make_axes_locatable(ax1); ca1 = dv1.append_axes("right", size = "3%", pad = 0.05)
    dv2 = make_axes_locatable(ax2); ca2 = dv2.append_axes("right", size = "3%", pad = 0.05)
    dv3 = make_axes_locatable(ax3); ca3 = dv3.append_axes("right", size = "3%", pad = 0.05)
    cb1 = plt.colorbar(pl1, cax = ca1, ticks = [])
    cb2 = plt.colorbar(pl2, cax = ca2, ticks = [])
    cb3 = plt.colorbar(pl3, cax = ca3); cb3.set_label("log flux")

    # plot best fitting
    ax4, ax5, ax6 = fig.add_subplot(5, 3, 4), fig.add_subplot(5, 3, 5), fig.add_subplot(5, 3, 6)
    pl4 = ax4.imshow(np.rot90(im450_fit), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet' )
    pl5 = ax5.imshow(np.rot90(im550_fit), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet' )
    pl6 = ax6.imshow(np.rot90(im700_fit), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet' )
    dv4 = make_axes_locatable(ax4); ca4 = dv4.append_axes("right", size = "3%", pad = 0.05)
    dv5 = make_axes_locatable(ax5); ca5 = dv5.append_axes("right", size = "3%", pad = 0.05)
    dv6 = make_axes_locatable(ax6); ca6 = dv6.append_axes("right", size = "3%", pad = 0.05)
    cb4 = plt.colorbar(pl4, cax = ca4, ticks = [])
    cb5 = plt.colorbar(pl5, cax = ca5, ticks = [])
    cb6 = plt.colorbar(pl6, cax = ca6); cb6.set_label("log flux")

    # plot sersic
    ax7, ax8, ax9 = fig.add_subplot(5, 3, 7), fig.add_subplot(5, 3, 8), fig.add_subplot(5, 3, 9)
    pl7 = ax7.imshow(np.rot90(im450_s), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet' )
    pl8 = ax8.imshow(np.rot90(im550_s), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet' )
    pl9 = ax9.imshow(np.rot90(im700_s), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet' )
    dv7 = make_axes_locatable(ax7); ca7 = dv7.append_axes("right", size = "3%", pad = 0.05)
    dv8 = make_axes_locatable(ax8); ca8 = dv8.append_axes("right", size = "3%", pad = 0.05)
    dv9 = make_axes_locatable(ax9); ca9 = dv9.append_axes("right", size = "3%", pad = 0.05)
    cb7 = plt.colorbar(pl7, cax = ca7, ticks = [])
    cb8 = plt.colorbar(pl8, cax = ca8, ticks = [])
    cb9 = plt.colorbar(pl9, cax = ca9); cb9.set_label("log flux")

    # plot exp
    ax10, ax11, ax12 = fig.add_subplot(5, 3, 10), fig.add_subplot(5, 3, 11), fig.add_subplot(5, 3, 12)
    pl10 = ax10.imshow(np.rot90(im450_e), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet' )
    pl11 = ax11.imshow(np.rot90(im550_e), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet' )
    pl12 = ax12.imshow(np.rot90(im700_e), extent = plt_ext,
        interpolation = 'nearest', vmin = vcmin, vmax = vcmax, cmap = 'jet' )
    dv10 = make_axes_locatable(ax10); ca10 = dv10.append_axes("right", size = "3%", pad = 0.05)
    dv11 = make_axes_locatable(ax11); ca11 = dv11.append_axes("right", size = "3%", pad = 0.05)
    dv12 = make_axes_locatable(ax12); ca12 = dv12.append_axes("right", size = "3%", pad = 0.05)
    cb10 = plt.colorbar(pl10, cax = ca10, ticks = [])
    cb11 = plt.colorbar(pl11, cax = ca11, ticks = [])
    cb12 = plt.colorbar(pl12, cax = ca12); cb12.set_label("log flux")

    # plot diff
    ax13 , ax14 , ax15  = fig.add_subplot(5, 3, 13), fig.add_subplot(5, 3, 14), fig.add_subplot(5, 3, 15)
    pl13  = ax13 .imshow(np.rot90(im450_dt), extent = plt_ext,
        interpolation = 'nearest', vmin = vdmin , vmax = vdmax , cmap = 'seismic')
    pl14  = ax14 .imshow(np.rot90(im550_dt), extent = plt_ext,
        interpolation = 'nearest', vmin = vdmin , vmax = vdmax , cmap = 'seismic')
    pl15  = ax15 .imshow(np.rot90(im700_dt), extent = plt_ext,
        interpolation = 'nearest', vmin = vdmin , vmax = vdmax , cmap = 'seismic')
    dv13  = make_axes_locatable(ax13 ); ca13  = dv13 .append_axes("right", size = "3%", pad = 0.05)
    dv14  = make_axes_locatable(ax14 ); ca14  = dv14 .append_axes("right", size = "3%", pad = 0.05)
    dv15  = make_axes_locatable(ax15 ); ca15  = dv15 .append_axes("right", size = "3%", pad = 0.05)
    cb13  = plt.colorbar(pl13 , cax = ca13, ticks = [] )
    cb14  = plt.colorbar(pl14 , cax = ca14, ticks = [])
    cb15  = plt.colorbar(pl15 , cax = ca15 ); cb15.set_label("$\Delta$ flux")

    # labels
    ax1.set_title("Rest-frame 4500A"); ax2.set_title("Rest-frame 5500A"); ax3.set_title("Rest-frame 7000A")
    ax1.set_ylabel("Observed"); ax4.set_ylabel("Best-fitting"); ax7.set_ylabel("Sersic"); ax10.set_ylabel("Exponential")
    ax13.set_ylabel("Residue"); ax13.set_xlabel('RA [arcsec]')
    ax14.set_xlabel('RA [arcsec]'); ax15.set_xlabel('RA [arcsec]')

    plt.savefig(basename + '-imgsm.pdf', bbox_inches = 'tight')
    plt.clf(), plt.close()

  # triangle plot
  if True:

    fig = corner.corner(chain, labels = par_names, truths = par_opt)
    fig.savefig(basename + '-prob.pdf', bbox_inches = 'tight')
