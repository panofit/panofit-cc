#!/usr/bin/python

'''
  Plot the result of dual-ssp disk + bulge fitting.
  YJ Qin, Jan 2016 @ Shanghai

  Parameters to optimize: sersic_n,   sersic_re,  sersic_Ie,  sersic_phi, sersic_q
                          sersic_c,   sersic_age, exp_Ic,     exp_h,      exp_phi,
                          exp_q,      exp_c,      exp_age,
'''

import sys
import numpy as np
import matplotlib.pyplot as plt
import corner

N_dim, N_walkers = 13, 26
chain_file = 'fsps-dssp-specfit-cube-NGC0001V500-PSF2.tx.npy'
# par_opt = (2.75, 7.5, 0.3, 0.3, 0.75, 2.25, 1.0, 2., 7.5, 0.85, 0.5, 2.0, 0.3)
par_names = ["Sersic n", "Bulge Re", "Bulge Ie", "Bulge phi", "Bulge b/a",
             "Bulge c", "Bulge log(Age)", "Disk I_c", "Disk R_h", "Disk phi",
             "Disk b/a", "Disk c", "Disk log(Age)"]

if __name__ == "__main__":

  chain_file = sys.argv[-1]
  chain = np.load(chain_file)

  # remove bad points
  chain = chain.reshape((-1, N_dim))
  chain = chain[chain[:,  4] <= 1., :]
  chain = chain[chain[:, 10] <= 1., :]

  # fig = corner.corner(chain, labels = par_names, truths = par_opt)
  fig = corner.corner(chain, labels = par_names)
  fig.savefig(chain_file + '.pdf')
