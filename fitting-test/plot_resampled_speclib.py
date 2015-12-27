#!/usr/bin/python
# DEBUG

'''
  Plot spectral library
'''

import numpy as np
import matplotlib.pyplot as plt

def _plot_spec_lib(fname):

  # read file
  arr = np.load(fname)

  # read dimension
  print 'array shape', arr.shape

  # plot a spectrum
  plt.plot(arr[-1, 6, :]), plt.show()

if __name__ == "__main__":
  _plot_spec_lib("./NGC0001/speclib-V500.dat.dbg.npy")
  _plot_spec_lib("./NGC0001/speclib-V1200.dat.dbg.npy")

# wo8yh0ZRmDJbKnKosUFC YJQIN DEC15
