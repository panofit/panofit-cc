#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

arr = np.fromfile(sys.argv[1], dtype = 'f8').reshape(int(sys.argv[2]), int(sys.argv[3]))
plt.imshow(np.rot90(arr), interpolation = 'nearest'), plt.colorbar(), plt.show()
