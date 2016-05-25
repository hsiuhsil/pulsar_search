import sys
import os
import os.path

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            ploting(filename)
        except (IOError, ValueError):
            print IOError

def ploting(filename):
    this_file = h5py.File(filename, "r")
    data = np.zeros((this_file['DATA_FOLDING'].shape[0], this_file['DATA_FOLDING'].shape[2]))
    for ii in range(len(data)):
        data[ii]=this_file['DATA_FOLDING'][ii][0][:].flatten()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.imshow(data, aspect='auto', cmap=cm.Greys)
    cbar = plt.colorbar(cax)
    plt.show()

if __name__ == '__main__':
    main()
