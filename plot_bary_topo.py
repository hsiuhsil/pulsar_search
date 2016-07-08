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
            plot_bary_topo(filename)
        except (IOError, ValueError):
            print IOError

def plot_bary_topo(filename):
    this_file = h5py.File(filename, "r")
    bary = this_file['BARY_TIME'][:]
    topo = this_file['TOPO_TIME'][:]
    diff = bary - topo
    plt.plot(diff, bary)
    plt.axis([diff[0],diff[-1], bary[0], bary[-1]], labelsize=20)
    plt.xlabel('BARY - TOPO', fontsize=20)
    plt.ylabel('BARY', fontsize=20)
    plt.show()

if __name__ == '__main__':
    main()

