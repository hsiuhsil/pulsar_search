import sys
import os
import os.path

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from folding import *

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            plot_bary_diff(filename)
        except (IOError, ValueError):
            print IOError

def plot_bary_diff(filename):
    this_file = h5py.File(filename, "r")

#    index = [[1,3], [4,6], [7,8], [11,12], [15,17], [18,19], [20,22], [23,27], [28,30], [33,38], [39,42], [44,49], [50,55], [61,62], [63,64], [65,66], [67,72], [73,76]]
#    max_phase = [89, 90, 90, 89, 89, 89, 60, 59, 59, 56, 59, 56, 59, 55, 57, 56, 47, 47] 

#    index = [[4,8], [140,143], [163,165], [493,495], [545,547], [578, 580]]
#    max_phase = [1, 93, 79, 46, 96, 56]

    index = [[4,8], [15, 19],[28, 32],[44, 49],[67, 83],[98,113],[114,117],[118,128],[146,152]]
    max_phase = [69, 69, 63, 65, 90, 97, 91, 95, 5]

    bary_diff = np.zeros(len(index))

    for ii in range(len(index)):
        bary_diff[ii] = ((this_file['BARY_TIME'][index[ii][0]]+this_file['BARY_TIME'][index[ii][1]])/2. -this_file['BARY_TIME'][0])*86400/60/60

    title = 'delta_t: '+str(delta_t)+' sec.'
    
    plt.plot(bary_diff, max_phase, 'bo')
#    plt.axis([diff[0],diff[-1], bary[0], bary[-1]], labelsize=20)
    plt.xlabel('Bary diff (hours)', fontsize=20)
    plt.ylabel('Max Phase Bins Number', fontsize=20)
    plt.title(title, fontsize=20)
    plt.show()

if __name__ == '__main__':
    main()

