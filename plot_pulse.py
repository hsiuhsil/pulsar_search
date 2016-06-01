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

def rebin_spec(input_data, rebin_factor):
    output_data = np.zeros((input_data.shape[0], input_data.shape[2]/rebin_factor))    
    for ii in range(len(output_data)):
        for jj in range(len(output_data[1])):
            output_data[ii,jj]=np.mean(input_data[ii,0,jj*rebin_factor:(jj+1)*rebin_factor])
    return output_data

def ploting(filename):
    this_file = h5py.File(filename, "r")
    rebin = True
    rebin_factor = 128
  
    if rebin == True:
        data = rebin_spec(this_file['DATA_FOLDING'], rebin_factor)
    else:
        data = np.zeros((this_file['DATA_FOLDING'].shape[0], this_file['DATA_FOLDING'].shape[2]))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.imshow(data, aspect='auto', cmap=cm.Greys)
    cbar = plt.colorbar(cax)
    plt.show()

if __name__ == '__main__':
    main()
