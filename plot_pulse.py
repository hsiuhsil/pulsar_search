import sys
import os
import os.path

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from folding import *

rebin = True
rebin_factor = 128
dm = 36.0
dedisperse = True

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            ploting(filename)
        except (IOError, ValueError):
            print IOError

def rebin_spec(input_data):
    output_data = np.zeros((input_data.shape[0], input_data.shape[2]/rebin_factor))    
    for ii in range(len(output_data)):
        for jj in range(len(output_data[1])):
            output_data[ii,jj]=np.mean(input_data[ii,0,jj*rebin_factor:(jj+1)*rebin_factor])
    return output_data

def dedisperse_index(freq1, freq2):
    delta_t = pulsar_period/phase_bins
    DM_CONST = 4148.808
    index = np.int64(np.around((DM_CONST * dm * (freq1**-2 - freq2**-2))/delta_t))
    return index

def dedisperse_spec(input_data):
    freq = np.arange(900., 700., -200./input_data.shape[1])    
    dedis_index_range = dedisperse_index(freq[:], np.amax(freq)) % phase_bins
    output_data = np.zeros((input_data.shape[0], input_data.shape[1]))
    scan_range = output_data.shape[0] - np.max(dedis_index_range)

    for jj in range(output_data.shape[1]):
        for ii in range(dedis_index_range[jj], dedis_index_range[jj]+scan_range):
            output_data[ii-dedis_index_range[jj],jj] = input_data[ii, jj]
    return output_data

def ploting(filename):
    this_file = h5py.File(filename, "r")
    '''Do rebin first and then do dedisperse'''
  
    if rebin == True and dedisperse == True:
        data_first = rebin_spec(this_file['DATA_FOLDING'])
        data = dedisperse_spec(data_first)
    elif rebin == True and dedisperse == False:
        data = rebin_spec(this_file['DATA_FOLDING'])
    elif rebin == False and dedisperse == True:
        data = dedisperse_spec(this_file['DATA_FOLDING'][:, 0, :, 0])
    else:
        data = this_file['DATA_FOLDING'][:, 0, :, 0]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.imshow(data, aspect='auto', cmap=cm.Greys)
    cbar = plt.colorbar(cax)
    plt.show()

if __name__ == '__main__':
    main()
