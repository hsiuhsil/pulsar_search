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
            ploting(filename)
        except (IOError, ValueError):
            print IOError

def rebin_spec(input_data):
    output_data = np.zeros((input_data.shape[0], input_data.shape[1]/rebin_factor))    
    for ii in range(len(output_data)):
        for jj in range(len(output_data[1])):
            output_data[ii,jj]=np.mean(input_data[ii, jj*rebin_factor:(jj+1)*rebin_factor])
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
        output_data[:,jj] = np.roll(input_data[:,jj], -dedis_index_range[jj], axis = 0)
    return output_data    

def ploting(filename):
    this_file = h5py.File(filename, "r")
    '''Do rebin first and then do dedisperse'''
  
    if rebin == True and dedisperse == True:
        data_first = dedisperse_spec(this_file['DATA_FOLDING'][:, 0, :, 0])
        data = rebin_spec(data_first).T               
    elif rebin == True and dedisperse == False:
        data = rebin_spec(this_file['DATA_FOLDING']).T
    elif rebin == False and dedisperse == True:
        data = dedisperse_spec(this_file['DATA_FOLDING'][:, 0, :, 0]).T
    else:
        data = this_file['DATA_FOLDING'][:, 0, :, 0].T

    '''get mean over phase_bins'''
    data2 = [0.]*len(this_file['DATA_FOLDING'][:, 0, :, 0])
    for ii in range(len(data2)):
        data2[ii] = np.mean(data_first[ii, :])  
    phase_bin_max_ind = np.argmax(data2)
    
    subtitle = 'max amp at phase bin: '+str(phase_bin_max_ind) 

    fig = plt.figure()
    fig.suptitle(subtitle, fontsize=20)
    ax1 = fig.add_subplot(211)
    ax1.set_ylabel('Freq(MHz)', fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    cax1 = ax1.imshow(data, extent=[0, 99, 700., 900.],aspect='auto', cmap=cm.Greys)
    cbar = plt.colorbar(cax1)
    ax2 = fig.add_subplot(212)
    ax2.set_ylabel('Mean Amp', fontsize=20)
    ax2.set_xlabel('Phase Bins', fontsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    cax2 = plt.plot(data2)
    fig.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
