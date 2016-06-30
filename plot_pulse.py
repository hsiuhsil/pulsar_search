import sys
import os
import os.path

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from folding import *

dedisperse = True
dm = 0.0

rebin = True
rebin_factor = 128


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
        data_third = dedisperse_spec(this_file['DATA_FOLDING_TOPO'][:, 0, :, 0])
        data3 = rebin_spec(data_third).T
    elif rebin == True and dedisperse == False:
        data_first = this_file['DATA_FOLDING'][:, 0, :, 0]
        data = rebin_spec(data_first).T
        data_third = this_file['DATA_FOLDING_TOPO'][:, 0, :, 0]
        data3 = rebin_spec(data_third).T
    elif rebin == False and dedisperse == True:
        data = dedisperse_spec(this_file['DATA_FOLDING'][:, 0, :, 0]).T
        data3 = dedisperse_spec(this_file['DATA_FOLDING_TOPO'][:, 0, :, 0]).T
    else:
        data = this_file['DATA_FOLDING'][:, 0, :, 0]
        data3 = this_file['DATA_FOLDING_TOPO'][:, 0, :, 0]

    '''get mean over phase_bins'''
    data2 = [0.]*data.shape[1]
    data4 = [0.]*data3.shape[1]

    for ii in range(len(data2)):
        data2[ii] = np.mean(data[:,ii])
    phase_bin_max_ind = np.argmax(data2)

    for ii in range(len(data4)):
        data4[ii] = np.mean(data3[:,ii])
    phase_bin_max_ind3 = np.argmax(data4)


    subtitle = 'max amp at phase bin for BARY: '+str(phase_bin_max_ind) +', TOPO: '+str(phase_bin_max_ind3) 

    fig = plt.figure()
    fig.suptitle(subtitle, fontsize=20)
    ax1 = fig.add_subplot(411)
    ax1.set_ylabel('Freq(MHz)', fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    cax1 = ax1.imshow(data, extent=[0, phase_bins-1, 700., 900.],aspect='auto', cmap=cm.Greys)
    cbar = plt.colorbar(cax1)
    ax2 = fig.add_subplot(412)
    ax2.set_ylabel('Mean Amp', fontsize=20)
    ax2.set_xlabel('Phase Bins', fontsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    cax2 = plt.plot(data2)
    ax3 = fig.add_subplot(413)
    ax3.set_ylabel('Freq(MHz)', fontsize=20)
    ax3.tick_params(axis='both', which='major', labelsize=20)
    cax3 = ax3.imshow(data3, extent=[0, phase_bins-1, 700., 900.],aspect='auto', cmap=cm.Greys)
    cbar = plt.colorbar(cax3)
    ax4 = fig.add_subplot(414)
    ax4.set_ylabel('Mean Amp', fontsize=20)
    ax4.set_xlabel('Phase Bins', fontsize=20)
    ax4.tick_params(axis='both', which='major', labelsize=20)
    cax4 = plt.plot(data4)
#    fig.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
