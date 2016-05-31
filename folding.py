import sys
import os
import os.path
sys.path.append('/home/p/pen/hsiuhsil/burst_search/')
from burst_search import preprocess

import h5py
import numpy as np
import math

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            folding(filename)
        except (IOError, ValueError):
            print IOError

def time_slope(input_data):
    ntime = 2048
    slope_mode = np.arange(ntime)
    slope_mode -= np.mean(slope_mode)
    slope_mode /= math.sqrt(np.sum(slope_mode**2))
    slope_amplitude = np.sum(input_data * slope_mode)
    input_data -= slope_amplitude * slope_mode
    return input_data

def preprocessing(input_data):

    sigma_threshold = 5
    remove_period = 64

    data = input_data[:,0,:,0].T
    m = np.mean(data[:],axis=1)
    m[m==0]=1 
    data = data / m[:,None] - 1
    preprocess.remove_noisy_freq(data, sigma_threshold)
    data = data-np.mean(data)
    data = time_slope(data)
    preprocess.remove_noisy_freq(data, sigma_threshold)
    preprocess.remove_periodic(data, remove_period)
    return data.T 

def folding(filename):
    global this_file
    this_file = h5py.File(filename, "r+")
    phase_bins = 100
    pulsar_period = 0.312470
    n_bins = 2048
    tbin = 0.001024
    do_preprocess = True   

    first_data = this_file['DATA'][0]
    data_folding = np.zeros((phase_bins,) + (first_data.shape[0], first_data.shape[2]))
    
    '''collecting data which satisfies the folding condition'''
    same_modulo_num = [0]*phase_bins
#    for ii in range(1):
    for ii in range(len(this_file['BARY_TIME'])-1):
        print 'ii = ' + str(ii)
        sample_BAT = this_file['BARY_TIME'][ii] + np.arange(-n_bins/2.0 + 0.5, n_bins/2.0 + 0.5)*tbin
        modulo_num = np.int64(np.around((sample_BAT % pulsar_period)/(pulsar_period/phase_bins)))
        for jj in range(len(modulo_num)):
            if modulo_num[jj] == phase_bins:
                modulo_num[jj] = 0
        for kk in range(len(same_modulo_num)):
            same_modulo_num[kk] += np.count_nonzero(modulo_num == kk)      

        if do_preprocess == True:
           this_record_data = preprocessing(this_file['DATA'][ii])
        else:
           this_record_data = this_file['DATA'][ii]

        for ll in range(len(modulo_num)):
            data_folding[modulo_num[ll],...] += this_record_data[ll]

    for mm in range(len(same_modulo_num)):
        if same_modulo_num[mm] != 0:
            data_folding[mm,...] = data_folding[mm,...]/same_modulo_num[mm]

    this_file.create_dataset('DATA_FOLDING', data_folding.shape, maxshape = data_folding.shape, dtype=data_folding.dtype, chunks=True)
    this_file['DATA_FOLDING'][...]=data_folding[...]

if __name__ == '__main__':
    main()
