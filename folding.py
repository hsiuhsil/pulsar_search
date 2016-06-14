import sys
import os
import os.path
sys.path.append('/home/p/pen/hsiuhsil/burst_search/')
from burst_search import preprocess

import h5py
import numpy as np
import math

'''Define variables'''
ntime = this_file['NCHAN'][0]
tbin = this_file['TBIN'][0]

do_preprocess = True
sigma_threshold = 5
remove_period = 64
phase_bins = 100
pulsar_period = 0.312470



def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            folding(filename)
        except (IOError, ValueError):
            print IOError

def time_slope(input_data):
    slope_mode = np.arange(ntime)
    slope_mode -= np.mean(slope_mode)
    slope_mode /= math.sqrt(np.sum(slope_mode**2))
    slope_amplitude = np.sum(input_data * slope_mode[:,None], 0)
    input_data -= slope_amplitude * slope_mode[:,None]
    return input_data

def preprocessing(input_data):
    output_data = np.zeros(input_data.shape)
    data = input_data[:,0,:,0].T
    m = np.mean(data[:],axis=1)
    m[m==0]=1 
    data = data / m[:,None] - 1
    preprocess.remove_noisy_freq(data, sigma_threshold)
    data = data-np.mean(data)
    data = time_slope(data)
    preprocess.remove_noisy_freq(data, sigma_threshold)
#    preprocess.remove_periodic(data, remove_period)
    output_data[:,0,:,0] = data.T
    output_data[:,1:4,:,0] = input_data[:,1:4,:,0]
    return output_data

def folding(filename):

    this_file = h5py.File(filename, "r+")   

    first_data = this_file['DATA'][0][0]
    data_folding = np.zeros((phase_bins,) + first_data.shape)
    
    '''collecting data which satisfies the folding condition'''
    same_modulo_num = [0]*phase_bins
#    for ii in range(1):
    for ii in range(len(this_file['BARY_TIME'])):
        print 'ii = ' + str(ii)
        sample_BAT = this_file['BARY_TIME'][ii] + np.arange(-ntime/2.0 + 0.5, ntime/2.0 + 0.5)*tbin
        modulo_num = np.int64(np.around((sample_BAT % pulsar_period)/(pulsar_period/phase_bins)))
        print 'modulo_num done'
        for jj in range(len(modulo_num)):
            if modulo_num[jj] == phase_bins:
                modulo_num[jj] = 0
        for kk in range(len(same_modulo_num)):
            same_modulo_num[kk] += np.count_nonzero(modulo_num == kk)      

        if do_preprocess == True:
           this_record_data = preprocessing(this_file['DATA'][ii])
           print 'preprocess done'
        else:
           this_record_data = this_file['DATA'][ii]

        for ll in range(len(modulo_num)):
            data_folding[modulo_num[ll],...] += this_record_data[ll]
        
    for mm in range(len(same_modulo_num)):
        if same_modulo_num[mm] != 0:
            data_folding[mm,...] = data_folding[mm,...]/same_modulo_num[mm]

    this_file.create_dataset('DATA_FOLDING', data_folding.shape, maxshape = data_folding.shape, dtype=data_folding.dtype, chunks=True)
    this_file['DATA_FOLDING'][...]=data_folding[...]

    '''Data folding for topocentric time'''

    data_folding_topo = np.zeros((phase_bins,) + first_data.shape)

    '''collecting data which satisfies the folding condition'''
    same_modulo_num_topo = [0]*phase_bins
#    for ii in range(1):
    for ii in range(len(this_file['TOPO_TIME'])):
        print 'ii = ' + str(ii)
        sample_BAT_topo = this_file['TOPO_TIME'][ii] + np.arange(-ntime/2.0 + 0.5, ntime/2.0 + 0.5)*tbin
        modulo_num_topo = np.int64(np.around((sample_BAT_topo % pulsar_period)/(pulsar_period/phase_bins)))
        print 'modulo_num done'
        for jj in range(len(modulo_num_topo)):
            if modulo_num_topo[jj] == phase_bins:
                modulo_num_topo[jj] = 0
        for kk in range(len(same_modulo_num_topo)):
            same_modulo_num_topo[kk] += np.count_nonzero(modulo_num_topo == kk)

        if do_preprocess == True:
           this_record_data_topo = preprocessing(this_file['DATA'][ii])
           print 'preprocess done'
        else:
           this_record_data_topo = this_file['DATA'][ii]

        for ll in range(len(modulo_num_topo)):
            data_folding_topo[modulo_num_topo[ll],...] += this_record_data[ll]

    for mm in range(len(same_modulo_num_topo)):
        if same_modulo_num_topo[mm] != 0:
            data_folding_topo[mm,...] = data_folding_topo[mm,...]/same_modulo_num_topo[mm]

    this_file.create_dataset('DATA_FOLDING_TOPO', data_folding_topo.shape, maxshape = data_folding_topo.shape, dtype=data_folding_topo.dtype, chunks=True)
    this_file['DATA_FOLDING_TOPO'][...]=data_folding_topo[...]


if __name__ == '__main__':
    main()
