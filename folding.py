import sys
import os
import os.path
sys.path.append('/home/p/pen/hsiuhsil/burst_search/')
from burst_search import preprocess

import h5py
import numpy as np
import math

'''Define variables'''

initial = 43
final = 44

do_preprocess = True
sigma_threshold = 5
remove_period = 64
phase_bins = 100

# J2139: delta_t = -5.8068471774736951e-06 -2.5301321417973068e-07 -4.2596495836139012e-07 + 2.9914132482230684e-07 
# J2139: pulsar_period = 0.312470 + delta_t

    #J0051
#delta_t = -3.2083428581595565e-07 +1.4583350248377532e-08 -7.0619503654410745e-09 +5.1170668232711475e-10 
#pulsar_period = 0.35473179890 + delta_t

# J1046: 
delta_t = 0.0
pulsar_period = 0.326271446035

# J1038: delta_t = 0.0
# J1038: pulsar_period = 0.02885155795131

# J0030: delta_t = 0.0
# J0030: pulsar_period = 0.0048654532073692


def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            folding(filename)
        except (IOError, ValueError):
            print IOError

def time_slope(input_data):
    print "start time_slope"
    slope_mode = np.arange(np.float(input_data.shape[1]))
    slope_mode -= np.mean(slope_mode)
    slope_mode /= math.sqrt(np.sum(slope_mode**2))
    slope_amplitude = np.sum(input_data * slope_mode[None,:], 1)
    input_data -= slope_amplitude[:,None] * slope_mode
    return input_data

def preprocessing(input_data):
    '''note: preprocess need data.shape = (nfreq, ntime)'''
    output_data = np.zeros(input_data.shape)
    data = input_data[:,0,:,0].T.astype(np.float64).copy()
    preprocess.remove_periodic(data, remove_period)
    m = np.mean(data[:],axis=1)
    m[m==0]=1
    data = data / m[:,None] - 1
    preprocess.remove_noisy_freq(data, sigma_threshold)
    data = data-np.mean(data)
    data = time_slope(data)
    preprocess.remove_bad_times(data, sigma_threshold)
    preprocess.remove_noisy_freq(data, sigma_threshold)
    output_data[:,0,:,0] = data.T
    output_data[:,1:4,:,0] = input_data[:,1:4,:,0]
    return output_data

def folding(filename):

    this_file = h5py.File(filename, "r+")   
    ntime = this_file['DATA'].shape[1]
    tbin = this_file['TBIN'][0]

    first_data = this_file['DATA'][0][0]
    data_folding = np.zeros((phase_bins,) + first_data.shape)
    
    '''collecting data which satisfies the folding condition'''
    same_modulo_num = np.zeros((phase_bins,), dtype=np.int)
    for ii in range(initial, final+1):
#    for ii in range(len(this_file['BARY_TIME'])):
        print 'ii = ' + str(ii)
        sample_BAT = this_file['BARY_TIME'][ii]*86400 + np.arange(-ntime/2.0 + 0.5, ntime/2.0 + 0.5)*tbin
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
    same_modulo_num_topo = np.zeros((phase_bins,), dtype=np.int)
    for ii in range(initial, final+1):
#    for ii in range(len(this_file['TOPO_TIME'])):
        print 'ii = ' + str(ii)
        sample_TOPO = this_file['TOPO_TIME'][ii]*86400 + np.arange(-ntime/2.0 + 0.5, ntime/2.0 + 0.5)*tbin
        modulo_num_topo = np.int64(np.around((sample_TOPO % pulsar_period)/(pulsar_period/phase_bins)))
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
            data_folding_topo[modulo_num_topo[ll],...] += this_record_data_topo[ll]

    for mm in range(len(same_modulo_num_topo)):
        if same_modulo_num_topo[mm] != 0:
            data_folding_topo[mm,...] = data_folding_topo[mm,...]/same_modulo_num_topo[mm]

    this_file.create_dataset('DATA_FOLDING_TOPO', data_folding_topo.shape, maxshape = data_folding_topo.shape, dtype=data_folding_topo.dtype, chunks=True)
    this_file['DATA_FOLDING_TOPO'][...]=data_folding_topo[...]



if __name__ == '__main__':
    main()
