import sys
import os
import os.path

import h5py
import numpy as np

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            folding(filename)
        except (IOError, ValueError):
            print IOError

def folding(filename):
    this_file = h5py.File(filename, "r+")
    phase_bins = 100
    pulsar_period = 0.312470
    n_bins = 2048
    tbin = 0.001024
    
    first_data = this_file['DATA'][0][0]
    data_folding = np.zeros((phase_bins,) + first_data.shape)

    '''collecting data which satisfies the folding condition'''
    same_modulo_num = [0]*phase_bins
    for ii in range(len(this_file['BARY_TIME'])-1):
#        print 'ii = ' + str(ii)
        sample_BAT = this_file['BARY_TIME'][ii] + np.arange(-n_bins/2.0 + 0.5, n_bins/2.0 + 0.5)*tbin
        modulo_num = np.int64(np.around((sample_BAT % pulsar_period)/(pulsar_period/phase_bins)))
        for jj in range(len(modulo_num)):
            if modulo_num[jj] == phase_bins:
                modulo_num[jj] = 0
        for kk in range(len(same_modulo_num)):
            same_modulo_num[kk] += np.count_nonzero(modulo_num == kk)      

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
