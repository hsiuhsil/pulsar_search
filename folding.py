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

    '''Create 'DATA_FOLDING' to record folding data'''
    
    first_data = this_file['DATA'][0][0]
    this_file.create_dataset('DATA_FOLDING', (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)

    '''Define phase_bins parameter'''

    phase_bins = 100
    pulsar_period = 0.312470
    for ii in range(0,phase_bins):
        modulo = pulsar_period/phase_bins*ii
        current_len = this_file['DATA_FOLDING'].shape[0]
        this_file['DATA_FOLDING'].resize(current_len + 1, 0)
        this_file['DATA_FOLDING'][ii] = np.float64(modulo)

    '''collecting data which satisfies the folding condition'''
    n_bins = 2048
    tbin = 0.001024
    print tbin
    same_modulo_num = [0]*phase_bins
#    for ii in range(len(f['BARY_TIME'])-1):
    for ii in range(1):
        for jj in range(n_bins):
            print jj
            sample_BAT = this_file['BARY_TIME'][ii] + ( jj - n_bins/2.0 + 0.5) * tbin
            modulo_num = np.int64(np.around((sample_BAT % pulsar_period)/(pulsar_period/phase_bins)))
            if modulo_num == phase_bins:
                modulo_num = 0
            print modulo_num
            same_modulo_num[modulo_num] += 1
            this_file['DATA_FOLDING'][modulo_num,...] += this_file['DATA'][ii][jj]
    return this_file['DATA_FOLDING']

    for modulo_num in range(len(same_modulo_num)):
        if same_modulo_num[modulo_num] != 0:
            this_file['DATA_FOLDING'][ii,...] = this_file['DATA_FOLDING'][ii,...]/same_modulo_num[ii]
    return this_file['DATA_FOLDING']

if __name__ == '__main__':
    main()
