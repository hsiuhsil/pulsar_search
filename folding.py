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
    file = h5py.File(filename, "r")
    new_data = h5py.File('folding_'+filename, "w")
    keys = file.keys()+['DATA_FOLDING']

    '''Create 'DATA_FOLDING' to record folding data'''
    
    for dataset_name in keys:
        if dataset_name == 'DATA_FOLDING':
            first_data = file['DATA'][0][0]
            new_data.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
        else:
            first_data = file[dataset_name][0]
            new_data.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)

    '''copy data to new file'''

    for ii in range(len(file['BARY_TIME'])):
        for dataset_name in file.keys():
            current_len = new_data[dataset_name].shape[0]
            new_data[dataset_name].resize(current_len + 1, 0)
            new_data[dataset_name][current_len-1,...] = file[dataset_name][ii]    

    '''Define phase_bins parameter'''

    phase_bins = 100
    pulsar_period = 0.312470
    for ii in range(0,phase_bins):
        modulo = pulsar_period/phase_bins*ii
        current_len = new_data['DATA_FOLDING'].shape[0]
        new_data['DATA_FOLDING'].resize(current_len + 1, 0)
        new_data['DATA_FOLDING'][ii] = np.float64(modulo)

    '''collecting data which satisfies the folding condition'''
    n_bins = 2048
    tbin = 0.001024
    print tbin
    same_modulo_num = [0]*phase_bins
#    for ii in range(len(f['BARY_TIME'])-1):
    for ii in range(1):
        for jj in range(n_bins):
            sample_BAT = new_data['BARY_TIME'][ii] + ( jj - n_bins/2.0 + 0.5) * tbin
            modulo_num = np.int32(np.around((sample_BAT % pulsar_period)/(pulsar_period/phase_bins)))
            if modulo_num == phase_bins:
                return modulo_num == 0
            current_len = new_data['DATA_FOLDING'].shape[0]
            new_data['DATA_FOLDING'].resize(current_len + 1, 0)
            new_data['DATA_FOLDING'][modulo_num,...] += new_data['DATA'][ii][jj]
            same_modulo_num[modulo_num] += 1

    for ii in range(len(new_data['DATA_FOLDING'])):
        new_data['DATA_FOLDING'][ii,...] = new_data['DATA_FOLDING'][ii,...]/same_modulo_num[ii]


if __name__ == '__main__':
    main()
