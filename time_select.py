import sys
import os
import os.path

import h5py
import numpy as np
from collections import Counter

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            time_select(filename)
        except (IOError, ValueError):
            print IOError

def time_select(filename):
    '''In this function, we select the data which has the maximal same MJD'''
    this_file = h5py.File(filename, "r")
    new_data = h5py.File('same_MJD_'+filename, "w")
    keys = this_file.keys()
    
    for dataset_name in keys:
        first_data = this_file[dataset_name][0]
        new_data.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)

    most_common_MJD = Counter(np.int64(this_file['BARY_TIME'][:])).most_common()[0][0]
    for ii in range(len(this_file['BARY_TIME'])):
        if np.int64(this_file['BARY_TIME'][ii]) == most_common_MJD:
            for dataset_name in keys:
                current_len = new_data[dataset_name].shape[0]
                new_data[dataset_name].resize(current_len + 1, 0)
                new_data[dataset_name][current_len-1,...] = this_file[dataset_name][ii]

if __name__ == '__main__':
    main()

