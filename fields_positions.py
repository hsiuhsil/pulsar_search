import sys
import os.path

import h5py
import pyfits
import numpy as np

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            hduls = pyfits.open(filename)
            search_positions(hduls, filename)
            print filename
        except (IOError, ValueError):
            print 'Skipped:'+ filename

def search_positions(hduls, filename):

    keys = ['POS']

    '''generate sets of RA_series and DEC_series from raw data'''
    positions = np.ndarray(shape=(len(hduls[1].data),2),dtype=float)
    positions[:,0] = hduls[1].data[:]['RA_SUB']
    positions[:,1] = hduls[1].data[:]['DEC_SUB']

    '''create dataset for each field position'''
    files = {} 
    if os.path.isfile('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/wz22hr/positionsh5') == True:
        files = h5py.File('positionsh5',"r+")
    else:
        this_file = h5py.File('positionsh5',"w")
        for dataset_name in keys:
            if dataset_name == 'POS':
                first_data = np.array([hduls[1].data[0]['RA_SUB'], hduls[1].data[0]['DEC_SUB']])
                this_file.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
        files = this_file


    '''searches fields positions, and records the information to h5df files'''
    for ii in xrange(len(positions)):
        for dataset_name in keys:
            current_len = files[dataset_name].shape[0]
            files[dataset_name].resize(current_len + 1, 0)
            if dataset_name == 'POS':
                files['POS'][current_len-1,...] = positions[ii]


if __name__ == '__main__':
    main()


