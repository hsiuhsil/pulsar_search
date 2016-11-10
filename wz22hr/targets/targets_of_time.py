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
            search_targets(hduls, filename)
        except (IOError, ValueError):
            print 'Skipped:'+ filename

def search_targets(hduls, filename):
    space = 0.15
    ra_scopes = np.arange(320, 330., space)
    dec_scopes = np.arange(-2., 2., space)
    targets = []
    for ii in xrange(len(ra_scopes)):
        for jj in xrange(len(dec_scopes)):
            filename = 'RA_'+str(ra_scopes[ii])+'_DEC_'+str(dec_scopes[jj])
            targets.append([filename, ra_scopes[ii], dec_scopes[jj]])

    """ targets[i][0] is a target name, targets[i][1] is its ra, and targets[i][2] is its dec."""

    keys = hduls[1].columns.names + ['ABS_TIME'] + ['TBIN'] + ['RA_sets'] + ['DEC_sets']    

    '''generate sets of RA_series and DEC_series from raw data'''
    RA_series = np.ndarray(shape=(len(hduls[1].data),3),dtype=float)
    RA_series[:,0] = np.append(np.nan, hduls[1].data[:-1]['RA_SUB'])
    RA_series[:,1] = hduls[1].data[:]['RA_SUB']
    RA_series[:,2] = np.append(hduls[1].data[1:]['RA_SUB'], np.nan)
    DEC_series = np.ndarray(shape=(len(hduls[1].data),3),dtype=float)
    DEC_series[:,0] = np.append(np.nan, hduls[1].data[:-1]['DEC_SUB'])
    DEC_series[:,1] = hduls[1].data[:]['DEC_SUB']
    DEC_series[:,2] = np.append(hduls[1].data[1:]['DEC_SUB'], np.nan)

    '''create dataset for each target location'''
    files = {} 
    if os.path.isfile('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/wz22hr/targets/RA_329.9_DEC_-0.05h5') == True:
        for i in xrange(len(targets)):
            files[targets[i][0]] = h5py.File(targets[i][0] + 'h5',"r+")
    else:
        for i in xrange(len(targets)):
            this_file = h5py.File(targets[i][0] +  'h5',"w")
            for dataset_name in keys:
                if dataset_name == 'ABS_TIME':
                    first_data = hduls[1].data[0]['OFFS_SUB']
                    this_file.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
                elif dataset_name == 'TBIN':
                    first_data = hduls[1].data[0]['TSUBINT']
                    this_file.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
                elif dataset_name == 'RA_sets':
                    first_data = np.array([hduls[1].data[0]['RA_SUB'], hduls[1].data[0]['RA_SUB'], hduls[1].data[0]['RA_SUB']])
                    this_file.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
                elif dataset_name == 'DEC_sets':
                    first_data = np.array([hduls[1].data[0]['DEC_SUB'], hduls[1].data[0]['DEC_SUB'], hduls[1].data[0]['DEC_SUB']])
                    this_file.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
                else:
                    first_data = hduls[1].data[0][dataset_name]
                    this_file.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
            files[targets[i][0]] = this_file


    '''search targets location, and copy the data to h5df file if condition matched'''

    for k in xrange(len(hduls[1].data)):
        for i in xrange(len(targets)):
            delta_ra = hduls[1].data[k]['RA_SUB'] - targets[i][1]
            delta_dec = hduls[1].data[k]['DEC_SUB'] - targets[i][2]
            scope = np.sqrt(delta_ra**2 + delta_dec**2)
            if scope <= 0.15:
                print filename
                print targets[i][0]
                for dataset_name in keys:
                    current_len = files[targets[i][0]][dataset_name].shape[0]
                    files[targets[i][0]][dataset_name].resize(current_len + 1, 0)
                    if dataset_name == 'ABS_TIME':
                        abs_time = hduls[0].header['STT_IMJD']*86400 + hduls[0].header['STT_SMJD'] + hduls[0].header['STT_OFFS'] + hduls[1].data[k]['OFFS_SUB']
                        files[targets[i][0]]['ABS_TIME'][current_len-1,...] = abs_time
                    elif dataset_name == 'TBIN':
                        files[targets[i][0]]['TBIN'][current_len-1,...] = hduls[1].header['TBIN']
                    elif dataset_name == 'RA_sets':
                        files[targets[i][0]]['RA_sets'][current_len-1,...] = RA_series[k]
                    elif dataset_name == 'DEC_sets':
                        files[targets[i][0]]['DEC_sets'][current_len-1,...] = DEC_series[k]
                    else:
                        files[targets[i][0]][dataset_name][current_len-1,...] = hduls[1].data[k][dataset_name]

if __name__ == '__main__':
    main()


