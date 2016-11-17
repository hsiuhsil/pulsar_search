import sys
import os.path

import h5py
import pyfits
import numpy as np
import time

def main():
    args = sys.argv[1:]
 
    space = 0.15
    ra_scopes = np.arange(324.5, 325.5, space)
    dec_scopes = np.arange(-2., 2., space)
    targets = []
    for ii in xrange(len(ra_scopes)):
        for jj in xrange(len(dec_scopes)):
            filename = 'RA_'+str(ra_scopes[ii])+'_DEC_'+str(dec_scopes[jj])
            targets.append([filename, ra_scopes[ii], dec_scopes[jj]])

    for filename in args:
        try:
            print filename
            hduls = pyfits.open(filename)
            search_targets(hduls, filename, targets, ra_scopes, dec_scopes)
        except (IOError, ValueError):
            print 'Skipped:'+ filename

def search_targets_index(hduls, ra_scopes, dec_scopes):
    targets_index = []    
    for k in xrange(len(hduls[1].data)):
        delta_ra = hduls[1].data[k]['RA_SUB'] - ra_scopes
        delta_dec = hduls[1].data[k]['DEC_SUB'] - dec_scopes
        distance = np.zeros(len(ra_scopes)*len(dec_scopes))
        for ii in xrange(len(ra_scopes)):
            distance[ii*len(dec_scopes):(ii+1)*len(dec_scopes)] = np.sqrt(delta_ra[ii]**2 + delta_dec**2)
        if len(np.where(distance <=0.15)[0]) == 0:
            targets_index.append([])
        elif len(np.where(distance <=0.15)[0]) != 0:    
            targets_index.append(np.where(distance <=0.15)[0])
    return targets_index


def search_targets(hduls, filename, targets, ra_scopes, dec_scopes):
    t0_time_line45 = time.time()
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
    if os.path.isfile('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/wz22hr/targets/RA_324.5_DEC_0.5h5') == True:
        for i in xrange(len(targets)):
            files[targets[i][0]] = h5py.File(targets[i][0] + 'h5',"r+")
    else:
        for i in xrange(len(targets)):
            this_file = h5py.File(targets[i][0] +  'h5',"w")
            for dataset_name in keys:
                if dataset_name == 'ABS_TIME':
                    first_data = hduls[1].data[0]['OFFS_SUB']
                    this_file.create_dataset(dataset_name, (1,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
                elif dataset_name == 'TBIN':
                    first_data = hduls[1].data[0]['TSUBINT']
                    this_file.create_dataset(dataset_name, (1,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
                elif dataset_name == 'RA_sets':
                    first_data = np.array([hduls[1].data[0]['RA_SUB'], hduls[1].data[0]['RA_SUB'], hduls[1].data[0]['RA_SUB']])
                    this_file.create_dataset(dataset_name, (1,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
                elif dataset_name == 'DEC_sets':
                    first_data = np.array([hduls[1].data[0]['DEC_SUB'], hduls[1].data[0]['DEC_SUB'], hduls[1].data[0]['DEC_SUB']])
                    this_file.create_dataset(dataset_name, (1,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
                else:
                    first_data = hduls[1].data[0][dataset_name]
                    this_file.create_dataset(dataset_name, (1,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
            files[targets[i][0]] = this_file

    print 'done create/write initial h5py files.'

    '''search targets location, and copy the data to h5df file if condition matched'''

    targets_index = search_targets_index(hduls, ra_scopes, dec_scopes)
    print 'get targets_index'
    print 'time of line 45-90: '+str(time.time() - t0_time_line45)+' sec'
    t0_time_file= time.time()
    for k in xrange(len(hduls[1].data)):
        print 'k= ' +str(k) +', targets: ', targets_index[k]
        t0_time_k = time.time()
        if len(targets_index[k]) == 0:
            print 'a nearby target is nan'
        else:
            record = hduls[1].data[k]
            for jj in xrange(len(targets_index[k])):
                for dataset_name in keys:
                    current_len = files[targets[targets_index[k][jj]][0]][dataset_name].shape[0]
                    files[targets[targets_index[k][jj]][0]][dataset_name].resize(current_len + 1, 0)
                    if dataset_name == 'ABS_TIME':
                        abs_time = hduls[0].header['STT_IMJD']*86400 + hduls[0].header['STT_SMJD'] + hduls[0].header['STT_OFFS'] + record['OFFS_SUB']
                        files[targets[targets_index[k][jj]][0]]['ABS_TIME'][current_len-1,...] = abs_time
                    elif dataset_name == 'TBIN':
                        files[targets[targets_index[k][jj]][0]]['TBIN'][current_len-1,...] = hduls[1].header['TBIN']
                    elif dataset_name == 'RA_sets':
                        files[targets[targets_index[k][jj]][0]]['RA_sets'][current_len-1,...] = RA_series[k]
                    elif dataset_name == 'DEC_sets':
                        files[targets[targets_index[k][jj]][0]]['DEC_sets'][current_len-1,...] = DEC_series[k]
                    else:
                        files[targets[targets_index[k][jj]][0]][dataset_name][current_len-1,...] = record[dataset_name]
            print 'time of this k: '+str(time.time() - t0_time_k)+' sec'
    print 'len. of this file:'+str(len(hduls[1].data))+', time of this file: '+str((time.time() - t0_time_file)/60.)+' min'

if __name__ == '__main__':
    main()


