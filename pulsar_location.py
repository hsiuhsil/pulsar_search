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
            search_pulsar(hduls, filename)
        except (IOError, ValueError):
            print 'Skipped:'+ filename

def search_pulsar(hduls, filename):

    pulsar = [['J0030+0451', 7.61428, 4.86103],
              ['J0051+0423', 12.87542, 4.38028],
              ['J1038+0032', 159.61222, 0.54544],
              ['J1046+0304', 161.68013, 3.06858],
              ['J1501-0046', 225.43732, -0.77320],
              ['J1518+0204A', 229.63882, 2.09099],
              ['J1518+0204B', 229.63107, 2.08763],
              ['J1518+0204C', 229.63662, 2.07995],
              ['J1518+0204D', 229.64167, 2.08278],
              ['J1518+0204E', 229.64167, 2.08278],
              ['J2139+00', 324.92500, 0.60000],
              ['J2222-0137', 335.52487, -1.62103]
             ]
    """ pulsar[i][0] is pulsar's name, pulsar[i][1] is its ra, and  pulsar[i][2] is its dec."""

    keys = hduls[1].columns.names + ['ABS_TIME'] + ['TBIN']    

    '''create dataset for each pulsar location'''
    files = {} 
    if os.path.isfile('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_search/J2222-0137h5') == True:
        for i in xrange(len(pulsar)):
            files[pulsar[i][0]] = h5py.File(pulsar[i][0] + 'h5',"r+")
    else:
        for i in xrange(len(pulsar)):
            this_file = h5py.File(pulsar[i][0] + 'h5',"w")
            for dataset_name in keys:
                if dataset_name == 'ABS_TIME':
                    first_data = hduls[1].data[0]['OFFS_SUB']
                    this_file.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
                elif dataset_name == 'TBIN':
                    first_data = hduls[1].data[0]['TSUBINT']
                    this_file.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
                else:
                    first_data = hduls[1].data[0][dataset_name]
                    this_file.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
            files[pulsar[i][0]] = this_file


    '''search pulsar location, and copy the data to h5df file if condition matched'''

    for k in xrange(len(hduls[1].data)):
        for i in xrange(len(pulsar)):
            delta_ra = hduls[1].data[k]['RA_SUB'] - pulsar[i][1]
            delta_dec = hduls[1].data[k]['DEC_SUB'] - pulsar[i][2]
            scope = np.sqrt(delta_ra**2 + delta_dec**2)
            if scope <= 0.15:
                print filename
                print pulsar[i][0]
                print k
                for dataset_name in keys:
                    current_len = files[pulsar[i][0]][dataset_name].shape[0]
                    files[pulsar[i][0]][dataset_name].resize(current_len + 1, 0)
                    if dataset_name == 'ABS_TIME':
                        abs_time = hduls[0].header['STT_IMJD']*86400 + hduls[0].header['STT_SMJD'] + hduls[0].header['STT_OFFS'] + hduls[1].data[k]['OFFS_SUB']
#                        print abs_time
                        files[pulsar[i][0]]['ABS_TIME'][current_len-1,...] = abs_time
                    elif dataset_name == 'TBIN':
                        files[pulsar[i][0]]['TBIN'][current_len-1,...] = hduls[1].header['TBIN']
#                        print hduls[1].header['TBIN']
                    else:
                        files[pulsar[i][0]][dataset_name][current_len-1,...] = hduls[1].data[k][dataset_name]

if __name__ == '__main__':
    main()


