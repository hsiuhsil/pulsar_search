import sys
import os
import os.path

import h5py
import numpy as np
from collections import Counter

mins = 10
time_threshold = mins*60./86400. # convertion to unit of day 

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
#            global filename
            print filename
            time_select(filename)
        except (IOError, ValueError):
            print IOError



def time_select(filename):
    '''In this function, we select the data which has the maximal same MJD'''
    this_file = h5py.File(filename, "r")
    keys = this_file.keys()

    '''set time selecting threshold'''
    time_threshold_count = np.zeros((len(this_file['BARY_TIME'])), dtype=int)    
    for ii in xrange(len(this_file['BARY_TIME'])):
        print ii
        for jj in xrange(ii+1, len(this_file['BARY_TIME'])):
            if np.abs(this_file['BARY_TIME'][ii]-this_file['BARY_TIME'][jj]) <= time_threshold:
                time_threshold_count[ii] +=1
    print time_threshold_count

#    most_common_MJD = Counter(np.int64(this_file['BARY_TIME'][:])).most_common()[0][0]

#    for ii in range(len(this_file['BARY_TIME'])):
#        print 'ii = '+str(ii)
#        if np.int64(this_file['BARY_TIME'][ii]) == most_common_MJD:
#            for dataset_name in keys:
#                current_len = new_data[dataset_name].shape[0]
#                new_data[dataset_name].resize(current_len + 1, 0)
#                new_data[dataset_name][current_len-1,...] = this_file[dataset_name][ii]

#    most_common_time_ind = np.where(time_threshold_count==np.amax(time_threshold_count))[0][0]
    most_common_time_ind = np.where(time_threshold_count==10)[0]
    print most_common_time_ind

#    for ii in range(most_common_time_ind,most_common_time_ind+time_threshold_count[most_common_time_ind]):
    for ii in range(len(most_common_time_ind)):
        print 'ii = '+str(ii)
        global index
        index = most_common_time_ind[ii]
        new_data = h5py.File('same_MJD_'+str(index)+'_'+filename, "w")
        for dataset_name in keys:
            first_data = this_file[dataset_name][0]
            new_data.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
        for jj in range(index, index + 10 +1):
            for dataset_name in keys:
                current_len = new_data[dataset_name].shape[0]
                new_data[dataset_name].resize(current_len + 1, 0)
                new_data[dataset_name][current_len-1,...] = this_file[dataset_name][jj]




if __name__ == '__main__':
    main()

