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

def preprocess(data):
    
    sigma_threshold = 5
    print sigma_threshold
    bad_chans = False
    var = np.empty(data.shape[3], dtype=np.float64)
    skew = np.empty(data.shape[3], dtype=np.float64)
    kurt = np.empty(data.shape[3], dtype=np.float64)
#    for ii in range(len(data)-1):
    for ii in range(1):
        data[ii] = data[ii]/np.mean(data[ii]) -1
        std = np.sqrt(np.var(data[ii]))
        for jj in range(data.shape[3]):
            var[jj] = np.var(data[ii,:,0,jj])
        print 'jj_done'
        bad_chans_var = abs(var - np.mean(var)) > sigma_threshold * np.std(var)
        bad_chans = np.logical_or(bad_chans, bad_chans_var)
        var[bad_chans] = np.mean(var)
        for kk in range(len(np.where(bad_chans==True)[0])):
            data[ii,:,0,np.where(bad_chans==True)[0][kk]]=0
        print 'kk_done'
        data[ii] = data[ii]-np.mean(data[ii])
        for ll in range(data.shape[3]):
            var[ll] = np.var(data[ii,:,0,ll])
            skew[ll] = np.mean((data[ii,:,0,ll] - np.mean(data[ii,:,0,ll])**3))
            kurt[ll] = np.mean((data[ii,:,0,ll] - np.mean(data[ii,:,0,ll])**4))
        print 'll_done'
        bad_chans_var = abs(var - np.mean(var)) > sigma_threshold * np.std(var)
        bad_chans_skew = abs(skew - np.mean(skew)) > sigma_threshold * np.std(skew)
        bad_chans_kurt = abs(kurt - np.mean(kurt)) > sigma_threshold * np.std(kurt)
        bad_chans = np.logical_or(bad_chans, bad_chans_var)
        bad_chans = np.logical_or(bad_chans, bad_chans_skew)
        bad_chans = np.logical_or(bad_chans, bad_chans_kurt)
        var[bad_chans] = np.mean(var)
        skew[bad_chans] = np.mean(skew)
        kurt[bad_chans] = np.mean(kurt)
        for mm in range(len(np.where(bad_chans==True)[0])):
            data[ii,:,0,np.where(bad_chans==True)[0][mm]]=0
        print 'mm_done'
#    return data
 

def folding(filename):
    this_file = h5py.File(filename, "r+")
    phase_bins = 100
    pulsar_period = 0.312470
    n_bins = 2048
    tbin = 0.001024
    preprocessing = True
 
    if preprocessing == True:
        data_preprocessing=preprocess(this_file['DATA'])    

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
