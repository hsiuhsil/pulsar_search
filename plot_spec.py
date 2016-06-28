import sys
import os
import os.path
sys.path.append('/home/p/pen/hsiuhsil/burst_search/')
from burst_search import preprocess

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math


interval = 1
dedisperse = False
dm = 151.082

rebin = True
rebin_time = 16
rebin_freq = 16

sigma_threshold = 5
remove_period = 128

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            plot_spec(filename)
        except (IOError, ValueError):
            print IOError

def rebin_spec(input_data):
    output_data = np.zeros((input_data.shape[0]/rebin_time, input_data.shape[1]/rebin_freq))
    for ii in range(len(output_data)):
        for jj in range(len(output_data[1])):
            output_data[ii,jj]=np.mean(input_data[ii*rebin_time:(ii+1)*rebin_time, jj*rebin_freq:(jj+1)*rebin_freq])
    return output_data

def dedisperse_time(freq1, freq2):
    DM_CONST = 4148.808
    time = DM_CONST * dm * (freq1**-2 - freq2**-2)
    return time

def dedisperse_index(freq1, freq2, tbin):
    delta_t = tbin
    DM_CONST = 4148.808
    index = np.int64(np.around((DM_CONST * dm * (freq1**-2 - freq2**-2))/delta_t))
    return index

def dedisperse_spec(input_data, tbin):
    '''note: input_data.shape should be (ntime, nfreq), and dedisperse with shape of (nfreq, ntime)'''
    freq_time = input_data.T
    freq =  np.arange(900., 700., -200./freq_time.shape[0])
    time_move = dedisperse_index(np.amin(freq), freq, tbin)
    out_data = np.zeros((freq_time.shape[0], freq_time.shape[1]+np.max(time_move)))
    for ii in range(freq_time.shape[0]):
        for jj in range(freq_time.shape[1]):
            out_data[ii, jj + time_move[ii]] = freq_time[ii, jj]
    return out_data.T
    '''note: the out_data.shape is (nfreq, ntime), and return as out_data.T with shape of (ntime, nfreq)'''

def time_slope(input_data):
    print "start time_slope"
    slope_mode = np.arange(np.float(input_data.shape[1]))
    slope_mode -= np.mean(slope_mode)
    slope_mode /= math.sqrt(np.sum(slope_mode**2))
    slope_amplitude = np.sum(input_data * slope_mode[None,:], 0)
    input_data -= slope_amplitude * slope_mode[None,:]
    return input_data

def preprocessing(input_data):
    '''note: preprocess need data.shape = (nfreq, ntime)'''
    output_data = np.zeros(input_data.shape)
    data = input_data[:,0,:,0].T
    m = np.mean(data[:],axis=1)
    m[m==0]=1
    data = data / m[:,None] - 1
    preprocess.remove_noisy_freq(data, sigma_threshold)
    data = data-np.mean(data)
    data = time_slope(data)
    preprocess.remove_noisy_freq(data, sigma_threshold)
    preprocess.remove_periodic(data, remove_period)
    output_data[:,0,:,0] = data.T
    output_data[:,1:4,:,0] = input_data[:,1:4,:,0]
    return output_data

def plot_spec(filename):
    this_file = h5py.File(filename, "r")
    ntime = this_file['DATA'].shape[1]
    tbin = this_file['TBIN'][0]
    nfreq = this_file['DATA'].shape[3]

    data_preprocessed = np.zeros((interval*ntime, 4, nfreq, 1))

    for ii in range(0,interval):
#    for ii in range(len(this_file['TOPO_TIME'])):
        print 'ii = ' + str(ii)
        this_record_data = preprocessing(this_file['DATA'][ii])
        data_preprocessed[ii*ntime:(ii+1)*ntime,:] = this_record_data

    '''data is for the plot with shape of (nfreq, ntime)'''

    if rebin == True and dedisperse == True:
        data_first = dedisperse_spec(data_preprocessed[:, 0, :, 0], tbin)
        data = rebin_spec(data_first).T
    elif rebin == True and dedisperse == False:
        data_first = data_preprocessed[:,0,:,0]
        data = rebin_spec(data_first).T
    elif rebin == False and dedisperse == True:
        data_first = dedisperse_spec(data_preprocessed[:, 0, :, 0], tbin)
        data = data_first.T
    else:
        data_first = data_preprocessed[:, 0, :, 0]
        data = data_first.T

    '''get averaged freq  amplitude to time '''
    data2 = [0.]*data.shape[1]

    for ii in range(len(data2)):
        data2[ii] = np.mean(data[:,ii])
    
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.set_ylabel('Freq(MHz)', fontsize=20)
    ax1.set_xlabel('time (sec)', fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    cax1 = ax1.imshow(data, extent=[0, data_first.shape[0]*tbin, 700., 900.],aspect='auto', cmap=cm.Greys)
    cbar = plt.colorbar(cax1)
    ax2 = fig.add_subplot(212)
    ax2.set_ylabel('Mean Amp', fontsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    cax2 = plt.plot(data2)   
    plt.xticks([0, 0.25*len(data2), 0.5*len(data2), 0.75*len(data2), len(data2)], [str(0), str(round(0.25*data_first.shape[0]*tbin, 4)), str(round(0.5*data_first.shape[0]*tbin, 4)), str(round(0.75*data_first.shape[0]*tbin,4)), str(round(data_first.shape[0]*tbin,4))])
    plt.show()

if __name__ == '__main__':
    main()
