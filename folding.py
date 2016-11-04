import sys
import os
import os.path
sys.path.append('/home/p/pen/hsiuhsil/burst_search/')
from burst_search import preprocess

import h5py
import numpy as np
import math

'''Define variables'''

initial, final = np.int(np.float(sys.argv[1])), np.int(np.float(sys.argv[2]))
assert initial <= final

#initial = 0
#final = 29
beam_model = True
RA_source = 161.68013
DEC_source = 3.06858

do_preprocess = True
sigma_threshold = 5
remove_period = 64
phase_bins = 100

#'''J2139 wz from folding'''
#delta_t = -5.8068471774736951e-06 -2.5301321417973068e-07 -4.2596495836139012e-07 + 2.9914132482230684e-07 
#pulsar_period = 0.312470 + delta_t
#p_dot = 4.8408891115682e-11
#p0_ind = np.int(0)

#'''J2139 57178 from folding'''
#delta_t = 6.2703317271681677e-06 - 5.7005305828936162e-07
#pulsar_period = 0.312470 -5.8068471774736951e-06 -2.5301321417973068e-07 -4.2596495836139012e-07 + 2.9914132482230684e-07 + delta_t

#'''J2139 57178 from tempo2'''
#delta_t = 0.0
#pulsar_period = 0.3124680322895416 + delta_t

#'''J0051'''
#delta_t = -3.2083428581595565e-07 +1.4583350248377532e-08 -7.0619503654410745e-09 +5.1170668232711475e-10 
#pulsar_period = 0.35473179890 + delta_t

# J1046: 
delta_t = 2e-08
pulsar_period = 0.326271446035 + delta_t

# J1038: 
#delta_t = 0.0
#pulsar_period = 0.02885155795131 + delta_t

# J0030: 
#delta_t = 0.0
#pulsar_period = 0.0048654532073692 + delta_t


def main():
    args = sys.argv[3:]
    for filename in args:
        try:
#            print filename
            folding(filename)
        except (IOError, ValueError):
            print IOError

def beam_weighting(input_data, ntime, nfreq, kk, total_kk, RA_sets, DEC_sets):
    '''note: beam weighting needs data.shape in (ntime, nfreq)'''
    output_data = np.zeros(input_data.shape)
    data = input_data[:,0,:,0].astype(np.float64).copy()
    freq = np.arange(900., 700., -200./nfreq)
    sigma_nu = 0.307 + (0.245-0.307)/(900-700)*(freq-700)
    tt = np.arange(-ntime/2.0 + 0.5, ntime/2.0 + 0.5)
    if kk ==0:
        if np.isnan(RA_sets[0]) == True and np.isnan(DEC_sets[0]) == True:
            RA_time = RA_sets[1] + (RA_sets[1]-RA_sets[2])/ntime*tt
            DEC_time = DEC_sets[1] + (DEC_sets[1]-DEC_sets[2])/ntime*tt
        else:
            print 'ERROR: RA_sets[0] or DEC_sets[0] is not np.NAN'
    elif kk == total_kk -1:
        if np.isnan(RA_sets[2]) == True and np.isnan(DEC_sets[2]) == True:
            RA_time = RA_sets[1] + (RA_sets[1]-RA_sets[0])/ntime*tt
            DEC_time = DEC_sets[1] + (DEC_sets[1]-DEC_sets[0])/ntime*tt
        else:
            print 'ERROR: RA_sets[2] or DEC_sets[2] is not np.NAN'
    else:
        RA_time_1 = RA_sets[1] + (RA_sets[1]-RA_sets[0])/ntime*tt[:ntime/2]
        RA_time_2 = RA_sets[1] + (RA_sets[1]-RA_sets[2])/ntime*tt[ntime/2:]
        RA_time = np.append(RA_time_1, RA_time_2)
        DEC_time_1 = DEC_sets[1] + (DEC_sets[1]-DEC_sets[0])/ntime*tt[:ntime/2]
        DEC_time_2 = DEC_sets[1] + (DEC_sets[1]-DEC_sets[2])/ntime*tt[ntime/2:]
        DEC_time = np.append(DEC_time_1, DEC_time_2)
    beam = np.zeros(data.shape, dtype=np.float64)
    for ii in range(len(beam)):
        beam[ii] = np.exp(-((RA_source - RA_time[ii])**2 + (DEC_source - DEC_time[ii])**2) / 2 / sigma_nu**2)
    data *= beam
    output_data[:,0,:,0] = data
    output_data[:,1:4,:,0] = input_data[:,1:4,:,0]
    return beam, output_data

def time_slope(input_data):
    print "start time_slope"
    slope_mode = np.arange(np.float(input_data.shape[1]))
    slope_mode -= np.mean(slope_mode)
    slope_mode /= math.sqrt(np.sum(slope_mode**2))
    slope_amplitude = np.sum(input_data * slope_mode[None,:], 1)
    input_data -= slope_amplitude[:,None] * slope_mode
    return input_data

def preprocessing(input_data):
    '''note: preprocess need data.shape = (nfreq, ntime)'''
    output_data = np.zeros(input_data.shape)
    data = input_data[:,0,:,0].T.astype(np.float64).copy()
    preprocess.remove_periodic(data, remove_period)
    m = np.mean(data[:],axis=1)
    m[m==0]=1
    data = data / m[:,None] - 1
    preprocess.remove_noisy_freq(data, sigma_threshold)
    data = data-np.mean(data)
    data = time_slope(data)
    preprocess.remove_bad_times(data, sigma_threshold)
    data -= np.mean(data, 0)[None,:]
    preprocess.remove_noisy_freq(data, sigma_threshold)
    output_data[:,0,:,0] = data.T
    output_data[:,1:4,:,0] = input_data[:,1:4,:,0]
    return output_data

def folding(filename):

    this_file = h5py.File(filename, "r+")   
    nfreq = this_file['DATA'].shape[3]
    ntime = this_file['DATA'].shape[1]
    tbin = this_file['TBIN'][0]

    first_data = this_file['DATA'][0][0]
    data_folding = np.zeros((phase_bins,) + first_data.shape)
    beam_ave = np.zeros(ntime, dtype=np.float64)
    beam_ave_phase = np.zeros(phase_bins, dtype=np.float64)
    
    '''collecting data which satisfies the folding condition'''
    same_modulo_num = np.zeros((phase_bins,), dtype=np.int)
    for ii in range(initial, final+1):
        print 'ii = ' + str(ii)  
        sample_BAT = this_file['BARY_TIME'][ii]*86400 + np.arange(-ntime/2.0 + 0.5, ntime/2.0 + 0.5)*tbin
        modulo_num = np.int64(np.around((sample_BAT % pulsar_period)/(pulsar_period / phase_bins)))
        print 'modulo_num done'
        for jj in range(len(modulo_num)):
            if modulo_num[jj] == phase_bins:
                modulo_num[jj] = 0
        for kk in range(len(same_modulo_num)):
            same_modulo_num[kk] += np.count_nonzero(modulo_num == kk)      
        if do_preprocess == True:
            this_file_data = preprocessing(this_file['DATA'][ii])
            print 'preprocess done'
        else:
            this_file_data = this_file['DATA'][ii]
        if beam_model == True:
            beam, this_record_data = beam_weighting(this_file_data, ntime, nfreq, ii, len(this_file['DATA']), this_file['RA_sets'][ii], this_file['DEC_sets'][ii])
            print 'beam model done'            
        else:
            this_record_data = this_file_data

        for ll in range(len(modulo_num)):
            data_folding[modulo_num[ll],...] += this_record_data[ll]
        print 'data_folding done'
        for ii in range(len(beam_ave)):
            beam_ave[ii] += np.average(beam[ii])
        print 'beam_ave done'        

    for ii in range(len(beam_ave)):
        beam_ave_phase[modulo_num[ii]] += beam_ave[ii]/np.sum(beam_ave)
    print 'beam_ave_phase done'

    for mm in range(len(same_modulo_num)):
        if same_modulo_num[mm] != 0:
            data_folding[mm,...] = data_folding[mm,...]*beam_ave_phase[mm]
#            data_folding[mm,...] = data_folding[mm,...]/same_modulo_num[mm]

    this_file.create_dataset('DATA_FOLDING', data_folding.shape, maxshape = data_folding.shape, dtype=data_folding.dtype, chunks=True)
    this_file['DATA_FOLDING'][...]=data_folding[...]


    '''Data folding for topocentric time'''
    data_folding_topo = np.zeros((phase_bins,) + first_data.shape)
    beam_ave_topo = np.zeros(ntime, dtype=np.float64)
    beam_ave_topo_phase = np.zeros(phase_bins, dtype=np.float64)

    '''collecting data which satisfies the folding condition'''
    same_modulo_num_topo = np.zeros((phase_bins,), dtype=np.int)
    for ii in range(initial, final+1):
        print 'ii = ' + str(ii)
        sample_TOPO = this_file['TOPO_TIME'][ii]*86400 + np.arange(-ntime/2.0 + 0.5, ntime/2.0 + 0.5)*tbin
        modulo_num_topo = np.int64(np.around((sample_TOPO % pulsar_period)/(pulsar_period/phase_bins)))
        print 'modulo_num done'
        for jj in range(len(modulo_num_topo)):
            if modulo_num_topo[jj] == phase_bins:
                modulo_num_topo[jj] = 0
        for kk in range(len(same_modulo_num_topo)):
            same_modulo_num_topo[kk] += np.count_nonzero(modulo_num_topo == kk)
        if do_preprocess == True:
            this_file_data = preprocessing(this_file['DATA'][ii])
            print 'preprocess done'
        else:
            this_file_data = this_file['DATA'][ii]
        if beam_model == True:
            beam_topo, this_record_data_topo = beam_weighting(this_file_data, ntime, nfreq, ii, len(this_file['DATA']), this_file['RA_sets'][ii], this_file['DEC_sets'][ii])
            print 'beam model done'
        else:
            this_record_data_topo = this_file_data

        for ll in range(len(modulo_num_topo)):
            data_folding_topo[modulo_num_topo[ll],...] += this_record_data_topo[ll]
        print 'data_folding done'
        for ii in range(len(beam_ave_topo)):
            beam_ave_topo[ii] += np.average(beam_topo[ii])
        print 'beam_ave done'

    for ii in range(len(beam_ave_topo)):
        beam_ave_topo_phase[modulo_num_topo[ii]] += beam_ave_topo[ii]/np.sum(beam_ave_topo)
    print 'beam_ave_phase done'

    for mm in range(len(same_modulo_num_topo)):
        if same_modulo_num_topo[mm] != 0:
            data_folding_topo[mm,...] = data_folding_topo[mm,...]*beam_ave_topo_phase[mm]

    this_file.create_dataset('DATA_FOLDING_TOPO', data_folding_topo.shape, maxshape = data_folding_topo.shape, dtype=data_folding_topo.dtype, chunks=True)
    this_file['DATA_FOLDING_TOPO'][...]=data_folding_topo[...]



if __name__ == '__main__':
    main()
