import numpy as np
import h5py
import bary_time

'''Physical parameters'''
AU = 149597870700.0      # m
C = 299792458.0    # m/s

'''Parameters of pulsar J0051+0423'''
#RA = bary_time.convHMS(str('00:51:30.1'))
#DEC = bary_time.convDMS(str('+04:22:49'))
#DM = 13.9
#T = 0.35473179890 


'''Parameters of pulsar J2139+00'''
RA = 324.8428583333333  # deg
DEC = 0.6959230555555556 # deg
DM = 31.7262
T = 0.312470

'''Parameters of fitting'''
NPHASEBIN_wz = 200
NPHASEBIN = NPHASEBIN_wz
NPHASEBIN_1hr = 800
SCALE = np.float(NPHASEBIN_wz) / np.float(NPHASEBIN_1hr)

PHASE_DIFF_wz_1hr = -0.37306754
PHASE_DIFF_wz_1hr_err = 0.02899884

#PHASE_DIFF_wz_1hr = -6.20587204
#PHASE_DIFF_wz_1hr_err = 0.02407328


TIME0 = 55707.   # MJD pivot

'''fit_pars is in the sequence of acceleration, period correction, phase offset, correction of RA, correction of DEC. The units are [bins/hr/hr, bins/hr, bins, radians, radians]'''

fit_pars = [9.64290817e-08, -3.34808812e+00, 3.24820522e+01, 2.13144311e-04, -4.08205113e-04]
fit_pars_err = [2.98682646e-09, 1.05550057e-04, 3.68840062e-01, 9.70029049e-07,   1.10959147e-06]
#fit_pars =  [  9.41132738e-08, -3.34801450e+00, 3.27048275e+01, 2.14716157e-04, -4.11270223e-04]
#fit_pars = [9.37426903e-08,  -3.34800232e+00,   3.26992424e+01,   2.14807618e-04,  -4.11457078e-04]
old_fit_pars =  [3.47609206e-07,  -3.37370520e+00,   4.74201487e+01,  -1.39689904e-05,  2.21425268e-05]
fit_time_start = -1500
fit_time_end = 36000

'''Files of pulsar J2139+00'''
this_file_wz = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/J2139+00_wzonlyh5', "r")
'''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
bin_number_wz = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/bin_number_2139_delta_new2.txt')
phase_amp_bin_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_amp_bin_wz_fft.npy')
phase_npy_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_wz.npy')
random_res_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/random_fitting_phase_fft_wz_.npy')


this_file_1hr = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/J2139+00_57178h5', "r")
'''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
bin_number_1hr = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/bin_number_2139_57178.txt')
#phase_amp_bin_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/phase_amp_bin_57178.npy')
phase_amp_bin_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_amp_bin_57178_fft.npy')
phase_npy_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_1hr.npy')
random_res_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/random_fitting_phase_fft_57178_.npy')

