import numpy as np
import h5py

'''Physical parameters'''
AU = 149597870700.0      # m
C = 299792458.0    # m/s

'''Parameters of pulsar J2139+00'''
RA = 324.8428583333333  # deg
DEC = 0.6959230555555556 # deg
T = 0.312470

'''Parameters of fitting'''
NPHASEBIN_wz = 200
NPHASEBIN = NPHASEBIN_wz
NPHASEBIN_1hr = 800
SCALE = np.float(NPHASEBIN_wz) / np.float(NPHASEBIN_1hr)

PHASE_DIFF_wz_1hr = 0.01735147
PHASE_DIFF_wz_1hr_err = 0.02428669

#PHASE_DIFF_wz_1hr = -6.20587204
#PHASE_DIFF_wz_1hr_err = 0.02407328


TIME0 = 55707.   # MJD pivot

fit_pars = [9.37426903e-08,  -3.34800232e+00,   3.26992424e+01,   2.14807618e-04,  -4.11457078e-04]
old_fit_pars =  [3.47609206e-07,  -3.37370520e+00,   4.74201487e+01,  -1.39689904e-05,  2.21425268e-05]
fit_time_start = -1500
fit_time_end = 36000

'''Files of pulsar J2139+00'''
this_file_wz = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_ANTF_delta_ra_dec_20170116/J2139+00_wzonlyh5', "r")
'''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
bin_number_wz = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_ANTF_delta_ra_dec_20170116/bin_number_2139_delta_new2.txt')
phase_amp_bin_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_ANTF_delta_ra_dec_20170116/phase_amp_bin_wz.npy')
phase_npy_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_ANTF_delta_ra_dec_20170116/phase_wz.npy')


this_file_1hr = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/J2139+00_57178h5', "r")
'''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
bin_number_1hr = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/bin_number_2139_57178.txt')
#phase_amp_bin_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/phase_amp_bin_57178.npy')
phase_amp_bin_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/phase_amp_bin_57178_fft.npy')
phase_npy_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/phase_1hr.npy')

