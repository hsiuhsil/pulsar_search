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

'''GBT parameters'''
Tsys = 25 # K
G = 2 # telescope gain, in the unit of K/Jy

'''Parameters of pulsar J2139+00'''
RA = 324.8428583333333  # deg
DEC = 0.6959230555555556 # deg
DM = 31.7262
T = 0.312470

'''Parameters of fitting'''
REFERNCE_FREQ = 900.
NPHASEBIN_wz = 200
NPHASEBIN = NPHASEBIN_wz
NPHASEBIN_1hr = 800
NCENTRALBINS = 100
NHARMONIC = 90
NCENTRALBINSMAIN = 160
SCALE = np.float(NPHASEBIN_wz) / np.float(NPHASEBIN_1hr)

len_5sec_stack = 40
len_5sec_all = 640

PHASE_DIFF_wz_1hr = -0.37306754
PHASE_DIFF_wz_1hr_err = 0.02899884

#PHASE_DIFF_wz_1hr = -6.20587204
#PHASE_DIFF_wz_1hr_err = 0.02407328


mid_mjd = 56436.5067837 # the average of the first and last MJD.
TIME0 = 55707   # MJD pivot


'''fit_pars is in the sequence of period derivative, period correction, phase offset, correction of RA, correction of DEC. The units are [bins/hr/hr, bins/hr, bins, radians, radians]'''

if TIME0 == 55707.:
    fit_time_start = -1500
    fit_time_end = 36000
    fit_pars = [  1.01403480e-07,  -3.34825452e+00,   3.20959712e+01,
         2.12231724e-04,  -4.05746692e-04]
    fit_pars_err = [  3.40364430e-09,   1.20394706e-04,   4.24574275e-01,
         1.13113130e-06,   1.28514882e-06]
elif TIME0 == mid_mjd:
    fit_time_start = -19500
    fit_time_end = 18000
    fit_pars = [  1.01407539e-07,  -3.34470374e+00,   4.13937006e+01,
         2.12230505e-04,  -4.05747569e-04]
    fit_pars_err = [  3.40364206e-09,   1.55647377e-06,   6.93423948e-01,
         1.13113145e-06,   1.28516943e-06]


'''fixed fit_pars[0:3] and get delta_ra and delta_dec for 2011 and 2015'''
'''array in the sequence of [delta_ra, delta_dec, delta_ra_err, delta_dec_err]'''
mean_mjd_11w = 55719.9836688
mean_mjd_15wp = 57143.1010341
delta_pos_11 = [0.000212243583, -0.000405771534, 4.31536460e-07, 9.78991105e-07]
delta_pos_15 = [0.000212200383 , -0.000405675702, 5.66098743e-07, 1.29119724e-06]

'''Files of pulsar J2139+00'''
this_file_wz = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/J2139+00_wzonlyh5', "r")
'''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
bin_number_wz = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/bin_number_2139_delta_new2.txt')
phase_amp_bin_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_amp_bin_wz_0_2_full_lik.npy')
#phase_amp_bin_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_amp_bin_wz_40epochs_6modes_90hars_lik.npy')
phase_npy_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_wz.npy')
random_res_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/random_fitting_phase_fft_wz_.npy')


this_file_1hr = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/J2139+00_57178h5', "r")
this_file_1hr_folding = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/J2139+00_1hr_foldingh5', "r")
'''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
bin_number_1hr = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/bin_number_2139_57178_40epochs.txt')
bin_number_1hr_5sec = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/bin_number_2139_57178_5sec.txt')
#phase_amp_bin_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/phase_amp_bin_57178.npy')
phase_amp_bin_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_amp_bin_57178_40epochs_6modes_90hars_lik.npy')
phase_amp_bin_1hr_5sec = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_amp_bin_57178_fft_5sec.npy')
phase_npy_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_57178_40epochs.npy')
phase_npy_1hr_5sec = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_57178_5sec.npy')
phase_npy_1hr_5sec_stack_40 = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_57178_5sec_stack_40.npy')
random_res_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/random_fitting_phase_fft_57178_.npy')
lik_57178_5sec_stack0 =  np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/lik_57178_5sec_stack0.npy')
sum_lik_57178_5sec_16 = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/sum_lik_57178_5sec_16.npy')


