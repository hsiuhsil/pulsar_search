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
REFERNCE_FREQ = 900.
NPHASEBIN_wz = 200
NPHASEBIN = NPHASEBIN_wz
NPHASEBIN_1hr = 800
NCENTRALBINS = 100
NHARMONIC =90
NCENTRALBINSMAIN = 160
SCALE = np.float(NPHASEBIN_wz) / np.float(NPHASEBIN_1hr)

PHASE_DIFF_wz_1hr = -0.37306754
PHASE_DIFF_wz_1hr_err = 0.02899884

#PHASE_DIFF_wz_1hr = -6.20587204
#PHASE_DIFF_wz_1hr_err = 0.02407328


TIME0 = 55707.   # MJD pivot

'''fit_pars is in the sequence of acceleration, period correction, phase offset, correction of RA, correction of DEC. The units are [bins/hr/hr, bins/hr, bins, radians, radians]'''

fit_pars = [  1.01403480e-07,  -3.34825452e+00,   3.20959712e+01,
         2.12231724e-04,  -4.05746692e-04]
fit_pars_err = [  3.40364430e-09,   1.20394706e-04,   4.24574275e-01,
         1.13113130e-06,   1.28514882e-06]


#fit_pars = [1.03109306e-07,-3.34831520e+00,3.23012605e+01,2.11843445e-04,-4.06378828e-04]
#fit_pars_err = [3.81231775e-09,1.34847651e-04,4.75612907e-01,1.26815054e-06,1.43954560e-06]

#fit_pars = [1.03138228e-07, -3.34831607e+00, 3.23090695e+01, 2.11789545e-04, -4.06362698e-04]
#fit_pars_err = [3.82134432e-09, 1.35167094e-04, 4.76714143e-01, 1.27063424e-06, 1.44293425e-06]

#fit_pars = [ 9.49517623e-08, -3.34803706e+00, 3.19152109e+01, 2.14696364e-04,-4.06802430e-04]
#fit_pars_err = [ 4.16686453e-09, 1.47268785e-04, 5.06014988e-01, 1.35249739e-06, 1.55025727e-06]

old_fit_pars = [9.99185321e-08, -3.34821079e+00, 3.24597760e+01, 2.13575158e-04, -4.09169152e-04]
#fit_pars_err =[2.75365010e-09, 9.73423789e-05, 3.31501320e-01, 8.87387431e-07,   1.03577324e-06]

#fit_pars = [9.92822906e-08, -3.34819045e+00, 3.25195830e+01, 2.12595220e-04, -4.06821843e-04]
#fit_pars_err = [3.83875405e-09, 1.35650109e-04, 4.76117106e-01, 1.23588577e-06, 1.41247969e-06]
#fit_pars = [9.64290817e-08, -3.34808812e+00, 3.24820522e+01, 2.13144311e-04, -4.08205113e-04]
#fit_pars_err = [2.98682646e-09, 1.05550057e-04, 3.68840062e-01, 9.70029049e-07,   1.10959147e-06]
#fit_pars =  [  9.41132738e-08, -3.34801450e+00, 3.27048275e+01, 2.14716157e-04, -4.11270223e-04]
#fit_pars = [9.37426903e-08,  -3.34800232e+00,   3.26992424e+01,   2.14807618e-04,  -4.11457078e-04]
#old_fit_pars =  [3.47609206e-07,  -3.37370520e+00,   4.74201487e+01,  -1.39689904e-05,  2.21425268e-05]
'''fixed fit_pars[0:3] and get delta_ra and delta_dec for 2011 and 2015'''
'''array in the sequence of [delta_ra, delta_dec, delta_ra_err, delta_dec_err]'''
mean_mjd_11w = 55719.9836688
mean_mjd_15wp = 57143.1010341
delta_pos_11 = [0.000212243583, -0.000405771534, 4.31536460e-07, 9.78991105e-07]
delta_pos_15 = [0.000212200383 , -0.000405675702, 5.66098743e-07, 1.29119724e-06]


fit_time_start = -1500
fit_time_end = 36000

'''Files of pulsar J2139+00'''
this_file_wz = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/J2139+00_wzonlyh5', "r")
'''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
bin_number_wz = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/bin_number_2139_delta_new2.txt')
phase_amp_bin_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_amp_bin_wz_0_2_full_lik.npy')
#phase_amp_bin_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_amp_bin_wz_40epochs_6modes_90hars_lik.npy')
phase_npy_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_wz.npy')
random_res_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/random_fitting_phase_fft_wz_.npy')


this_file_1hr = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/J2139+00_57178h5', "r")
'''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
bin_number_1hr = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/bin_number_2139_57178_40epochs.txt')
bin_number_1hr_5sec = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/bin_number_2139_57178_5sec.txt')
#phase_amp_bin_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/phase_amp_bin_57178.npy')
phase_amp_bin_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_amp_bin_57178_40epochs_6modes_90hars_lik.npy')
phase_amp_bin_1hr_5sec = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_amp_bin_57178_fft_5sec.npy')
phase_npy_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_57178_40epochs.npy')
phase_npy_1hr_5sec = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_57178_5sec.npy')
random_res_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/random_fitting_phase_fft_57178_.npy')
