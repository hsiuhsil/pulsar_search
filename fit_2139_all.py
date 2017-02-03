import sys
import os
import os.path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import leastsq
from astropy.time import Time

def main():
    try:
        plot_bary_diff()
    except None:
        print IOError


RA = 324.8428583333333  # deg
DEC = 0.6959230555555556 # deg
AU = 149597870700.0      # m
C = 299792458.0    # m/s
NPHASEBIN_wz = 200
NPHASEBIN = NPHASEBIN_wz
NPHASEBIN_1hr = 800
SCALE = np.float(NPHASEBIN_wz) / np.float(NPHASEBIN_1hr)
T = 0.312470
PHASE_DIFF_wz_1hr = 6.20587204
PHASE_DIFF_wz_1hr_err = 0.02407328

TIME0 = 55707.   # MJD pivot


def transform_time(time_mjd):
    return (time_mjd - TIME0) * 24

def untrans_time(delta_time_hrs):
    return delta_time_hrs / 24. + TIME0

def timing_model_1(parameters, time_mjd, dBATdra, dBATddec, wz_range):
    time = transform_time(time_mjd)

    out1 = parameters[0] * time**2
    out1 += parameters[1] * time
    out1 += parameters[2]
    for ii in xrange(len(time)):
        if ii < wz_range:
            NPHASEBIN = NPHASEBIN_wz
#            print 'NPHASEBIN, <236', NPHASEBIN
            out1[ii] +=  (NPHASEBIN / T) * (dBATdra[ii] * 86400 * 180 / np.pi * parameters[3] + dBATddec[ii] * 86400 * 180 / np.pi * parameters[4])
            out1[ii] = out1[ii] % NPHASEBIN
#            print 'out1 < wz_range', out1
        elif ii >= wz_range:
#            print 'NPHASEBIN, >=236', NPHASEBIN
            NPHASEBIN = NPHASEBIN_1hr
            out1[ii] +=  (NPHASEBIN / T) * (dBATdra[ii] * 86400 * 180 / np.pi * parameters[3] + dBATddec[ii] * 86400 * 180 / np.pi * parameters[4])
            out1[ii] = out1[ii] % NPHASEBIN
            out1[ii] = out1[ii] * SCALE + PHASE_DIFF_wz_1hr
#            print 'out1 >= wz_range', out1 

    return out1

def residuals_1(parameters, time_mjd, dBATdra, dBATddec, phase_data, phase_data_err, wz_range):
    model_1 = timing_model_1(parameters, time_mjd, dBATdra, dBATddec, wz_range)
    
    res_1 = phase_data - model_1
    res_1 = (res_1 + NPHASEBIN / 2.) % NPHASEBIN - NPHASEBIN / 2.
#    print 'NPHASEBIN in res', NPHASEBIN
    res_1 = res_1 / phase_data_err
    return res_1

def models_plot(time_mjd, dBATdra, dBATddec, filename, time_range=None):

    if time_range is None:
        time_range = (time_mjd[0], time_mjd[-1])

    equinox_date = ['2010-03-20T17:32:00','2011-03-20T23:21:00','2012-03-20T05:14:00','2013-03-20T11:02:00','2014-03-20T16:57:00','2015-03-20T22:45:00','2016-03-20T04:30:00','2017-03-20T10:28:00']
    t = Time(equinox_date, format='isot', scale='utc')
    equinox_mjd = t.mjd
    theta = np.zeros(len(time_mjd))
    for ii in range(len(theta)):
        theta[ii] = (time_mjd[ii] - equinox_mjd[np.argmin(np.absolute(time_mjd[ii] - equinox_mjd))]) /  365.259636*2*np.pi

    old_RA = (NPHASEBIN / T) * (-AU / C * np.sin(-RA*np.pi/180 + theta))
    new_RA = (NPHASEBIN / T) * dBATdra * 86400 * 180 / np.pi # unit is bins/rad
    new_DEC = (NPHASEBIN / T) * dBATddec * 86400 * 180 / np.pi # unit is bins/rad

    plt.subplot(2,1,1)
    plt.plot(transform_time(time_mjd), old_RA, 'r--')
    plt.plot(transform_time(time_mjd), new_RA, 'bo')
    plt.xlabel('Bary diff (hours)', fontsize=14)
    plt.ylabel('Quantities', fontsize=14)
    plt.xlim(transform_time(time_range[0]), transform_time(time_range[1]))

    plt.subplot(2,1,2)
    plt.plot(transform_time(time_mjd), old_RA, 'r--')
    plt.plot(transform_time(time_mjd), new_DEC, 'g^')
    plt.xlabel('Bary diff (hours)', fontsize=14)
    plt.ylabel('Quantities', fontsize=14)
    plt.xlim(transform_time(time_range[0]), transform_time(time_range[1]))

    plt.savefig(filename)

def make_save_plot(parameters, model, res, time, dBATdra, dBATddec, data, data_err, wz_range, filename, time_range=None):
    
    if time_range is None:
        time_range = (time[0], time[-1])

    num_points = len(time)
    model_time = np.linspace(time_range[0], time_range[1], num_points)

    model = timing_model_1(parameters, model_time, dBATdra, dBATddec, wz_range)
    res = residuals_1(parameters, time, dBATdra, dBATddec, data, data_err, wz_range)

    plt.subplot(2,1,1)
    plt.plot(transform_time(time), data, 'bo')
    plt.plot(transform_time(model_time), model, 'r--')
    plt.xlabel('Bary diff (hours)', fontsize=14)
    plt.ylabel('Max Phase Bins Number', fontsize=14)
    plt.xlim(transform_time(time_range[0]), transform_time(time_range[1]))

    plt.subplot(2,1,2)
    plt.plot(transform_time(time), res, 'bo')
    plt.xlim(transform_time(time_range[0]), transform_time(time_range[1]))
    if len(time_range) > 2:
        plt.ylim(time_range[2], time_range[3])
    plt.xlabel('Bary diff (hours)', fontsize=14)
    plt.ylabel('Phase bin residuals', fontsize=14)

    plt.savefig(filename)

def time_pattern(this_file, bin_number, phase_amp_bin, NPHASEBIN):
    time_mjd = np.zeros(len(bin_number))
    dBATdra = np.zeros(len(bin_number))
    dBATddec = np.zeros(len(bin_number))
    for ii in range(len(bin_number)):
        time_mjd[ii] = (this_file['BARY_TIME'][bin_number[ii][0]] + this_file['BARY_TIME'][bin_number[ii][1]])/2.
        dBATdra[ii] = (this_file['dBATdra'][bin_number[ii][0]] + this_file['dBATdra'][bin_number[ii][1]])/2.
        dBATddec[ii] = (this_file['dBATddec'][bin_number[ii][0]] + this_file['dBATddec'][bin_number[ii][1]])/2.

    if NPHASEBIN == NPHASEBIN_wz:
        phase_data = phase_amp_bin[:,1]
        phase_data_err = phase_amp_bin[:,3]
    elif NPHASEBIN == NPHASEBIN_1hr:
        'As WZ to be the template, need to rescale the phase bin and then add the difference between two templates.'
        phase_data = phase_amp_bin[:,1] * SCALE + PHASE_DIFF_wz_1hr
        phase_data_err = np.sqrt(phase_amp_bin[:,3]**2 * SCALE**2 + PHASE_DIFF_wz_1hr_err**2)
    else:
        print 'NPHASEBIN error'


    return time_mjd, dBATdra, dBATddec, phase_data, phase_data_err


def plot_bary_diff():

    this_file_wz = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_ANTF_delta_ra_dec_20170116/J2139+00_wzonlyh5', "r")
    '''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
    bin_number_wz = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_ANTF_delta_ra_dec_20170116/bin_number_2139_delta_new2.txt')
    phase_amp_bin_wz = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_ANTF_delta_ra_dec_20170116/phase_amp_bin_wz.npy')


    this_file_1hr = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/J2139+00_57178h5', "r")
    '''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
    bin_number_1hr = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/bin_number_2139_57178.txt')
    phase_amp_bin_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/phase_amp_bin_57178.npy')

    time_mjd_wz, dBATdra_wz, dBATddec_wz, phase_data_wz, phase_data_err_wz = time_pattern(this_file_wz, bin_number_wz, phase_amp_bin_wz, NPHASEBIN_wz)
    time_mjd_1hr, dBATdra_1hr, dBATddec_1hr, phase_data_1hr, phase_data_err_1hr = time_pattern(this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, NPHASEBIN_1hr)

#    print 'phase_data_1hr',phase_data_1hr
#    print 'phase_data_err_1hr',phase_data_err_1hr

    time_mjd = np.concatenate((time_mjd_wz, time_mjd_1hr))
    dBATdra = np.concatenate((dBATdra_wz, dBATdra_1hr))
    dBATddec = np.concatenate((dBATddec_wz, dBATddec_1hr))
    phase_data = np.concatenate((phase_data_wz, phase_data_1hr))
    phase_data_err = np.concatenate((phase_data_err_wz, phase_data_err_1hr))


    # What data to fit.
    fit_range = (untrans_time(-1500), untrans_time(36000))
    sl = np.logical_and(time_mjd > fit_range[0], time_mjd < fit_range[1])
    n = 5
    pars_init_1 = [3.93934677e-07,  -3.37523445e+00 + n*3e-3,   4.82981320e+01, 1e-05, 1e-05 ]
#    pars_init_1 = [2.e-07  , -3.353   ,3.4e+01  , 1e-05, 1e-05]
    p = [ -7.35e-7, -1.63079989e+00,  -1.13336394e+02,   7.82201828e-04] 

    print 'timing_model_1', timing_model_1(pars_init_1, time_mjd, dBATdra, dBATddec, len(time_mjd_wz))


    '''parameters for new RA correction'''
    fit_pars_1, pcov_1, infodict_1, errmsg_1, success_1 = leastsq(residuals_1, pars_init_1,
            args=(time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl], len(time_mjd_wz)), xtol = 1e-6, ftol=1e-6, full_output=1)
    print 'n value: ', n
    print "Fit parameters: ", fit_pars_1
    print "sucess?:", success_1
    print "Chi-squared: ", np.sum(residuals_1(fit_pars_1, time_mjd[sl], dBATdra[sl], dBATddec[sl],
                         phase_data[sl], phase_data_err[sl], len(time_mjd_wz))**2), "DOF: ", len(phase_data[sl])-len(pars_init_1)


    if (len(phase_data) > len(pars_init_1)) and pcov_1 is not None:
        s_sq_1 = (residuals_1(fit_pars_1, time_mjd[sl], dBATdra[sl], dBATddec[sl],
            phase_data[sl], phase_data_err[sl], len(time_mjd_wz))**2).sum()/(len(phase_data[sl])-len(pars_init_1))
        pcov_1 = pcov_1 * s_sq_1
    else:
        pcov_1 = np.inf

    error_1 = []
    for i in range(len(fit_pars_1)):
        try:
            error_1.append(np.absolute(pcov_1[i][i])**0.5)
        except:
            error_1.append( 0.00 )
    pfit_leastsq_1 = fit_pars_1
    perr_leastsq_1 = np.array(error_1)

#    print("\nFit paramters and parameter errors from lestsq method :")
#    print("pfit = ", pfit_leastsq_1)
#    print("perr = ", perr_leastsq_1)


    models_plot(time_mjd, dBATdra, dBATddec,
                   'models1.png', (untrans_time(-500), untrans_time(1800)))
    plt.figure()
    models_plot(time_mjd, dBATdra, dBATddec,
                   'models2.png', (untrans_time(33300), untrans_time(34500)))

    plt.figure()
    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl], len(time_mjd_wz), 'new_phase_fit_J2139_init.png')
    plt.figure()
    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl], len(time_mjd_wz), 'new_phase_fit_J2139.png')
    plt.figure()
    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl], len(time_mjd_wz), 'new_phase_fit_J2139_zoom.png',
            (untrans_time(-500), untrans_time(1500)))
    plt.figure()
    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl], len(time_mjd_wz), 'new_phase_fit_J2139_zoom2.png',
            (untrans_time(33400), untrans_time(34400)))
    plt.figure()
    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd, dBATdra, dBATddec,  phase_data, phase_data_err[sl], len(time_mjd_wz), 'new_phase_fit_J2139_all.png')


if __name__ == '__main__':
    main()
