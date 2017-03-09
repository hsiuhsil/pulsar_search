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
import pars

def main():
    try:
        plot_bary_diff()
    except None:
        print IOError


RA = pars.RA
DEC = pars.DEC
AU = pars.AU
C = pars.C
NPHASEBIN_wz = pars.NPHASEBIN_wz
NPHASEBIN = pars.NPHASEBIN
NPHASEBIN_1hr = pars.NPHASEBIN_1hr
SCALE = pars.SCALE
T = pars.T
#PHASE_DIFF_wz_1hr = pars.PHASE_DIFF_wz_1hr
#PHASE_DIFF_wz_1hr_err = pars.PHASE_DIFF_wz_1hr_err

TIME0 = pars.TIME0

def transform_time(time_mjd):
    return (time_mjd - TIME0) * 24

def untrans_time(delta_time_hrs):
    return delta_time_hrs / 24. + TIME0

def timing_model_1(parameters, time_mjd, dBATdra, dBATddec, NPHASEBIN=None, RESCALE=None):

    time = transform_time(time_mjd)

    out1 = parameters[0] * time**2
    out1 += parameters[1] * time
    out1 += parameters[2]

    if (NPHASEBIN == None) or (NPHASEBIN == NPHASEBIN_wz):
        NPHASEBIN = NPHASEBIN_wz
    elif NPHASEBIN == NPHASEBIN_1hr:
        NPHASEBIN = NPHASEBIN_wz
    
    out1 += (NPHASEBIN / T) * (-1) *(dBATdra * 86400 * 180 / np.pi * parameters[3] + dBATddec * 86400 * 180 / np.pi * parameters[4])

    out1 = out1 % NPHASEBIN 

    if (time_mjd[0] > 57178) and RESCALE == None :
        out1 = out1 / SCALE

    for ii in xrange(len(out1)):
        if out1[ii] == NPHASEBIN:
            out1[ii] == 0

    return out1

def residuals_1(parameters, time_mjd, dBATdra, dBATddec, phase_data, phase_data_err, NPHASEBIN=None, RESCALE=None):

    if (NPHASEBIN == None) or (NPHASEBIN == NPHASEBIN_wz):
        NPHASEBIN = NPHASEBIN_wz
    elif NPHASEBIN == NPHASEBIN_1hr:
        NPHASEBIN = NPHASEBIN_1hr
    model_1 = timing_model_1(parameters, time_mjd, dBATdra, dBATddec, NPHASEBIN, RESCALE)
#    print 'timing_model_1',model_1
    res_1 = phase_data - model_1
    res_1 = (res_1 + NPHASEBIN / 2.) % NPHASEBIN - NPHASEBIN / 2.
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

def make_save_plot(parameters, model, res, time, dBATdra, dBATddec, data, data_err, random_res_std, wz_range, filename, time_range=None):
    
    if time_range is None:
        time_range = (time[0], time[-1])

    num_points = len(time)
    model_time = np.linspace(time_range[0], time_range[1], num_points)

    model = timing_model_1(parameters, model_time, dBATdra, dBATddec)
    res = residuals_1(parameters, time, dBATdra, dBATddec, data, data_err) 
    '''plot res, rather than res/errors'''
    res *= data_err 
    '''plot in phase, rather than phase bin'''
    res /= NPHASEBIN
#    print 'res= ', res
    yerr = data_err / NPHASEBIN

    markersize = 2.0

#    plt.subplot(2,1,1)
#    plt.plot(transform_time(time), data, 'bo', markersize=markersize)
#    plt.errorbar(transform_time(time), data, yerr= data_err, fmt=None, color='b')
#    plt.plot(transform_time(model_time), model, 'r--')
#    plt.xlabel('Bary diff (hours)', fontsize=14)
#    plt.ylabel('Max Phase Bins Number', fontsize=14)
#    plt.xlim(transform_time(time_range[0]), transform_time(time_range[1]))

#    plt.subplot(2,1,2)
    plt.plot(time, res, 'bo', markersize=markersize)
    plt.errorbar(time, res , yerr= yerr, fmt=None, ecolor='b')
#    plt.plot(transform_time(time), random_res_std, 'ks', markersize=markersize)
#    plt.errorbar(transform_time(time), res , yerr= random_res_std, fmt=None, ecolor='r')
    plt.xlim(time_range[0], time_range[1])
    if len(time_range) > 2:
        plt.ylim(time_range[2], time_range[3])
    plt.xlabel('MJD ', fontsize=14)
    plt.ylabel('Phase residuals', fontsize=14)
    plt.ticklabel_format(useOffset=False)

    plt.savefig(filename)

def make_save_plot_pointing(parameters, model, res, time, dBATdra, dBATddec, data, data_err, random_res_std, wz_range, filename, time_range=None):

    time_range = np.linspace(time[wz_range], time[-1], 64)
#    print time_range
    timing_model_1_pointing = timing_model_1(parameters, time, dBATdra, dBATddec)[236:]

    phase_data_pointing =  data[236:]
    phase_data_pointing_err =  data_err[236:]
    residuals_1_pointing = residuals_1(parameters, time, dBATdra, dBATddec,  data, data_err)[236:]
    random_res_std = random_res_std[236:]

    '''plot res, rather than res/errors'''
    residuals_1_pointing *= phase_data_pointing_err
    '''plot in phase, rather than phase bin'''
    residuals_1_pointing /= NPHASEBIN
    phase_data_pointing_err /= NPHASEBIN

#    print 'phase_data_pointing_err',phase_data_pointing_err
    markersize = 3.0

#    plt.subplot(2,1,1)
#    plt.plot(time_range, phase_data_pointing, 'bo', markersize=markersize)
#    plt.errorbar(time_range, phase_data_pointing, yerr = phase_data_pointing_err)
#    plt.plot(time_range, timing_model_1_pointing, 'r--')
#    plt.xlabel('Bary diff (hours)', fontsize=14)
#    plt.ylabel('Max Phase Bins Number', fontsize=14)
#    plt.subplot(2,1,2)

    plt.plot(time_range, residuals_1_pointing, 'bo', markersize=markersize)
    plt.errorbar(time_range, residuals_1_pointing , yerr= phase_data_pointing_err, fmt=None, ecolor='b')
#    plt.plot(time_range, random_res_std, 'ks', markersize=markersize)
#    plt.errorbar(time_range, residuals_1_pointing , yerr= random_res_std, fmt=None, ecolor='r')
    plt.xlabel('MJD ', fontsize=10)
    plt.ylabel('Phase residuals', fontsize=10)
    plt.ticklabel_format(useOffset=False)

    plt.savefig(filename)

def time_pattern(this_file, bin_number, phase_amp_bin, NPHASEBIN = None):
    time_mjd = np.zeros(len(bin_number))
    dBATdra = np.zeros(len(bin_number))
    dBATddec = np.zeros(len(bin_number))
    for ii in range(len(bin_number)):
        time_mjd[ii] = (this_file['BARY_TIME'][bin_number[ii][0]] + this_file['BARY_TIME'][bin_number[ii][1]])/2.
        dBATdra[ii] = (this_file['dBATdra'][bin_number[ii][0]] + this_file['dBATdra'][bin_number[ii][1]])/2.
        dBATddec[ii] = (this_file['dBATddec'][bin_number[ii][0]] + this_file['dBATddec'][bin_number[ii][1]])/2.

    '''As pointed data to be the template, need to add the difference between two templates. Note: The phase_data of 1hr has been rescaled to 200, rather than 800 in origin.'''

    if (NPHASEBIN == None) or (NPHASEBIN == NPHASEBIN_wz):
        phase_data = phase_amp_bin[:,0] 
        phase_data_err = phase_amp_bin[:,phase_amp_bin.shape[1]/2]
    elif NPHASEBIN == NPHASEBIN_1hr:
        phase_data = phase_amp_bin[:,0] 
        phase_data_err = phase_amp_bin[:,phase_amp_bin.shape[1]/2]
    else:
        print 'NPHASEBIN error'

    return time_mjd, dBATdra, dBATddec, phase_data, phase_data_err

def fitting(pars_init_1, time_mjd, dBATdra, dBATddec, phase_data, phase_data_err, random_res):

    len_time_mjd_wz = len(pars.bin_number_wz) 
    fit_range = (untrans_time(pars.fit_time_start), untrans_time(pars.fit_time_end))
    sl = np.logical_and(time_mjd > fit_range[0], time_mjd < fit_range[1])

    '''parameters for new RA correction'''
    fit_pars_1, pcov_1, infodict_1, errmsg_1, success_1 = leastsq(residuals_1, pars_init_1,
            args=(time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl]), xtol = 1e-6, ftol=1e-6, full_output=1)

    print "Fit parameters: ", fit_pars_1
    print "sucess?:", success_1
    print "Chi-squared: ", np.sum(residuals_1(fit_pars_1, time_mjd[sl], dBATdra[sl], dBATddec[sl],
                         phase_data[sl], phase_data_err[sl], NPHASEBIN=None)**2), "DOF: ", len(phase_data[sl])-len(pars_init_1)

    if (len(phase_data) > len(pars_init_1)) and pcov_1 is not None:
        s_sq_1 = (residuals_1(fit_pars_1, time_mjd[sl], dBATdra[sl], dBATddec[sl],
            phase_data[sl], phase_data_err[sl], NPHASEBIN=None)**2).sum()/(len(phase_data[sl])-len(pars_init_1))
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

    print("\nFit paramters and parameter errors from lestsq method :")
    print("pfit = ", pfit_leastsq_1)
    print("perr = ", perr_leastsq_1)

    '''Std of random res'''
    random_res_std = np.zeros(len(random_res))
    for ii in xrange(len(random_res_std)):
        random_res_std[ii] = np.std(random_res[ii, :, 0])

    plt.figure()
    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl], random_res_std[sl], len_time_mjd_wz, 'new_phase_fit_J2139.png')
    plt.figure()
    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl], random_res_std[sl], len_time_mjd_wz, 'new_phase_fit_J2139_zoom.png',
           (untrans_time(-500), untrans_time(1500)))
    plt.figure()
    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl], random_res_std[sl], len_time_mjd_wz, 'new_phase_fit_J2139_zoom2.png',
           (untrans_time(33400), untrans_time(35314)))
    plt.figure()
    make_save_plot_pointing(fit_pars_1, 'model_1', 'res_1', time_mjd, dBATdra, dBATddec, phase_data, phase_data_err, random_res_std, len_time_mjd_wz, 'new_phase_fit_J2139_zoom3.png')
#            (untrans_time(35312), untrans_time(35314)))
    plt.figure()
    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd, dBATdra, dBATddec,  phase_data, phase_data_err, random_res_std, len_time_mjd_wz, 'new_phase_fit_J2139_all.png')
#    plt.figure()
#    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl], random_res_std[sl], len_time_mjd_wz, 'new_phase_fit_J2139_zoom_600.png',
#            (untrans_time(570), untrans_time(600)))


if __name__ == '__main__':
    main()
