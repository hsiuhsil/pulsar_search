import sys
import os
import os.path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
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
#    out1 = pars.fit_pars[0] * time**2
#    out1 += pars.fit_pars[1] * time
#    out1 += pars.fit_pars[2]

    if (NPHASEBIN == None) or (NPHASEBIN == NPHASEBIN_wz):
        NPHASEBIN = NPHASEBIN_wz
    elif NPHASEBIN == NPHASEBIN_1hr:
        NPHASEBIN = NPHASEBIN_wz
    
    out1 += (NPHASEBIN / T) * (-1) *(dBATdra * 86400 * 180 / np.pi * parameters[3] + dBATddec * 86400 * 180 / np.pi * parameters[4])
#    out1 += (NPHASEBIN / T) * (-1) *(dBATdra * 86400 * 180 / np.pi * parameters[0] + dBATddec * 86400 * 180 / np.pi * parameters[1])

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
    yerr = data_err / NPHASEBIN
    '''Plot in time unit of ms, rather than phase'''
    res *= (T*1000)
    yerr *= (T*1000)

    zeros_line = np.zeros(len(res))

    markersize = 2.0
    fontsize = 16

    plt.plot(time, res, 'bo', markersize=markersize)
    plt.plot(time, zeros_line,'r--')
    plt.errorbar(time, res , yerr= yerr, fmt='none', ecolor='b')
    plt.xlim(time_range[0], time_range[1])
    if len(time_range) > 2:
        plt.ylim(time_range[2], time_range[3])
    plt.xlabel('MJD ', fontsize=fontsize)
    plt.ylabel('Timing residuals (ms)', fontsize=fontsize)
#    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.ticklabel_format(useOffset=False)

    plt.savefig(filename, bbox_inches='tight')

def make_save_plot_pointing(parameters, model, res, time, dBATdra, dBATddec, data, data_err, random_res_std, wz_range, filename, time_range=None):

    time_range = np.linspace(time[wz_range], time[-1], 40)
    timing_model_1_pointing = timing_model_1(parameters, time, dBATdra, dBATddec)[231:]

    phase_data_pointing =  data[231:]
    phase_data_pointing_err =  data_err[231:]

    residuals_1_pointing = residuals_1(parameters, time, dBATdra, dBATddec,  data, data_err)[231:]
    random_res_std = random_res_std[231:]

    phase_data_old = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_amp_bin_57178_40epochs_6modes_90hars.npy')
    phase_data_pointing_old = phase_data_old[:,0]
    phase_data_pointing_err_old = phase_data_old[:,phase_data_old.shape[1]/2]
    data_old = np.concatenate((data[0:231], phase_data_pointing_old))
    data_err_old = np.concatenate((data_err[0:231], phase_data_pointing_err_old))
    residuals_old_pointing = residuals_1(parameters, time, dBATdra, dBATddec, data_old, data_err_old)[231:]


    '''plot res, rather than res/errors'''
    residuals_1_pointing *= phase_data_pointing_err
    residuals_old_pointing *= phase_data_pointing_err_old
    '''plot in phase, rather than phase bin'''
    residuals_1_pointing /= NPHASEBIN
    phase_data_pointing_err = data_err[231:] / NPHASEBIN

    residuals_old_pointing /= NPHASEBIN
    phase_data_pointing_err_old /= NPHASEBIN

    '''Plot in time unit of ms, rather than phase'''
    residuals_1_pointing *= (T*1000)
    phase_data_pointing_err *= (T*1000)

    residuals_old_pointing *= (T*1000)
    phase_data_pointing_err_old *= (T*1000)


    zeros_line = np.zeros(len(residuals_1_pointing))
    markersize = 2.0
    fontsize = 16

    plt.close('all')
    plt.plot(time_range, zeros_line,'r--')
    plt.plot(time_range, residuals_1_pointing, 'bo', markersize=markersize)
    plt.errorbar(time_range, residuals_1_pointing , yerr= phase_data_pointing_err, fmt='none', ecolor='b')
    plt.plot(time_range, residuals_old_pointing, 'gs', markersize=markersize)
    plt.errorbar(time_range, residuals_old_pointing , yerr= phase_data_pointing_err_old, fmt='none', ecolor='g')
    plt.xlabel('MJD ', fontsize=fontsize)
    plt.ylabel('Timing residuals (ms)', fontsize=fontsize)
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.ticklabel_format(useOffset=False)
   
    plt.savefig(filename, bbox_inches='tight')

    plt.close('all')
    plt.hist(residuals_1_pointing, bins='auto')
    plt.xlabel('Timing residuals (ms)', fontsize=fontsize)
    plt.ylabel('Number', fontsize=fontsize)
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.savefig('hist_'+filename, bbox_inches='tight')

def make_save_plot_panel(parameters, model, res, time, dBATdra, dBATddec, data, data_err, random_res_std, wz_range, filename, time_range=None):

    if time_range is None:
        time_range = (time[0], time[-1])

    num_points = len(time)
    model_time = np.linspace(time_range[0], time_range[1], num_points)

    model = timing_model_1(parameters, model_time, dBATdra, dBATddec)
    res = residuals_1(parameters, time, dBATdra, dBATddec, data, data_err)

    '''plot res, rather than res/errors'''
    res *= data_err
    yerr = data_err
    '''plot in phase, rather than phase bin'''
    res /= NPHASEBIN
    yerr = data_err / NPHASEBIN
    '''Plot in time unit of ms, rather than phase'''
    res *= (T*1000)
    yerr *= (T*1000)

    print 'res and err in ms'
    print 'res/yerr >3.0 index:', np.where( res/yerr >3)
    print (res/yerr)[np.where( res/yerr >3)]
    print 'res[2]/yerr[2]',res[2]/yerr[2]

    zeros_line = np.zeros(len(res))

    # outlier cases parameters from full timing solution.
    time_2 = 55694.6342468 
    res_2 = 8.36257644547 
    yerr_2 = 0.890570241555
    time_46 = 55704.5156269 
    res_46 = 1.18155430711 
    yerr_46 = 0.3088385541
    time_197 = 57103.6670504 
    res_197 = 3.13366921351 
    yerr_197 = 0.942867422513

#    print 'wz_2: ', time[2], res[2], yerr[2]
#    print 'wz_46: ', time[46], res[46], yerr[46]
#    print 'wz_197: ', time[197], res[197], yerr[197]

    markersize = 2.0
    fontsize = 16

    time_11wz = time[0:185]
    res_11wz = res[0:185]
    yerr_11wz = yerr[0:185]
    zeros_line_11wz = zeros_line[0:185]

    time_15wz = time[185:231]
    res_15wz = res[185:231]
    yerr_15wz = yerr[185:231]
    zeros_line_15wz = zeros_line[185:231]

    time_1hr = time[231:]
    res_1hr = res[231:]
    yerr_1hr = yerr[231:]
    zeros_line_1hr = zeros_line[231:]

    plt.figure(0, figsize=(16,9))
    ax1 = plt.subplot2grid((2,3), (0,0), colspan=3)
    ax2 = plt.subplot2grid((2,3), (1, 0))
    ax3 = plt.subplot2grid((2,3), (1, 1), sharey=ax2)
    ax4 = plt.subplot2grid((2,3), (1, 2))

    ax1.plot(time, res, 'bo', markersize=markersize)
    ax1.plot(time, zeros_line,'r--')
    ax1.errorbar(time, res , yerr= yerr, fmt='none', ecolor='b')
    ax1.plot(time_2, res_2, 'gs', markersize=markersize*2)
    ax1.errorbar(time_2, res_2, yerr=yerr_2, fmt='none', ecolor='g', elinewidth=2.0)
#   ax1.plot(time_46, res_46, 'gs', markersize=markersize*2)
#   ax1.errorbar(time_46, res_46, yerr=yerr_46, fmt='none', ecolor='g')
#   ax1.plot(time_197, res_197, 'gs', markersize=markersize*2)
#   ax1.errorbar(time_197, res_197, yerr=yerr_197, fmt='none', ecolor='g')
    ax1.set_xlim([np.amin(time),np.amax(time)])
    ax1.tick_params(axis='both', which='major', labelsize=fontsize)
    ax1.set_xlabel('MJD', fontsize=fontsize)
    ax1.set_ylabel('Timing residuals (ms)', fontsize=fontsize)

    ax2.plot(time_11wz, res_11wz, 'bo', markersize=markersize)
    ax2.plot(time_11wz, zeros_line_11wz,'r--')
    ax2.errorbar(time_11wz, res_11wz , yerr= yerr_11wz, fmt='none', ecolor='b')
    ax2.plot(time_2, res_2, 'gs', markersize=markersize*2)
    ax2.errorbar(time_2, res_2, yerr=yerr_2, fmt='none', ecolor='g', elinewidth=2.0)
#   ax2.plot(time_46, res_46, 'gs', markersize=markersize*2)
#   ax2.errorbar(time_46, res_46, yerr=yerr_46, fmt='none', ecolor='g', elinewidth=2.0)
#    ax2.xaxis.set_ticks([np.percentile(time_11wz, 25), np.percentile(time_11wz, 75), np.percentile(time_11wz, 90)])
    ax2.xaxis.set_ticks([55700, 55730, 55760])
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax2.yaxis.set_ticks([-8, -4, 0, 4, 8, 10])
    ax2.set_ylim([-8,10])
#    ax2.get_xaxis().get_major_formatter().set_useOffset(False)
    ax2.tick_params(axis='both', which='major', labelsize=fontsize)
    ax2.set_xlabel('MJD', fontsize=fontsize)
    ax2.set_ylabel('Timing residuals (ms)', fontsize=fontsize)

    ax3.plot(time_15wz, res_15wz, 'bo', markersize=markersize)
    ax3.plot(time_15wz, zeros_line_15wz,'r--')
    ax3.errorbar(time_15wz, res_15wz , yerr= yerr_15wz, fmt='none', ecolor='b')
#   ax3.plot(time_197, res_197, 'gs', markersize=markersize*2)
#   ax3.errorbar(time_197, res_197, yerr=yerr_197, fmt='none', ecolor='g')
    ax3.set_xlim([np.amin(time_15wz),np.amax(time_15wz)])
#    ax3.xaxis.set_ticks([np.percentile(time_15wz, 25), np.percentile(time_15wz, 75), np.percentile(time_15wz, 90)])
    ax3.xaxis.set_ticks([57100, 57110, 57120, 57130])
    ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax3.yaxis.set_ticks([-8, -4, 0, 4, 8, 10])
    ax3.set_ylim([-8,10])
#    ax3.get_xaxis().get_major_formatter().set_useOffset(False)
    ax3.tick_params(axis='both', which='major', labelsize=fontsize)
    ax3.set_xlabel('MJD', fontsize=fontsize)

    ax4.plot(time_1hr, res_1hr, 'bo', markersize=markersize)
    ax4.plot(time_1hr, zeros_line_1hr,'r--')
    ax4.errorbar(time_1hr, res_1hr , yerr= yerr_1hr, fmt='none', ecolor='b')
#    ax4.set_xlim([np.amin(time_1hr),np.amax(time_1hr)])
    ax4.xaxis.set_ticks([np.percentile(time_1hr, 25), np.percentile(time_1hr, 75)])
    ax4.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#    ax4.yaxis.set_ticks([-1.2, 0, 0.6])
#    ax4.set_ylim([-1.2,0.6])
#    ax4.get_xaxis().get_major_formatter().set_useOffset(False)
    ax4.tick_params(axis='both', which='major', labelsize=fontsize)
    ax4.set_xlabel('MJD', fontsize=fontsize)


#    plt.ticklabel_format(useOffset=False)
    plt.savefig(filename, bbox_inches='tight')

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

    len_time_mjd_wz = 231
    fit_range = (untrans_time(pars.fit_time_start), untrans_time(pars.fit_time_end))
    sl = np.logical_and(time_mjd > fit_range[0], time_mjd < fit_range[1])

    '''parameters for new RA correction'''
    fit_pars_1, pcov_1, infodict_1, errmsg_1, success_1 = leastsq(residuals_1, pars_init_1,
            args=(time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl]), xtol = 1e-6, ftol=1e-6, full_output=1)

    print "Fit parameters: ", fit_pars_1
    print "sucess?:", success_1
    print "Chi-squared: ", np.sum(residuals_1(fit_pars_1, time_mjd[sl], dBATdra[sl], dBATddec[sl],
                         phase_data[sl], phase_data_err[sl], NPHASEBIN=None)**2), "DOF: ", len(phase_data[sl])-len(pars_init_1)

    res = residuals_1(fit_pars_1, time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl], NPHASEBIN=None) * phase_data_err[sl]
#    print res
    print 'max res: ', np.amax(np.abs(res))
    print 'max res index: ', np.argmax(np.abs(res))

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
#    print("pfit_0 = ", "%.12f" % pfit_leastsq_1[0])
#    print("pfit_1 = ", "%.12f" % pfit_leastsq_1[1])
#    print("perr_0 = ", "%.12f" % perr_leastsq_1[0])
#    print("perr_1 = ", "%.12f" % perr_leastsq_1[1])

    '''Std of random res'''
    random_res_std = np.zeros(len(random_res))
    for ii in xrange(len(random_res_std)):
        random_res_std[ii] = np.std(random_res[ii, :, 0])

    plt.figure()
    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl], random_res_std[sl], len_time_mjd_wz, 'new_phase_fit_J2139_zoom1.png',
           (untrans_time(-300), untrans_time(1500)))
    plt.figure()
    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], phase_data_err[sl], random_res_std[sl], len_time_mjd_wz, 'new_phase_fit_J2139_zoom2.png',
           (untrans_time(33400), untrans_time(34300)))
    plt.figure()
    plt.close('all')
    make_save_plot_pointing(fit_pars_1, 'model_1', 'res_1', time_mjd, dBATdra, dBATddec, phase_data, phase_data_err, random_res_std, len_time_mjd_wz, 'new_phase_fit_J2139_zoom3.png',        (untrans_time(35312.4), untrans_time(35313.4)))
    plt.close('all')
    plt.figure()
    make_save_plot(fit_pars_1, 'model_1', 'res_1', time_mjd, dBATdra, dBATddec,  phase_data, phase_data_err, random_res_std, len_time_mjd_wz, 'new_phase_fit_J2139_all.png')
    plt.close('all')
    plt.figure()
    make_save_plot_panel(fit_pars_1, 'model_1', 'res_1', time_mjd, dBATdra, dBATddec,  phase_data, phase_data_err, random_res_std, len_time_mjd_wz, 'new_phase_fit_J2139_all_panel.png')

def plot_wz_hist(bin_number_wz):
    # the histogram of the data
    t = (bin_number_wz[:,1] - bin_number_wz[:,0] +1)
    print 'total count', np.count_nonzero(t)
    t *= (2048*0.001024)
    plt.hist(t, bins=np.arange(0, 32, 2.097))
#    plt.xticks(np.arange(0,34, 2.097))
    plt.xlim((0,30))
    plt.xlabel('Target time (sec)')
    plt.ylabel('Counts')
#    plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
    plt.savefig('wz_hist.png', bbox_inches='tight')

if __name__ == '__main__':
    main()
