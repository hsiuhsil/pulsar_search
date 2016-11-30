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
    args = sys.argv[1:]
    for filename in args:
        try:
            plot_bary_diff(filename)
        except None:
            print IOError


RA = 324.94166666666666  # deg
DEC = 0.6 # deg
AU = 149597870700.0      # m
C = 299792458.0    # m/s
NPHASEBIN = 100
T = 0.312470

TIME0 = 55707.   # MJD pivot


def transform_time(time_mjd):
    return (time_mjd - TIME0) * 24

def untrans_time(delta_time_hrs):
    return delta_time_hrs / 24. + TIME0


def timing_model_0(parameters, time_mjd, dBATdra, dBATddec):

    equinox_date = ['2010-03-20T17:32:00','2011-03-20T23:21:00','2012-03-20T05:14:00','2013-03-20T11:02:00','2014-03-20T16:57:00','2015-03-20T22:45:00','2016-03-20T04:30:00','2017-03-20T10:28:00']
    t = Time(equinox_date, format='isot', scale='utc')
    equinox_mjd = t.mjd
    theta = np.zeros(len(time_mjd))
    for ii in range(len(theta)):
        theta[ii] = (time_mjd[ii] - equinox_mjd[np.argmin(np.absolute(time_mjd[ii] - equinox_mjd))]) /  365.259636*2*np.pi

    time = transform_time(time_mjd)

#    out0 = -1.4e-6 * time**2
#    out0 += parameters[0] * time
#    out0 += parameters[1]
#    out0 += (-1*AU*parameters[2]/C*np.sin(RA*np.pi/180 + theta)) * (NPHASEBIN / T)
#    out0 = out0 % NPHASEBIN

    out1 = parameters[0] * time**2
    out1 += parameters[1] * time
    out1 += parameters[2]
    d_model_dBAT = 2 * parameters[0] * time + parameters[1]
    out1 +=  d_model_dBAT * (dBATdra * parameters[3] + dBATddec * parameters[4])
    out1 = out1 % NPHASEBIN

    return out1


def residuals_0(parameters, time_mjd, dBATdra, dBATddec, phase_data):
    model_0 = timing_model_0(parameters, time_mjd, dBATdra, dBATddec)

    res_0 = phase_data - model_0
    res_0 = (res_0 + NPHASEBIN / 2.) % NPHASEBIN - NPHASEBIN / 2.


    return res_0


def make_save_plot(parameters, time, dBATdra, dBATddec, data, filename, time_range=None):
    
    if time_range is None:
        time_range = (time[0], time[-1])

    num_points = 2500
    model_time = np.linspace(time_range[0], time_range[1], num_points)
    model_0 = timing_model_0(parameters, model_time, dBATdra, dBATddec)
    res_0 = residuals_0(parameters, time, dBATdra, dBATddec, data)

    plt.subplot(2,1,1)
    plt.plot(transform_time(time), data, 'bo')
    plt.plot(transform_time(model_time), model_0, 'r--')
    plt.xlabel('Bary diff (hours)', fontsize=14)
    plt.ylabel('Max Phase Bins Number', fontsize=14)
    plt.xlim(transform_time(time_range[0]), transform_time(time_range[1]))

    plt.subplot(2,1,2)
    plt.plot(transform_time(time), res_0, 'bo')
    plt.xlim(transform_time(time_range[0]), transform_time(time_range[1]))
    if len(time_range) > 2:
        plt.ylim(time_range[2], time_range[3])
    plt.xlabel('Bary diff (hours)', fontsize=14)
    plt.ylabel('Phase bin residuals', fontsize=14)

    plt.savefig(filename)


def plot_bary_diff(filename):
    this_file = h5py.File(filename, "r")
 
    '''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
    bin_number = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_RA_+21:39:46_DEC_+00:36:00/bin_number_2139.txt')
    

    time_mjd = np.zeros(len(bin_number))
    dBATdra = np.zeros(len(bin_number))
    dBATddec = np.zeros(len(bin_number))
    for ii in range(len(bin_number)):
        time_mjd[ii] = (this_file['BARY_TIME'][bin_number[ii][0]] + this_file['BARY_TIME'][bin_number[ii][1]])/2.
        dBATdra[ii] = (this_file['dBATdra'][bin_number[ii][0]] + this_file['dBATdra'][bin_number[ii][1]])/2.
        dBATddec[ii] = (this_file['dBATddec'][bin_number[ii][0]] + this_file['dBATddec'][bin_number[ii][1]])/2.
    phase_data = bin_number[:,2]


    # What data to fit.
    fit_range = (untrans_time(-1200), untrans_time(34000))
    #fit_range = (untrans_time(-1200), untrans_time(1800))
    sl = np.logical_and(time_mjd > fit_range[0], time_mjd < fit_range[1])

 
    #pars_init = (1e-06,  20.,   10.,   1e-04)
    #pars_init = [  8.35502582e-03,  -2.68046236e+01,  -3.85750551e+04,  -4.83510862e-01]
    #pars_init =  [  1.59790741e-03,   1.29654764e+01,  -6.58168024e+03, -8.29025454e-02]
    #pars_init = [  1.00413585e-03,   1.58649897e+01,  -4.24531396e+03,  -5.36469976e-02]
    #pars_init = [  3.41365396e-04,   1.93681348e+01,  -1.42306322e+03,  -1.83065288e-02]
    #pars_init = [  9.39658241e-06,   2.12734451e+01,   1.08429609e+02,   8.70694138e-04]
    #pars_init = [  1.33857708e-05,   2.12558751e+01,   9.41132616e+01,   6.91371967e-04]
#    pars_init = [  2.10776559e+01,  -5.13073206e+01,  -1.12983346e-03]
#    pars_init = [  2.11259893e+01,  -1.29592052e+01,  -6.54054280e-04]
#    pars_init = [  2.11573398e+01,   1.19145923e+01,  -3.45448859e-04]
#    pars_init = [  2.11886901e+01,   3.67880992e+01,  -3.68471222e-05]
    pars_init = [ -1.4e-06,  2.12200405e+01,   6.16618627e+01,  1e-04, 1e-04]
    pars_init = [ -1.4e-06,  2.12513909e+01,   8.65355252e+01,  -1e-04, -1e-04]
    #pars_init = [  2.12388508e+01,   7.65860533e+01,   4.56920031e-04]
    fit_pars, pcov, infodict, errmsg, success = leastsq(residuals_0, pars_init,
            args=(time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl]), full_output=1)
    print "Fit parameters: ", fit_pars
    print "sucess?:", success
    print "Chi-squared: ", np.sum(residuals_0(fit_pars, time_mjd[sl], dBATdra[sl], dBATddec[sl],
                         phase_data[sl])**2), "DOF: ", len(phase_data[sl])-len(pars_init)


    if (len(phase_data) > len(pars_init)) and pcov is not None:
        s_sq = (residuals_0(fit_pars, time_mjd[sl], dBATdra[sl], dBATddec[sl],
            phase_data[sl])**2).sum()/(len(phase_data[sl])-len(pars_init))
        pcov = pcov * s_sq
    else:
        pcov = np.inf

    error = [] 
    for i in range(len(fit_pars)):
        try:
          error.append(np.absolute(pcov[i][i])**0.5)
        except:
          error.append( 0.00 )
    pfit_leastsq = fit_pars
    perr_leastsq = np.array(error) 

    make_save_plot(pars_init, time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl],
            'phase_fit_J2139_init.png')
    plt.figure()
    make_save_plot(fit_pars, time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], 'phase_fit_J2139.png')
    plt.figure()
    make_save_plot(fit_pars, time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl], 'phase_fit_J2139_zoom.png',
            (untrans_time(-500), untrans_time(1500)))
    plt.figure()
    make_save_plot(pars_init, time_mjd[sl], dBATdra[sl], dBATddec[sl], phase_data[sl],
            'phase_fit_J2139_zoom2.png',
            (untrans_time(33500), untrans_time(33600)))
    plt.figure()
    make_save_plot(fit_pars, time_mjd, dBATdra, dBATddec,  phase_data, 'phase_fit_J2139_all.png')



if __name__ == '__main__':
    main()
