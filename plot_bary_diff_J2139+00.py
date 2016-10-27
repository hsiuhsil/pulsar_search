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

def plot_bary_diff(filename):
    this_file = h5py.File(filename, "r")
 
#    index = [[4,8], [15, 19],[28, 32],[44, 49],[67, 83],[98,113],[114,117],[118,128],[146,152]]
#    max_phase = [69, 69, 63+100, 65+100, 90+200, 97+200, 91+300, 95+300, 5+400]

    '''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''

#    '''wigglez data'''
#    bin_number = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/bin_number_2139_wz.txt')
    '''wa data with delta_ra'''
    bin_number = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_RA_324.86337029308453/bin_number_2139_wz_dra.txt')

    '''one hour pointing data'''
#    bin_number = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/bin_number_2139_57178_20.txt')

    bary_diff = np.zeros(len(bin_number))
    mjd_ave = np.zeros(len(bin_number))

    for ii in range(len(bin_number)):
        bary_diff[ii] = ((this_file['BARY_TIME'][bin_number[ii][0]]+this_file['BARY_TIME'][bin_number[ii][1]])/2. -this_file['BARY_TIME'][0])*24
        mjd_ave[ii] = ((this_file['TOPO_TIME'][bin_number[ii][0]]+this_file['TOPO_TIME'][bin_number[ii][1]])/2.)

    '''Try to fit a parabolic curve'''

    def qua_func(x, a, b, c):
        return a*x**2 + b*x + c

    n_phase_bin = 100
    period = 0.31246381331597484
    sl = np.logical_and(bary_diff > 0, bary_diff < 1100)
    data_i = bin_number[sl,2]
    time_i = bary_diff[sl]
    mjd_i = mjd_ave[sl]
    t0 = time_i[0]
    time_i += 0


    '''Find quadratic parameters by leastsq'''
#    funcQuad=lambda tpl,time_i,data_i : (((data_i - ( tpl[0]*time_i**2 + tpl[1]*time_i + tpl[2]) + n_phase_bin/2) % n_phase_bin - n_phase_bin/2))
#    func=funcQuad
#    ErrorFunc=lambda tpl,time_i,data_i : func(tpl,time_i,data_i)
##   tplInitial=(0.031, 4.95, 15)
#    tplInitial = (1.89674411e-05, 2.11809594e+01, -1.85784448e+03)
#    tplFinal,success=leastsq(funcQuad,tplInitial[:],args=(time_i,data_i))
#    print "quadratic fit: " ,tplFinal
#    print "sucess?:", success
#    print np.sum(funcQuad(tplFinal, time_i, data_i)**2)

    '''Fit for delta_RA'''
    equinox_date = ['2010-03-20T17:32:00','2011-03-20T23:21:00','2012-03-20T05:14:00','2013-03-20T11:02:00','2014-03-20T16:57:00','2015-03-20T22:45:00','2016-03-20T04:30:00','2017-03-20T10:28:00']
    t = Time(equinox_date, format='isot', scale='utc')
    equinox_mjd = t.mjd
    theta_i = np.zeros(len(mjd_i))
    for ii in range(len(theta_i)):
        theta_i[ii] = (mjd_i[ii] - equinox_mjd[np.argmin(np.absolute(mjd_i[ii] - equinox_mjd))]) /  365.259636*2*np.pi
        '''RA, theta_i and delta_RA(tpl[3]) are in degree, AU in m, c in m/s'''
#    RA = 324.92817 '''(without dealt_ra)'''    
    RA = 324.86337029308453
    AU = 149597870700.0 
    c = 299792458.0
#    funcQuad=lambda tpl,time_i,data_i, theta_i : (((data_i - ( (tpl[0]*time_i**2 + tpl[1]*time_i + tpl[2])*(n_phase_bin/period/3600) + (1*AU*tpl[3]/c*np.sin(RA*np.pi/180 - theta_i))*(period/n_phase_bin)) + n_phase_bin/2) % n_phase_bin - n_phase_bin/2))
    funcQuad=lambda tpl,time_i,data_i, theta_i : (((data_i - ( (tpl[0]*time_i**2 + tpl[1]*time_i + tpl[2]) + (-1*AU*tpl[3]/c*np.sin(RA*np.pi/180 + theta_i))*(n_phase_bin/period)) + n_phase_bin/2) % n_phase_bin - n_phase_bin/2))
    func=funcQuad
    rrorFunc=lambda tpl,time_i,data_i, theta_i : func(tpl,time_i,data_i, theta_i)
    tplInitial = (6.47369211e-04,   1.69470826e+01,  -3.25790318e+03,  -3.85431541e-02)
    tplFinal,success=leastsq(funcQuad,tplInitial[:],args=(time_i,data_i, theta_i))
    print "quadratic fit: " ,tplFinal
    print "sucess?:", success
    print np.sum(funcQuad(tplFinal, time_i, data_i, theta_i)**2)

    num_points = 50
    theta_fit = np.zeros(num_points)
    mjd_range = np.linspace(np.amin(mjd_i), np.amax(mjd_i), num_points)
    for ii in range(len(theta_fit)):
        theta_fit[ii] = (mjd_range[ii] - equinox_mjd[np.argmin(np.absolute(mjd_range[ii] - equinox_mjd))]) /  365.259636*2*np.pi

    x_axes = np.linspace(np.amin(time_i), np.amax(time_i), num_points)
    y = (tplFinal[0]*x_axes**2 + tplFinal[1]*x_axes + tplFinal[2] + 1*AU*tplFinal[3]/c*np.sin(RA*np.pi/180 - theta_fit)) % n_phase_bin
    
    plt.subplot(2,1,1)
    plt.plot(time_i, data_i, 'bo')
    plt.plot(x_axes, y, 'r--')
    plt.xlabel('Bary diff (hours)', fontsize=14)
    plt.ylabel('Max Phase Bins Number', fontsize=14)

    plt.subplot(2,1,2)
    plt.plot(time_i, funcQuad(tplFinal, time_i, data_i, theta_i), 'bo')
    plt.xlabel('Bary diff (hours)', fontsize=14)
    plt.ylabel('Phase bin residuals', fontsize=14)

#    plt.show()
    plt.savefig('phase_fit_test.png')



if __name__ == '__main__':
    main()

