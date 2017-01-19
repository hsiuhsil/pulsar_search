import sys
import os
import os.path

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import leastsq

#from __future__ import division, print_function

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            ploting(filename)
        except (IOError, ValueError):
            print IOError

RA = 324.8428583333333  # deg
DEC = 0.6959230555555556 # deg
AU = 149597870700.0      # m
C = 299792458.0    # m/s
NPHASEBIN = 200
T = 0.312470

TIME0 = 55707.   # MJD pivot

fit_pars =  [3.47609206e-07,  -3.37370520e+00,   4.74201487e+01,  -1.39689904e-05,  2.21425268e-05]


def transform_time(time_mjd):
    return (time_mjd - TIME0) * 24

def timing_model_1(parameters, time_mjd, dBATdra, dBATddec):
    time = transform_time(time_mjd)

    out1 = parameters[0] * time**2
    out1 += parameters[1] * time
    out1 += parameters[2]
    out1 +=  (NPHASEBIN / T) * (dBATdra * 86400 * 180 / np.pi * parameters[3] + dBATddec * 86400 * 180 / np.pi * parameters[4])
    out1 = out1 % NPHASEBIN
    out1 = np.around(out1).astype(int)
    for ii in xrange(len(out1)):
        if out1[ii] == NPHASEBIN:
            out1[ii] == 0

    return out1

def phase_curve(parameters, V):
    # parameters = [amp, phase_bin, offset]
    phase_curve = parameters[0] * (np.roll(V[0], int(parameters[1]))) + parameters[2]
    return phase_curve

def residuals(parameters, V,  data):
    model = phase_curve(parameters, V)
    res = data - model
    return res

def ploting(filename):

    this_file = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_ANTF_delta_ra_dec_20170116/J2139+00_wzonlyh5', "r")

    '''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
    bin_number = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_ANTF_delta_ra_dec_20170116/bin_number_2139_delta_new2.txt')
 
    

    time_mjd = np.zeros(len(bin_number))
    dBATdra = np.zeros(len(bin_number))
    dBATddec = np.zeros(len(bin_number))
    for ii in range(len(bin_number)):
        time_mjd[ii] = (this_file['BARY_TIME'][bin_number[ii][0]] + this_file['BARY_TIME'][bin_number[ii][1]])/2.
        dBATdra[ii] = (this_file['dBATdra'][bin_number[ii][0]] + this_file['dBATdra'][bin_number[ii][1]])/2.
        dBATddec[ii] = (this_file['dBATddec'][bin_number[ii][0]] + this_file['dBATddec'][bin_number[ii][1]])/2.
    phase_data = bin_number[:,2]


    phase_matrix_origin = np.load(filename)
    phase_model = timing_model_1(fit_pars, time_mjd, dBATdra, dBATddec) 

    phase_matrix_new = np.zeros(phase_matrix_origin.shape)
    for ii in xrange(len(phase_matrix_new)):
        phase_matrix_new[ii] = np.roll(phase_matrix_origin[ii], -1 * phase_model[ii] )
    
    U, s, V = np.linalg.svd(phase_matrix_new, full_matrices=True)

#    plt.figure()
#    plt.plot(np.arange(200), s, 'ro-')
#    plt.xlabel('phase bin number')
#    plt.ylabel('s values')
#    plt.savefig('phase_s.png')

#    plt.figure()
#    plt.plot(np.arange(-100, 100), np.roll(V[0]     , -100), 'r-',linewidth=2.5)
#    plt.plot(np.arange(-100, 100), np.roll(V[1] -0.1, -100), 'b-',linewidth=2.5)
#    plt.plot(np.arange(-100, 100), np.roll(V[2] -0.2, -100), 'g-',linewidth=2.5)
#    plt.plot(np.arange(-100, 100), np.roll(V[3] -0.3, -100), 'k-',linewidth=2.5)
#    plt.plot(np.arange(-100, 100), np.roll(V[4] -0.4, -100), 'y-',linewidth=2.5)
#    plt.plot(np.arange(-100, 100), np.roll(V[5] -0.5, -100), color = '0.9',linewidth=1.5)
#    plt.plot(np.arange(-100, 100), np.roll(V[6] -0.6, -100), color = '0.7',linewidth=1.5)
#    plt.plot(np.arange(-100, 100), np.roll(V[7] -0.7, -100), color = '0.5',linewidth=1.5)
#    plt.plot(np.arange(-100, 100), np.roll(V[8] -0.8, -100), color = '0.3',linewidth=1.5)
#    plt.plot(np.arange(-100, 100), np.roll(V[9] -0.9, -100), color = '0.1',linewidth=1.5)
#    plt.xlabel('phase bin number')
#    plt.ylabel('V values')
#    plt.savefig('phase_V.png')

    for ii in xrange(0,3):
        pars_init = [ -0.1, 5, 0.2]

        fit_pars_phase, pcov, infodict, errmsg, success = leastsq(residuals, pars_init,
               args=(V, phase_matrix_origin[ii]), full_output=1)
#                args=(V, phase_matrix_origin[ii]), xtol = 1e-6, ftol=1e-6, full_output=1)
        print 'ii = '+str(ii)
        print "Fit parameters: ", fit_pars_phase
        print "sucess?:", success
        print "Chi-squared: ", np.sum(residuals(fit_pars_phase, V, phase_matrix_origin[ii])**2), "DOF: ", len(phase_matrix_origin[ii])-len(pars_init)


# points in time series
    n= len(V)
# final time (initial time is 0)
    tfin= 200

# *end of changeable parameters*

# stepsize
    dt= tfin/(n-1)
# sample count
    s= np.arange(n)
# signal; somewhat arbitrary
    y= phase_matrix_origin[1]
# DFT
    fy= np.fft.fft(y)
# frequency spectrum in rad/sample
    wps= np.linspace(0,2*np.pi,n+1)[:-1]

# basis for DFT
# see, e.g., http://en.wikipedia.org/wiki/Discrete_Fourier_transform#equation_Eq.2
# and section "Properties -> Orthogonality"; the columns of 'basis' are the u_k vectors
# described there
    basis= 1.0/n*np.exp(1.0j * V * wps * s[:,np.newaxis]) 

# reconstruct signal from DFT coeffs and basis
    recon_y= np.dot(basis,fy)

# expect yerr to be "small"
    yerr= np.max(np.abs(y-recon_y))
    print('yerr:',yerr)

# find coefficients by fitting to basis
    lin_fy= np.linalg.solve(basis,y)

# fyerr should also be "small"
    fyerr= np.max(np.abs(fy-lin_fy))
    print('fyerr',fyerr)


   
if __name__ == '__main__':
    main()

