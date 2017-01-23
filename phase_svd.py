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

def fft(file):
    profile_fft = np.fft.fft(file)
    profile_fft[0] = 0
    return profile_fft

def ifft(file):
    profile_ifft = np.fft.ifft(file)
    return profile_ifft

def fft_phase_curve(parameters, profile_fft):
    freq = np.fft.fftfreq(len(profile_fft))
    n= len(profile_fft)
    fft_model = parameters[0] * np.exp(1.0j * 2 * np.pi * freq * (parameters[1] / n)) * profile_fft
    return fft_model

def residuals(parameters, model_fft, data_fft):
    model = fft_phase_curve(parameters, model_fft)
    res_Re = (data_fft - model).real
    res_Im = (data_fft - model).imag
    res = np.concatenate((res_Re, res_Im))
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

    if np.abs(np.amax(V[0])) < np.abs(np.amin(V[0])):
        V[0] = -V[0]


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

#    pars_init = [0.29, phase_model[1]]
    pars_init = [0.29, 3.0e+4]
    pars_init2 = [0.29, 3.1e+4]
    pars_init3 = [0.29, 3.2e+4]
    model_fft = fft(V[0])
    data_fft = fft(phase_matrix_origin[1])   

    print('V0_max', np.amax(V[0]), 'index', np.argmax(V[0]))
    print('data1_max', np.amax(phase_matrix_origin[1]), 'index', np.argmax(phase_matrix_origin[1]))
    print('estimate amp', np.amax(phase_matrix_origin[1])/np.amax(V[0]))
    print('estimate pars[1]', np.sqrt(sum(phase_matrix_origin[1]**2) / sum(V[0]**2)))
    fit_pars_phase, pcov, infodict, errmsg, success = leastsq(residuals, pars_init,
                   args=(model_fft, data_fft), full_output=1)

    print "Fit parameters: ", fit_pars_phase
    print "sucess?:", success
    '''DOF -1, since we set fft_file[0] = 0.'''
    res = residuals(fit_pars_phase, model_fft, data_fft)
    print "Chi-squared: ", np.sum(np.abs(res)**2), "DOF: ", len(data_fft)-len(pars_init)-1


    '''functions in Fourier space'''
    '''Real part'''
    model_fft_real = fft_phase_curve(fit_pars_phase, model_fft).real
    data_fft_real = data_fft.real
    init_fft_real = fft_phase_curve(pars_init, model_fft).real
    '''Imag part'''
    model_fft_imag = fft_phase_curve(fit_pars_phase, model_fft).imag
    data_fft_imag = data_fft.imag
    init_fft_imag = fft_phase_curve(pars_init, model_fft).imag

    '''functions in Real space (ifft)'''
    model_ifft = np.fft.ifft(fft_phase_curve(fit_pars_phase, model_fft))
    data_ifft = np.fft.ifft(data_fft)
    init_ifft = np.fft.ifft(fft_phase_curve(pars_init, model_fft))
    init_ifft2 = np.fft.ifft(fft_phase_curve(pars_init2, model_fft))    
    init_ifft3 = np.fft.ifft(fft_phase_curve(pars_init3, model_fft))

    freq_range = np.linspace(np.amin(np.fft.fftfreq(len(model_fft))), np.amax(np.fft.fftfreq(len(model_fft))), num = len(model_fft), endpoint=True)
    freq_min = np.amin(freq_range)
    freq_max = np.amax(freq_range)

    phase_range = np.arange(-100,100)

    plt.figure
    plt.plot(freq_range, np.roll(model_fft_real, -100),'r-')
    plt.plot(freq_range, np.roll(data_fft_real, -100),'b-')
    plt.plot(freq_range, np.roll(init_fft_real, -100),'k--')
    plt.xlabel('Frequency')
    plt.xlim((freq_min,freq_max))
    plt.savefig('phase_fft_real.png')

    plt.figure()
    plt.plot(freq_range, np.roll(model_fft_imag, -100),'r-')
    plt.plot(freq_range, np.roll(data_fft_imag, -100),'b-')
    plt.plot(freq_range, np.roll(init_fft_imag, -100),'k--')
    plt.xlabel('Frequency')
    plt.xlim((freq_min,freq_max))
    plt.savefig('phase_fft_imag.png')

    plt.figure()
    plt.plot(phase_range, np.roll(model_ifft, -100),'r-')
    plt.plot(phase_range, np.roll(data_ifft, -100),'b-')
    plt.plot(phase_range, np.roll(init_ifft, -100),'k--')
#    plt.plot(phase_range, np.roll(init_ifft2, -100),'y--')
#    plt.plot(phase_range, np.roll(init_ifft3, -100),'g-.')
    plt.xlabel('Phase bin number')
    plt.savefig('phase_ifft.png')


   
if __name__ == '__main__':
    main()

