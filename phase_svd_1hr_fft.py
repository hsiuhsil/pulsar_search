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
NPHASEBIN = 800
SCALE = 200/800.
T = 0.312470
TIME0 = 55707.   # MJD pivot

fit_pars =  [3.47609206e-07,  -3.37370520e+00,   4.74201487e+01,  -1.39689904e-05,  2.21425268e-05]

save_fit_pars = False

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
    fft_model = parameters[0] * np.exp(1.0j * 2 * np.pi * freq * ( n - parameters[1])) * profile_fft
    return fft_model

def fft_phase_curve_inverse(parameters, profile_fft):
    '''inverse phase for chaning 1.0j to -1 1.0j'''
    freq = np.fft.fftfreq(len(profile_fft))
    n= len(profile_fft)
    fft_model = parameters[0] * np.exp(-1.0j * 2 * np.pi * freq * ( n - parameters[1])) * profile_fft
    return fft_model


def residuals(parameters, model_fft, data_fft):
    '''Only use positive frequencies for residuals'''
    model = fft_phase_curve(parameters, model_fft)
    residuals_complex = (data_fft - model)[:len(model)/2]
    res_Re = residuals_complex.real
    res_Im = residuals_complex.imag
    res = np.concatenate((res_Re, res_Im))
    return res

def phase_fit(index, phase_matrix_origin, V, phase_model):

    pars_init = [np.amax(phase_matrix_origin[index])/np.amax(V[0]) , phase_model[index] % NPHASEBIN]
    pars_init2 = [0.29,phase_model[index]  + 0.8]
    pars_init3 = [0.29,phase_model[index]  + 0.1]
    model_fft = fft(V[0])
    data_fft = fft(phase_matrix_origin[index])

#    print('V0_max', np.amax(V[0]), 'index', np.argmax(V[0]))
#    print('data1_max', np.amax(phase_matrix_origin[1]), 'index', np.argmax(phase_matrix_origin[1]))
#    print('estimate amp', np.amax(phase_matrix_origin[1])/np.amax(V[0]))
#    print('estimate pars[1]', np.sqrt(sum(phase_matrix_origin[1]**2) / sum(V[0]**2)))
    fit_pars_phase, pcov, infodict, errmsg, success = leastsq(residuals, pars_init,
                   args=( model_fft, data_fft), full_output=1)

#    print "Fit parameters: ", fit_pars_phase
    print "sucess?:", success
    '''DOF -1, since we set fft_file[0] = 0.'''
    res = residuals(fit_pars_phase, model_fft, data_fft)
    print "Chi-squared: ", np.sum(np.abs(res)**2), "DOF: ", len(data_fft)-len(pars_init)-1

    if (len(data_fft) > len(pars_init)) and pcov is not None:
        s_sq = (residuals(fit_pars_phase, model_fft, data_fft)**2).sum()/(len(data_fft)-len(pars_init)-1)
        pcov = pcov * s_sq
    else:
        pcov = np.inf

    error = []
    for i in range(len(fit_pars_phase)):
        try:
          error.append(np.absolute(pcov[i][i])**0.5)
        except:
          error.append( 0.00 )
    pfit_leastsq = fit_pars_phase *SCALE
    perr_leastsq = np.array(error) *SCALE

    print("\nFit paramters and parameter errors from lestsq method :")
    print("pfit = ", pfit_leastsq)
    print("perr = ", perr_leastsq)

    '''save the fitting amp and bin as [amp, bin, amp_err, bin_err]'''
    npy_file = 'phase_amp_bin_57178_fft.npy'
    phase_amp_bin = np.concatenate((pfit_leastsq, perr_leastsq ))
    print 'phase_amp_bin',phase_amp_bin
    if save_fit_pars == True:
        if os.path.exists(npy_file):
            sequence = np.load(npy_file)
            np.save(npy_file, np.vstack((sequence, phase_amp_bin)))
        else:
            np.save(npy_file, phase_amp_bin)

    '''functions in Fourier space'''
    '''Real part'''
    model_fft_real = fft_phase_curve(fit_pars_phase, model_fft).real
    data_fft_real = data_fft.real
    init_fft_real = fft_phase_curve(pars_init, model_fft).real
    res_fft_real = np.concatenate((res[:(len(res)/2)], res[:(len(res)/2)][::-1]))
    '''Imag part'''
    model_fft_imag = fft_phase_curve(fit_pars_phase, model_fft).imag
    data_fft_imag = data_fft.imag
    init_fft_imag = fft_phase_curve(pars_init, model_fft).imag
    res_fft_imag = np.concatenate((res[(len(res)/2):], -res[(len(res)/2):][::-1]))

    '''functions in Real space (ifft)'''
    model_ifft = np.fft.ifft(fft_phase_curve(fit_pars_phase,  model_fft)).real
    data_ifft = np.fft.ifft(data_fft).real
    res_ifft = np.fft.ifft(fft_phase_curve(fit_pars_phase, model_fft) - data_fft).real
    init_ifft = np.fft.ifft(fft_phase_curve(pars_init, model_fft)).real

    freq_range = np.linspace(np.amin(np.fft.fftfreq(len(model_fft))), np.amax(np.fft.fftfreq(len(model_fft))), num = len(model_fft), endpoint=True)
    freq_min = np.amin(freq_range)
    freq_max = np.amax(freq_range)

    phase_range = np.arange(-400,400)

    plot_title = 'rescaled phase_bin = ' + str("%.3f" % phase_amp_bin[1]) + ' +/- ' + str("%.3f" % phase_amp_bin[3])
    plot_name = 'phase_57178_new_' + str(index) + '_'

    '''Plot for real part in the Fourier space'''
    plt.figure()
    plt.subplot(2,1,1)
    plt.title(plot_title)
    plt.plot(freq_range, np.roll(model_fft_real, -400),'r-')
    plt.plot(freq_range, np.roll(data_fft_real, -400),'b-')
    plt.plot(freq_range, np.roll(init_fft_real, -400),'k--')
    plt.xlabel('Frequency')
    plt.xlim((freq_min,freq_max))

    plt.subplot(2,1,2)
    plt.plot(freq_range, np.roll(res_fft_real, -400),'bo')
    plt.xlabel('Frequency')
    plt.ylabel('Residuals')
    plt.xlim((freq_min,freq_max))
    plt.savefig(plot_name + 'fft_real.png')

    '''Plot for imag part in the Fourier space'''
    plt.figure()
    plt.subplot(2,1,1)
    plt.title(plot_title)
    plt.plot(freq_range, np.roll(model_fft_imag, -400),'r-')
    plt.plot(freq_range, np.roll(data_fft_imag, -400),'b-')
    plt.plot(freq_range, np.roll(init_fft_imag, -400),'k--')
    plt.xlabel('Frequency')
    plt.xlim((freq_min,freq_max))

    plt.subplot(2,1,2)
    plt.plot(freq_range, np.roll(res_fft_imag, -400),'bo')
    plt.xlabel('Frequency')
    plt.ylabel('Residuals')
    plt.xlim((freq_min,freq_max))
    plt.savefig(plot_name + 'fft_imag.png')

    '''Plot for real part in real space'''
    plt.figure()
    plt.subplot(2,1,1)
    plt.title(plot_title)
    plt.plot(phase_range, np.roll(model_ifft, -400),'r-')
    plt.plot(phase_range, np.roll(data_ifft, -400),'b-')
    plt.plot(phase_range, np.roll(init_ifft, -400),'k--')
    plt.xlabel('Phase bin number')

    plt.subplot(2,1,2)
    plt.plot(phase_range, np.roll(res_ifft, -400),'bo')
    plt.xlabel('Phase bin number')
    plt.ylabel('Residuals')
    plt.savefig(plot_name + 'ifft.png')

def ploting(filename):

    this_file = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/J2139+00_57178h5', "r")

    '''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''
    bin_number = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/bin_number_2139_57178.txt')
 
    

    time_mjd = np.zeros(len(bin_number))
    dBATdra = np.zeros(len(bin_number))
    dBATddec = np.zeros(len(bin_number))
    for ii in range(len(bin_number)):
        time_mjd[ii] = (this_file['BARY_TIME'][bin_number[ii][0]] + this_file['BARY_TIME'][bin_number[ii][1]])/2.
        dBATdra[ii] = (this_file['dBATdra'][bin_number[ii][0]] + this_file['dBATdra'][bin_number[ii][1]])/2.
        dBATddec[ii] = (this_file['dBATddec'][bin_number[ii][0]] + this_file['dBATddec'][bin_number[ii][1]])/2.
    phase_data = bin_number[:,2]


    phase_matrix_origin = np.load(filename)
#    print 'phase_matrix_origin', phase_matrix_origin
    phase_model = timing_model_1(fit_pars, time_mjd, dBATdra, dBATddec) 
#    print 'phase_model', phase_model
    phase_matrix_new = np.zeros(phase_matrix_origin.shape)
    freq = np.fft.fftfreq(len(phase_matrix_origin[0]))
    n = len(fft(phase_matrix_origin[0]))
    for ii in xrange(len(phase_matrix_new)):
#        phase_matrix_new[ii] = np.roll(phase_matrix_origin[ii], -1 * phase_model[ii] )
#        profile_aligned[ii] = np.roll(profile[ii], -int(round(timing_model[ii])))
        phase_matrix_new[ii] = ifft(fft_phase_curve_inverse([1, phase_model[ii]], fft(phase_matrix_origin[ii])))
#        phase_matrix_new[ii] = ifft(fft(phase_matrix_origin[ii]) * np.exp(-1.0j * 2 * np.pi * freq * (n - phase_model[ii])))
    
    print 'finish phase_matrix_new' 
    U, s, V = np.linalg.svd(phase_matrix_new, full_matrices=True)

    if np.abs(np.amax(V[0])) < np.abs(np.amin(V[0])):
        V[0] = -V[0]

    print 'V[0][0:10]:', V[0][0:10]

    phase_fit(1, phase_matrix_origin, V, phase_model)

#    for ii in xrange(len(phase_model)):
#        print 'ii= '+str(ii)
#        phase_fit(ii, phase_matrix_origin, V, phase_model)

#    print 'phase_matrix_new.shape', phase_matrix_new.shape
#    print 's.shape', s.shape

    plt.figure()
    plt.plot(np.arange(64), s, 'ro-')
    plt.xlabel('phase bin number')
    plt.ylabel('s values')
    plt.savefig('phase_57178_new_s.png')

    n_step = -0.3
    plt.figure()
    plt.plot(np.arange(-400, 400), np.roll(V[0] + 0 *n_step, -400), 'r-',linewidth=1.0)
    plt.plot(np.arange(-400, 400), np.roll(V[1] + 1 *n_step, -400), 'b-',linewidth=1.0)
    plt.plot(np.arange(-400, 400), np.roll(V[2] + 2 *n_step, -400), 'g-',linewidth=1.0)
    plt.plot(np.arange(-400, 400), np.roll(V[3] + 3 *n_step, -400), 'k-',linewidth=1.0)
    plt.plot(np.arange(-400, 400), np.roll(V[4] + 4 *n_step, -400), 'y-',linewidth=1.0)
    plt.plot(np.arange(-400, 400), np.roll(V[5] + 5 *n_step, -400), color = '0.9',linewidth=1.0)
    plt.plot(np.arange(-400, 400), np.roll(V[6] + 6 *n_step, -400), color = '0.7',linewidth=1.0)
    plt.plot(np.arange(-400, 400), np.roll(V[7] + 7 *n_step, -400), color = '0.5',linewidth=1.0)
    plt.plot(np.arange(-400, 400), np.roll(V[8] + 8 *n_step, -400), color = '0.3',linewidth=1.0)
    plt.plot(np.arange(-400, 400), np.roll(V[9] + 9 *n_step, -400), color = '0.1',linewidth=1.0)
    plt.xlabel('phase bin number')
    plt.ylabel('V values')
    plt.savefig('phase_57178_new_V.png')


  
if __name__ == '__main__':
    main()

