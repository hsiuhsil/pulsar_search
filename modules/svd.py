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

import pars
import fit_timing


def main():
    try:
        ploting()
    except (IOError, ValueError):
        print IOError

aligned_1bin = False 
aligned_fft = True
save_fit_pars = True


RA = pars.RA
DEC = pars.DEC
AU = pars.AU
C = pars.C
NPHASEBIN_wz = pars.NPHASEBIN_wz
NPHASEBIN = pars.NPHASEBIN
NPHASEBIN_1hr = pars.NPHASEBIN_1hr
SCALE = pars.SCALE
T = pars.T
PHASE_DIFF_wz_1hr = pars.PHASE_DIFF_wz_1hr
PHASE_DIFF_wz_1hr_err = pars.PHASE_DIFF_wz_1hr_err

TIME0 = pars.TIME0

fit_pars =  pars.fit_pars

def fft(file):
    profile_fft = np.fft.fft(file)
    profile_fft[0] = 0
    return profile_fft

def ifft(file):
    profile_ifft = np.fft.ifft(file)
    return profile_ifft

def fft_phase_curve(parameters, profile, freq=None):
    if freq is None:
        freq = np.fft.fftfreq(len(profile))
    else:
        freq = freq

    profile_fft = 0.
    for ii in xrange(0,len(parameters)-1):
        if ii == 0:
            profile_fft += parameters[ii + 1] * fft(profile[ii,:])
        elif ii >= 1:
            profile_ii_fft = fft(profile[ii,:])
            for jj in xrange(len(profile[ii,:])):
                if (np.abs(np.fft.fftfreq(len(profile[ii,:]), 1./len(profile[ii,:]))) > 100)[jj] ==True:
                    profile_ii_fft[jj] = 0
            profile_fft += parameters[ii + 1] * profile_ii_fft

        
#    profile_fft = fft(model)
    n= len(freq)
    fft_model = np.exp(1.0j * 2 * np.pi * freq * ( n - parameters[0])) * parameters[1] * np.concatenate((profile_fft[:len(freq)/2], profile_fft[-len(freq)/2:]))

    return fft_model

def fft_phase_curve_inverse(parameters, profile_fft):
    '''inverse phase for chaning 1.0j to -1 1.0j'''
    freq = np.fft.fftfreq(len(profile_fft))
    n= len(profile_fft)
    fft_model = parameters[1] * np.exp(-1.0j * 2 * np.pi * freq * ( n - parameters[0])) * profile_fft
    return fft_model

def residuals(parameters, model, data_fft, freq=None):
    '''Only use positive frequencies for residuals'''
    if freq is None:
        freq = None
    else:
        freq = freq  
 
    model = fft_phase_curve(parameters, model, freq)
    residuals_complex = (data_fft - model)[:len(model)/2]
    res_Re = residuals_complex.real
    res_Im = residuals_complex.imag
    res = np.concatenate((res_Re, res_Im))
    return res

def phase_fit(index, phase_matrix_origin, V, plot_name, NPHASEBIN=None, RESCALE=None):

    if (NPHASEBIN == None) or (NPHASEBIN == NPHASEBIN_wz):
        NPHASEBIN = NPHASEBIN_wz
    elif NPHASEBIN == NPHASEBIN_1hr:
        NPHASEBIN = NPHASEBIN_1hr

    '''pars_init = [phase_bin, amp_V0, amp_V1, amp_V2]'''
    '''initial value of phase_bin'''
    time_mjd_wz, dBATdra_wz, dBATddec_wz, phase_data_wz, phase_data_err_wz = fit_timing.time_pattern(pars.this_file_wz, pars.bin_number_wz, pars.phase_amp_bin_wz, pars.NPHASEBIN_wz)

    model_phase_bin_wz = fit_timing.timing_model_1(pars.fit_pars, time_mjd_wz, dBATdra_wz, dBATddec_wz, pars.NPHASEBIN_wz, RESCALE=None)

    time_mjd_1hr, dBATdra_1hr, dBATddec_1hr, phase_data_1hr, phase_data_err_1hr = fit_timing.time_pattern(pars.this_file_1hr, pars.bin_number_1hr, pars.phase_amp_bin_1hr, pars.NPHASEBIN_1hr)

    model_phase_bin_1hr = fit_timing.timing_model_1(pars.fit_pars, time_mjd_1hr, dBATdra_1hr, dBATddec_1hr, pars.NPHASEBIN_1hr, RESCALE=None)

    model = V
    data_fft = fft(phase_matrix_origin[index])
    freq = np.fft.fftfreq(len(data_fft))    

    if len(model) > len(data_fft):
        amp_V0 = np.amax(phase_matrix_origin[index])/np.amax(V[0]) * 0.5
    else:
        amp_V0 = np.amax(phase_matrix_origin[index])/np.amax(V[0])
    pars_init = [ model_phase_bin_wz[index], amp_V0, amp_V0/10, amp_V0/5000, amp_V0/500000]

#    pars_init = [ np.argmax(phase_matrix_origin[index]) % NPHASEBIN, amp_V0, amp_V0/10, amp_V0/5000, amp_V0/500000]
    print 'pars_init', pars_init

#    model_fft = fft(V[0])

    fit_pars_phase, pcov, infodict, errmsg, success = leastsq(residuals, pars_init,
                   args=( model, data_fft, freq), full_output=1)

    print "Fit parameters: ", fit_pars_phase
    print "sucess?:", success
    '''DOF -1, since we set fft_file[0] = 0.'''
    res = residuals(fit_pars_phase, model, data_fft, freq)
    print "Chi-squared: ", np.sum(np.abs(res)**2), "DOF: ", len(data_fft)-len(pars_init)-1

    if (len(data_fft) > len(pars_init)) and pcov is not None:
        s_sq = (residuals(fit_pars_phase, model, data_fft, freq)**2).sum()/(len(data_fft)-len(pars_init)-1)
        pcov = pcov * s_sq
    else:
        pcov = np.inf

    error = []
    for i in range(len(fit_pars_phase)):
        try:
          error.append(np.absolute(pcov[i][i])**0.5)
        except:
          error.append( 0.00 )
    pfit_leastsq = fit_pars_phase 
    perr_leastsq = np.array(error)

    print("\nFit paramters and parameter errors from lestsq method :")
    print("pfit = ", pfit_leastsq)
    print("perr = ", perr_leastsq)

    '''save the fitting amp and bin as [amp, bin, amp_err, bin_err]'''
    npy_file = 'phase_amp_bin_wz_fft.npy'

    if (phase_matrix_origin.shape[1] == pars.phase_npy_1hr.shape[1]):
        phase_amp_bin = np.concatenate(([pfit_leastsq[0]*SCALE], pfit_leastsq[1:], [perr_leastsq[0]*SCALE], perr_leastsq[1:] ))
    else:
        phase_amp_bin = np.concatenate((pfit_leastsq, perr_leastsq ))

    if save_fit_pars == True:
        if os.path.exists(npy_file):
            sequence = np.load(npy_file)
            np.save(npy_file, np.vstack((sequence, phase_amp_bin)))
        else:
            np.save(npy_file, phase_amp_bin)
    '''functions in Fourier space'''
    model_fft = fft_phase_curve(fit_pars_phase, model, freq)  
    init_fft = fft_phase_curve(pars_init, model, freq)

    '''Real part'''
    model_fft_real = model_fft.real
    data_fft_real = data_fft.real
    init_fft_real = init_fft.real
    res_fft_real = np.concatenate((res[:(len(res)/2)], res[:(len(res)/2)][::-1]))

    '''Imag part'''
    model_fft_imag = model_fft.imag
    data_fft_imag = data_fft.imag
    init_fft_imag = init_fft.imag
    res_fft_imag = np.concatenate((res[(len(res)/2):], -res[(len(res)/2):][::-1]))

    '''functions in Real space (ifft)'''
    model_ifft = np.fft.ifft(model_fft).real
    data_ifft = np.fft.ifft(data_fft).real
    res_ifft = np.fft.ifft(model_fft - data_fft).real
    init_ifft = np.fft.ifft(init_fft).real

    freq_range = np.linspace(np.amin(np.fft.fftfreq(len(data_fft))), np.amax(np.fft.fftfreq(len(data_fft))), num = len(data_fft), endpoint=True)
    freq_min = np.amin(freq_range)
    freq_max = np.amax(freq_range)

    phase_range = np.arange(-int(NPHASEBIN/2), int(NPHASEBIN/2))

    plot_title = 'rescaled phase_bin = ' + str("%.3f" % phase_amp_bin[0]) + ' +/- ' + str("%.3f" % phase_amp_bin[len(phase_amp_bin)/2])
    plot_name += str(index) + '_'

    '''Plot for real part in the Fourier space'''
    plt.figure()
    plt.subplot(2,1,1)
    plt.title(plot_title)
    plt.plot(freq_range, np.roll(model_fft_real, -int(NPHASEBIN/2)),'r-')
    plt.plot(freq_range, np.roll(data_fft_real, -int(NPHASEBIN/2)),'b-')
    plt.plot(freq_range, np.roll(init_fft_real, -int(NPHASEBIN/2)),'k--')
    plt.xlabel('Frequency')
    plt.xlim((freq_min,freq_max))
    plt.subplot(2,1,2)
    plt.plot(freq_range, np.roll(res_fft_real, -int(NPHASEBIN/2)),'bo')
    plt.xlabel('Frequency')
    plt.ylabel('Residuals')
    plt.xlim((freq_min,freq_max))
    plt.savefig(plot_name + 'fft_real.png')
  
    '''Plot for imag part in the Fourier space'''
    plt.figure()
    plt.subplot(2,1,1)
    plt.title(plot_title)
    plt.plot(freq_range, np.roll(model_fft_imag, -int(NPHASEBIN/2)),'r-')
    plt.plot(freq_range, np.roll(data_fft_imag, -int(NPHASEBIN/2)),'b-')
    plt.plot(freq_range, np.roll(init_fft_imag, -int(NPHASEBIN/2)),'k--')
    plt.xlabel('Frequency')
    plt.xlim((freq_min,freq_max))

    plt.subplot(2,1,2)
    plt.plot(freq_range, np.roll(res_fft_imag, -int(NPHASEBIN/2)),'bo')
    plt.xlabel('Frequency')
    plt.ylabel('Residuals')
    plt.xlim((freq_min,freq_max))
    plt.savefig(plot_name + 'fft_imag.png')

    '''Plot for real part in real space'''
    plt.figure()
    plt.subplot(2,1,1)
    plt.title(plot_title)
    plt.plot(phase_range, np.roll(model_ifft, -int(NPHASEBIN/2)),'r-')
    plt.plot(phase_range, np.roll(data_ifft, -int(NPHASEBIN/2)),'b-')
    plt.plot(phase_range, np.roll(init_ifft, -int(NPHASEBIN/2)),'k--')
    plt.xlabel('Phase bin number')

    plt.subplot(2,1,2)
    plt.plot(phase_range, np.roll(res_ifft, -100),'bo')
    plt.xlabel('Phase bin number')
    plt.ylabel('Residuals')
    plt.savefig(plot_name + 'ifft.png')


def scale_matrix(old_matrix, SCALE):
    matrix = np.zeros((old_matrix.shape[0], int((old_matrix.shape[1])*SCALE)))
    for ii in xrange(len(matrix)):
        for jj in xrange(len(matrix[0])):
            matrix[ii, jj] = np.average(old_matrix[ii,int(jj/SCALE):int((jj+1)/SCALE)])
    return matrix


def svd(this_file, bin_number, phase_amp_bin, phase_npy, NPHASEBIN, RESCALE):

    if (NPHASEBIN == None) or (NPHASEBIN == NPHASEBIN_wz):
        NPHASEBIN = NPHASEBIN_wz
    elif NPHASEBIN == NPHASEBIN_1hr:
        NPHASEBIN = NPHASEBIN_1hr

    time_mjd, dBATdra, dBATddec, _ , _ = fit_timing.time_pattern(this_file, bin_number, phase_amp_bin, NPHASEBIN)
    
    phase_model = fit_timing.timing_model_1(fit_pars, time_mjd, dBATdra, dBATddec, NPHASEBIN, RESCALE)
#    print 'phase_model', phase_model

    if (phase_npy.shape[1] == pars.phase_npy_1hr.shape[1]) and (RESCALE == True):
        phase_matrix_origin = scale_matrix(phase_npy, pars.SCALE)
        print 'phase_matrix_origin.shape',phase_matrix_origin.shape
    else:
        phase_matrix_origin = phase_npy
    phase_matrix_new = np.zeros(phase_matrix_origin.shape)
#    print'phase_matrix_new.shape', phase_matrix_new.shape

    for ii in xrange(len(phase_matrix_new)):
        if aligned_1bin == True and aligned_fft == False:
            phase_matrix_new[ii] = np.roll(phase_matrix_origin[ii], -1 * int(round(phase_model[ii])))
        elif aligned_1bin == False and aligned_fft == True:
            phase_matrix_new[ii] = ifft(fft_phase_curve_inverse([phase_model[ii], 1], fft(phase_matrix_origin[ii]))).real
        else:
            print 'conflict setting of aligned_1bin and aligned_fft'
#    print 'phase_matrix_new', phase_matrix_new

#    print 'finish phase_matrix_new'
    U, s, V = np.linalg.svd(phase_matrix_new, full_matrices=True)

    if np.abs(np.amax(V[0])) < np.abs(np.amin(V[0])):
        V[0] = -V[0]

#    if (phase_npy.shape[1] == pars.phase_npy_1hr.shape[1]) and (RESCALE == True):
#        V_new = scale_matrix(V, pars.SCALE)
#        V = V_new


    return U, s, V, phase_model

def plot_two_temps(this_file_wz, bin_number_wz, phase_amp_bin_wz, phase_npy_wz, this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, phase_npy_1hr, RESCALE=None):

    U_wz, s_wz, V_wz, phase_model_wz = svd(this_file_wz, bin_number_wz, phase_amp_bin_wz, phase_npy_wz, pars.NPHASEBIN_wz, RESCALE)
    U_1hr_origin, s_1hr_origin, V_1hr_origin, phase_model_1hr = svd(this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, phase_npy_1hr, pars.NPHASEBIN_1hr, RESCALE)
    
    V0_wz = V_wz[0]
    '''resacle V0_1hr'''
    V_1hr = np.zeros((V_1hr_origin.shape[0], int((V_1hr_origin.shape[1])*SCALE)))  
    for ii in xrange(len(V_1hr)):
        for jj in xrange(len(V_1hr[0])):
            V_1hr[ii, jj] = np.average(V_1hr_origin[ii,int(jj/SCALE):int(jj/SCALE)+4])    

    phase_fit(0, V_wz, V_1hr, 'phase_wz_1hr_')
    phase_fit(0, V_1hr, V_wz, 'phase_1hr_wz_')

    plt.figure()
    plt.plot(np.arange(-100, 100), np.roll(V_wz[0]     , -100), 'r-',linewidth=2.5)
    plt.plot(np.arange(-100, 100), np.roll(V_1hr[0] -0.5, -100), 'b-',linewidth=2.5)
    plt.xlabel('phase bin number')
    plt.ylabel('amp')
    plt.savefig('phase_wz_and_1hr.png')

def plot_svd(this_file, bin_number, phase_amp_bin, phase_npy, plot_name, NPHASEBIN=None, RESCALE=None):

    U, s, V, _ = svd(this_file, bin_number, phase_amp_bin, phase_npy, NPHASEBIN, RESCALE)

    print 'len(V[0])', len(V[0])
    print 's.shape', s.shape

    plt.figure()
    x_range = np.arange(0, len(s))
    plt.plot(x_range, s, 'ro-')
    plt.xlabel('phase bin number')
    plt.ylabel('s values')
    plot_name_s = plot_name + '_s.png'
    plt.savefig(plot_name_s)

    plt.figure()
    n_step = -0.3
    x_range = np.arange(-len(V[0])/2, len(V[0])/2)
    color = ['r', 'b', 'g', 'k', 'y', '0.9', '0.7', '0.5', '0.3', '0.1']
    for ii in xrange(len(color)):
        plt.plot(x_range, np.roll(V[ii] + ii *n_step, -len(V[0])/2), color[ii], linewidth=1.0)
    plt.xlabel('phase bin number')
    plt.ylabel('V values')
    plot_name_V = plot_name + '_V.png'
    plt.savefig(plot_name_V)




  
if __name__ == '__main__':
    main()

