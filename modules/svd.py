import sys
import os
import os.path

import h5py
import numpy as np
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors as mcolors
from scipy import fftpack, optimize, interpolate, linalg, integrate
from scipy.optimize import curve_fit, minimize, leastsq

import pars
import fit_timing


def main():
    try:
        ploting()
    except (IOError, ValueError):
        print IOError

aligned_1bin = False 
aligned_fft = True
save_fit_pars = False
random_res = False
save_random_res = False
save_phase_plot = True
phase_fit_lik = True
save_lik_pars = False

RA = pars.RA
DEC = pars.DEC
AU = pars.AU
C = pars.C
NPHASEBIN_wz = pars.NPHASEBIN_wz
NPHASEBIN = pars.NPHASEBIN
NPHASEBIN_1hr = pars.NPHASEBIN_1hr
NCENTRALBINS = pars.NCENTRALBINS
NHARMONIC = pars.NHARMONIC
NCENTRALBINSMAIN = pars.NCENTRALBINSMAIN
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
#            cut = abs(freq) > freq[NHARMONIC]
#            profile_ii_fft[cut] = 0
            for jj in xrange(len(profile[ii,:])):
                if (np.abs(np.fft.fftfreq(len(profile[ii,:]), 1./len(profile[ii,:]))) > 90)[jj] ==True:
                    profile_ii_fft[jj] = 0
            profile_fft += parameters[ii + 1] * profile_ii_fft
        
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
    res_Re = residuals_complex.real[1:NHARMONIC]
    res_Im = residuals_complex.imag[1:NHARMONIC]
    res = np.concatenate((res_Re, res_Im))
    return res

def random_residuals_data(parameters, model, data_fft, freq=None):
    '''Only use positive frequencies for residuals'''
    if freq is None:
        freq = None
    else:
        freq = freq

    fit_res = data_fft - fft_phase_curve(parameters, model, freq)
    random_res = (fit_res * np.exp(1.0j * 2 * np.pi * np.random.random()))[:len(fit_res)/2]
    random_res_Re = random_res.real
    random_res_Im = random_res.imag
    random_res = np.concatenate((random_res_Re, random_res_Im))
    # Take care here that the negative frequencies remain the conjugate of the positive ones.
    random_res_data = random_res + fft_phase_curve(parameters, model, freq)

    return random_res_data

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
    n_phase_bins = V.shape[0]
    V[1:,NCENTRALBINS//2:-NCENTRALBINS//2] = 0
    if False:
        spline = interpolate.splrep(
                np.arange(n_phase_bins - NCENTRALBINSMAIN),
                V[0,NCENTRALBINSMAIN//2:-NCENTRALBINSMAIN//2],
                s=0.0015,
                )
        V[0,NCENTRALBINSMAIN//2:-NCENTRALBINSMAIN//2] = interpolate.splev(
                np.arange(n_phase_bins - NCENTRALBINSMAIN),
                spline,
                )
#        plt.plot(np.roll(V[0], n_phase_bins//2))

    data_fft = fft(phase_matrix_origin[index])
    freq = np.fft.fftfreq(len(data_fft))    
#    freq[abs(freq) > freq[NHARMONIC]] = 0

    if len(model) > len(data_fft):
        amp_V0 = np.amax(phase_matrix_origin[index])/np.amax(V[0]) * 0.5
    else:
        amp_V0 = np.amax(phase_matrix_origin[index])/np.amax(V[0])

    if phase_matrix_origin.shape[1] == pars.NPHASEBIN_wz:
        model_phase_bin = model_phase_bin_wz
    elif phase_matrix_origin.shape[1] == pars.NPHASEBIN_1hr:   
        model_phase_bin = model_phase_bin_1hr

    pars_init = [ model_phase_bin[index], amp_V0, amp_V0/10, amp_V0/50000, amp_V0/500000, amp_V0/500000, amp_V0/500000]

#    pars_init = [ np.argmax(phase_matrix_origin[index]) % NPHASEBIN, amp_V0, amp_V0/10, amp_V0/5000, amp_V0/500000]
    print 'pars_init', pars_init

#    model_fft = fft(V[0])
    fit_pars_phase, pcov, infodict, errmsg, success = leastsq(residuals, pars_init,
                   args=( model, data_fft, freq), full_output=1)


    print "Fit parameters: ", fit_pars_phase
    print "sucess?:", success
    '''DOF -1, since we set fft_file[0] = 0.'''
    res = residuals(fit_pars_phase, model, data_fft, freq)
    print 'len(res)',len(res)
    chi2_value = np.sum(np.abs(res)**2)
    dof = len(res)-len(pars_init)
    red_chi2 = chi2_value / dof
    print "Chi-squared: ", chi2_value, ", DOF: ", dof,' ,red_chi2:', red_chi2

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

    '''save the fitting amp and bin as [bin, amps, bin_err, amp_errs]'''
    npy_file = 'phase_amp_bin_wz_40epochs_6modes_90hars.npy'

    if (phase_matrix_origin.shape[1] == pars.phase_npy_1hr.shape[1]):
        phase_amp_bin = np.concatenate(([pfit_leastsq[0]*SCALE], pfit_leastsq[1:], [perr_leastsq[0]*SCALE], perr_leastsq[1:] ))
    else:
        phase_amp_bin = np.concatenate((pfit_leastsq, perr_leastsq ))

    if save_fit_pars== True:
        if os.path.exists(npy_file):
            sequence = np.load(npy_file)
            np.save(npy_file, np.vstack((sequence, phase_amp_bin)))
        else:
            np.save(npy_file, phase_amp_bin)

    '''Simulation test on random residuals phase'''
    if random_res == True:
        random_times = 100
        '''save the random res fitting amp and bin as [bin, amps, bin_err, amp_errs]'''
        random_npy_file = 'random_'+ plot_name + '.npy'
        if os.path.exists(random_npy_file) == False:
            random_phase_amp_bin = np.zeros((phase_matrix_origin.shape[0], random_times, len(fit_pars_phase)))
        else:
            random_phase_amp_bin = np.load(random_npy_file)

        for ii in xrange(random_times):
            print 'ii', ii
            random_res_data = random_residuals_data(fit_pars_phase, model, data_fft, freq)
            random_fit_pars_phase, pcov, infodict, errmsg, success = leastsq(residuals, 
                       pars_init, args=( model, random_res_data, freq), full_output=1)
            print 'random_fit_pars_phase', random_fit_pars_phase

            if (phase_matrix_origin.shape[1] == pars.phase_npy_1hr.shape[1]):
                random_phase_amp_bin[index, ii, ...] = np.concatenate(([random_fit_pars_phase[0]*SCALE], random_fit_pars_phase[1:] ))
            else:
                random_phase_amp_bin[index, ii, ...] = random_fit_pars_phase

        if save_random_res == True:
            np.save(random_npy_file, random_phase_amp_bin)
            
    if save_phase_plot == True:
        '''functions in Fourier space'''
        model_fft = fft_phase_curve(fit_pars_phase, model, freq)  
        init_fft = fft_phase_curve(pars_init, model, freq)

        '''Real part'''
        model_fft_real = model_fft.real
        data_fft_real = data_fft.real
        init_fft_real = init_fft.real
        res_fft_real = data_fft_real - model_fft_real
#        res_fft_real = np.concatenate((res[:(len(res)/2)], res[:(len(res)/2)][::-1]))

        '''Imag part'''
        model_fft_imag = model_fft.imag
        data_fft_imag = data_fft.imag
        init_fft_imag = init_fft.imag
        res_fft_imag = data_fft_imag - model_fft_imag
#        res_fft_imag = np.concatenate((res[(len(res)/2):], -res[(len(res)/2):][::-1]))

        '''functions in Real space (ifft)'''
        model_ifft = np.fft.ifft(model_fft).real
        data_ifft = np.fft.ifft(data_fft).real
        res_ifft = np.fft.ifft(data_fft - model_fft).real
        init_ifft = np.fft.ifft(init_fft).real

        freq_range = np.linspace(np.amin(np.fft.fftfreq(len(data_fft))), np.amax(np.fft.fftfreq(len(data_fft))), num=len(data_fft), endpoint=True)
        freq_min = np.amin(freq_range)
        freq_max = np.amax(freq_range)

#        phase_range = np.arange(-int(NPHASEBIN/2), int(NPHASEBIN/2)) / np.float(NPHASEBIN)
        phase_range = np.linspace(-0.5, 0.5, num=NPHASEBIN, endpoint=True)
    
        title_phase = phase_amp_bin[0]/np.float(NPHASEBIN_wz)
        if title_phase >= 0.5:
            title_phase -= 1
        title_phase_err = phase_amp_bin[len(phase_amp_bin)/2]/np.float(NPHASEBIN_wz)

        plot_title = 'Phase: ' + str("%.5f" % (title_phase)) + ' +/- ' + str("%.5f" % (title_phase_err))
        plot_name += str(index) + '_'
        fontsize = 16
        markersize = 5.0
        markersize_dot = 8.0

        '''Plot for real and imag parts in the Fourier space.'''
        plt.close('all')
        f, ((ax1, ax2, ax3, ax4)) = plt.subplots(4, 1, sharex='col', figsize=(9, 16))
        #sharex='col', sharey='row'
        f.subplots_adjust(wspace=0.00, hspace=0.00)
#        plt.figure(0, figsize=(9,16))
#        ax1 = plt.subplot2grid((4,1), (0,0))
#        ax2 = plt.subplot2grid((4,1), (1,0), sharex=ax1, hspace = 0.05)
#        ax3 = plt.subplot2grid((4,1), (2,0), sharex=ax1)
#        ax4 = plt.subplot2grid((4,1), (3,0), sharex=ax1, hspace = 0.05)

        mode_range = np.linspace(-len(freq)/2, len(freq)/2, num=len(freq), endpoint=True)
        xmax = NHARMONIC  #np.amax(mode_range)
        xmin = 1 #np.amin(mode_range)
        ax1.plot(mode_range, np.roll(model_fft_real, -int(NPHASEBIN/2)),'r-')
        ax1.plot(mode_range, np.roll(data_fft_real, -int(NPHASEBIN/2)),'b.', markersize=markersize_dot)
        ax1.set_title('Real', size=fontsize)
#        ax1.legend(['Real'], fontsize=fontsize)
        ax1.set_ylabel('T/Tsys', fontsize=fontsize)
        ax1.set_xlim([xmin,xmax])
        ax1.set_ylim([-0.35, 0.2])
        ax1.tick_params(axis='both', which='major', labelsize=fontsize)

        ax2.plot(mode_range, np.roll(res_fft_real, -int(NPHASEBIN/2)),'gs', markersize=markersize)
        ax2.set_xlabel('Harmonic modes', fontsize=fontsize)
        ax2.set_ylabel('Residuals', fontsize=fontsize)
        ax2.set_xlim([xmin,xmax])
        ax2.set_ylim([-0.05, 0.04])
        ax2.tick_params(axis='both', which='major', labelsize=fontsize)

        ax3.plot(mode_range, np.roll(model_fft_imag, -int(NPHASEBIN/2)),'r-')
        ax3.plot(mode_range, np.roll(data_fft_imag, -int(NPHASEBIN/2)),'b.', markersize=markersize_dot)
        ax3.set_title('Imag', size=fontsize)
#        ax3.legend(['Imag'], fontsize=fontsize)
        ax3.set_ylabel('T/Tsys', fontsize=fontsize)
        ax3.set_xlim([xmin,xmax])
        ax3.set_ylim([-0.35, 0.4])
        ax3.tick_params(axis='both', which='major', labelsize=fontsize)

        ax4.plot(mode_range, np.roll(res_fft_imag, -int(NPHASEBIN/2)),'gs', markersize=markersize)
        ax4.set_xlabel('Harmonic modes', fontsize=fontsize)
        ax4.set_ylabel('Residuals', fontsize=fontsize)
        ax4.set_xlim([xmin,xmax])
        ax4.tick_params(axis='both', which='major', labelsize=fontsize)
        
        plt.savefig(plot_name + 'fft.png', bbox_inches='tight')

        '''Plot for real part in real space'''
        plt.close('all')
        f, ((ax1, ax2)) = plt.subplots(2, 1, sharex='col', figsize=(8,9))
        f.subplots_adjust(hspace=0.07)
        xmax = np.amax(phase_range)
        xmin = np.amin(phase_range)
        ax1.plot(phase_range, np.roll(model_ifft, -int(NPHASEBIN/2)),'r-')
        ax1.plot(phase_range, np.roll(data_ifft, -int(NPHASEBIN/2)),'b.', markersize=markersize_dot)
        ax1.set_title(plot_title, size=fontsize)
        ax1.set_xlim([xmin,xmax])
        ax1.set_ylabel('T/Tsys', fontsize=fontsize)
        ax1.tick_params(axis='both', which='major', labelsize=fontsize)

        ax2.plot(phase_range, np.roll(res_ifft, -int(NPHASEBIN/2)),'gs', markersize=markersize)
        ax2.set_xlim([xmin,xmax])
        ax2.set_xlabel('Phase ', fontsize=fontsize)
        ax2.set_ylabel('Residuals (T/Tsys)', fontsize=fontsize)
        ax2.tick_params(axis='both', which='major', labelsize=fontsize)

        plt.savefig(plot_name + 'ifft.png', bbox_inches='tight')

    if phase_fit_lik == True:
        # fitting phase by likeli`hood
#        phase_diff_samples = np.arange(-50, 50, 0.3) * perr_leastsq[0]
#        phase_diff_samples = np.arange(-50, 50, 1)
        phase_diff_samples = np.arange(-20, 20, 0.02)
        chi2_samples = []
        for p in phase_diff_samples:
            this_phase = p + pars_init[0]
            if True:
                P = shift_trunc_modes(this_phase, model, NPHASEBIN)
#                print 'P',P
                d = pick_harmonics(data_fft)
#                print 'd',d
                N = linalg.inv(np.dot(P, P.T))
                this_pars_fit = np.dot(N, np.sum(P * d, 1))
#                print 'this_pars_fit', this_pars_fit
#            print '[this_phase] + list(this_pars_fit)', [this_phase] + list(this_pars_fit)
            chi2_sample = chi2(
                        [this_phase] + list(this_pars_fit),
                        model, data_fft, NPHASEBIN,
                        1. / red_chi2,
                        )
#            print 'chi2_sample', chi2_sample
            chi2_samples.append(chi2_sample)
        chi2_samples = np.asarray(chi2_samples)


#        print 'chi2_samples', chi2_samples
#        phase_diff_samples = np.array(phase_diff_samples[np.where(chi2_samples< 20*np.amin(chi2_samples))])
#        print 'phase_diff_samples', phase_diff_samples
#        chi2_samples = chi2_samples[np.where(chi2_samples< 20*np.amin(chi2_samples))]
#        print 'chi2_samples', chi2_samples      

#        plot_name += str(index) + '_'
#        plt.close('all')
#        plt.plot(phase_diff_samples, chi2_samples)
#        plt.xlabel('Phase_diff')
#        plt.ylabel('Chi2')
#        plt.savefig(plot_name+'phase_chi2.png', bbox_inches='tight')

        # Integrate the full liklihood, taking first and second moments to
        # get mean phase and variance.
        likelihood = np.exp(-chi2_samples / 2) 
        
        print 'likelihood',likelihood
#        norm1 = np.sum(likelihood)
        norm = simpson(likelihood, 0, len(likelihood)-1)
#        norm = integrate.simps(likelihood)
#        print 'norm1', norm1
        print 'norm',norm
#        mean = np.sum(phase_diff_samples * likelihood) / norm
        mean = simpson(phase_diff_samples * likelihood, 0, len(likelihood)-1) / norm 
#        mean = integrate.simps(phase_diff_samples * likelihood) / norm
        print 'mean', mean
#        var = np.sum(phase_diff_samples**2 * likelihood) / norm - mean**2
        var = simpson(np.array(phase_diff_samples)**2 * likelihood, 0, len(likelihood)-1) / norm - mean**2
#        var = integrate.simps(phase_diff_samples**2 * likelihood) / norm - mean**2
        std = np.sqrt(var)
        print 'std',std
        print "Integrated Liklihood:", pars_init[0] + mean, std
        phases_lik = []
        phase_errors_lik = []
        phases_lik.append(pars_init[0] + mean)
        phase_errors_lik.append(std)

#        print np.where((likelihood / norm)>np.amax(likelihood / norm) * 10**-4)
#        print np.where((likelihood / norm)>np.amax(likelihood / norm) * 10**-4)[0][0]
#        print np.where((likelihood / norm)>np.amax(likelihood / norm) * 10**-4)[0][-1]
        print 'max:', np.amax(likelihood / norm / (0.02/NPHASEBIN))

        plot_name += str(index) + '_'
        plt.close('all')
        phase_diff_range = np.linspace(np.amin(phase_diff_samples)/NPHASEBIN, np.amax(phase_diff_samples)/NPHASEBIN, num=len(phase_diff_samples), endpoint=True)
        plt.semilogy(phase_diff_range, likelihood / norm  / (0.02/NPHASEBIN))
        plt.xlabel('Phase', fontsize=fontsize)
        plt.ylabel('log(Likelihood)', fontsize=fontsize)
        plt.xlim((phase_diff_range[np.where((likelihood / norm  / (0.02/NPHASEBIN))>np.amax(likelihood / norm  / (0.02/NPHASEBIN)) * 10**-4)[0][0]],phase_diff_range[np.where((likelihood / norm  / (0.02/NPHASEBIN))>np.amax(likelihood / norm  / (0.02/NPHASEBIN)) * 10**-4)[0][-1]]))
        plt.ylim((np.amax(likelihood / norm / (0.02/NPHASEBIN)) * 10**-4, np.amax(likelihood / norm / (0.02/NPHASEBIN))*4.5))
        plt.tick_params(axis='both', which='major', labelsize=fontsize)
        plt.savefig(plot_name+'phase_chi2.png', bbox_inches='tight')



        '''save the fitting amp and bin as [bin, amps, bin_err, amp_errs]'''
        npy_lik_file = 'phase_amp_bin_57178_40epochs_6modes_90hars_lik.npy'

        if (phase_matrix_origin.shape[1] == pars.phase_npy_1hr.shape[1]):
            phase_amp_bin_lik = np.concatenate(([phases_lik[0]*SCALE], pfit_leastsq[1:], [phase_errors_lik[0]*SCALE], perr_leastsq[1:]))
        else:
            phase_amp_bin_lik = np.concatenate(([phases_lik[0]], pfit_leastsq[1:], [phase_errors_lik[0]], perr_leastsq[1:]))

        if save_lik_pars== True:
            if os.path.exists(npy_lik_file):
                sequence = np.load(npy_lik_file)
                np.save(npy_lik_file, np.vstack((sequence, phase_amp_bin_lik)))
            else:
                np.save(npy_lik_file, phase_amp_bin_lik)

def simpson(func, a, b):
    '''a is init, b is final, and n is number of steps.'''

    s = func[a] + func[b]
    n = len(func)
    h = np.float((b - a)) / n
    s += (2*np.sum(func[::2]) + 4*np.sum(func[1::2]))

    return s * h / 3

def shift_trunc_modes(phase_shift, model, NPHASEBIN):
    model_fft = fftpack.fft(model, axis=1)
#    model_fft = fft(model)
#    print 'model_fft[0][0]', model_fft[0][0]
    model_shift = apply_phase_shift(model_fft, phase_shift, NPHASEBIN)
    V_harmonics = pick_harmonics(model_shift)
    return V_harmonics[:6]


def chi2(parameters, model, data_fft, NPHASEBIN, norm=1):
    return np.sum(residuals_lik(parameters, model, data_fft, NPHASEBIN)**2) * norm

def model_lik(parameters, V, NPHASEBIN):
    phase = parameters[0]
    amplitudes = np.array(parameters[1:])
    shifted_modes = shift_trunc_modes(phase, V, NPHASEBIN)
    return np.sum(amplitudes[:,None] * shifted_modes, 0)

def residuals_lik(parameters, model, profile_fft, NPHASEBIN):
    return pick_harmonics(profile_fft) - model_lik(parameters, model, NPHASEBIN)

def pick_harmonics(profile_fft):
    harmonics = profile_fft[..., 1:NHARMONIC]
    harmonics = np.concatenate((harmonics.real, harmonics.imag), -1)
    return harmonics

def apply_phase_shift(profile_fft, phase_shift, NPHASEBIN):
    "Parameter *phase_shift* takes values [0 to 1)."

    n = profile_fft.shape[-1]
#    print 'n', n
    freq = fftpack.fftfreq(n, 1./n)
    phase = np.exp(-2j * np.pi * phase_shift /NPHASEBIN  * freq)
#    print 'phase',phase
    return profile_fft * phase 

def scale_array(old_array, SCALE):
    array = np.zeros(int(len(old_array)*SCALE))
    for ii in xrange(len(array)):
        array[ii] = np.average(old_array[int(ii/SCALE):int((ii+1)/SCALE)])
    return array

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
#    phase_model_old = fit_timing.timing_model_1(pars.old_fit_pars, time_mjd, dBATdra, dBATddec, NPHASEBIN, RESCALE)

#    phase_diff = phase_model - phase_model_old
#    print 'diff max',np.amax(phase_diff)
#    print 'diff min',np.amin(phase_diff)

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
    print 'phase_matrix_new.shape', phase_matrix_new.shape

    if False:
        fontsize = 14
        plt.figure()
        plt.close('all')
        plt.figure()
        n_step = -0.2
        x_range = np.arange(-len(phase_matrix_new[0])/2 , len(phase_matrix_new[0])/2) / np.float(NPHASEBIN)
        color = ['r', 'g', 'b', 'y', 'c', '0.0', '0.2', '0.4', '0.6', '0.8']
#    color = ['r', 'g', 'b']
        for ii in xrange(len(color)):
            plt.plot(x_range, np.roll(phase_matrix_new[ii] + ii *n_step, -len(phase_matrix_new[0])/2), color[ii], linewidth=1.0)
            plt.xlim((-0.5, 0.5))
            plt.xlabel('Phase', fontsize=fontsize)
            plt.ylabel('amp', fontsize=fontsize)
            plt.tick_params(axis='both', which='major', labelsize=fontsize)
            plt.savefig('phase_npy_1hr_5sec.png', bbox_inches='tight')

#    print 'finish phase_matrix_new'
    U, s, V = np.linalg.svd(phase_matrix_new, full_matrices=True)

    if np.abs(np.amax(V[0])) < np.abs(np.amin(V[0])):
        V[0] = -V[0]

#    if (phase_npy.shape[1] == pars.phase_npy_1hr.shape[1]) and (RESCALE == True):
#        V_new = scale_matrix(V, pars.SCALE)
#        V = V_new

    if len(V[0]) == NPHASEBIN_1hr:
#        profile_stack = 16
#        nprof = len(V)
#        nprof -= nprof % profile_stack
#        V = V[:nprof].reshape(nprof // profile_stack, profile_stack, NPHASEBIN_1hr)
#        V = np.mean(V, 1)
        n_phase_bins = V.shape[0]
        V[1:,NCENTRALBINS//2:-NCENTRALBINS//2] = 0
        if True:
            spline = interpolate.splrep(
                    np.arange(n_phase_bins - NCENTRALBINSMAIN),
                    V[0,NCENTRALBINSMAIN//2:-NCENTRALBINSMAIN//2],
                    s=0.0015,
                    )
            V[0,NCENTRALBINSMAIN//2:-NCENTRALBINSMAIN//2] = interpolate.splev(
                    np.arange(n_phase_bins - NCENTRALBINSMAIN),
                    spline,
                    )

    np.save('V_1hr_5sec_liksol.npy', V)
    return U, s, V, time_mjd, phase_model

def plot_two_temps(this_file_wz, bin_number_wz, phase_amp_bin_wz, phase_npy_wz, this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, phase_npy_1hr, RESCALE=None):

    U_wz, s_wz, V_wz, phase_model_wz = svd(this_file_wz, bin_number_wz, phase_amp_bin_wz, phase_npy_wz, pars.NPHASEBIN_wz, RESCALE)
    U_1hr_origin, s_1hr_origin, V_1hr_origin, phase_model_1hr = svd(this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, phase_npy_1hr, pars.NPHASEBIN_1hr, RESCALE)
    
#    V0_wz = V_wz[0]
#    '''resacle V0_1hr'''
#    V_1hr = np.zeros((V_1hr_origin.shape[0], int((V_1hr_origin.shape[1])*SCALE)))  
#    for ii in xrange(len(V_1hr)):
#        for jj in xrange(len(V_1hr[0])):
#            V_1hr[ii, jj] = np.average(V_1hr_origin[ii,int(jj/SCALE):int(jj/SCALE)+4])    

    phase_fit(0, V_wz, V_1hr_origin, 'phase_wz_1hr_')
#    phase_fit(0, V_1hr_origin, V_wz, 'phase_1hr_wz_')


    plt.figure()
    plt.plot(np.arange(-100, 100), np.roll(V_wz[0]     , -100), 'r-',linewidth=2.5)
    plt.plot(np.arange(-100, 100), np.roll(V_1hr[0] -0.5, -100), 'b-',linewidth=2.5)
    plt.xlabel('phase bin number')
    plt.ylabel('amp')
    plt.savefig('phase_wz_and_1hr.png')

def plot_svd(this_file, bin_number, phase_amp_bin, phase_npy, plot_name, NPHASEBIN=None, RESCALE=None):

    U, s, V, time_mjd, _ = svd(this_file, bin_number, phase_amp_bin, phase_npy, NPHASEBIN, RESCALE)

    print 'len(V[0])', len(V[0])
    print 's.shape', s.shape

    fontsize = 16

    plt.figure()
    x_range = np.arange(0, len(s))
    plt.plot(x_range, s, 'ro-')
    plt.xlabel('Phase')
    plt.ylabel('s values')
    plot_name_s = plot_name + '_s.png'
    plt.savefig(plot_name_s)

    plt.close('all')
    plt.figure()
    n_step = -0.3
    x_range = np.arange(-len(V[0])/2 , len(V[0])/2) / np.float(NPHASEBIN)
    color = ['r', 'g', 'b', 'y', 'c', '0.0', '0.2', '0.4', '0.6', '0.8']
#    color = ['r', 'g', 'b']
    for ii in xrange(len(color)):
        plt.plot(x_range, np.roll(V[ii] + ii *n_step, -len(V[0])/2), color[ii], linewidth=1.0)
    plt.xlim((-0.2, 0.2))
    plt.xlabel('Phase', fontsize=fontsize)
    plt.ylabel('V values', fontsize=fontsize)
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plot_name_V = plot_name + '_V.png'
    plt.savefig(plot_name_V, bbox_inches='tight')

    plt.close('all')
    plt.figure()
    n_step = -0.2
    x_range = np.arange(-len(V[0])/2 , len(V[0])/2) / np.float(NPHASEBIN)
#    color = ['r', 'g', 'b', 'y', 'c', '0.0', '0.2', '0.4', '0.6', '0.8']
    color = ['r', 'g', 'b']
    for ii in xrange(len(color)):
        plt.plot(x_range, np.roll(V[ii] + ii *n_step, -len(V[0])/2), color[ii], linewidth=1.0)
    plt.xlim((-0.5, 0.5))
    plt.xlabel('Phase', fontsize=fontsize)
    plt.ylabel('V values', fontsize=fontsize)
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plot_name_V = plot_name + '_V_zoom.png'
    plt.savefig(plot_name_V, bbox_inches='tight')

    '''Plot U function'''
    ut = U.T
#    np.save('1hr_5sec_ut.npy', ut)
    print 'ut.shape', ut.shape
    print 'len(time_mjd)', len(time_mjd)
    plt.close('all')
    plt.figure()
    n_step = -0.2
#    x_range = np.arange(-len(ut[0])/2 , len(ut[0])/2) / np.float(len(ut[0]))
    x_range = time_mjd
    color = ['r', 'g', 'b', 'y', 'c', '0.0', '0.2', '0.4', '0.6', '0.8']
#    color = ['r', 'g', 'b']
    for ii in xrange(len(color)):
        plt.plot(x_range, np.roll(ut[ii] + ii *n_step, 0), color[ii], linewidth=1.0)
    plt.xlabel('Time (MJD)', fontsize=fontsize)
    plt.ylabel('U values', fontsize=fontsize)
    plt.xlim((time_mjd[0], time_mjd[639]))
    plt.xticks([time_mjd[0], time_mjd[320], time_mjd[639]])
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plot_name_U = plot_name + '_U.png'
    plt.ticklabel_format(useOffset=False)
    plt.savefig(plot_name_U, bbox_inches='tight')


def fft_plot(npy_file):
    data_fft = fft(npy_file)
    data_fft_re = data_fft.real 
    data_fft_im = data_fft.imag

    SCALE_factor = T/(8.192e-5*4096*320)
    scale_data_fft_re = scale_array(data_fft_re, SCALE_factor)
    scale_data_fft_im = scale_array(data_fft_im, SCALE_factor)

    freq = np.fft.fftfreq(len(data_fft))
    freq_scale = np.fft.fftfreq(len(scale_data_fft_re))

    freq_range = np.linspace(np.amin(np.fft.fftfreq(len(data_fft))), np.amax(np.fft.fftfreq(len(data_fft))), num=len(data_fft), endpoint=True)
    freq_min = np.amin(freq_range)
    freq_max = np.amax(freq_range)

    freq_range_scale = np.linspace(np.amin(np.fft.fftfreq(len(scale_data_fft_re))), np.amax(np.fft.fftfreq(len(scale_data_fft_re))), num=len(scale_data_fft_re), endpoint=True)
    freq_min_scale = np.amin(freq_range_scale)
    freq_max_scale = np.amax(freq_range_scale)

    mode_range = np.linspace(-len(freq)/2, len(freq)/2, num=len(freq), endpoint=True)
    mode_range_scale = np.linspace(-len(freq_scale)/2, len(freq_scale)/2, num=len(freq_scale), endpoint=True)

    xmax = np.amax(mode_range)
    xmin = np.amin(mode_range)
    xmax_scale = np.amax(mode_range_scale)
    xmin_scale = np.amin(mode_range_scale)

    '''Plotting part'''
    plt.close('all')
    fontsize = 16
    f, ((ax1, ax2)) = plt.subplots(2, 1, sharex='col', figsize=(8,9))
    f.subplots_adjust(hspace=0.07)

#    ax1.plot(mode_range, np.roll(data_fft_re, -int(len(freq)/2)),'r-')
    ax1.plot(mode_range_scale, np.roll(scale_data_fft_re, -int(len(freq_scale)/2)),'b-', linewidth=1)
    ax1.set_title('Real', size=fontsize)
    ax1.set_xlim([xmin_scale,xmax_scale])
    ax1.tick_params(axis='both', which='major', labelsize=fontsize)

#    ax2.plot(mode_range, np.roll(data_fft_im, -int(len(freq)/2)),'r-')
    ax2.plot(mode_range_scale, np.roll(scale_data_fft_im, -int(len(freq_scale)/2)),'b-', linewidth=1)
    ax2.set_title('Imag', size=fontsize)
    ax2.set_xlim([xmin_scale,xmax_scale])
    ax2.set_xlabel('Harmonic Modes ', fontsize=fontsize)
#    ax2.set_ylabel('Residuals (T/Tsys)', fontsize=fontsize)
    ax2.tick_params(axis='both', which='major', labelsize=fontsize)

    plt.savefig('pointed_0_319_fft_scale.png', bbox_inches='tight')


  
if __name__ == '__main__':
    main()

