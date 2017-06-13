import numpy as np
import h5py
import svd
import fit_timing
import pars
import psr_svd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import fftpack, optimize, interpolate, linalg, integrate

REFERNCE_FREQ = pars.REFERNCE_FREQ
T = pars.T

NPHASEBIN_1hr = pars.NPHASEBIN_1hr
NPHASEBIN = NPHASEBIN_1hr #for pointed data
this_file = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/J2139+00_1hr_foldingh5','r')

series = ['DATA_FOLDING_0_255', 'DATA_FOLDING_256_511', 'DATA_FOLDING_512_767', 'DATA_FOLDING_768_1023', 'DATA_FOLDING_1024_1279', 'DATA_FOLDING_1280_1535', 'DATA_FOLDING_1536_1791', 'DATA_FOLDING_1792_2047', 'DATA_FOLDING_2048_2303', 'DATA_FOLDING_2304_2559', 'DATA_FOLDING_2561_2815', 'DATA_FOLDING_2816_3071', 'DATA_FOLDING_3072_3327', 'DATA_FOLDING_3328_3583', 'DATA_FOLDING_3584_3839', 'DATA_FOLDING_3840_4095', 'DATA_FOLDING_4096_4351', 'DATA_FOLDING_4352_4607', 'DATA_FOLDING_4608_4863', 'DATA_FOLDING_4864_5119', 'DATA_FOLDING_5120_5375', 'DATA_FOLDING_5376_5631', u'DATA_FOLDING_5632_5887', 'DATA_FOLDING_5888_6143', 'DATA_FOLDING_6144_6399', 'DATA_FOLDING_6400_6655', 'DATA_FOLDING_6656_6911', 'DATA_FOLDING_6912_7167', 'DATA_FOLDING_7168_7423', 'DATA_FOLDING_7424_7679', 'DATA_FOLDING_7680_7935', 'DATA_FOLDING_7936_8191', 'DATA_FOLDING_8192_8447', 'DATA_FOLDING_8448_8703', 'DATA_FOLDING_8704_8959', 'DATA_FOLDING_8960_9215', 'DATA_FOLDING_9216_9471', 'DATA_FOLDING_9472_9727', 'DATA_FOLDING_9728_9983', 'DATA_FOLDING_9984_10239']

model = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/pointed_phase_model.npy')

dis_stack = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/shift_profiles_stack.npy')
V_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/V_1hr_5sec_liksol.npy')
V0 = V_1hr[0]

plt.close('all')
plt.plot(dis_stack[1023])
plt.savefig('dis_stack1023.png')

plt.close('all')
plt.imshow(dis_stack)
plt.savefig('dis_stack0_all.png')

def main():
    if False:
        shift_profiles = np.zeros((len(series), this_file['DATA'].shape[3], NPHASEBIN_1hr))
        for ii in xrange(len(series)):
            shift_profiles[ii] = shift_phase_bin(model[ii], this_file[series[ii]][:,0,:,0].T)
        #plt.close('all')
        #plt.imshow(psr_svd.rebin_spec(shift_profiles[0], 32, 1))
        #plt.savefig('profile_0_shift.png')

    #stack profiles
    if False:
        profile_stack = 1
        nprof = len(shift_profiles)
        nprof -= nprof % profile_stack
        shift_profiles = shift_profiles[:nprof].reshape(nprof // profile_stack, profile_stack, this_file['DATA'].shape[3], NPHASEBIN_1hr)
        shift_profiles = np.mean(shift_profiles, 1)
        shift_profiles = np.mean(shift_profiles, 0)
        plt.imshow(psr_svd.rebin_spec(shift_profiles, 32, 1))
        plt.savefig('profiles_all_shifted.png')
        #np.save('shift_profiles.npy', shift_profiles)

    # fit DM
    V_fft = fftpack.fft(V0)
    # dis_stack is the stacking dispersion file, which is in the shpae of (frequencies, phase_bins)
#    for ii in xrange(len(dis_stack)):
    for ii in xrange(1023, 1024):
        freq = this_file['DAT_FREQ'][0][ii]
        print 'freq', freq
        profile_fft = fftpack.fft(dis_stack[ii])
        # pars_init: [Amplitude, factor of power law, random phase, DM]
#        pars_init = [np.amax(profile_fft.real)/np.amax(V_fft.real), -2., np.argmax(dis_stack[ii])/NPHASEBIN_1hr, 31.726]
        pars_init = [np.amax(profile_fft.real)/np.amax(V_fft.real), -2., 0.80, 31.726]
        pars_fit, cov, infodict, mesg, ier = optimize.leastsq(
                residuals,
                pars_init,
                (freq, profile_fft, V_fft),
                full_output=1,
                )
        fit_res = residuals(pars_fit, freq, profile_fft, V_fft)
        chi2_fit = chi2(pars_fit, freq, profile_fft, V_fft)
        dof = len(fit_res) - len(pars_init)
        red_chi2 = chi2_fit / dof
        print "chi1, dof, chi2/dof:", chi2_fit, dof, red_chi2

        cov_norm = cov * red_chi2

        errs = np.sqrt(cov_norm.flat[::len(pars_init) + 1])
        corr = cov_norm / errs[None,:] / errs[:,None]

        print "amp, factor, phase, DM:"
        print pars_fit
        print errs
        print "correlations:"
        print corr

        if True:
            plot_name = 'DM_fit_'+str(ii)+'_'
            fit_plot(pars_fit, pars_init, freq, profile_fft, V_fft, plot_name)

def chi2(parameters, freq, profile_fft, V_fft, norm=1):
    return np.sum(residuals(parameters, freq, profile_fft, V_fft)**2) * norm

def residuals(parameters, freq, profile_fft, V_fft):  
    res = pick_harmonics(profile_fft) - pick_harmonics(model(parameters, freq, V_fft))
    return res

def pick_harmonics(profile_fft):
#    harmonics = profile_fft[..., 1:NHARMONIC]
    harmonics = np.concatenate((profile_fft.real, profile_fft.imag), -1)
    return harmonics

def model(parameters, freq, V_fft):
    amp = parameters[0]
    factor = parameters[1]
    phi_0 = parameters[2]
    DM = parameters [3]
    phi_f = (phi_0 + dedisperse_time(DM, freq)/T) 
    template = amp * (freq/800.)**factor * shift_phase_bin_fft(phi_f * NPHASEBIN_1hr, V_fft)
    return template

def dedisperse_time(DM, freq):
    DM_CONST = 4148.808
    time = DM_CONST * DM * (freq**-2 - REFERNCE_FREQ**-2)
    return time

def shift_phase_bin_fft(phase_bin_shift, profile_fft):
    # profile_fft is an 1D array with fft already
    phase_bin_shift = phase_bin_shift % NPHASEBIN_1hr
    profile_shift_fft = np.zeros((profile_fft.shape), dtype=complex)
    n = profile_shift_fft.shape[-1]
    freq = fftpack.fftfreq(n, 1./n)
    phase = np.exp(-2j * np.pi * (-phase_bin_shift) / NPHASEBIN_1hr * freq)
    profile_shift_fft = profile_fft * phase
    return profile_shift_fft

def shift_phase_bin(phase_bin_shift, profile):
    # the input profile is required in the shape of (frequencies, phase_bins), without fft yet 
    profile_shift = np.zeros(profile.shape)
    profile_fft = np.zeros((profile.shape), dtype=complex)
    for ii in xrange(profile.shape[0]):   
        profile_fft[ii] = fftpack.fft(profile[ii])
        profile_shift[ii] = fftpack.ifft(shift_phase_bin_fft(phase_bin_shift, profile_fft[ii])).real
    return profile_shift 

def fit_plot(pars_fit, pars_init, freq, profile_fft, V_fft, plot_name):
    '''functions in Fourier space'''
    model_fft = model(pars_fit, freq, V_fft)
    data_fft = profile_fft
    init_fft = model(pars_init, freq, V_fft)

    '''Real part'''
    model_fft_real = model_fft.real
    data_fft_real = data_fft.real
    init_fft_real = init_fft.real
    res_fft_real = data_fft_real - model_fft_real
#    res_fft_real = np.concatenate((res[:(len(res)/2)], res[:(len(res)/2)][::-1]))

    '''Imag part'''
    model_fft_imag = model_fft.imag
    data_fft_imag = data_fft.imag
    init_fft_imag = init_fft.imag
    res_fft_imag = data_fft_imag - model_fft_imag
#    res_fft_imag = np.concatenate((res[(len(res)/2):], -res[(len(res)/2):][::-1]))

    '''functions in Real space (ifft)'''
    model_ifft = np.fft.ifft(model_fft).real
    data_ifft = np.fft.ifft(data_fft).real
    res_ifft = np.fft.ifft(data_fft - model_fft).real
    init_ifft = np.fft.ifft(init_fft).real

    freq_range = np.linspace(np.amin(np.fft.fftfreq(len(data_fft))), np.amax(np.fft.fftfreq(len(data_fft))), num=len(data_fft), endpoint=True)
    freq_min = np.amin(freq_range)
    freq_max = np.amax(freq_range)

#        phase_range = np.arange(-int(NPHASEBIN/2), int(NPHASEBIN/2)) / np.float(NPHASEBIN)
    phase_range = np.linspace(-0.5, 0.5, num=NPHASEBIN_1hr, endpoint=True)
    
#    title_phase = phase_amp_bin[0]/np.float(NPHASEBIN_wz)
#    if title_phase >= 0.5:
#        title_phase -= 1
#        title_phase_err = phase_amp_bin[len(phase_amp_bin)/2]/np.float(NPHASEBIN_wz)

#        plot_title = 'Phase: ' + str("%.5f" % (title_phase)) + ' +/- ' + str("%.5f" % (title_phase_err))
#    plot_name += str(index) + '_'
    fontsize = 16

    '''Plot for real and imag parts in the Fourier space.'''
    plt.close('all')
    f, ((ax1, ax2, ax3, ax4)) = plt.subplots(4, 1, sharex='col', figsize=(9, 16))
    #sharex='col', sharey='row'
    f.subplots_adjust(wspace=0.00, hspace=0.00)
    mode_freq = np.fft.fftfreq(len(data_fft))
    mode_range = np.linspace(-len(mode_freq)/2, len(mode_freq)/2, num=len(mode_freq), endpoint=True)
    xmax = np.amax(mode_range) #NHARMONIC
    xmin = 1 #np.amin(mode_range)
    ax1.plot(mode_range, np.roll(model_fft_real, -int(NPHASEBIN/2)),'r-')
    ax1.plot(mode_range, np.roll(data_fft_real, -int(NPHASEBIN/2)),'b-')
    ax1.set_title('Real', size=fontsize)
#        ax1.legend(['Real'], fontsize=fontsize)
    ax1.set_xlim([xmin,xmax])
#    ax1.set_ylim([-0.35, 0.2])
    ax1.tick_params(axis='both', which='major', labelsize=fontsize)

    ax2.plot(mode_range, np.roll(res_fft_real, -int(NPHASEBIN/2)),'bo')
    ax2.set_xlabel('Harmonic modes', fontsize=fontsize)
    ax2.set_ylabel('Residuals (T/Tsys)', fontsize=fontsize)
    ax2.set_xlim([xmin,xmax])
#        ax2.set_ylim([-0.05, 0.04])
    ax2.tick_params(axis='both', which='major', labelsize=fontsize)

    ax3.plot(mode_range, np.roll(model_fft_imag, -int(NPHASEBIN/2)),'r-')
    ax3.plot(mode_range, np.roll(data_fft_imag, -int(NPHASEBIN/2)),'b-')
    ax3.set_title('Imag', size=fontsize)
#        ax3.legend(['Imag'], fontsize=fontsize)
    ax3.set_xlim([xmin,xmax])
#        ax3.set_ylim([-0.35, 0.4])
    ax3.tick_params(axis='both', which='major', labelsize=fontsize)

    ax4.plot(mode_range, np.roll(res_fft_imag, -int(NPHASEBIN/2)),'bo')
    ax4.set_xlabel('Harmonic modes', fontsize=fontsize)
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
    ax1.plot(phase_range, np.roll(data_ifft, -int(NPHASEBIN/2)),'b-')
    ax1.plot(phase_range, np.roll(init_ifft, -int(NPHASEBIN/2)),'k-')
#    ax1.set_title(plot_title, size=fontsize)
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylabel('T/Tsys', fontsize=fontsize)
    ax1.tick_params(axis='both', which='major', labelsize=fontsize)

    ax2.plot(phase_range, np.roll(res_ifft, -int(NPHASEBIN/2)),'bo')
    ax2.set_xlim([xmin,xmax])
    ax2.set_xlabel('Phase ', fontsize=fontsize)
    ax2.set_ylabel('Residuals (T/Tsys)', fontsize=fontsize)
    ax2.tick_params(axis='both', which='major', labelsize=fontsize)

    plt.savefig(plot_name + 'ifft.png', bbox_inches='tight')



if __name__ == '__main__':
    main()
