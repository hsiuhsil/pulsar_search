import sys
import os
import os.path

import numpy as np
import svd
import fit_timing
import pars
import psr_svd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from scipy import optimize
from scipy.optimize import curve_fit, minimize, leastsq


def main():
    try:
        compare_integration_lik()        
#        generate_funcitons()
#        justify_svd()
#        two_pulses(32, 25)
    except (IOError, ValueError):
        print IOError

this_file_wz = pars.this_file_wz
bin_number_wz = pars.bin_number_wz
phase_amp_bin_wz = pars.phase_amp_bin_wz
phase_npy_wz = pars.phase_npy_wz

this_file_1hr = pars.this_file_1hr
bin_number_1hr = pars.bin_number_1hr
bin_number_1hr_5sec = pars.bin_number_1hr_5sec
phase_amp_bin_1hr = pars.phase_amp_bin_1hr
phase_amp_bin_1hr_5sec = pars.phase_amp_bin_1hr_5sec
phase_npy_1hr = pars.phase_npy_1hr
phase_npy_1hr_5sec = pars.phase_npy_1hr_5sec
random_res_1hr = pars.random_res_1hr

#U_1hr, s_1hr, V_1hr, phase_model_1hr = svd.svd(this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, phase_npy_1hr, pars.NPHASEBIN_1hr, RESCALE=None)

V_1hr = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/V_1hr_5sec_liksol.npy')

#print 'U_1hr.shape', U_1hr.shape
#print 's_1hr.shape', s_1hr.shape
#print 'V_1hr.shape', V_1hr.shape


def U_analysis():
    '''import U.T of pointed observation data'''
    u_data = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/1hr_5sec_ut.npy')
    du_0 = np.gradient(u_data[0])
    du_1 = np.gradient(u_data[1])
    template = np.vstack((du_0, du_1))
    psr_svd.phase_fitting(u_data, template)

def generate_funcitons():

#    U_wz, s_wz, V_wz, phase_model_wz = svd.svd(this_file_wz, bin_number_wz, phase_npy_wz, NPHASEBIN_wz)
#    U_1hr_origin, s_1hr_origin, V_1hr_origin, phase_model_1hr = svd.svd(this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, pars.phase_npy_1hr, pars.NPHASEBIN_1hr, RESCALE=None)
    

#    svd.plot_svd(this_file_wz, bin_number_wz, phase_amp_bin_wz, phase_npy_wz, 'phase_wz',pars.NPHASEBIN_wz, RESCALE=None)
    svd.plot_svd(this_file_1hr, bin_number_1hr_5sec, phase_amp_bin_1hr_5sec, phase_npy_1hr_5sec, 'phase_1hr_5sec', pars.NPHASEBIN_1hr, RESCALE=None)
#    svd.plot_two_temps(this_file_wz, bin_number_wz, phase_amp_bin_wz, phase_npy_wz, this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, phase_npy_1hr, RESCALE=None)

#    np.save('V_pointed_5sec.npy', V_1hr)
  
#    V_1hr_scaled = svd.scale_matrix(V_1hr, pars.SCALE)

    '''fitting phase for pointed data'''
#    phase_npy_57178_40epochs = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/phase_57178_40epochs.npy')
#    svd.phase_fit(0, phase_npy_1hr, V_1hr, 'fitting_phase_fft_57178_', pars.NPHASEBIN_1hr)
#    svd.phase_fit(39, phase_npy_1hr, V_1hr, 'fitting_phase_fft_57178_', pars.NPHASEBIN_1hr)
#    for ii in xrange(40):
#    for ii in xrange(len(s_1hr)):
#        print 'ii=' + str(ii)
#        svd.phase_fit(ii, phase_npy_1hr, V_1hr, 'fitting_phase_fft_57178_40epochs_', pars.NPHASEBIN_1hr)

    '''fitting phase for WZ data'''
#    svd.phase_fit(2, phase_npy_wz, V_1hr, 'fitting_phase_fft_wz_', pars.NPHASEBIN_wz)
#    svd.phase_fit(9, phase_npy_wz, V_1hr, 'fitting_phase_fft_wz_', pars.NPHASEBIN_wz)
#    svd.phase_fit(197, phase_npy_wz, V_1hr, 'fitting_phase_fft_wz_', pars.NPHASEBIN_wz)
#    for ii in xrange(0, 236):
#        print 'ii= '+str(ii)
#        svd.phase_fit(ii, phase_npy_wz, V_1hr, 'fitting_phase_fft_wz_40epochs_6modes_90hars_', pars.NPHASEBIN_wz)

#    pointed_0_319_npy = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/pointed_0_319.npy')
   
#    svd.fft_plot(pointed_0_319_npy)

def justify_svd():
    '''Justify SVD is a better way'''
    all_profiles = pars.phase_npy_1hr_5sec # in the shape of (640, 800), (# of profiles, bin numbers)
    stack_profiles = pars.phase_npy_1hr_5sec_stack_40 # in the shape of (40, 800), (# of profiles, bin numbers)
    for ii in xrange(40):
        print 'ii: ', ii
        init = ii*16 
        final = (ii+1)*16

        if False:
            '''Firstly, average 16 5sec profiles and fit'''
            profile_stack = 16
            nprof = 640
            nprof -= nprof % profile_stack
            profile = profiles[int(init) : int(final)].reshape(nprof // profile_stack, profile_stack, pars.NPHASEBIN_1hr)
            profile = np.mean(profile, 1)
            svd.phase_fit(ii, stack_profiles, V_1hr, 'fitting_phase_fft_57178_stack_', pars.NPHASEBIN_1hr)
            '''Secondly, individually fit the 16 profiles, and add the chi-squared curves.'''
            for jj in xrange(init, final):
                svd.phase_fit(jj, all_profiles, V_1hr, 'fitting_phase_fft_57178_5sec_', pars.NPHASEBIN_1hr)
        '''lik1 is the first way, and lik2 is the second way'''
        lik_npy_file = 'lik_57178_5sec_stack_' + str(ii) + '.npy'
        sum_lik_file = 'sum_lik_57178_5sec_seg_' + str(ii) + '.npy'
        lik1 = np.load(lik_npy_file) #pars.lik_57178_5sec_stack0
        lik2_origin = np.load(sum_lik_file) #pars.sum_lik_57178_5sec_16 
        lik2 = np.ones(lik2_origin.shape[1], dtype=np.float128)
        for kk in xrange(len(lik2_origin)):
            lik2 *= lik2_origin[kk]

        print 'lik2', lik2
        plot_name = 'two_lik_phase_chi2_' + str(ii) + '.png'
        plot_two_likelihoods(lik1, lik2, pars.NPHASEBIN_1hr, plot_name)

def plot_two_likelihoods(lik1, lik2, NPHASEBIN, plot_name):

    phase_diff_samples = np.arange(-20, 20, 0.02)
    norm1, mean1, std1 = svd.lik_norm_mean_std(lik1, phase_diff_samples)
    norm2, mean2, std2 = svd.lik_norm_mean_std(lik2, phase_diff_samples)

    '''save [norm1, mean1, std1, norm2, mean2, std2]. Note these values are in 800 Phase bins. The first way is for the average profile, and the second way is for the individual summing up profile.'''
    two_lik_mean_std_file = 'two_lik_mean_std.npy'
    two_lik_mean_std = np.concatenate(([[norm1], [mean1], [std1], [norm2], [mean2], [std2]])).astype(np.float64)
    print 'two_lik_mean_std: ', two_lik_mean_std

    if False:
        if os.path.exists(two_lik_mean_std_file):
            sequence = np.load(two_lik_mean_std_file)
            np.save(two_lik_mean_std_file, np.vstack((sequence, two_lik_mean_std)))
        else:
            np.save(two_lik_mean_std_file, two_lik_mean_std)


    fontsize = 16
    linewidth = 3.0
    plt.close('all')
    phase_diff_range = np.linspace(np.amin(phase_diff_samples)/NPHASEBIN, np.amax(phase_diff_samples)/NPHASEBIN, num=len(phase_diff_samples), endpoint=True)
    plt.semilogy(phase_diff_range, lik1 / norm1  / (0.02/NPHASEBIN), c=colors.cnames['darkblue'], linestyle='-', linewidth=linewidth, label='Likelihood 1')
    plt.semilogy(phase_diff_range, lik2 / norm2  / (0.02/NPHASEBIN), c=colors.cnames['orange'], linestyle='--', linewidth=linewidth, label='Likelihood 2')
    plt.xlabel('Phase', fontsize=fontsize)
    plt.ylabel('log(Likelihood)', fontsize=fontsize) 
    phase_min = min(phase_diff_range[np.where((lik1 / norm1  / (0.02/NPHASEBIN))>np.amax(lik1 / norm1  / (0.02/NPHASEBIN)) * 10**-4)[0][0]], phase_diff_range[np.where((lik2 / norm2  / (0.02/NPHASEBIN))>np.amax(lik2 / norm2 / (0.02/NPHASEBIN)) * 10**-4)[0][0]])
    phase_max = max(phase_diff_range[np.where((lik1 / norm1  / (0.02/NPHASEBIN))>np.amax(lik1 / norm1  / (0.02/NPHASEBIN)) * 10**-4)[0][-1]], phase_diff_range[np.where((lik2 / norm2  / (0.02/NPHASEBIN))>np.amax(lik2 / norm2 / (0.02/NPHASEBIN)) * 10**-4)[0][-1]])
#    plt.xlim((phase_diff_range[np.where((lik1 / norm1  / (0.02/NPHASEBIN))>np.amax(lik1 / norm1  / (0.02/NPHASEBIN)) * 10**-4)[0][0]],phase_diff_range[np.where((lik2 / norm2  / (0.02/NPHASEBIN))>np.amax(lik2 / norm2  / (0.02/NPHASEBIN)) * 10**-4)[0][-1]]))
    plt.xlim((phase_min, phase_max))
    plt.ylim((np.amax(lik1 / norm1 / (0.02/NPHASEBIN)) * 10**-4, np.amax(lik1 / norm1 / (0.02/NPHASEBIN))*4.5))
    plt.legend(loc='upper left', fontsize=fontsize)
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.savefig(plot_name, bbox_inches='tight')

def compare_integration_lik():

    '''Note the two_lik_mean_std_file is saved in the sequence of [norm1, mean1, std1, norm2, mean2, std2]. Note these values are in 800 Phase bins. The first way is for the average profile, and the second way is for the individual summing up profile.'''
    two_lik_mean_std_file = np.load('two_lik_mean_std.npy')
    lik1_mean = two_lik_mean_std_file[:,1]
    lik1_std = two_lik_mean_std_file[:,2]
    lik2_mean = two_lik_mean_std_file[:,4]
    lik2_std = two_lik_mean_std_file[:,5] 

    time_mjd_1hr, dBATdra_1hr, dBATddec_1hr, phase_data_1hr, phase_data_err_1hr = fit_timing.time_pattern(this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, pars.NPHASEBIN_1hr)

    '''change the unit from phase bins to phase'''
    lik1_mean_phase = two_lik_mean_std_file[:,1] / pars.NPHASEBIN_1hr
    lik1_std_phase = two_lik_mean_std_file[:,2] / pars.NPHASEBIN_1hr
    lik2_mean_phase = two_lik_mean_std_file[:,4] / pars.NPHASEBIN_1hr
    lik2_std_phase = two_lik_mean_std_file[:,5] / pars.NPHASEBIN_1hr

    '''Least squares polynomial fit'''
    poly1 = np.polyfit(time_mjd_1hr, lik1_mean_phase, 1, w=(1/lik1_std_phase), cov=True) 
    poly2 = np.polyfit(time_mjd_1hr, lik2_mean_phase, 1, w=(1/lik2_std_phase), cov=True)
    poly_fit1 = np.poly1d(poly1[0])
    poly_fit2 = np.poly1d(poly2[0])  
    lik1_fit = poly_fit1(time_mjd_1hr)
    lik2_fit = poly_fit2(time_mjd_1hr) 
    chi_squared1 = np.sum((np.polyval(poly1[0], time_mjd_1hr) - lik1_mean_phase) ** 2)
    chi_squared2 = np.sum((np.polyval(poly2[0], time_mjd_1hr) - lik2_mean_phase) ** 2)

    print 'chi_squared1: ',chi_squared1
    print 'chi_squared2: ',chi_squared2


    '''Print the phase rate of change. Change the unit from 1/MJD to 1/sec'''
    print 'the phase rate of change of lik1 is:' + str(poly1[0][0]/86400) + '\+-' + str(poly1[1][0][0]/86400) + '(1/sec)'
    print 'the phase rate of change of lik2 is:' + str(poly2[0][0]/86400) + '\+-' + str(poly2[1][0][0]/86400) + '(1/sec)'
    

    zeros_line = np.zeros(len(lik1_mean))
    markersize = 2.0
    fontsize = 16

    plt.close('all')
    plt.plot(time_mjd_1hr, lik1_fit,'r--')
    plt.plot(time_mjd_1hr, lik1_mean_phase, 'bo', markersize=markersize)
    plt.errorbar(time_mjd_1hr, lik1_mean_phase, yerr= lik1_std_phase, fmt='none', ecolor='b')
    plt.plot(time_mjd_1hr, lik2_fit,'k')
    plt.plot(time_mjd_1hr, lik2_mean_phase, 'gs', markersize=markersize)
    plt.errorbar(time_mjd_1hr, lik2_mean_phase, yerr= lik2_std_phase, fmt='none', ecolor='g')
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.savefig('two_lik_mean_std.png',bbox_inches='tight')

   
def pulses(index):

    freq = np.fft.fftfreq(phase_npy_wz.shape[1])

    pulse_raw = phase_npy_wz[index]

    '''model with 1st to 4th V mode '''    
    paras_mode_1 = np.zeros(phase_amp_bin_wz.shape[1]/2)
    paras_mode_1[0:2] = phase_amp_bin_wz[index][0:2]
    paras_mode_1_to_4 = phase_amp_bin_wz[index][0:4]    

    pulse_mode1 = np.fft.ifft(svd.fft_phase_curve(paras_mode_1, V_1hr, freq)).real    
    pulse_mode1_to_4 = np.fft.ifft(svd.fft_phase_curve(paras_mode_1_to_4, V_1hr, freq)).real

    phase_fitting = phase_amp_bin_wz[index][0]

    return pulse_raw, pulse_mode1, pulse_mode1_to_4, phase_fitting

def two_pulses(index_1, index_2):

    pulse_raw_1, pulse_mode1_1, pulse_mode1_to_4_1, phase_fitting_1 = pulses(index_1)
    pulse_raw_2, pulse_mode1_2, pulse_mode1_to_4_2, phase_fitting_2 = pulses(index_2)    

    phase_range = np.linspace(-0.5, 0.5, num=pars.NPHASEBIN, endpoint=True)
    xmax = 0.2
    xmin = -0.2
#    xmax = np.amax(phase_range)
#    xmin = np.amin(phase_range)
    fontsize = 16
    linewidth = 3.0

    '''Plot pulse raw data'''

    plt.close('all')
    plt.plot(phase_range, np.roll(pulse_raw_1, -int(np.round(phase_fitting_1)-100)),c=colors.cnames['darkgreen'], linestyle='-', linewidth=linewidth, label='Pulse 1')
    plt.plot(phase_range, np.roll(pulse_raw_2, -int(np.round(phase_fitting_2)-100)),c=colors.cnames['darkred'], linestyle='--', linewidth=linewidth, label='Pulse 2')
    plt.xlim([xmin,xmax])
    plt.ylabel('T/Tsys', fontsize=fontsize)
    plt.xlabel('Phase', fontsize=fontsize)
    plt.legend(loc='upper right', fontsize=fontsize)
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.savefig('two_pulses_raw.png', bbox_inches='tight')
#    plt.show()

    '''Plot pulse profiles in real space'''
    plt.close('all')
    f, ((ax1, ax2)) = plt.subplots(2, 1, sharex='col', figsize=(8,9))
    f.subplots_adjust(hspace=0.07)

    ax1.plot(phase_range, np.roll(pulse_raw_1, -int(np.round(phase_fitting_1)-100)),c=colors.cnames['darkred'], linestyle='--', linewidth=1.5, label='Pulse 1')
    ax1.plot(phase_range, np.roll(pulse_raw_2, -int(np.round(phase_fitting_2)-100)),c=colors.cnames['darkgreen'], linestyle='--', linewidth=1.5, label='Pulse 2')
    ax1.plot(phase_range, np.roll(pulse_mode1_1, -int(np.round(phase_fitting_1)-100)),'r-', linewidth=1, label='Fitting pulse 1')
    ax1.plot(phase_range, np.roll(pulse_mode1_2, -int(np.round(phase_fitting_2)-100)),'g-', linewidth=1, label='Fitting pulse 2')
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylabel('T/Tsys', fontsize=fontsize)
    ax1.set_title('The 1st V mode')
    legend = ax1.legend(loc='upper right')

    ax2.plot(phase_range, np.roll(pulse_raw_1, -int(np.round(phase_fitting_1)-100)),c=colors.cnames['darkred'], linestyle='--', linewidth=1.5, label='Pulse 1')
    ax2.plot(phase_range, np.roll(pulse_raw_2, -int(np.round(phase_fitting_2)-100)),c=colors.cnames['darkgreen'], linestyle='--', linewidth=1.5, label='Pulse 2')    
    ax2.plot(phase_range, np.roll(pulse_mode1_to_4_1, -int(np.round(phase_fitting_1)-100)),'r-', label='Fitting pulse 1')
    ax2.plot(phase_range, np.roll(pulse_mode1_to_4_2, -int(np.round(phase_fitting_2)-100)),'g-', label='Fitting pulse 2')
    ax2.set_xlim([xmin,xmax])
    ax2.set_xlabel('Phase ', fontsize=fontsize)
    ax2.set_ylabel('T/Tsys', fontsize=fontsize)
    ax2.set_title('The 1st to 4th V mode')
    legend = ax2.legend(loc='upper right')

    plt.savefig('two_pulses.png')
  
    '''Plot pulse profiles separately  in real space'''
    plt.close('all')
    f, ((ax1, ax2)) = plt.subplots(2, 1, sharex='col', figsize=(8,9))
    f.subplots_adjust(hspace=0.07)

    ax1.plot(phase_range, np.roll(pulse_raw_1, -int(np.round(phase_fitting_1)-100)),c=colors.cnames['darkred'], linestyle='--', linewidth=1.5, label='Data')
    ax1.plot(phase_range, np.roll(pulse_mode1_1, -int(np.round(phase_fitting_1)-100)),c=colors.cnames['darkgreen'], linestyle='-', linewidth=1, label='1st V mode')
    ax1.plot(phase_range, np.roll(pulse_mode1_to_4_1, -int(np.round(phase_fitting_1)-100)),'b-', linewidth=1, label='1st to 4th V mode')
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylabel('T/Tsys', fontsize=fontsize)
    ax1.set_title('Pulse 1')
    legend = ax1.legend(loc='upper right')


    ax2.plot(phase_range, np.roll(pulse_raw_2, -int(np.round(phase_fitting_2)-100)),c=colors.cnames['darkred'], linestyle='--', linewidth=1.5, label='Data')
    ax2.plot(phase_range, np.roll(pulse_mode1_2, -int(np.round(phase_fitting_2)-100)),c=colors.cnames['darkgreen'], linestyle='-', linewidth=1, label='1st V mode')
    ax2.plot(phase_range, np.roll(pulse_mode1_to_4_2, -int(np.round(phase_fitting_2)-100)),'b-', linewidth=1, label='1st to 4th V mode')
    ax2.set_xlim([xmin,xmax])
    ax2.set_xlabel('Phase ', fontsize=fontsize)
    ax2.set_ylabel('T/Tsys', fontsize=fontsize)
    ax2.set_title('Pulse 2')
    legend = ax2.legend(loc='upper right')

    plt.savefig('two_pulses_separate.png')

 
if __name__ == '__main__':
    main()

