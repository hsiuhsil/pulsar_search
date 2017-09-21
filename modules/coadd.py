import os
import sys
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm

import subprocess
import fit_timing
import pars
import bary_time
import convert
import svd
import plot_pulse

TIME0 = pars.TIME0
T = pars.T
NPHASEBIN = pars.NPHASEBIN_wz

def main():
    '''Get timing profiles'''

    try:
        this_file = h5py.File(sys.argv[1],'r')
        bin_number = np.loadtxt(sys.argv[2])
    except(RuntimeError, TypeError, NameError, ValueError, IndexError):
        pass

    if pars.PSR == 'J2139+0040_ANTF':
        fit_pars_init = pars.fit_pars     
    else:
        fit_pars_init = [0.]*5
        fit_pars_init[0] = convert.pdot_to_fit_pars0(pars.T1)

    print 'fit_pars_init', fit_pars_init 

    # generate predicted phase bins
    time_mjd, dBATdra, dBATddec = fit_timing.time_pattern(this_file, bin_number)
    print 'time_mjd: ', time_mjd
    print '[pars.T1 , pars.T, pars.PEPOCH, pars.NPHASEBIN]: ',[pars.T1 , pars.T, pars.PEPOCH, pars.NPHASEBIN]
    predicted_phase_bins = (pars.T1 / (2*pars.T**2) * ((time_mjd - pars.PEPOCH)*86400) **2) 
    predicted_phase_bins = (predicted_phase_bins * pars.NPHASEBIN) % pars.NPHASEBIN
    print 'predicted_phase_bins: ', predicted_phase_bins

    model_phase_bins = fit_timing.timing_model_1(fit_pars_init, time_mjd, dBATdra, dBATddec, NPHASEBIN, RESCALE=None) # output of timing_model is phase bins 
    print 'predicted model_phase_bin', model_phase_bins

    co_add_folding(this_file, bin_number, predicted_phase_bins, NPHASEBIN)

def co_add_folding(this_file, bin_number, model_phase_bins, NPHASEBIN):
    '''phase transform in FFT, back to ifft, and co-add folding profiles'''
    coadd_folding = np.zeros((this_file['DAT_FREQ'].shape[1], NPHASEBIN))

#    for ii in [0,1,2,3,4,5]:
    for ii in xrange(len(model_phase_bins)):
        print 'ii: ',ii
        initial, final = np.int(bin_number[ii][0]), np.int(bin_number[ii][1])
        fold_origin = plot_pulse.dedisperse_spec(this_file['DATA_FOLDING'+'_'+str(initial)+'_'+str(final)][:, 0, :, 0], this_file).T #fold_origin.shape = (freq, phase bins)
        profile_fft = svd.fft(fold_origin)
        phase_bin_shifts = -model_phase_bins[ii] # minus sign for inverse the phase
        fold_trans = svd.apply_phase_shift(profile_fft, phase_bin_shifts, NPHASEBIN)
        coadd_folding += svd.ifft(fold_trans).real

    plot_folding(this_file, coadd_folding)    

def plot_folding(this_file, data):

    # data in the shape of (freq, phase bins)
    # find the average over frequencies
    data_rebin = (plot_pulse.rebin_spec(data.T)).T
    data_ave = data.mean(axis=0)

    fontsize = 16

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.set_ylabel('Freq(MHz)', fontsize=fontsize)
    ax1.tick_params(axis='both', which='major', labelsize=fontsize)
    cax1 = ax1.imshow(data_rebin, extent=[0, NPHASEBIN-1, this_file['DAT_FREQ'][0,-1], this_file['DAT_FREQ'][0,0]],aspect='auto', cmap=cm.Greys)
    cbar = plt.colorbar(cax1)
    ax2 = fig.add_subplot(212)
    ax2.set_ylabel('Mean Amp', fontsize=12)
    ax2.set_xlabel('Phase Bins', fontsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    cax2 = plt.plot(data_ave)
    plot_name = 'coadd_'+str(sys.argv[1])+'.png'
#    plot_name = 'coadd_J2139_1hr.png'
    plt.savefig(plot_name, dpi = 300, bbox_inches='tight')


if __name__ == '__main__':
    main()

