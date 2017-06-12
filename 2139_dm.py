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

NPHASEBIN_1hr = pars.NPHASEBIN_1hr
this_file = h5py.File('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/J2139+00_1hr_foldingh5','r')

series = ['DATA_FOLDING_0_255']
#series = ['DATA_FOLDING_0_255', 'DATA_FOLDING_256_511', 'DATA_FOLDING_512_767', 'DATA_FOLDING_768_1023', 'DATA_FOLDING_1024_1279', 'DATA_FOLDING_1280_1535', 'DATA_FOLDING_1536_1791', 'DATA_FOLDING_1792_2047', 'DATA_FOLDING_2048_2303', 'DATA_FOLDING_2304_2559', 'DATA_FOLDING_2561_2815', 'DATA_FOLDING_2816_3071', 'DATA_FOLDING_3072_3327', 'DATA_FOLDING_3328_3583', 'DATA_FOLDING_3584_3839', 'DATA_FOLDING_3840_4095', 'DATA_FOLDING_4096_4351', 'DATA_FOLDING_4352_4607', 'DATA_FOLDING_4608_4863', 'DATA_FOLDING_4864_5119', 'DATA_FOLDING_5120_5375', 'DATA_FOLDING_5376_5631', u'DATA_FOLDING_5632_5887', 'DATA_FOLDING_5888_6143', 'DATA_FOLDING_6144_6399', 'DATA_FOLDING_6400_6655', 'DATA_FOLDING_6656_6911', 'DATA_FOLDING_6912_7167', 'DATA_FOLDING_7168_7423', 'DATA_FOLDING_7424_7679', 'DATA_FOLDING_7680_7935', 'DATA_FOLDING_7936_8191', 'DATA_FOLDING_8192_8447', 'DATA_FOLDING_8448_8703', 'DATA_FOLDING_8704_8959', 'DATA_FOLDING_8960_9215', 'DATA_FOLDING_9216_9471', 'DATA_FOLDING_9472_9727', 'DATA_FOLDING_9728_9983', 'DATA_FOLDING_9984_10239']

model = np.load('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_1hr/pointed_phase_model.npy')

plt.close('all')
plt.imshow(psr_svd.rebin_spec(this_file[series[0]][:,0,:,0].T, 32, 1))
plt.savefig('profile_0.png')

def shift_phase_bin(phase_bin_shift, profile):
    # the input profile is required in the shape of (frequencies, phase_bins) 
    profile_shift = np.zeros(profile.shape)
    profile_fft = np.zeros((profile.shape), dtype=complex)
    for ii in xrange(profile.shape[0]):   
        profile_fft[ii] = fftpack.fft(profile[ii])
        n = profile_fft.shape[-1]
        freq = fftpack.fftfreq(n, 1./n)
        phase = np.exp(-2j * np.pi * phase_bin_shift / NPHASEBIN_1hr * freq)
        profile_shift[ii] = fftpack.ifft(profile_fft[ii] * phase).real
    return profile_shift 


shift_profiles = np.zeros((len(series), this_file['DATA'].shape[3], NPHASEBIN_1hr))
for ii in xrange(len(series)):
    shift_profiles[ii] = shift_phase_bin(model[ii], this_file[series[ii]][:,0,:,0].T)
#    shift_profiles[ii] = shift_phase_bin(0, (this_file[series[ii]][:,0,:,0].T))

plt.close('all')
plt.imshow(psr_svd.rebin_spec(shift_profiles[0], 32, 1))
plt.savefig('profile_0_shift.png')


#stack profiles
if True:
    profile_stack = 1
    nprof = len(shift_profiles)
    nprof -= nprof % profile_stack
    shift_profiles = shift_profiles[:nprof].reshape(nprof // profile_stack, profile_stack, this_file['DATA'].shape[3], NPHASEBIN_1hr)
    shift_profiles = np.mean(shift_profiles, 1)
    shift_profiles = np.mean(shift_profiles, 0)
    plt.imshow(psr_svd.rebin_spec(shift_profiles, 32, 1))
    plt.savefig('profiles_all_shifted.png')
