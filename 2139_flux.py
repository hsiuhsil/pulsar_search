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

T = pars.T
Tsys = pars.Tsys
G = pars.G
NPHASEBIN_1hr = pars.NPHASEBIN_1hr

# pars.this_file_1hr_folding['DAT_FREQ'][0][1229] gives the freq at 799.98047 MHz
# profile = pars.this_file_1hr_folding['DATA_FOLDING_dedispersed_0_255'][1229]
profile = pars.phase_npy_1hr[0]

# Flux density = Tpeak * W / (G*P) 
amp_max = np.amax(profile)
Tpeak = amp_max * Tsys

# calculate width of pulse
w_i = np.where (profile > np.amax(profile)/2)[0][0]
w_f = np.where (profile > np.amax(profile)/2)[0][-1]
#print 'w_i', w_i
#print 'w_f', w_f
W = np.float(w_f - w_i) / NPHASEBIN_1hr * T # in the unit of sec

print "Tpeak", Tpeak
print "W", W
print 'G',G
print 'T',T

flux = Tpeak * W / (G*T)
print 'flux:', flux

