import numpy as np
import svd
import pars

def main():
    try:
        ploting()
    except (IOError, ValueError):
        print IOError

def ploting():

    this_file_wz = pars.this_file_wz
    bin_number_wz = pars.bin_number_wz
    phase_amp_bin_wz = pars.phase_amp_bin_wz
    phase_npy_wz = pars.phase_npy_wz

    this_file_1hr = pars.this_file_1hr
    bin_number_1hr = pars.bin_number_1hr
    phase_amp_bin_1hr = pars.phase_amp_bin_1hr
    phase_npy_1hr = pars.phase_npy_1hr 

#    U_wz, s_wz, , V_wz, phase_model_wz = svd.svd(this_file_wz, bin_number_wz, phase_npy_wz, NPHASEBIN_wz)
#    U_1hr_origin, s_1hr_origin, V_1hr_origin, phase_model_1hr = svd.svd(this_file_1hr, bin_number_1hr, phase_npy_1hr, NPHASEBIN_1hr)

    svd.plot_svd(this_file_wz, bin_number_wz, phase_amp_bin_wz, phase_npy_wz, 'phase_wz',pars.NPHASEBIN_wz, RESCALE=None)
    svd.plot_svd(this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, phase_npy_1hr, 'phase_1hr', pars.NPHASEBIN_1hr, RESCALE=None)

#    U_1hr, s_1hr, V_1hr, phase_model_1hr = svd.svd(this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, phase_npy_1hr, pars.NPHASEBIN_1hr, RESCALE=None)

#    svd.phase_fit(1, phase_npy_1hr, V_1hr, 'fitting_phase_fft_57178_', pars.NPHASEBIN_1hr)

#    svd.plot_two_temps(this_file_wz, bin_number_wz, phase_amp_bin_wz, phase_npy_wz, this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, phase_npy_1hr, RESCALE=None)


#    for ii in xrange(len(s_1hr)):
#        print 'ii= '+str(ii)
#        svd.phase_fit(ii, phase_npy_1hr, V_1hr, 'fitting_phase_fft_57178_', pars.NPHASEBIN_1hr)  
   
   


  
if __name__ == '__main__':
    main()

