import numpy as np
import h5py
import fit_timing
import pars

def main():
    try:
        plot_bary_diff()
    except None:
        print IOError


def plot_bary_diff():

    '''Get timing profiles'''

    this_file_wz = pars.this_file_wz
    bin_number_wz = pars.bin_number_wz
    phase_amp_bin_wz = pars.phase_amp_bin_wz
    random_res_wz = pars.random_res_wz

    this_file_1hr = pars.this_file_1hr
    bin_number_1hr = pars.bin_number_1hr
    phase_amp_bin_1hr = pars.phase_amp_bin_1hr
    random_res_1hr = pars.random_res_1hr
 
    time_mjd_wz, dBATdra_wz, dBATddec_wz, phase_data_wz, phase_data_err_wz = fit_timing.time_pattern(this_file_wz, bin_number_wz, phase_amp_bin_wz)
    time_mjd_1hr, dBATdra_1hr, dBATddec_1hr, phase_data_1hr, phase_data_err_1hr = fit_timing.time_pattern(this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, pars.NPHASEBIN_1hr)

    time_mjd = np.concatenate((time_mjd_wz, time_mjd_1hr))
    dBATdra = np.concatenate((dBATdra_wz, dBATdra_1hr))
    dBATddec = np.concatenate((dBATddec_wz, dBATddec_1hr))
    phase_data = np.concatenate((phase_data_wz, phase_data_1hr))
    phase_data_err = np.concatenate((phase_data_err_wz, phase_data_err_1hr))
    random_res = np.concatenate((random_res_wz, random_res_1hr))

    '''Find timing solution'''
    fit_timing.fitting(pars.fit_pars, time_mjd, dBATdra, dBATddec, phase_data, phase_data_err, random_res)





if __name__ == '__main__':
    main()
