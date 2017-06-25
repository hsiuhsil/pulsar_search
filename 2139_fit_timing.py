import numpy as np
import h5py
import subprocess
import fit_timing
import pars
import bary_time

TIME0 = pars.TIME0
T = pars.T

def main():
    '''Get timing profiles'''

    this_file_wz = pars.this_file_wz
    bin_number_wz = pars.bin_number_wz
    phase_amp_bin_wz = pars.phase_amp_bin_wz
    random_res_wz = pars.random_res_wz
    bin_number_wz_hist = np.delete(bin_number_wz, [89,90,91,92], axis=0)
#    fit_timing.plot_wz_hist(bin_number_wz_hist)

    this_file_1hr = pars.this_file_1hr
    bin_number_1hr = pars.bin_number_1hr
    phase_amp_bin_1hr = pars.phase_amp_bin_1hr
    random_res_1hr = pars.random_res_1hr
 
    time_mjd_wz, dBATdra_wz, dBATddec_wz, phase_data_wz, phase_data_err_wz = fit_timing.time_pattern(this_file_wz, bin_number_wz, phase_amp_bin_wz)
    time_mjd_1hr, dBATdra_1hr, dBATddec_1hr, phase_data_1hr, phase_data_err_1hr = fit_timing.time_pattern(this_file_1hr, bin_number_1hr, phase_amp_bin_1hr, pars.NPHASEBIN_1hr)

    model_phase_bin_wz = fit_timing.timing_model_1(pars.fit_pars, time_mjd_wz, dBATdra_wz, dBATddec_wz, pars.NPHASEBIN_wz, RESCALE=None)
    model_phase_bin_1hr = fit_timing.timing_model_1(pars.fit_pars, time_mjd_1hr, dBATdra_1hr, dBATddec_1hr, pars.NPHASEBIN_1hr, RESCALE=None)

    model_phase_wz = model_phase_bin_wz / pars.NPHASEBIN_wz
    model_phase_1hr = model_phase_bin_1hr / pars.NPHASEBIN_1hr

    b = np.loadtxt('./data_file/bin_number_2139_57178_5sec.txt')
    time_mjd_1hr_stack = np.zeros(40)
    dBATdra_1hr_stack = np.zeros(40)
    dBATddec_1hr_stack = np.zeros(40)

#    for ii in xrange(len(time_mjd_1hr_stack)):
#        a = []
#        r = []
#        d = []
#        for jj in xrange(16):
#            a.append(np.mean((this_file_1hr['BARY_TIME'][b[jj+16*ii][0]],this_file_1hr['BARY_TIME'][b[jj+16*ii][1]])))
#            r.append(np.mean((this_file_1hr['dBATdra'][b[jj+16*ii][0]],this_file_1hr['dBATdra'][b[jj+16*ii][1]])))
#            d.append(np.mean((this_file_1hr['dBATddec'][b[jj+16*ii][0]],this_file_1hr['dBATddec'][b[jj+16*ii][1]])))
#        time_mjd_1hr_stack[ii] = np.mean(a)
#        dBATdra_1hr_stack[ii] = np.mean(r)
#        dBATddec_1hr_stack[ii] = np.mean(d)

    # note that wz[2] is a outlier case to be excluded
    # note that wz[89:93] is the same as [82:86]

##   time_mjd = np.concatenate((time_mjd_wz[0:2], time_mjd_wz[3:89], time_mjd_wz[93:], time_mjd_1hr))
##    dBATdra = np.concatenate((dBATdra_wz[0:2], dBATdra_wz[3:89], dBATdra_wz[93:], dBATdra_1hr))
##    dBATddec = np.concatenate((dBATddec_wz[0:2], dBATddec_wz[3:89], dBATddec_wz[93:], dBATddec_1hr))
##    phase_data = np.concatenate((phase_data_wz[0:2], phase_data_wz[3:89], phase_data_wz[93:], phase_data_1hr))
##    phase_data_err = np.concatenate((phase_data_err_wz[0:2], phase_data_err_wz[3:89], phase_data_err_wz[93:], phase_data_err_1hr))
#    random_res = np.concatenate((random_res_wz, random_res_1hr[:32]))
##    random_res = np.zeros((271,100,5))


    time_mjd = np.concatenate((time_mjd_wz, time_mjd_1hr))
    dBATdra = np.concatenate((dBATdra_wz, dBATdra_1hr))
    dBATddec = np.concatenate((dBATddec_wz, dBATddec_1hr))
    phase_data = np.concatenate((phase_data_wz, phase_data_1hr))
    phase_data_err = np.concatenate((phase_data_err_wz, phase_data_err_1hr))
    model_phase = np.concatenate((model_phase_wz, model_phase_1hr))


    rev_index = [2, 89, 90, 91, 92] # 2 is a outlier, 89-92 are the same as 82-85.

    '''data of 2011 wz'''
    time_mjd_11w = np.delete(time_mjd[0:190], rev_index)
    print 'mean of mjd_11w', np.mean(time_mjd_11w)
    dBATdra_11w = np.delete(dBATdra[0:190], rev_index)
    dBATddec_11w = np.delete(dBATddec[0:190], rev_index)
    phase_data_11w = np.delete(phase_data[0:190], rev_index)
    phase_data_err_11w = np.delete(phase_data_err[0:190], rev_index)
    random_res_11w = np.zeros((len(time_mjd),100,5))[0:185]

    '''data of 2015 wz'''
    time_mjd_15wp = time_mjd[190:276]
    print 'mean of mjd_15wp', np.mean(time_mjd_15wp)
    dBATdra_15wp = dBATdra[190:276]
    dBATddec_15wp = dBATddec[190:276]
    phase_data_15wp = phase_data[190:276]
    phase_data_err_15wp = phase_data_err[190:276]
    random_res_15wp = np.zeros((len(time_mjd),100,5))[190:276]

    '''data of 11w, 15w, and pointed observation'''
    time_mjd = np.delete(time_mjd, rev_index)
    dBATdra = np.delete(dBATdra, rev_index)
    dBATddec = np.delete(dBATddec, rev_index)
    phase_data = np.delete(phase_data, rev_index)
    phase_data_err = np.delete(phase_data_err, rev_index)
    random_res = np.zeros((len(time_mjd),100,5))
    model_phase = np.delete(model_phase, rev_index)




#    np.save('mjd.npy', time_mjd)
    print 'time_mjd.shape',time_mjd.shape
    print 'dBATdra.shape',dBATdra.shape
    print 'dBATddec.shape',dBATddec.shape
    print 'phase_data.shape',phase_data.shape
    print 'phase_data_err.shape',phase_data_err.shape
    print 'random_res.shape', random_res.shape
    print 'model_phase.shape', model_phase.shape 
    print 'model_phase', model_phase

    '''generates TOAs and errors'''
#    TOAs = time_mjd + pars.T * phase_data / 86400 
#    TOAs_err = pars.T * phase_data_err / 86400
#    print 'TOAs: ', TOAs
#    print 'TOAs_err ', TOAs_err

    '''Find timing solution'''
    if False:
        '''For 11w data'''
        print 'Timing solution of 11w, with fixed period derivative, period, and offset'
        pars_11w = [1e-04, 1e-04]
        fit_timing.fitting(pars_11w, time_mjd_11w, dBATdra_11w, dBATddec_11w, phase_data_11w, phase_data_err_11w, random_res_11w)
        '''For 15wp  ata'''
        print 'Timing solution of 15w and 15 pointed, with fixed period derivative, period, and offset'
        pars_15wp = [1e-04, 1e-04]
        fit_timing.fitting(pars_15wp, time_mjd_15wp, dBATdra_15wp, dBATddec_15wp, phase_data_15wp, phase_data_err_15wp, random_res_15wp)

    if False:
        '''For all data'''
        fit_timing.fitting(pars.fit_pars, time_mjd, dBATdra, dBATddec, phase_data, phase_data_err, random_res)
        print 'pars.fit_pars', pars.fit_pars

    '''Check timing model'''
#    model = fit_timing.timing_model_1(pars.fit_pars, time_mjd, dBATdra, dBATddec)
#    np.save('pointed_phase_model.npy',model[231:]*4)
#    print len(model)
    
    '''Check timing solution with TEMPO'''
    #TOA = reference_time + [(time - reference_time) // T + phase] * T
    TOAs_topo = np.zeros(len(time_mjd))
    TOAs_bary = np.zeros(len(time_mjd))
    RA = bary_time.deg_to_HMS(pars.RA)
    DEC = bary_time.deg_to_DMS(pars.DEC)

    for ii in xrange(len(TOAs_bary)):
        print 'ii:',ii
        TOAs_bary[ii] =  TIME0 + ((time_mjd[ii] - TIME0) // pars.T + model_phase[ii]) * pars.T
        bary_in = repr(TOAs_bary[ii])
        topo_guess = bary_in
        error = 1.
        while np.abs(error) *3600 * 24 > 1e-4:
            p = subprocess.Popen(["bary", "GBT", RA, DEC], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            bary_guess = np.float64(p.communicate(input=topo_guess)[0].split()[1])
            error = bary_guess - np.float64(bary_in)
            topo_guess = repr(np.float64(topo_guess) - error)
        TOAs_topo[ii] = topo_guess
   

    diff = np.zeros(len(time_mjd))
    for ii in xrange(len(diff)):
        topo_time = repr(TOAs_topo[ii])
        p = subprocess.Popen(["bary", "GBT", RA, DEC], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bary_TOAs_topo = np.float64(p.communicate(input=topo_time)[0].split()[1])
        diff[ii] = (bary_TOAs_topo - TOAs_bary[ii])*86400*1000 #change unit from day to ms

    print 'diff', diff
    print 'diff_min', np.amin(diff)
    print 'diff_max', np.amax(diff)
if __name__ == '__main__':
    main()

