import sys
import os
import os.path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import leastsq

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            plot_bary_diff(filename)
        except None:
            print IOError

def plot_bary_diff(filename):
    this_file = h5py.File(filename, "r")
 
#    index = [[4,8], [15, 19],[28, 32],[44, 49],[67, 83],[98,113],[114,117],[118,128],[146,152]]
#    max_phase = [69, 69, 63+100, 65+100, 90+200, 97+200, 91+300, 95+300, 5+400]

    '''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''

    '''wigglez data'''
    bin_number = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/bin_number_2139_wz.txt')
    '''one hour pointing data'''
#    bin_number = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/bin_number_2139_57178_20.txt')

    bary_diff = np.zeros(len(bin_number))

    for ii in range(len(bin_number)):
        bary_diff[ii] = ((this_file['BARY_TIME'][bin_number[ii][0]]+this_file['BARY_TIME'][bin_number[ii][1]])/2. -this_file['BARY_TIME'][0])*24

#   title = 'delta_t: '+str(delta_t)+' sec.'

    '''Try to fit'''

    def qua_func(x, a, b, c):
        return a*x**2 + b*x + c

#    '''Find quadratic parameters by minimizing chi-squared'''
    n_phase_bin = 100
#    data_i = bin_number[191:237,2]
#    time_i = bary_diff[191:237]
#    index_j = np.arange(0, len(time_i),1)

    sl = np.logical_and(bary_diff > 00, bary_diff < 800)
    data_i = bin_number[sl,2]
    time_i = bary_diff[sl]
    t0 = time_i[0]
    time_i -= 250


#    chi_squ = lambda x:  np.sum(((data_i - (x[0]*bary_diff[:]**2 + x[1]*bary_diff[:] + x[2]) + n_phase_bin/2) % n_phase_bin - n_phase_bin/2)**2)   
#    res = minimize(chi_squ, ([100, popt[1], popt[2]]), tol=1e3)
#    print "best [a, b, c]: "+str(res.x)

    '''Find quadratic parameters by leastsq'''
    funcQuad=lambda tpl,time_i,data_i : (((data_i - (0*time_i**3 + tpl[0]*time_i**2 + tpl[1]*time_i + tpl[2]) + n_phase_bin/2) % n_phase_bin - n_phase_bin/2))
    func=funcQuad
    ErrorFunc=lambda tpl,time_i,data_i : func(tpl,time_i,data_i)
    tplInitial=(0.0226,10.52, 12.8)
    tplFinal,success=leastsq(funcQuad,tplInitial[:],args=(time_i,data_i))
    print "quadratic fit: " ,tplFinal
    print "sucess?:", success
    print np.sum(funcQuad(tplFinal, time_i, data_i)**2)

    '''Find quadratic parameters by complex method'''

    data_im = np.exp(1j * data_i / n_phase_bin * 2 * np.pi)
    model_im = np.exp(1j * (tplFinal[0]*time_i**2 + tplFinal[1]*time_i + tplFinal[2]) / n_phase_bin * 2 * np.pi)
    funcQuad=lambda tpl_im,time_i,data_i : (np.sum(np.absolute(np.exp(1j * data_i / n_phase_bin * 2 * np.pi) -np.exp(1j * (tpl_im[0]*time_i**2 + tpl_im[1]*time_i + tpl_im[2]) / n_phase_bin * 2 * np.pi)))**2)
    func=funcQuad
    ErrorFunc=lambda tpl_im,time_i,data_i : func(tpl_im,time_i,data_i)
    tpl_imInitial=(0.0226,10.52, 12.8)
    tpl_imFinal,success=leastsq(funcQuad,tpl_imInitial[:],args=(time_i,data_i))
    print "quadratic fit: " ,tpl_imFinal
    print "sucess?:", success
    print np.sum(funcQuad(tpl_imFinal, time_i, data_i)**2)



    

    x_axes = np.linspace(time_i[0], time_i[-1],5000)
    y = (tplFinal[0]*x_axes**2 + tplFinal[1]*x_axes + tplFinal[2]) % n_phase_bin
    y_im = (tpl_imFinal[0]*x_axes**2 + tpl_imFinal[1]*x_axes + tpl_imFinal[2]) % n_phase_bin

#   title = 'P(t)= '+str(np.round(tplFinal[0],3))+'*t**2 +'+str(np.round(tplFinal[1],3))+'*t +'+ str(np.round(tplFinal[2],3))
    
    plt.subplot(4,1,1)
    plt.plot(time_i, data_i, 'bo')
    plt.plot(x_axes, y, 'r--')
    plt.xlabel('Bary diff (hours)', fontsize=14)
    plt.ylabel('Max Phase Bins Number', fontsize=14)

    plt.subplot(4,1,2)
    plt.plot(time_i, funcQuad(tplFinal, time_i, data_i), 'bo')
    plt.xlabel('Bary diff (hours)', fontsize=14)
    plt.ylabel('Phase bin residuals', fontsize=14)

    plt.subplot(4,1,3)
    plt.plot(time_i, data_im, 'bo')
    plt.plot(x_axes, y_im, 'r--')
    plt.xlabel('Bary diff (hours)', fontsize=14)
    plt.ylabel('Max Phase Bins Number', fontsize=14)

    plt.subplot(4,1,4)
    plt.plot(time_i, funcQuad(tpl_imFinal, time_i, data_im), 'bo')
    plt.xlabel('Bary diff (hours)', fontsize=14)
    plt.ylabel('Phase bin residuals', fontsize=14)

    plt.savefig('phase_fit_test.png')



if __name__ == '__main__':
    main()

