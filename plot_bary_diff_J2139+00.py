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

    bin_number = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/bin_number_2139_wz_pdot.txt')

    bary_diff = np.zeros(len(bin_number))

    for ii in range(len(bin_number)):
        bary_diff[ii] = ((this_file['BARY_TIME'][bin_number[ii][0]]+this_file['BARY_TIME'][bin_number[ii][1]])/2. -this_file['BARY_TIME'][0])*24

#   title = 'delta_t: '+str(delta_t)+' sec.'

    '''Try to fit'''

    def qua_func(x, a, b, c):
        return a*x**2 + b*x + c

    #popt, pcov = curve_fit(qua_func, bary_diff, bin_number[:,2])
    #print 'popt: ' + str(popt)
    #print 'pcov: ' + str(pcov)

#    '''Find quadratic parameters by minimizing chi-squared'''
    n_phase_bin = 100
    data_i = bin_number[:,2]
    time_i = bary_diff[:]
#    chi_squ = lambda x:  np.sum(((data_i - (x[0]*bary_diff[:]**2 + x[1]*bary_diff[:] + x[2]) + n_phase_bin/2) % n_phase_bin - n_phase_bin/2)**2)   
#    res = minimize(chi_squ, ([100, popt[1], popt[2]]), tol=1e3)
#    print "best [a, b, c]: "+str(res.x)


    '''Find quadratic parameters by leastsq'''
    funcQuad=lambda tpl,time_i,data_i : (((data_i - (0*time_i**3 + tpl[0]*time_i**2 + tpl[1]*time_i + tpl[2]) + n_phase_bin/2) % n_phase_bin - n_phase_bin/2))
    func=funcQuad
    ErrorFunc=lambda tpl,time_i,data_i : func(tpl,time_i,data_i)
    tplInitial=(0.03,-0.05,70)
    tplFinal,success=leastsq(funcQuad,tplInitial[:],args=(time_i,data_i))
    print "quadratic fit: " ,tplFinal
    print "sucess?:", success
    print np.sum(funcQuad(tplFinal, time_i, data_i)**2)

    

    x_axes = np.linspace(0, bary_diff[-1],300)
    y = (tplFinal[0]*x_axes**2 + tplFinal[1]*x_axes + tplFinal[2]) % n_phase_bin
    
    plt.subplot(2,1,1)
    plt.plot(time_i, data_i, 'bo')
    plt.plot(x_axes, y, 'r--')
    plt.xlabel('Bary diff (hours)', fontsize=20)
    plt.ylabel('Max Phase Bins Number', fontsize=20)
#    plt.title(title, fontsize=20)

    plt.subplot(2,1,2)
    plt.plot(time_i, funcQuad(tplFinal, time_i, data_i), 'bo')
    plt.xlabel('Bary diff (hours)', fontsize=20)
    plt.ylabel('Phase bin residuals', fontsize=20)
    #plt.show()
    plt.savefig('phase_fit.png')

if __name__ == '__main__':
    main()

