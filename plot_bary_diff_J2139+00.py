import sys
import os
import os.path

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit
from scipy.optimize import minimize

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            plot_bary_diff(filename)
        except (IOError, ValueError):
            print IOError

def plot_bary_diff(filename):
    this_file = h5py.File(filename, "r")
 
#    index = [[4,8], [15, 19],[28, 32],[44, 49],[67, 83],[98,113],[114,117],[118,128],[146,152]]
#    max_phase = [69, 69, 63+100, 65+100, 90+200, 97+200, 91+300, 95+300, 5+400]

    '''bin_number[0,1,2] = [initial, final, maximal phase bin number]'''

    bin_number = np.loadtxt('./bin_number_2139_wz_pdot.txt')

    bary_diff = np.zeros(len(bin_number))

    for ii in range(len(bin_number)):
        bary_diff[ii] = ((this_file['BARY_TIME'][bin_number[ii][0]]+this_file['BARY_TIME'][bin_number[ii][1]])/2. -this_file['BARY_TIME'][0])*86400

#   title = 'delta_t: '+str(delta_t)+' sec.'

    '''Try to fit'''

    def func(x, a, b, c):
        return a*x**2 + b*x + c

    popt, pcov = curve_fit(func, bary_diff, bin_number[:,2])
    print 'popt: ' + str(popt)
    print 'pcov: ' + str(pcov)

    n_phase_bin = 100
    data_i = bin_number[:,2]
    chi_squ = lambda x:  np.sum(((data_i - (x[0]*bary_diff[:]**2 + x[1]*bary_diff[:] + x[2]) + n_phase_bin/2) % n_phase_bin - n_phase_bin/2)**2)   
    res = minimize(chi_squ, popt, tol=1e0)
    print "best [a, b, c]: "+str(res.x)

    x_axes = np.linspace(0, bary_diff[-1],50)
    y = popt[0]*x_axes**2 + popt[1]*x_axes + popt[2]

    
#   plt.plot(bary_diff, bin_number[:,2].tolist(), 'bo')
    plt.plot(bary_diff[:], bin_number[:,2].tolist(), 'bo')
#    plt.axis([diff[0],diff[-1], bary[0], bary[-1]], labelsize=20)
#   plt.plot(x_axes, y, 'r--')
    plt.xlabel('Bary diff (hours)', fontsize=20)
    plt.ylabel('Max Phase Bins Number', fontsize=20)
#    plt.title(title, fontsize=20)
    plt.show()

if __name__ == '__main__':
    main()

