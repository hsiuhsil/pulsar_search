import sys
import os
import os.path

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            plot_bary_topo(filename)
        except (IOError, ValueError):
            print IOError

def plot_bary_topo(filename):
    this_file = h5py.File(filename, "r")
    p = np.loadtxt('/scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/J2139+00_all/data_file/bin_number_2139_57178_5sec.txt')

    bary = np.zeros(len(p)) 
    topo = np.zeros(len(p))
    diff = np.zeros(len(p))

    for ii in xrange(len(p)):
        initial, final = int(p[ii][0]), int(p[ii][1])
        bary[ii] = (this_file['BARY_TIME'][initial]+this_file['BARY_TIME'][final])/2
        topo[ii] = (this_file['TOPO_TIME'][initial]+this_file['TOPO_TIME'][final])/2
        diff[ii] = (this_file['BARY_TIME'][initial]+this_file['BARY_TIME'][final])/2 - (this_file['TOPO_TIME'][initial]+this_file['TOPO_TIME'][final])/2

    '''change unit of diff from day to second'''
    diff *= 86400

    fontsize = 16

    plt.close('all')
    plt.plot(topo, bary)
    plt.axis([topo[0],topo[-1], bary[0], bary[-1]], labelsize=fontsize)
    plt.xlabel('TOPO', fontsize=fontsize)
    plt.ylabel('BARY', fontsize=fontsize)
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.savefig('57178_bary_topo.png', bbox_inches='tight')

    plt.close('all')
    plt.plot(diff, bary)
    plt.axis([diff[0],diff[-1], bary[0], bary[-1]], labelsize=fontsize)
    plt.xlabel('BARY - TOPO (sec)', fontsize=fontsize)
    plt.ylabel('BARY', fontsize=fontsize)
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.savefig('57178_bary_diff.png', bbox_inches='tight') 

if __name__ == '__main__':
    main()

