import os
import sys
import numpy as np

filename = sys.argv[1]
index_file = sys.argv[2]
p = np.loadtxt(index_file)
 
list = []
for ii in xrange(len(p)):
    initial, final = int(p[ii][0]), int(p[ii][1])
    os.system('python /scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/modules/folding.py ' + str(initial) + ' ' + str(final) + ' ' + filename)
    os.system('python /scratch2/p/pen/hsiuhsil/gbt_data/pulsar_folding/pulsar_search/modules/plot_pulse.py ' + str(initial) + ' ' + str(final) + ' ' + filename + ' >> bin_number_' + sys.argv[1]+'.txt')
