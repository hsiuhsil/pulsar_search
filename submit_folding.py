import os
import sys
import numpy as np

#p = np.loadtxt('./test.txt', delimiter=',')
p = np.loadtxt('./index_2139_wz.txt')
filename = sys.argv[1]
list = []
for ii in xrange(len(p)):
    initial, final = int(p[ii][0]), int(p[ii][1])
    os.system('cp ' + filename + ' temph5')
    os.system('python folding.py ' + str(initial) + ' ' + str(final) + ' temph5')
    os.system('python plot_pulse.py ' + str(initial) +' ' + str(final) + ' temph5 >> bin_number_2139_wz.txt')
    os.system('rm temph5')
