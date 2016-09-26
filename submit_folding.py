import os
import sys
import numpy as np

#p = np.loadtxt('./test.txt', delimiter=',')
p = np.loadtxt('./test.txt')
filename = sys.argv[1]
for i in xrange(len(p)):
    initial, final = int(p[i][0]), int(p[i][1])
    os.system('cp ' + filename + ' temph5')
    os.system('python folding.py ' + str(initial) + ' ' + str(final) + ' temph5')
    os.system('python plot_pulse.py ' + str(initial) +' ' + str(final) + ' temph5')
    os.system('rm temph5')
