import h5py
import numpy as np
import matplotlib.pyplot as plt


def plot_histogram():

f = h5py.File('/home/p/pen/hsiuhsil/pulsar_search/J2139+00h5', 'r')

time = []
time_interval = 1000
bin = []
for ii in xrange(len(f['ABS_TIME'])):
    time += f['ABS_TIME'][ii]

for ii in xrange((max(time)-min(time))//100+1):
    bin += [min(time)+ii*100]

plt.hist(time, bins=bin)
plt.show()
