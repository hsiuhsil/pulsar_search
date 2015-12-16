import h5py
import numpy as np
import matplotlib.pyplot as plt


f = h5py.File('/home/p/pen/hsiuhsil/pulsar_search/J2139+00h5', 'r')


time = []
for ii in xrange(len(f['ABS_TIME'])-1):
    time.append(f['ABS_TIME'][ii])


print len(time) 

plt.hist(time, bins=20)
plt.title('J2139+00', fontsize=20)
plt.xlabel('time(s)', fontsize=20)
plt.ylabel('number of data', fontsize=20)
plt.show()
