import h5py
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

this_file = h5py.File('positionsh5', "r")
x = this_file['POS'][:-1,0]
y = this_file['POS'][:-1,1]
plt.scatter(x, y)
plt.xlabel('RA (degs)', fontsize=14)
plt.ylabel('DEC (degs)', fontsize=14)
plt.savefig('wz_pos.png')
print 'save png'
#print ' len(np.arange(np.amin(this_file['POS'][:-1,0]), np.amax(this_file['POS'][:,0]), 0.15))'
