import sys
import os
import os.path

import h5py
import numpy as np
from numpy import array

mins = 5.0

'''new[ii] would be index kk in this_file[kk]'''
all = np.arange(0,1002)
'''index 379 to 439 is repeated, and [9,10,...917] is for strong RFI'''
index_removed = np.append(np.arange(379,439),np.array([9, 10, 540,764, 856, 867, 889, 917]))
new = np.delete(all, index_removed)


#all = np.arange(0,100)
#index_removed = np.arange(20,50)
#new = np.delete(all, index_removed)


def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            index_select(filename)
        except (IOError, ValueError):
            print IOError

def index_select(filename):
    this_file = h5py.File(filename, "r")
    previous_initial = np.zeros((len(new),), dtype=np.int)
    previous_final = np.zeros((len(new),), dtype=np.int)
    for ii in range(len(new)):
        kk = new[ii]
        previous_initial[ii] = np.int(kk)
        '''convert to mins'''
        diff = 24*60*(this_file['BARY_TIME'][:]-this_file['BARY_TIME'][kk]) 
        index = np.where(np.logical_and(0<diff, diff<mins))
        if len(index[0]) == 0:
            previous_final[ii] = np.int(kk)
        else:
            previous_final[ii] = np.int(index[0][-1])

    list = []
    list.append([previous_initial[0],previous_final[0]])
    for ii in range(len(previous_final)-1):
        if previous_final[ii] != previous_final[ii+1]:
            list.append([previous_initial[ii+1],previous_final[ii+1]])
    scan = np.array(list) 

    np.savetxt('test.txt', scan, fmt='%4d')    

if __name__ == '__main__':
    main()




