import sys
import os
import os.path

import h5py
import pyfits
import numpy as np
import subprocess

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            find_topo_bary(filename)
        except (IOError, ValueError):
            print IOError
            

def find_topo_bary(filename):

    this_file = h5py.File(filename, "r+")
    first_data = this_file['ABS_TIME'][0]
    this_file.create_dataset('BARY_TIME', (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
    this_file.create_dataset('TOPO_TIME', (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)

    for ii in range(len(this_file['ABS_TIME'])-1):
        print ii
#        RA = str(this_file['RA_SUB'][ii])
#        DEC = str(this_file['DEC_SUB'][ii])
#        RA = str(324.92500)
#        DEC = str(0.60000)
        RA = str(324.92817)
        DEC = str(0.60222)
        topo_time = repr(this_file['ABS_TIME'][ii]/86400)
        p = subprocess.Popen(["bary", "GBT", RA, DEC], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bary_time = p.communicate(input=topo_time)[0].split()[1]
        current_len_topo = this_file['TOPO_TIME'].shape[0]
        current_len_bary = this_file['BARY_TIME'].shape[0]
        this_file['TOPO_TIME'].resize(current_len_topo + 1, 0)
        this_file['TOPO_TIME'][ii] = np.float64(topo_time)
        this_file['BARY_TIME'].resize(current_len_bary + 1, 0)
        this_file['BARY_TIME'][ii] = np.float64(bary_time)
  
if __name__ == '__main__':
    main()
