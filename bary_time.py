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
            find_bary(filename)
        except (IOError, ValueError):
            print IOError
            

def find_bary(filename):

    this_file = h5py.File(filename, "r+")
    hduls = pyfits.open('/scratch2/p/pen/hsiuhsil/gbt_data/FRB110523/guppi_55704_wigglez22hrst_0258_0001.fits')
    first_data = hduls[1].data[0]['OFFS_SUB']
    this_file.create_dataset('BARY_TIME', (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
    os.system('bary')
#    print len(this_file['ABS_TIME'])
    for ii in range(len(this_file['ABS_TIME'])-1):
        RA = str(this_file['RA_SUB'][ii])
        DEC = str(this_file['DEC_SUB'][ii])
        topo_time = str(this_file['ABS_TIME'][ii]/86400)
        p = subprocess.Popen(["bary", "GBT", RA, DEC], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bary_time = p.communicate(input=topo_time)[0].split()[1]
        current_len = this_file['BARY_TIME'].shape[0]
        this_file['BARY_TIME'].resize(current_len + 1, 0)
        this_file['BARY_TIME'][ii] = np.float64(bary_time)
  
  
  if __name__ == '__main__':
    main()
