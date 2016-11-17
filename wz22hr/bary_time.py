import sys
import os
import os.path

import h5py
import pyfits
import numpy as np
import re
import subprocess

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            print filename
            find_topo_bary(filename)
        except (IOError, ValueError):
            print IOError

def convHMS(ra):
   try :
      sep1 = ra.find(':')
      hh=int(ra[0:sep1])
      sep2 = ra[sep1+1:].find(':')
      mm=int(ra[sep1+1:sep1+sep2+1])
      ss=float(ra[sep1+sep2+2:])
   except:
      raise
   else:
      pass  
   return(hh*15.+mm/4.+ss/240.)

def convDMS(dec):
   Csign=dec[0]
   if Csign=='-':
      sign=-1.
      off = 1
   elif Csign=='+':
      sign= 1.
      off = 1
   else:
      sign= 1.
      off = 0
   try :
      sep1 = dec.find(':')
      deg=int(dec[off:sep1])
      sep2 = dec[sep1+1:].find(':')
      arcmin=int(dec[sep1+1:sep1+sep2+1])
      arcsec=float(dec[sep1+sep2+2:])
   except:
      raise
   else:
      pass
   return(sign*(deg+(arcmin*5./3.+arcsec*5./180.)/100.))            

def deg_to_HMS( RA ):
    if(RA<0):
        sign = -1
        ra   = -RA
    else:
        sign = 1
        ra   = RA
    h = int( ra/15. )
    ra -= h*15.
    m = int( ra*4.)
    ra -= m/4.
    s = ra*240.
    if(sign == -1):
        out = '-%02d:%02d:%06.3f'%(h,m,s)
    else:
        out = '+%02d:%02d:%06.3f'%(h,m,s)
    return out

def deg_to_DMS( Dec ):
    if(Dec<0):
        sign = -1
        dec  = -Dec
    else:
        sign = 1
        dec  = Dec
    d = int( dec )
    dec -= d
    dec *= 100.
    m = int( dec*3./5. )
    dec -= m*5./3.
    s = dec*180./5.
    if(sign == -1):
        out = '-%02d:%02d:%06.3f'%(d,m,s)
    else:
        out = '+%02d:%02d:%06.3f'%(d,m,s)
    return out


def find_topo_bary(filename):

    this_file = h5py.File(filename, "r+")
    first_data = this_file['ABS_TIME'][0]
    this_file.create_dataset('BARY_TIME', (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
    this_file.create_dataset('TOPO_TIME', (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)

    '''Make sure RA and DEC are in the unit of hh:mm:ss'''
    RA = deg_to_HMS(float(re.split('_',filename)[1]))
    DEC = deg_to_DMS(float(re.split('_',filename)[3]))

    for ii in range(len(this_file['ABS_TIME'])-1):
        print ii
#        RA = str(this_file['RA_SUB'][ii])
#        DEC = str(this_file['DEC_SUB'][ii])
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
