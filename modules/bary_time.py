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
    this_file.create_dataset('dBATdra', (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)
    this_file.create_dataset('dBATddec', (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype, chunks=True)

    for ii in range(len(this_file['ABS_TIME'])-1):
        print ii
#        RA = str(this_file['RA_SUB'][ii])
#        DEC = str(this_file['DEC_SUB'][ii])
        '''J2139+00, tempo + delta_ra'''
#        RA = deg_to_HMS(324.86337029308453)
#        DEC = deg_to_DMS(0.60222)
        '''J2139+00, tempo'''
#        RA = deg_to_HMS(324.92817)
#        DEC = deg_to_DMS(0.60222)
        '''J2139+00, atnf'''
#        RA = str("+21:39:46")
#        DEC = str("+00:36:00")
        '''J2139+00, atnf + delta(ra, dec)'''
        RA = convHMS(str("+21:39:46"))
        '''test new location with factor of -1'''
        RA += -1 * 1.72450646e-03 * 180 / np.pi # convert rad to deg
        RA = deg_to_HMS(RA)
        DEC = convDMS(str("+00:36:00"))
        '''test new location with factor of -1'''
        DEC += -1 * -1.67417143e-03 * 180 / np.pi # convert rad to deg
        DEC = deg_to_DMS(DEC)
        '''J0030+0451'''
#        RA = str('00:30:27.42823')
#        DEC = str('+04:51:39.7112')
        '''J0051+0423'''   
#        RA = str('00:51:30.1')
#        DEC = str('+04:22:49')
        '''J1038+0032'''
#        RA = str('10:38:26.933')
#        DEC = str('+00:32:43.6') 
        '''J1046+0304'''
#        RA = str('10:46:43.23')
#        DEC = str('+03:04:06.9')

        RA1 = deg_to_HMS(convHMS(RA)+0.1)
        DEC1 = deg_to_DMS(convDMS(DEC)+0.1)
  
        topo_time = repr(this_file['ABS_TIME'][ii]/86400)
        p = subprocess.Popen(["bary", "GBT", RA, DEC], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        BAT = p.communicate(input=topo_time)[0].split()[1]

        p_ra1 = subprocess.Popen(["bary", "GBT", RA1, DEC], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        BAT_ra1 = p_ra1.communicate(input=topo_time)[0].split()[1]
        dBATdra = str((np.float64(BAT_ra1) - np.float64(BAT)) / 0.1)

        p_dec1 = subprocess.Popen(["bary", "GBT", RA, DEC1], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        BAT_dec1 = p_dec1.communicate(input=topo_time)[0].split()[1]
        dBATddec = str((np.float64(BAT_dec1) - np.float64(BAT)) / 0.1)

        current_len_topo = this_file['TOPO_TIME'].shape[0]
        current_len_bary = this_file['BARY_TIME'].shape[0]
        current_len_dBATdra = this_file['dBATdra'].shape[0]
        current_len_dBATddec = this_file['dBATddec'].shape[0]

        this_file['TOPO_TIME'].resize(current_len_topo + 1, 0)
        this_file['TOPO_TIME'][ii] = np.float64(topo_time)
        this_file['BARY_TIME'].resize(current_len_bary + 1, 0)
        this_file['BARY_TIME'][ii] = np.float64(BAT)
        this_file['dBATdra'].resize(current_len_dBATdra + 1, 0)
        this_file['dBATdra'][ii] = np.float64(dBATdra)
        this_file['dBATddec'].resize(current_len_dBATddec + 1, 0)
        this_file['dBATddec'][ii] = np.float64(dBATddec)

if __name__ == '__main__':
    main()
