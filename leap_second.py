import sys
import os
import os.path

import h5py
import pyfits
import numpy as np
import subprocess

MJD = 56109.00 #54832         
RA = str(324.92817)
DEC = str(0.60222)

def check_leap_second():
    topo_before = repr(MJD -1 + 86399/86400.)
    topo_middle = repr(MJD)
    p = subprocess.Popen(["bary", "GBT", RA, DEC], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bary_before = p.communicate(input=topo_before)[0].split()[1]
    p = subprocess.Popen(["bary", "GBT", RA, DEC], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bary_middle = p.communicate(input=topo_middle)[0].split()[1]
    bary_diff_standard = (np.float(bary_middle) - np.float(bary_before))*86400
    print 'bary_diff_standard: '+str(bary_diff_standard)
    for ii in range(86400):
        topo_1 = repr(MJD + ii/86400.)
        topo_2 = repr(MJD +(ii+1)/86400.)
        p = subprocess.Popen(["bary", "GBT", RA, DEC], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bary_1 = p.communicate(input=topo_1)[0].split()[1]
        p = subprocess.Popen(["bary", "GBT", RA, DEC], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bary_2 = p.communicate(input=topo_2)[0].split()[1]
        bary_diff = (np.float(bary_2) - np.float(bary_1))*86400
        if bary_diff_standard != bary_diff:
            print 'at ii:'+str(ii)+', bary_diff for one sec:' +str(bary_diff)

check_leap_second()
