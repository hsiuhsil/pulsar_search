import h5py
import pyfits
import numpy as np


J0030 = h5py.File('J0030.hdf5','a')
J0030_TSUBINT = J0030.create_dataset('TSUBINT', (10,1), maxshape=(None,1), dtype=np.float64)
J0030_OFFS_SUB = J0030.create_dataset('OFFS_SUB', (10,1), maxshape=(None,1), dtype=np.float64)
J0030_LST_SUB = J0030.create_dataset('LST_SUB', (10,1), maxshape=(None,1), dtype=np.float64)
J0030_RA_SUB = J0030.create_dataset('RA_SUB', (10,1), maxshape=(None,1), dtype=np.float64)
J0030_DEC_SUB = J0030.create_dataset('DEC_SUB', (10,1), maxshape=(None,1), dtype=np.float64)
J0030_GLON_SUB = J0030.create_dataset('GLON_SUB', (10,1), maxshape=(None,1), dtype=np.float64)
J0030_GLAT_SUB = J0030.create_dataset('GLAT_SUB', (10,1), maxshape=(None,1), dtype=np.float64)
J0030_FD_ANG = J0030.create_dataset('FD_ANG', (10,1), maxshape=(None,1), dtype=np.float32)
J0030_POS_ANG = J0030.create_dataset('POS_ANG', (10,1), maxshape=(None,1), dtype=np.float32)
J0030_PAR_ANG = J0030.create_dataset('PAR_ANG', (10,1), maxshape=(None,1), dtype=np.float32)
J0030_TEL_AZ = J0030.create_dataset('TEL_AZ', (10,1), maxshape=(None,1), dtype=np.float32)
J0030_TEL_ZEN = J0030.create_dataset('TEL_SEN', (10,1), maxshape=(None,1), dtype=np.float32)
J0030_DAT_FREQ = J0030.create_dataset('DAT_FREQ', (10,4096), maxshape=(None,4096), dtype=np.float32)
J0030_DAT_WTS = J0030.create_dataset('DAT_WTS', (10,4096), maxshape=(None,4096), dtype=np.float32)
J0030_DAT_OFFS = J0030.create_dataset('DAT_OFFS', (10,16384), maxshape=(None,16384), dtype=np.float32)
J0030_DAT_SCL = J0030.create_dataset('DAT_SCL', (10,16384), maxshape=(None,16384), dtype=np.float32)
J0030_DATA = J0030.create_dataset('DATA', (10,2048,4,4096,1), maxshape=(None,2048,4,4096,1))

J0051 = h5py.File('J0051.hdf5','a')
J0051_TSUBINT = J0051.create_dataset('TSUBINT', (10,1), maxshape=(None,1), dtype=np.float64)
J0051_OFFS_SUB = J0051.create_dataset('OFFS_SUB', (10,1), maxshape=(None,1), dtype=np.float64)
J0051_LST_SUB = J0051.create_dataset('LST_SUB', (10,1), maxshape=(None,1), dtype=np.float64)
J0051_RA_SUB = J0051.create_dataset('RA_SUB', (10,1), maxshape=(None,1), dtype=np.float64)
J0051_DEC_SUB = J0051.create_dataset('DEC_SUB', (10,1), maxshape=(None,1), dtype=np.float64)
J0051_GLON_SUB = J0051.create_dataset('GLON_SUB', (10,1), maxshape=(None,1), dtype=np.float64)
J0051_GLAT_SUB = J0051.create_dataset('GLAT_SUB', (10,1), maxshape=(None,1), dtype=np.float64)
J0051_FD_ANG = J0051.create_dataset('FD_ANG', (10,1), maxshape=(None,1), dtype=np.float32)
J0051_POS_ANG = J0051.create_dataset('POS_ANG', (10,1), maxshape=(None,1), dtype=np.float32)
J0051_PAR_ANG = J0051.create_dataset('PAR_ANG', (10,1), maxshape=(None,1), dtype=np.float32)
J0051_TEL_AZ = J0051.create_dataset('TEL_AZ', (10,1), maxshape=(None,1), dtype=np.float32)
J0051_TEL_ZEN = J0051.create_dataset('TEL_SEN', (10,1), maxshape=(None,1), dtype=np.float32)
J0051_DAT_FREQ = J0051.create_dataset('DAT_FREQ', (10,4096), maxshape=(None,4096), dtype=np.float32)
J0051_DAT_WTS = J0051.create_dataset('DAT_WTS', (10,4096), maxshape=(None,4096), dtype=np.float32)
J0051_DAT_OFFS = J0051.create_dataset('DAT_OFFS', (10,16384), maxshape=(None,16384), dtype=np.float32)
J0051_DAT_SCL = J0051.create_dataset('DAT_SCL', (10,16384), maxshape=(None,16384), dtype=np.float32)
J0051_DATA = J0051.create_dataset('DATA', (10,2048,4,4096,1), maxshape=(None,2048,4,4096,1))




class FileSource(object):
    def __init__(self, fitsfile):
        self._filename = fitsfile

    def get_records(self):
        hduls = pyfits.open(self._filename, 'readonly')
        data = read_records(hduls)
        hduls.close()
        return data


def search_pulsar(hduls):

    pulsar = [['wigglez1hr_centre',
               ['J0030+0451', 7.61428, 4.86103],
               ['J0051+0423', 12.87542, 4.38028],
              ],
              ['wigglez11hr_centre',
               ['J1038+0032', 159.61222, 0.54544],
               ['J1046+0304', 161.68013, 3.06858],
              ],
              ['wigglez15hr_centre',
               ['J1501-0046', 225.43732, -0.77320],
               ['J1518+0204A', 229.63882, 2.09099],
               ['J1518+0204B', 229.63107, 2.08763],
               ['J1518+0204C', 229.63662, 2.07995],
               ['J1518+0204D', 229.64167, 2.08278],
               ['J1518+0204E', 229.64167, 2.08278],
              ],
              ['wigglez22hr_centre',
               ['J2139+00', 324.92500, 0.60000],
               ['J2222-0137', 335.52487, -1.62103]
              ]
             ]
    """ pulsar[i][0] is wigglez field, pulsar[i][j][1] is its ra, and  pulsar[i][j][2] is its dec."""

    keys = hduls[1].columns.names
    ii = 0
    for k in xrange(len(hduls[1].data)):
        for i in xrange(len(pulsar)):
            if hduls[0].header['SRC_NAME'] == pulsar[i]:
                for j in xrange(1,len(pulsar[i])):
                    delta_ra = hduls[1].data[k]['RA_SUB'] - pulsar[i][j][1]
                    delta_dec = hduls[1].data[k]['DEC_SUB'] - pulsar[i][j][2]
                    scope = np.sqrt(delta_ra**2 + delta_dec**2)
                    if scope <= 0.25:
                        if pulsar[i][j][0] == 'J0030+0451':
                            for jj in xrange(len(keys)):
                                J0030[keys[jj]][ii,:] = hduls[1].data[k][keys[jj]]
                                ii += 1
                        elif pulsar[i][j][0] == 'J0051+0423':
                            for jj in xrange(len(keys)):
                                J0051[keys[jj]][ii,:] = hduls[1].data[k][keys[jj]]
                                ii += 1
