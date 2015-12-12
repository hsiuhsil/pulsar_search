import h5py
import pyfits
import numpy as np

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
    
    '''create dataset for each pulsar location'''
    
    files = {}
    for i in xrange(len(pulsar)):
        for j in xrange(1,len(pulsar[i])):
            this_file = h5py.File(pulsar[i][j][0] + 'h5','a')
            for dataset_name in keys:
                first_data = hduls[1].data[0][dataset_name]
                this_file.create_dataset(dataset_name, (0,) + first_data.shape, maxshape = (None,) +first_data.shape, dtype=first_data.dtype)
            files[pulsar[i][j][0]] = this_file


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
