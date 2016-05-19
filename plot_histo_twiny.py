import sys
import os.path

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import pylab
from pylab import *

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            plot_time(filename)
        except (IOError, ValueError):
            print 'Skipped:'+ filename


def plot_time(filename):

    f = h5py.File(filename, 'r')

    time_mjd = []
    time_days = []
#    for ii in xrange(len(f['ABS_TIME'])-1):
    for ii in xrange(300):
        mjd = f['ABS_TIME'][ii]/86400
        time_mjd.append(mjd)
        days = (f['ABS_TIME'][ii]-f['ABS_TIME'][0])/86400
        time_days.append(days)

    xmax = time_days[-1]
    xmin = time_days[0]
        
    total_number = len(time_days)
    file_name = f.filename[48:]
    text = str(file_name) + ', total number: ' + str(total_number)
    offset = time_mjd[0]  

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()

    ax1.hist(time_days, bins=2000)
    ax1.set_xlabel('days,'+' offset of MJD: '+str(int(offset)), fontsize=20)
    ax1.set_ylabel('number of data', fontsize=20)
    ax = plt.gca()
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.set_xlim((xmin,xmax))

    def firstday_years(mjd):
        if 55562 <= mjd < 55926:
            year = 2011
        elif 55927 <= mjd < 56292:
            year = 2012
        elif 56293 <= mjd < 56657:
            year = 2013
        elif 56658 <= mjd < 57022:
            year = 2014
        elif 57023 <= mjd :
            year = 2015
        return year

    year_label = []
    first_year = min(firstday_years(time_mjd[0])+1,firstday_years(time_mjd[-1]))
    last_year = max(firstday_years(time_mjd[0])+1,firstday_years(time_mjd[-1]))
    space = last_year - first_year
    for ii in xrange(first_year, last_year + 1):
        year_label.append(ii)
 
    year_label_position = [] 
    year_day_diff = [55562 - time_mjd[0],
                     55927 - time_mjd[0],
                     56293 - time_mjd[0],
                     56658 - time_mjd[0],
                     57023 - time_mjd[0]
                    ] 

    min_diff = []
    for ii in xrange(len(year_day_diff)):
        if year_day_diff[ii] >= 0:
            min_diff.append(year_day_diff[ii])
            min_year = min(min_diff)
    for ii in xrange(0,space+1):
        loc = min_year + 365*ii                 
        year_label_position.append(loc) 

    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(year_label_position)
    ax2.set_xticklabels(year_label)
    ax2.set_xlabel('Years', fontsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=20)

    title = ax1.set_title(text, fontsize=20)
    title.set_y(1.1)
    fig.subplots_adjust(top=0.85)

    plt.show()


if __name__ == '__main__':
    main()
