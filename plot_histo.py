import sys
import os.path

import h5py
import numpy as np
import matplotlib.pyplot as plt

def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            plot_time(filename)
        except (IOError, ValueError):
            print 'Skipped:'+ filename

def plot_time(filename):
    f = h5py.File(filename, 'r')
    time_days = []
    time_months = []
    for ii in xrange(len(f['ABS_TIME'])-1):
        time_days.append(f['ABS_TIME'][ii]/86400)
        time_months.append(f['ABS_TIME'][ii]/86400/30)

    total_number = len(time_days)
    file_name = f.filename[35:]
    text = str(file_name) + ', total number: ' + str(total_number)

    plt.figure(1)
    plt.subplot(211)
    plt.hist(time_days, bins=200)
    plt.title(text, fontsize=20)
    plt.xlabel('MJD', fontsize=20)
    plt.ylabel('number of data', fontsize=20)
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    plt.tick_params(axis='both', which='major', labelsize=20)

    plt.subplot(212)
    plt.hist(time_months, bins=200)
    plt.xlabel('time(months)', fontsize=20)
    plt.ylabel('number of data', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()

