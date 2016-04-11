"""Calculate plane-of-sky stats from moment maps

Create a table of time series.  The maps should be previously
generated using ppv_to_moments.py

"""
from __future__ import print_function
import sys
import glob
import numpy as np
from astropy.io import fits
from astropy.table import Table


def stats_from_momfiles(fn):
    """Return a dict of plane-of-sky average and standard dev of the
line-of-sight mean velocity and rms widths

    """
    # ridiculous hack to extract time from filename
    itime = int(fn.split('-')[0].split('_')[-1].replace('vsum', ''))
    # Read the maps
    vsum = fits.open(fn)[0].data
    vmean = fits.open(fn.replace('vsum', 'vmean'))[0].data
    vsig = fits.open(fn.replace('vsum', 'vsig'))[0].data
    # Find POS mean and s.d of each
    vm_mean = np.nansum(vsum*vmean)/np.nansum(vsum)
    vm_sd = np.sqrt(np.nansum(vsum*(vmean - vm_mean)**2)/np.nansum(vsum))
    vs_mean = np.nansum(vsum*vsig)/np.nansum(vsum)
    vs_sd = np.sqrt(np.nansum(vsum*(vsig - vs_mean)**2)/np.nansum(vsum))
    return {
        'time': itime,
        'vmean': vm_mean,
        'd vmean': vm_sd,
        'vsig': vs_mean,
        'd vsig': vs_sd,
    }


COLUMNS = ['time', 'vmean', 'd vmean', 'vsig', 'd vsig']
if __name__ == '__main__':
    try:
        prefix = sys.argv[1]   # e.g., 04052012_4
        suffix = sys.argv[2]    # e.g., xp-N26584
    except:
        print('Usage: {} PREFIX SUFFIX'.format(sys.argv[0]))

    pattern = '{}_*vsum-{}.fits'.format(prefix, suffix)
    filenames = glob.glob(pattern)
    filenames.sort()
    data = []
    for fn in filenames:
        data.append(stats_from_momfiles(fn))
    tab = Table(rows=data, names=COLUMNS)
    tabname = '{}-POS-stats-{}.tab'.format(prefix, suffix)
    tab.write(tabname, format='ascii.tab')
    
