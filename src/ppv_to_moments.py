"""Calculate velocity moments from PPV cube

Assumes that first axis (in FITS sense) is the velocity and that it
has 1 km/s pixels

"""
from __future__ import print_function
import sys
import numpy as np
from astropy.io import fits

# 06 May 2016: Based on the Fortran implementation in
# /Users/will/Work/BobKPNO/src/newlinemod.f90
def find_fwhm(f, v, frac=0.5):
    ipeak = np.argmax(f)
    fpeak = f[ipeak]
    ileft = np.argmin(f[f >= frac*fpeak])
    iright = np.argmax(f[f >= frac*fpeak])
    if ileft <= 0:
        uleft = v[0]
    elif ileft >= len(f):
        uleft = v[-1]
    else:
        uleft = (
            v[ileft] -
            (v[ileft] - v[ileft-1])
            * (f[ileft] - frac*fpeak)
            / (f[ileft] - f[ileft-1])
        )
    if iright < 0:
        uright = v[0]
    elif iright >= len(f):
        uright = v[-1]
    else:
        uright = (
            v[iright] +
            (v[iright+1] - v[iright])
            * (f[iright] - frac*fpeak)
            / (f[iright] - f[iright-1])
        )
    return uright - uleft


def find_peak(f, v):
    raise NotImplementedError


if __name__ == '__main__':
    try:
        fn = sys.argv[1]
    except:
        print('Usage {} FILENAME'.format(sys.argv[0]))
        
    uleft, uright
    
    hdu, = fits.open(fn)
    # Velocity is first FITS axis; last python axis
    ny, nx, nv = hdu.data.shape
    assert(nv == 81)
    vels = np.linspace(-40.0, 40.0, nv).reshape((1, 1, nv))

    # Eliminate negative values
    hdu.data[hdu.data < 0.0] = 0.0
    
    # zeroth moment
    m0 = np.sum(hdu.data, axis=-1)
    # first moment
    m1 = np.sum(hdu.data*vels, axis=-1)
    vmean = m1/m0
    # centered 2nd moment
    m2 = np.sum(hdu.data*(vels - vmean[:, :, None])**2, axis=-1)
    sigma = np.sqrt(m2/m0)

    fits.PrimaryHDU(data=m0).writeto(fn.replace('vcube', 'vsum'), clobber=True)
    fits.PrimaryHDU(data=vmean).writeto(fn.replace('vcube', 'vmean'), clobber=True)
    fits.PrimaryHDU(data=sigma).writeto(fn.replace('vcube', 'vsig'), clobber=True)

    
