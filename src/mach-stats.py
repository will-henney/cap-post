import sys
import numpy as np
from astropy.io import fits

def find_mach_stats(prefix='04052012_4_0030'):
    u = fits.open(prefix + 'u.fits')[0].data 
    v = fits.open(prefix + 'v.fits')[0].data 
    w = fits.open(prefix + 'w.fits')[0].data 
    eo3 = fits.open(prefix + 'e-O35007.fits')[0].data
    en2 = fits.open(prefix + 'e-N26584.fits')[0].data

    for label, e in ('[O III]', eo3), ('[N II]', en2):
        u0 = np.average(u, weights=e)
        v0 = np.average(v, weights=e)
        w0 = np.average(w, weights=e)

        vv = (u - u0)**2 + (v - v0)**2 + (w - w0)**2
        vrms = np.sqrt(np.average(vv, weights=e))
        print(label, 'mean = [{:.2f}, {:.2f}, {:.2f}]'.format(u0, v0, w0))
        print(label, 'RMS = {:.2f}'.format(vrms))

if __name__ == '__main__':
    try:
        # Example prefix: 04052012_4_0030
        prefix = sys.argv[1]
    except IndexError:
        sys.exit('Usage: {} PREFIX'.format(sys.argv[0]))
    find_mach_stats(prefix)
