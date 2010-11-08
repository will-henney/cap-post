"""
Calculate various statistics of a FITS file
"""
import argparse
import pyfits
import scipy
from scipy.ndimage.measurements import extrema

parser = argparse.ArgumentParser(description="Calculate various statistics of a FITS file")
parser.add_argument('file', help="FITS file name")
args = parser.parse_args()

data = pyfits.open(args.file)[0].data
dmin, dmax, dminloc, dmaxloc = extrema(data)

print "Minimum: %s at %s" % (dmin, dminloc)
print "Maximum: %s at %s" % (dmax, dmaxloc)




