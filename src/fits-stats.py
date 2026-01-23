# /// script
# requires-python = ">=3.8"
# dependencies = [
#     "astropy>=4.0",
#     "scipy>=1.5",
# ]
# ///
"""
Calculate various statistics of a FITS file
"""

import argparse
from astropy.io import fits as pyfits
import scipy
from scipy.ndimage import extrema

parser = argparse.ArgumentParser(
    description="Calculate various statistics of a FITS file"
)
parser.add_argument("file", help="FITS file name")
args = parser.parse_args()

data = pyfits.open(args.file)[0].data
dmin, dmax, dminloc, dmaxloc = extrema(data)

print("Minimum: %s at %s" % (dmin, dminloc))
print("Maximum: %s at %s" % (dmax, dmaxloc))
