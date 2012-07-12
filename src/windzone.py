"""
Find the zone that is dominated by the stellar wind ram pressure
"""

import pyfits
import numpy as np
from scipy import ndimage
import argparse

def parseargs:
    parser = argparse.ArgumentParser(
        description="Find the zone that is dominated by the stellar wind ram pressure")
    parser.add_argument('id', help="File prefix", default="04052012_4_0030")
    return parser.parse_args()

cmd_args = parseargs()

p, = pyfits.open(cmd_args.id + "p.fits") # gas pressure
# cubewind must have been run first to generate these files
pr, = pyfits.open(cmd_args.id + "pr.fits") # HII region inward ram pressure
pw, = pyfits.open(cmd_args.id + "pw.fits") # Wind ram pressure


##
## First approximation - assume that shell is momentum driven
##
# logical array where wind ram pressure beats all else
iswind = pw.data > (pr.data + p.data)

# divide iswind into mutually connected "features"
structure = np.ones((3, 3, 3))  # include diagonal connections to boost percolation
labeled_array, num_features = ndimage.label(iswind, structure=structure)

# we are only interested in the central "feature" around the star
nz, ny, nx = p.shape
ilabel = labeled_array[nz/2, ny/2, nx/2] # find label of central feature
iswind = labeled_array == ilabel

##
## Note that this is not really realistic since if it were really
## momentum driven, then there would be radial shadowing
##

# X-ray emissivity
em = iswind.astype(np.float)    # 1.0 or 0.0 where iswind is True or False

hdu = pyfits.PrimaryHDU(em)
hdu.writeto(cmd_args.id + "e-X00mom.fits")

