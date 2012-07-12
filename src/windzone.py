"""
Find the zone that is dominated by the stellar wind ram pressure
"""

import pyfits
import numpy as np
from scipy import ndimage
import argparse

def parseargs():
    parser = argparse.ArgumentParser(
        description="Find the zone that is dominated by the stellar wind ram pressure")
    parser.add_argument('id', type=str, help="File prefix", default="04052012_4_0030")
    parser.add_argument('--verbose', '-v', action="store_true", help="Give verbose info about progress")
    parser.add_argument('--boost', type=float, default=1.0, 
                        help="Boost factor for wind momentum (wrt M/1e-7 V/1e3 = 4.2)")
    parser.add_argument('--connection', type=str, default="all", choices=("all", "ortho"), 
                        help="Allowed connection for a contiguous region")
    return parser.parse_args()

cmd_args = parseargs()

if cmd_args.verbose: print "Reading pressure cubes..."
p, = pyfits.open(cmd_args.id + "p.fits") # gas pressure
# cubewind must have been run first to generate these files
pr, = pyfits.open(cmd_args.id + "pr.fits") # HII region inward ram pressure
pw, = pyfits.open(cmd_args.id + "pw.fits") # Wind ram pressure


##
## First approximation - assume that shell is momentum driven
##
# logical array where wind ram pressure beats all else
if cmd_args.verbose: print "Comparing pressures..."
iswind = cmd_args.boost*pw.data > (pr.data + p.data)

# divide iswind into mutually connected "features"
if cmd_args.verbose: print "Labeling features..."

if cmd_args.connection == "all":
    # include diagonal connections to boost percolation
    structure = np.ones((3, 3, 3))  
else:
    # include only orthogonal connections
    structure = np.array([
        [[0, 0, 0], [0, 1, 0], [0, 0, 0]],
        [[0, 1, 0], [1, 1, 1], [0, 1, 0]],
        [[0, 0, 0], [0, 1, 0], [0, 0, 0]]
        ])

labeled_array, num_features = ndimage.label(iswind, structure=structure)

# we are only interested in the central "feature" around the star
nz, ny, nx = p.data.shape
ilabel = labeled_array[nz/2, ny/2, nx/2] # find label of central feature
if cmd_args.verbose: print "Isolating central feature..."
iswind = labeled_array == ilabel

# Find volume of wind region
ncells = iswind.sum()
print "Volume fraction occupied by wind: ", float(ncells) / (nx*ny*nz)

PARSEC = 3.085677582e18
YEAR = 3.15576e7
MSUN = 1.989e33
KM = 1.0e5
dx = (4.0 / nx) * PARSEC
t = (10000.0*YEAR)*float(cmd_args.id.split("_")[-1])

windvolume = ncells*dx**3
# canonical wind parameters from GAHA2001 sec 5.2 
Mdot = (1.e-6*MSUN/YEAR) * 0.35 * cmd_args.boost
Vwind = (1000.0*KM) * 1.2
# Total wind energy
Ewind = 0.5 * Mdot * Vwind**2
# Wind pressure: (gamma - 1) E / V
Pwind = (2./3.)*Ewind/windvolume

print "Shocked wind pressure in energy-driven case would be: ", Pwind

##
## Note that this is not really realistic since if it were really
## momentum driven, then there would be radial shadowing
##

# X-ray emissivity
if cmd_args.verbose: print "Calculating x-ray emissivity..."
em = iswind.astype(np.float32)    # 1.0 or 0.0 where iswind is True or False

if cmd_args.verbose: print "Writing emissivity cube..."
hdu = pyfits.PrimaryHDU(em)
hdu.writeto(cmd_args.id + "e-X{:02}mom.fits".format(int(cmd_args.boost)), clobber=True)

if cmd_args.verbose: print "Removing hole from Ha emissivity cube..."
ha, = pyfits.open(cmd_args.id + "e-Halpha.fits")
ha.data *= 1.0 - em             # remove X-ray hole from H alpha
ha.writeto(cmd_args.id + "e-HaX{:02}m.fits".format(int(cmd_args.boost)), clobber=True)
