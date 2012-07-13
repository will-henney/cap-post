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
    parser.add_argument('--energy', action="store_true", 
                        help="Calculate energy driven instead of momentum driven case")
    parser.add_argument('--stats', action="store_true", 
                        help="Calculate statistics for the wind zone (which may be slow if it is large)")
    parser.add_argument('--slowmethod', action="store_true", 
                        help="DEPRECATED Use a method that can find the inner wind zone even if it doesn't include the grid center")
    parser.add_argument('--efrac', type=float, default=5.0/11.0, 
                        help="Fraction of wind mechanical luminosity that goes into thermal energy of bubble")
    parser.add_argument('--pguess', type=float, default=3.0e-10, 
                        help="Fraction of wind mechanical luminosity that goes into thermal energy of bubble")
    parser.add_argument('--connection', type=str, default="ortho", choices=("all", "ortho"), 
                        help="Allowed connection for a contiguous region")
    parser.add_argument('--onlyinfo', action="store_true", help="Only give info; write no emissivity cubes")
    return parser.parse_args()

cmd_args = parseargs()

if cmd_args.verbose: print "Reading pressure cubes..."
p, = pyfits.open(cmd_args.id + "p.fits") # gas pressure
# cubewind must have been run first to generate these files
pr, = pyfits.open(cmd_args.id + "pr.fits") # HII region inward ram pressure
pw, = pyfits.open(cmd_args.id + "pw.fits") # Wind ram pressure

ptot = pr.data + p.data

##
## First approximation - assume that shell is momentum driven
##
# logical array where wind ram pressure beats all else
if cmd_args.verbose: print "Comparing pressures..."
if cmd_args.energy:
    # energy-driven case
    # but allow wind ram pressure to help too
    iswind = (cmd_args.pguess + cmd_args.boost*pw.data) > ptot
else:
    # momentum-driven case
    iswind = cmd_args.boost*pw.data > ptot

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

if cmd_args.slowmethod:
    # This should not be necessary any more
    # array of radii from central source
    if cmd_args.verbose: print "Constructing radius array..."
    x = np.arange(nx) - float(nx-1)/2
    y = x[:,None]
    z = y[:,:,None]
    r = np.sqrt(x**2 + y**2 + z**2)

    # find point closest to center of grid that satisfies the wind zone condition
    icenter = r[iswind].argmin()

    ilabel = labeled_array[iswind][icenter] # find label of central feature
    print "Grid center total pressure = ", ptot[nz/2, ny/2, nx/2] 
    print "Central feature total pressure = ", ptot[iswind][icenter], " at distance from center of ", r[iswind][icenter]
else:
    ilabel = labeled_array[nz/2, ny/2, nx/2]

print "Central feature label = ", ilabel
if ilabel == 0: ilabel = 1

if cmd_args.verbose: print "Isolating central feature..."
iswind = labeled_array == ilabel

# Find volume of wind region
ncells = iswind.sum()
print "Volume fraction occupied by wind: ", float(ncells) / (nx*ny*nz)

if cmd_args.stats:
    wpw = pw.data[iswind]
    wpt = pr.data[iswind] + p.data[iswind]
    print "Wind ram pressure stats in wind-dominated volume:"
    a, b, c = np.min(wpw), float(np.median(wpw)), np.max(wpw)
    print "Min/Median/Max = {:.2e}/{:.2e}/{:.2e}".format(1.0*a, 1.0*b, 1.0*c)
    print "Nebula ram+thermal pressure stats in wind-dominated volume:"
    print "Min/Median/Max = {:.2e}/{:.2e}/{:.2e}".format(1.0*np.min(wpt), 1.0*np.median(wpt), 1.0*np.max(wpt))

PARSEC = 3.085677582e18
YEAR = 3.15576e7
MSUN = 1.989e33
KM = 1.0e5
dx = (4.0 / nx) * PARSEC
t = (10000.0*YEAR)*float(cmd_args.id.split("_")[-1])

windvolume = ncells*dx**3
# canonical wind parameters from GAHA2001 sec 5.2 
Mdot = (1.e-6*MSUN/YEAR) * 0.35
Vwind = (1000.0*KM) * 1.2
# Total wind energy with efficiency factor of 5/11 (Freyer et al 2003, eq 10)
Ewind = cmd_args.efrac * 0.5 * Mdot * Vwind**2 * t
# Wind pressure: (gamma - 1.0) E / V
Pwind = (2./3.)*Ewind/windvolume

print "Shocked wind pressure in energy-driven case would be: ", Pwind
if cmd_args.energy:
    print "Guess for pressure balance contour: ", cmd_args.pguess
    print "This would work if bubble energy efficiency were {:.4f} times the canonical value".format(cmd_args.pguess/Pwind)

##
## Note that this is not really realistic since if it were really
## momentum driven, then there would be radial shadowing
##

if not cmd_args.onlyinfo:
    if cmd_args.energy:
        Xid = "X{:02}ene".format(int(1.e11*cmd_args.pguess))
    else:
        Xid = "X{:02}mom".format(int(cmd_args.boost))

    # X-ray emissivity
    if cmd_args.verbose: print "Calculating x-ray emissivity..."
    em = iswind.astype(np.float32)    # 1.0 or 0.0 where iswind is True or False

    if cmd_args.verbose: print "Writing emissivity cube..."
    hdu = pyfits.PrimaryHDU(em)
    hdu.writeto(cmd_args.id + "e-{}.fits".format(Xid), clobber=True)

    if cmd_args.verbose: print "Removing hole from Ha emissivity cube..."
    ha, = pyfits.open(cmd_args.id + "e-Halpha.fits")
    ha.data *= 1.0 - em             # remove X-ray hole from H alpha
    ha.writeto(cmd_args.id + "e-Ha{}.fits".format(Xid[:-2]), clobber=True)
