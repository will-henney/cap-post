"""
Make a movie of a rotating view of a datacube at a fixed time
"""
import subprocess, StringIO, os, sys
import myfitsutils as F
from PIL import Image


class Camera(object):
    """
    A camera with (theta, phi) orientation
    """
    theta = 0.0
    phi = 0.0
    def set_angles(self, theta, phi):
	self.theta = theta
	self.phi = phi
	self.fold_angles()
    def set_steps(self, dth, dph):
	self.dth = dth
	self.dph = dph
    def __init__(self, theta=0.0, dth=1.0, phi=0.0, dph=1.0):
	# Initial angles
	self.set_angles(theta, phi)
	# Atomic steps for angles
	self.set_steps(dth, dph)
    def advance(self, incr=1):
	"Advance angles by increment"
	self.theta += incr*self.dth
	self.phi += incr*self.dph
	self.fold_angles()

    def fold_angles(self):
	"""Reduce angles to the canonical range [0:360] 

	In reality, theta should be in [0:180] but my first stab at
	that didn't work, so I am leaving it
	"""
# THIS DOESN'T WORK
# 	if self.theta > 180.0:
# 	    self.theta = 360.0 - self.theta
# 	    self.phi = self.phi + 180.0
# 	if self.theta < 0.0:
# 	    self.theta = -self.theta
# 	    self.phi = self.phi + 180.0
	self.theta = self.theta % 360.0
	self.phi = self.phi % 360.0
	    

class Movie(object):
    runid = None
    time = 100
    emtypes = ["N26584", "Halpha", "O35007"]
    camera = Camera()
    datadir = '.'
    boxsize = 1.0 		# z extent in parsecs
    brightmax = 1.e7		# typical maximum brightness
    brightmaxfunc = None	# if set, used to dynamically update brightmax

    # relative scaling of the emission bands
    bandscales = [0.28, 1.0, 0.2]
    gamma = 2.0			# default but can be overwritten by makeRGBmap

    makerotmap = "makerotmap" 	# may be swapped to makerotmapsmall

#     simtype = None
#     simtypes = ['Garrelt', 'Cap3D', 'Fabio']
#     denid = {'Garrelt': 'd', 'Fabio': 'dd'}
#     temid = {'Garrelt': 't', 'Fabio': 'te'}
    def __init__(self, runid, movieid, 
		 time=100, dtime=0, nframes=100, frame=1, 
		 verbose=0, force=0):
	self.runid = runid
	self.movieid = movieid
	self.time = time
	self.dtime = dtime
	self.nframes = nframes
	self.verbose = verbose
	self.force = force 	# Whether to regenerate all files
	self.frame = frame	# Allow to start from other than frame 1
	self.execdir = os.path.split(sys.argv[0])[0]
	self.imageprefix = "rgb-NHO-%s-%s" % (self.runid, self.movieid)
	self.imagelist = []
# 	for simtype in self.simtypes:
# 	    if os.getcwd().find(simtype) >= 0:
# 		self.simtype = simtype
# 	if self.simtype is None:
# 	    raise ValueError, "Cannot determine simulation type"

    def mapprefix(self, emtype):
	"""
	Returns the prefix string for a given map

	E.g., glob12-128xmap-rot+060-015-Halpha0060
	"""
	return "%smap-rot%+3.3i%+3.3i-%s%4.4i" % (
	    self.runid,
	    self.camera.theta,
	    self.camera.phi,
	    emtype,
	    self.time,
	    )

    def emcubefilename(self, emtype):
	"""
	Returns the file name of the given emissivity cube

	E.g., Bstar-e-Halpha0060
	"""
	return "%s-e-%s%4.4i.fits" % (
	    self.runid,
	    emtype,
	    self.time,
	    )

    def makeRGBimage(self, gamma=None):
	"""
	Make a single RGB emision image

	Returns filename of image
	"""
	imagename = self.imageprefix + "-%4.4i.png" % self.frame
	if self.force > 0 or (self.force > -1 and not os.path.isfile(imagename)):
	    if self.verbose > 0: print "Creating image file %s" % imagename
	    try:
		# first check that all em files exist
		for emtype in self.emtypes:
		    if not os.path.isfile(os.path.join(self.datadir, self.emcubefilename(emtype))):
			print "Emissivity cube does not exist: %s" % (self.emcubefilename(emtype))
			raise IOError
		# now make the maps
		rfits, gfits, bfits = [self.makefitsmap(emtype) 
				       for emtype in self.emtypes]
                # WJH 25 May 2010 - added ability to have brightmax be band-dependent
                # This is necessary if we are using brightmaxfunc and we want a different time dependence in the different bands
                try: 
                    # case where brightmax is a triplet
                    rmax, gmax, bmax = [scale*thisbrightmax for scale, thisbrightmax in zip(self.bandscales, self.brightmax)]
                except TypeError:
                    # case where brightmax is a scalar
                    rmax, gmax, bmax = [scale*self.brightmax for scale in self.bandscales]
		if gamma is None: gamma = self.gamma
		if self.verbose > 0: 
		    print "R/G/B max: %.2e/%.2e/%.2e, gamma: %.2f: " % (rmax, gmax, bmax, gamma)
		image = F.RGB3Image(rfits, gfits, bfits, 
				    0, 0, 0, rmax, gmax, bmax, gamma)
		image = image.transpose(Image.FLIP_TOP_BOTTOM)
		image.save(imagename)
	    except IOError:
		print "Could not generate file %s - skipping to next" % imagename
	else:
	    if self.verbose > 0: print "Image file %s already exists" % imagename
	return imagename


    def makeimages(self):
	framerange = range(self.frame, self.frame+self.nframes)
	self.imagelist = []
	for self.frame in framerange:
	    if self.brightmaxfunc: # make sure brightmax is up to date
		self.brightmax = self.brightmaxfunc(self.time)
	    self.imagelist.append(self.makeRGBimage())
	    self.camera.advance(1)
	    self.time += self.dtime

    def makemovie(self):
	if len(self.imagelist) != self.nframes:
	    self.makeimages()
	encode_exec = "mencoder"
	encode_args = "-ovc lavc -lavcopts vbitrate=5000:vcodec=wmv2 -mf type=png:fps=15"
	file_args = r"-o %s.avi mf://%s\*.png" % (self.imageprefix, self.imageprefix)
	cmd = "%s %s %s" % (encode_exec, encode_args, file_args)
	subprocess.Popen(cmd, shell=True).wait()
	print "Written %s.avi" % (self.imageprefix)

    def makefitsmap(self, emtype):
	"""
	Make a single FITS emission map

	Returns fits file name
	"""
	fitsfilename = self.mapprefix(emtype) + ".fits"
	if self.force > 1 or (self.force > -1 and 
			      not os.path.isfile(os.path.join(self.datadir, fitsfilename))):
	    # try and make the fits file if it is not there, or if force flag is set
	    if self.verbose > 0: print "Creating FITS file %s" % fitsfilename
	    if self.verbose > 1:
		output = None 		# print out all shell output
	    else:
		output = subprocess.PIPE # we never read from this
	    cmd = os.path.join(self.execdir, self.makerotmap)
	    p = subprocess.Popen([cmd], stdin=subprocess.PIPE, stdout=output, 
				 cwd=self.datadir)
	    p.stdin.writelines(
		[str(x) + '\n' for x in [
			self.runid,
			self.time,
			emtype,
			self.boxsize,
			self.camera.theta,
			self.camera.phi,
			]]
		)
	    p.stdin.close()
	    p.wait()
	else:
	    if self.verbose > 0: print "FITS file %s already exists" % fitsfilename

	return fitsfilename


if __name__ == '__main__': 
    # Test the machinery

    try: 
	movie = Movie("glob12-128x", "testmovie100", 
		      time=100, nframes=70, verbose=2, force=0)
	movie.brightmax = 5.e7
	movie.camera.set_steps(5.0, 5.0)
	movie.makemovie()
    except:
	itime = 50
	movie = Movie("06102007_x", "testmovie%i" % itime, 
		      time=itime, nframes=2, verbose=2, force=2)
	movie.brightmax = 2.6e7*(1.0+(float(itime)/30)**2)/(1.0+(float(itime)/50)**3)
	movie.boxsize = 4.0
	movie.makerotmap = "makerotmapsmall"
	movie.camera.set_steps(5.0, 5.0)
	movie.makeimages()
