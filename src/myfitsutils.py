"""
Classes for reading fits files into image arrays

These start off pretty general, then get specific to our
velocity map data

WJH 01 Feb 2005

Changes:

18 Apr 2005: Add classes for doing HSV color images using the moments

"""

import Image
import pyfits
import numpy.numarray as numarray


class MyImage(Image.Image):
    """
    Helper class to create empty image of desired shape
    
    This was a bit tricky since the Image class does not
    seem to be designed for easy sub-classing
    """

    def __init__(self, mode='L', size=None):
        Image.Image.__init__(self)
        self.mode = mode
        self.size = size
        self.im = Image.core.new(self.mode, self.size)

    def show(self):
        "Correct the orientation before showing"
        Image.Image.show(self.transpose(Image.FLIP_TOP_BOTTOM))


class FitsData:
    """
    A class that contains the data from a (possible subimage of a) fits files
    """
    def __init__(self, fitsfile, bb=[None, None, None, None], 
		 icut=None, cutaxis='x'):
        self.bb = bb
	if not icut is None: 
	    # we have a data cube, from which we will extract one plane
	    if cutaxis=='x':
		fitsfulldata = pyfits.open(fitsfile)[0].data[:,:,icut]
	    elif cutaxis=='y':
		fitsfulldata = pyfits.open(fitsfile)[0].data[:,icut,:]
	    elif cutaxis=='z':
		fitsfulldata = pyfits.open(fitsfile)[0].data[icut,:,:]
	else:
	    # note that fitsdata array has y as the first dimension
	    fitsfulldata = pyfits.open(fitsfile)[0].data
        # crop bounding box on top-right
	# if self.bb[3] and self.bb[3] < numarray.shape(fitsfulldata)[0]:
	#     self.bb[3] = numarray.shape(fitsfulldata)[0]
	# if self.bb[2] and self.bb[2] < numarray.shape(fitsfulldata)[1]:
	#     self.bb[2] = numarray.shape(fitsfulldata)[1]
        self.fitsdata = fitsfulldata[self.bb[1]:self.bb[3],self.bb[0]:self.bb[2]]
        self.shape = numarray.shape(self.fitsdata)
        self.size = (self.shape[1],self.shape[0])
        
class FitsImage(MyImage, FitsData):
    """
    A fits file with associated grayscale image, with gamma scaling
    between fmin and fmax.
    """

    def __init__(self, fitsfile, fmin, fmax, gamma=1.0, bb=[0,0,-1,-1], 
		 icut=None, cutaxis=1, takelog=0):
        FitsData.__init__(self, fitsfile, bb, icut, cutaxis)
        MyImage.__init__(self, 'L', self.size)
        # 06 Dec 2010: add new attributes: min, max to allow outside inspection
        # 05 Dec 2010: If fmin or fmax is None, then the respective limit is found from the data
        self.min = fmin if not fmin is None else self.fitsdata.min()
        self.max = fmax if not fmax is None else self.fitsdata.max()
	if takelog:
	    self.flatnormaldata = (
		numarray.ravel(numarray.log10(self.fitsdata)) 
		- numarray.log10(self.min) ) / ( 
		numarray.log10(self.max) - numarray.log10(self.min))
	else:
	    self.flatnormaldata = (numarray.ravel(self.fitsdata)-self.min)/(self.max-self.min)
	# 09 Sep 2005: How on earth did this ever work before? It is necessary
	# to clip fitsdata to [fmin,fmax] before applying gamma correction, or
	# else we get sqrt of negative numbers....
	# Aha, there is a smart way to do this!
	self.flatnormaldata = numarray.clip(self.flatnormaldata, 0.0, 1.0)

	if gamma != 1.0:
	    self.flatnormaldata = self.flatnormaldata**(1./gamma)
	self.putdata(self.flatnormaldata, scale=255.0)


class RGB3Image(FitsImage):
    """
    A color image in which each of the 3 RGB channels comes
    from a different fits file.
    """

    def __init__(self,
                 redfile, greenfile, bluefile,
                 redmin, greenmin, bluemin,
                 redmax, greenmax, bluemax,
                 gamma=1.0, bb=[0,0,-1,-1]):
        self.red = FitsImage(redfile, redmin, redmax, gamma, bb)
        self.green = FitsImage(greenfile, greenmin, greenmax, gamma, bb)
        self.blue = FitsImage(bluefile, bluemin, bluemax, gamma, bb)
        MyImage.__init__(self, 'RGB', self.red.size)
        self.im.putband(self.red.im, 0)
        self.im.putband(self.green.im, 1)
        self.im.putband(self.blue.im, 2)

class HSV3Image(FitsData, MyImage):
    """
    A color image in which each of the 3 HSV channels comes
    from a different fits file.
    """

    def __init__(self,
                 huefile, saturationfile, valuefile,
                 huemin, saturationmin, valuemin,
                 huemax, saturationmax, valuemax,
                 gamma=1.0, bb=[0,0,-1,-1]):
        self.hdata = FitsData(huefile, bb)
        self.sdata = FitsData(saturationfile, bb)
        self.vdata = FitsData(valuefile, bb)

        # use the original data to calculate the RGB values
        hue = (numarray.ravel(self.hue.hdata)-huemin)/(huemax-huemin)
        hue = (numarray.ravel(self.hue.hdata)-huemin)/(huemax-huemin)
        hue = (numarray.ravel(self.hue.hdata)-huemin)/(huemax-huemin)

        self.red = MyImage.__init__(self, 'L', self.hdata.size)
        self.green = MyImage.__init__(self, 'L', self.hdata.size)
        self.blue = MyImage.__init__(self, 'L', self.hdata.size)
        
        MyImage.__init__(self, 'RGB', self.hdata.size)
        self.im.putband(self.hue.im, 0)
        self.im.putband(self.saturation.im, 1)
        self.im.putband(self.value.im, 2)


class RatioImage(MyImage, FitsData):
    """
    A grayscale image of the ratio of two fits files
    """

    # optional function to use instead of the simple ratio
    ratiofunc = None
    
    def __init__(self, file1, file2, fmin, fmax, gamma=1.0, bb=[0,0,-1,-1], minquality=1.0):
        self.numerator = FitsData(file1, bb)
        self.denominator = FitsData(file2, bb)
        self.size = self.numerator.size
        if self.size != self.denominator.size:
            raise IndexError, "Incompatible sizes in RatioImage!"
        MyImage.__init__(self, 'L', self.size)
        if self.ratiofunc == None:
            self.fitsdata = self.numerator.fitsdata / self.denominator.fitsdata
        else:
            self.fitsdata = self.ratiofunc( self.numerator.fitsdata, self.denominator.fitsdata )
	self.flatnormaldata = (numarray.ravel(self.fitsdata)-fmin)/(fmax-fmin)
	self.flatnormaldata = numarray.clip(self.flatnormaldata, 0.0, 1.0)
	if gamma != 1.0:
	    self.flatnormaldata = self.flatnormaldata**(1./gamma)
	self.putdata(self.flatnormaldata, scale=255.0)
        self.quality = numarray.where( self.denominator.fitsdata > minquality, 1, 0 )
        self.quality = numarray.where( self.numerator.fitsdata > minquality, self.quality, 0 )
        self.alpha = Image.new('L', self.size)
        self.alpha.putdata(numarray.ravel(self.quality), scale=255 )
        self.im = Image.composite( self, Image.new('L', self.size, 255 ), self.alpha ).im
        
class IsoVel:
    """
    Parent class for isovelocity maps
    """

    ##
    ## Class variables
    ##
    
    # channel width (km/s)
    deltav = 4.0

    # smoothing suffix
    suffix = '01'

    # middle part of filename
    midfix = 'wisomom-sum'

    # root directory
    root = './'

    # pixel size in arcsec
    pixel = 1.24

    # position of theta 1c in pixels
    starposdict = {
        # new values for the smooth2d images
        'siil': (172, 377),
        'siis': (172, 377),
        'oi': (172, 377),
        'siii': (172, 377),
##        'siil': (82, 189), # found from eden-tere.f90
##        'siis': (82, 189),
##        'oi': (178, 373),   # found from orionvimages.py
##        'siii': (178, 373),
        }

    # bounding box of view area in arcsec wrt star
    viewbox = (-120, -200, 80, 120)

    # possible shift of image wrt bottom left corner of window
    # (in arcsec)
    xshift = 0
    yshift = 0

    ##
    ## More class variables that should be overwritten
    ## by subclasses
    ##

    # emission line ID 
    emline = 'siil'

    def __init__(self):
        pass

    def getvrange(self, v):
        "Returns velocity range as string given central velocity"
        return "%+3.3i%+3.3i" % (v-0.5*self.deltav, v+0.5*self.deltav)

    def getfilename(self, v):
        "Returns filename corresponding to central velocity v"
        return self.root + self.emline + '_' + self.getvrange(v) \
               + '.' + self.midfix + '-' + self.suffix + '.fits'

    def getpixelcoords(self, x, y):
        "Returns pixel coordinates for a given offset in arcsec"
        try:
            px0, py0 = self.starposdict[self.emline]
        except KeyError:
            print "Warning: no star position found for %s" % self.emline
            px0, py0 = (162, 377)
        px = px0 + int(x/self.pixel)
        py = py0 + int(y/self.pixel)
        return px, py

    def getarcseccoords(self, px, py):
        try:
            px0, py0 = self.starposdict[self.emline]
        except KeyError:
            print "Warning: no star position found for %s" % self.emline
            px0, py0 = (0, 0)
        x = (px-px0)*self.pixel
        y = (py-py0)*self.pixel
        return x, y
    
    def getbb(self):
        "Returns pixel bounding box"
        if self.viewbox == None:
            bb = [0,0,-1,-1]
        else:
            bbllx, bblly = self.getpixelcoords(self.viewbox[0], self.viewbox[1])
            bburx, bbury = self.getpixelcoords(self.viewbox[2], self.viewbox[3])
            # we can only use these values if they are within the size of the fits image
            # we can check they aren't negative
            # also, if they are negative then we would need to shift the image
            if bbllx < 0:
                self.xshift = -bbllx*self.pixel
                bbllx = 0
            if bblly < 0:
                self.yshift = -bblly*self.pixel
                bblly = 0
            # but we can't check that they aren't > fits dimension
            # However, this seems to not matter too much
            bb = [bbllx, bblly, bburx, bbury]
            print bb
        return bb


class RatioMap(IsoVel, RatioImage):
    """
    A grayscale ratio map
    """

    def getfilename(self, emline, v):
        "Returns filename corresponding to central velocity v"
        return self.root + emline + '_' + self.getvrange(v) \
               + '.' + self.midfix + '-' + self.suffix + '.fits'

    def __init__(self, emlines, v, min, max, gamma=1.0, viewbox=None, minquality=1.0, ratiofunc=None):
        self.emlines = emlines
        self.emline = emlines[0]
        self.deltav = 3*IsoVel.deltav
        self.filenames = [ self.getfilename( emline, v) for emline in emlines ]
        if viewbox == None:
            self.viewbox = IsoVel.viewbox
        else:
            self.viewbox = viewbox
        self.bb = self.getbb()
        self.ratiofunc = ratiofunc
        RatioImage.__init__(self, self.filenames[0], self.filenames[1], min, max, gamma, self.bb, minquality)
        
class RGB3Map(RGB3Image, IsoVel):
    """
    A 3-color isovelocity map
    """

    def __init__(self, emline, v, maxlist, gamma=1.0, viewbox=None):
        self.emline = emline
        self.vels = [ v-self.deltav, v, v+self.deltav ]
        self.maxlist = maxlist
        self.filenames = [ self.getfilename(v) for v in self.vels ]
        if viewbox == None:
            self.viewbox = IsoVel.viewbox
        else:
            self.viewbox = viewbox
        self.bb = self.getbb()
        # most positive velocity first, since that is red
        RGB3Image.__init__(self,
                           self.filenames[2], self.filenames[1], self.filenames[0],
                           0, 0, 0, 
                           maxlist[2], maxlist[1], maxlist[0],
                           gamma, self.bb)
        
class RGB3MultiMap(RGB3Image, IsoVel):
    """
    A 3-color map in 3 different lines
    """

    def __init__(self, emlines, v, maxlist, gamma=1.0, viewbox=None):
        self.maxlist = maxlist
	self.emlines = emlines
	self.filenames = []
	for emline in self.emlines:
	    self.emline = emline
	    self.filenames.append(self.getfilename(v))
        if viewbox == None:
            self.viewbox = IsoVel.viewbox
        else:
            self.viewbox = viewbox
        self.bb = self.getbb()
        # most positive velocity first, since that is red
        RGB3Image.__init__(self,
                           self.filenames[0], self.filenames[1], self.filenames[2],
                           0, 0, 0, 
                           maxlist[0], maxlist[1], maxlist[2],
                           gamma, self.bb)
        
