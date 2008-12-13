""" image processing

$Header$

"""


#---------------------------------------------------------------------------
class AIT(object):
    """ Manage a full-sky image of a SkyProjection or SkyFunction, wrapping SkyImage
     """
    
    def __init__(self, skyfun, pixelsize=0.5, center=None, galactic=True, fitsfile='', proj='AIT', size=180, earth=False):
        """
        skyfun SkyProjection or SkyFunction object
        pixelsize [0.5] size, in degrees, of pixels
        galactic [True] galactic or equatorial coordinates
        fitsfile [''] if set, write the projection to a FITS file
        proj ['AIT'] could be 'CAR' for carree or 'ZEA': used by wcslib
        center [None] if default center at (0,0) in coord system
        size [180] make less for restricted size
        earth [False]  looking down at Earth

        """
        from skymaps import SkyImage, SkyDir
        from numpy import array, isnan, ma
        from pylab import normalize
        
        self.skyfun = skyfun
        self.galactic = galactic
        self.pixelsize = pixelsize
        self.size = size
        # set up, then create a SkyImage object to perform the projection to a grid
        if center is None:
            center = SkyDir(0,0, SkyDir.GALACTIC if galactic else SkyDir.EQUATORIAL)
        self.skyimage = SkyImage(center, fitsfile, pixelsize, size, 1, proj, galactic, earth)
        self.skyimage.fill(skyfun)
        
        # now extract stuff for the pylab image, creating a masked array to deal with the NaN values
        self.nx, self.ny = self.skyimage.naxis1(), self.skyimage.naxis2()
        self.image = array(self.skyimage.image()).reshape((self.ny, self.nx))
        self.mask = isnan(self.image)
        self.masked_image = ma.array( self.image, mask=self.mask)
        if not earth:
            self.extent = (180,-180, -90, 90) if size==180 else (size, -size, -size, size)
        else:
            self.extent = (-180,180, -90, 90) if size==180 else (-size, size, -size, size)
       
        self.vmin ,self.vmax = self.skyimage.minimum(), self.skyimage.maximum()

        # we want access to the projection object, to allow interactive display via pix2sph function
        self.proj = self.skyimage.projector()
        self.x = self.y = 100 # initial def

    def grid(self, fig=None, labels=True, symbol="0.5"):
	"""Draws gridlines and labels for map."""
        from pylab import figure
        from numpy import hstack
        def ait(l, b):
            " convert lon, lat to car "
            return self.proj.sph2pix(l, b)

        self.figure = figure(1, figsize = (10, 5), facecolor = 'w') if fig is None else fig
        self.axes = self.figure.add_subplot(111)

        from numpy import arange
        from pylab import plot, text, connect, show
        import matplotlib
        matplotlib.interactive(False)
        self.axes.set_autoscale_on(False)
        self.axes.set_xlim(0, 360/self.pixelsize)
        self.axes.set_ylim(0, 180/self.pixelsize)
        self.axes.set_axis_off()
        self.axes.set_aspect('equal')
        self.extent= (ait(180,0)[0],ait(180.001,0)[0], ait(0,-90)[1], ait(0,90)[1])
        bs = arange(-90, 91, 5)
        for l in hstack((arange(0, 360, 45),[180.01])):
            plot([ait(l,b)[0] for b in bs], [ait(l,b)[1] for b in bs], '--', color='gray')
            if labels:
                x,y = ait(l, 45) # Need to find right syntax for annotation position.
                text(x,y, '%3.0f'%l ,size='smaller')#, weight = 'bold')

        ls = hstack((arange(180, 0, -5), arange(355, 180,-5), [180.01]))
        for b in arange(-60, 61, 30):
            plot([ait(l,b)[0] for l in ls], [ait(l,b)[1] for l in ls], '--', color='gray')
            if labels:
                x,y = ait(-160, b)                                               
                text(x,y, '%+3.0f'%b, size='smaller')#, weight = 'bold')
        matplotlib.interactive(True)
        show()
                  
    def imshow(self,  title=None, scale='linear',  **kwargs):
        'run imshow'
        from numpy import ma
        import pylab
        # change defaults
        if 'origin'        not in kwargs: kwargs['origin']='lower'
        if 'interpolation' not in kwargs: kwargs['interpolation']='nearest'
        if 'extent'        not in kwargs: kwargs['extent']=self.extent
        
        if   scale=='linear':  pylab.imshow(self.masked_image,   **kwargs)
        elif scale=='log':     pylab.imshow(ma.log10(self.masked_image), **kwargs)
        else: raise Exception('bad scale: %s'%scale)
                                        
        pylab.colorbar(orientation='horizontal', shrink=1.0 if self.size==180 else 0.6)
        #pylab.axes().set_axis_off()
        if self.galactic:
            pylab.xlabel('glon'); pylab.ylabel('glat')
        else:
            pylab.xlabel('ra'); pylab.ylabel('dec')
        self.title(title)

        # for interactive formatting of the coordinates when hovering
        ##pylab.gca().format_coord = self.format_coord # replace the function on the fly!

    def pcolor(self,  title=None, scale='linear',  **kwargs):
        'run pcolor'
        from numpy import ma, array
        import pylab
        if self.galactic:
            xvalues=array([self.skydir(i,0).l() for i in range(self.nx+1)])
            yvalues=array([self.skydir(0,i).b() for i in range(self.ny+1)])
            pylab.xlabel('glon'); pylab.ylabel('glat')
        else:
             xvalues=array([self.skydir(i,0).ra() for i in range(self.nx+1)])
             yvalues=array([self.skydir(0,i).dec() for i in range(self.ny+1)])
             pylab.xlabel('ra'); pylab.ylabel('dec')

        if   scale=='linear':  pylab.pcolor(self.masked_image,   **kwargs)
        elif scale=='log':     pylab.pcolor(ma.log10(self.masked_image), **kwargs)
        else: raise Exception('bad scale: %s'%scale)
                                        
        pylab.colorbar(orientation='horizontal', shrink=1.0 if self.size==180 else 0.6)

        self.title(title)

    def axes(self, color='black',  **kwargs):
        ' overplot axis lines'
        import pylab
        pylab.axvline(0, color=color, **kwargs)
        pylab.axhline(0, color=color, **kwargs)
        pylab.axis(self.extent)
 
    def title(self, text=None, **kwargs):
        ' plot a title, default the name of the SkySpectrum'
        import pylab
        try:
            pylab.title( text if text else self.skyfun.name(), **kwargs)
        except AttributeError: #no name?
            pass

    def skydir(self, x, y):
        " from pixel coordinates to sky "
        from pointlike import SkyDir
        xpixel = (180-x)*float(self.nx)/360.
        ypixel = (y+90)*float(self.ny)/180.
        if self.proj.testpix2sph(xpixel,ypixel) !=0: return None #outside valid region
        sdir = SkyDir(x, y, self.proj)
        return sdir

    def format_coord(self, x, y):
        " replacement for Axes.format_coord"
        sdir = self.skydir(x,y)
        val  = self.skyfun(sdir)

        return 'ra,dec: (%7.2f,%6.2f); l,b: (%7.2f,%6.2f), value:%6.3g' %\
            ( sdir.ra(), sdir.dec(), sdir.l(), sdir.b(), val)
                
