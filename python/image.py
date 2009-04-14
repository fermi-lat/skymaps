""" image processing

$Header$

"""

import pylab
from numpy import arange, hstack
from skymaps import SkyImage, SkyDir, double2, SkyProj

def draw_grid(ait, labels=True, color='gray', pixelsize=0.5, textsize=8):
        label_offset = 5/pixelsize
        my_axes = pylab.axes() #creates figure and axes if not set

        pylab.matplotlib.interactive(False)
        my_axes.set_autoscale_on(False)
        my_axes.set_xlim(0, 360/pixelsize)
        my_axes.set_ylim(0, 180/pixelsize)
        my_axes.set_axis_off()
        my_axes.set_aspect('equal')
        #? extent= (ait(180,0)[0],ait(180.001,0)[0], ait(0,-90)[1], ait(0,90)[1])


        bs = arange(-90, 91, 5)
        for l in hstack((arange(0, 360, 45),[180.01])):
            lstyle = '-' if int(l)==180 or int(l)==0 else '--' 
            pylab.plot([ait(l,b)[0] for b in bs], [ait(l,b)[1] for b in bs], lstyle, color=color)
            if labels:
                x,y = ait(l, 45) 
                pylab.text(x,y, '%3.0f'%l ,size=textsize, ha='center')

        ls = hstack((arange(180, 0, -5), arange(355, 180,-5), [180.01]))
        for b in arange(-60, 61, 30):
            lstyle = '-' if int(b)==0 else '--'
            pylab.plot([ait(l,b)[0] for l in ls], [ait(l,b)[1] for l in ls], lstyle, color=color)
            if labels:
                x,y = ait(180.1, b)                                               
                pylab.text(x+label_offset,y+b/60*label_offset, '%+3.0f'%b, size=textsize, ha='center',va='center')
        if labels:
            for b in [90,-90]:
                x,y = ait(0,b)
                pylab.text(x,y+b/90*label_offset,'%+3.0f'%b, size=textsize, ha='center',va='center') 


class AIT_grid():

    def __init__(self, labels=True, color='gray', pixelsize=0.5, textsize=8):
	"""Draws gridlines and labels for map.
        
        """

        self.pixelsize = pixelsize

        xsize,ysize = 325,162
        crpix = double2(xsize/pixelsize/2., ysize/pixelsize/2.)
        crval = double2(0,0)
        cdelt = double2(-pixelsize, pixelsize)
        self.proj = SkyProj('AIT', crpix, crval, cdelt, 0, True)
        
        self.axes = pylab.axes() #creates figure and axes if not set

        pylab.matplotlib.interactive(False)
        self.axes.set_autoscale_on(False)
        self.axes.set_xlim(0, 360/self.pixelsize)
        self.axes.set_ylim(0, 180/self.pixelsize)
        self.axes.set_axis_off()
        self.axes.set_aspect('equal')
        self.extent= (self.ait(180,0)[0],self.ait(180.001,0)[0], self.ait(0,-90)[1], self.ait(0,90)[1])
        label_offset = 5/self.pixelsize
        bs = arange(-90, 91, 5)
        for l in hstack((arange(0, 360, 45),[180.01])):
            lstyle = '-' if int(l)==180 or int(l)==0 else '--' 
            pylab.plot([self.ait(l,b)[0] for b in bs], [self.ait(l,b)[1] for b in bs], lstyle, color=color)
            if labels:
                x,y = self.ait(l, 45) 
                pylab.text(x,y, '%3.0f'%l ,size=textsize, ha='center')

        ls = hstack((arange(180, 0, -5), arange(355, 180,-5), [180.01]))
        for b in arange(-60, 61, 30):
            lstyle = '-' if int(b)==0 else '--'
            pylab.plot([self.ait(l,b)[0] for l in ls], [self.ait(l,b)[1] for l in ls], lstyle, color=color)
            if labels:
                x,y = self.ait(180.1, b)                                               
                pylab.text(x+label_offset,y+b/60*label_offset, '%+3.0f'%b, size=textsize, ha='center',va='center')#, weight = 'bold')
        if labels:
            for b in [90,-90]:
                x,y = self.ait(0,b)
                pylab.text(x,y+b/90*label_offset,'%+3.0f'%b, size=textsize, ha='center',va='center') 
        pylab.matplotlib.interactive(True)
        pylab.show()

    def ait(self, l, b):
        " convert lon, lat to car "
        return self.proj.sph2pix(l, b)

    def plot(self, sources, symbol='+', text=None, fontsize=8, **kwargs):
        """ plot symbols at points
        text: optional text strings (same lenght as soruces)
        """
        X=[]
        Y=[]
        for i,s in enumerate(sources):
            x,y = self.ait(s.l(),s.b())
            X.append(x)
            Y.append(y)
            if text is not None:
                pylab.text(x,y,text[i],fontsize=fontsize)
        pylab.plot(X,Y, symbol,  **kwargs)

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
        from numpy import array, isnan, ma
        
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

    def grid(self, fig=None, labels=True, color='gray'):
	"""Draws gridlines and labels for map."""
        self.axes = pylab.axes() #creates figure and axes if not set

        pylab.matplotlib.interactive(False)
        self.axes.set_autoscale_on(False)
        self.axes.set_xlim(0, 360/self.pixelsize)
        self.axes.set_ylim(0, 180/self.pixelsize)
        self.axes.set_axis_off()
        self.axes.set_aspect('equal')
        ait = self.proj.sph2pix
        self.extent= (ait(180,0)[0],ait(180.001,0)[0], ait(0,-90)[1], ait(0,90)[1])

        draw_grid(self.proj.sph2pix, labels=labels, color=color, pixelsize=self.pixelsize) 
        pylab.show()

    def plot(self, sources, symbol='+',  **kwargs):
        " plot symbols at points"
        
        X=[]
        Y=[]
        for s in sources:
            x,y = self.proj.sph2pix(s.l(),s.b())
            X.append(x)
            Y.append(y)
        pylab.plot(X,Y, symbol,  **kwargs)

    def on_move(self, event):
        """Reports mouse's position in galactic coordinates."""
        from numpy import fabs
        if event.xdata == None or event.ydata == None:
            pass 
        else:
            try:
                coords = self.proj.pix2sph(event.xdata, event.ydata)
                self.poslabel.set_text("long=%1.2f\n  lat=%1.2f" %(coords[0],coords[1]))
            except:
                self.poslabel.set_text("")
	self.figure.canvas.draw()
                  
    def imshow(self,  title=None, scale='linear', factor=1.0, **kwargs):
        'run imshow'
        from numpy import ma
        # change defaults
        if 'origin'        not in kwargs: kwargs['origin']='lower'
        if 'interpolation' not in kwargs: kwargs['interpolation']='nearest'
        if 'extent'        not in kwargs: kwargs['extent']=self.extent
        
        if self.size==180: pylab.axes().set_axis_off()
        if   scale=='linear':  pylab.imshow(self.masked_image*factor,   **kwargs)
        elif scale=='log':     pylab.imshow(ma.log10(self.masked_image), **kwargs)
        else: raise Exception('bad scale: %s, expect either "linear" or "log"'%scale)
                                        
        self.colorbar =pylab.colorbar(orientation='horizontal', shrink=1.0 if self.size==180 else 0.6)
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
            pylab.title( text if text is not None else self.skyfun.name(), **kwargs)
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
                
