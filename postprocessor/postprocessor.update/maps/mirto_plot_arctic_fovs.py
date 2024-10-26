from __future__ import division

import georaster
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import glob
import matplotlib.pyplot as plt
#import pylab


from matplotlib.patches import Ellipse
from mpl_toolkits.basemap import Basemap
from osgeo import gdal, osr
from netCDF4 import Dataset


from matplotlib.patches import Polygon
from mpl_toolkits.basemap import pyproj

import sys
sys.path.insert(0,'/home/amethyst/amethyst/postprocessor')

# Fix Matplotlib KeyError with ProjLib
#os.environ['PROJ_LIB'] = "/home/oper/anaconda3/envs/gdal/bin/proj"

class Basemap(Basemap):
    def ellipse(self, x0, y0, a, b, n, ax=None, **kwargs):
        """
        Draws a polygon centered at ``x0, y0``. The polygon approximates an
        ellipse on the surface of the Earth with semi-major-axis ``a`` and
        semi-minor axis ``b`` degrees longitude and latitude, made up of
        ``n`` vertices.

        For a description of the properties of ellipsis, please refer to [1].

        The polygon is based upon code written do plot Tissot's indicatrix
        found on the matplotlib mailing list at [2].

        Extra keyword ``ax`` can be used to override the default axis instance.

        Other \**kwargs passed on to matplotlib.patches.Polygon

        RETURNS
            poly : a maptplotlib.patches.Polygon object.

        REFERENCES
            [1] : http://en.wikipedia.org/wiki/Ellipse

        """
        ax = kwargs.pop('ax', None) or self._check_ax()
        g = pyproj.Geod(a=self.rmajor, b=self.rminor)
        # Gets forward and back azimuths, plus distances between initial
        # points (x0, y0)
        azf, azb, dist = g.inv([x0, x0], [y0, y0], [x0+a, x0], [y0, y0+b])
        tsid = dist[0] * dist[1] # a * b

        # Initializes list of segments, calculates \del azimuth, and goes on
        # for every vertex
        seg = [self(x0+a, y0)]
        AZ = numpy.linspace(azf[0], 360. + azf[0], n)
        for i, az in enumerate(AZ):
            # Skips segments along equator (Geod can't handle equatorial arcs).
            if numpy.allclose(0., y0) and (numpy.allclose(90., az) or
                numpy.allclose(270., az)):
                continue

            # In polar coordinates, with the origin at the center of the
            # ellipse and with the angular coordinate ``az`` measured from the
            # major axis, the ellipse's equation  is [1]:
            #
            #                           a * b
            # r(az) = ------------------------------------------
            #         ((b * cos(az))**2 + (a * sin(az))**2)**0.5
            #
            # Azymuth angle in radial coordinates and corrected for reference
            # angle.
            azr = 2. * numpy.pi / 360. * (az + 90.)
            A = dist[0] * numpy.sin(azr)
            B = dist[1] * numpy.cos(azr)
            r = tsid / (B**2. + A**2.)**0.5
            lon, lat, azb = g.fwd(x0, y0, az, r)
            x, y = self(lon, lat)

            # Add segment if it is in the map projection region.
            if x < 1e20 and y < 1e20:
                seg.append((x, y))

        poly = Polygon(seg, **kwargs)
        ax.add_patch(poly)

        # Set axes limits to fit map region.
        self.set_axes_limits(ax=ax)

        return poly

def prepare_parser():
    from argparse import ArgumentParser

    v_levels = ['debug', 'info', 'warning']

    parser = ArgumentParser()
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='The level of verbosity of the software')
    parser.add_argument('--instrument','-i', default='cris', type=str,
                        help='Instrument, iasi or [cris]')
    parser.add_argument('--overpass','-o',default = None,type=str,
                        help='Overpass date and time YYYYMMDD_HHMMSS')

    return parser.parse_args()

def convertXY(xy_source, inproj, outproj):
    # function to convert coordinates
    
    shape = xy_source[0,:,:].shape
    size = xy_source[0,:,:].size

    # the ct object takes and returns pairs of x,y, not 2d grids
    # so the the grid needs to be reshaped (flattened) and back.
    ct = osr.CoordinateTransformation(inproj, outproj)
    xy_target = np.array(ct.TransformPoints(xy_source.reshape(2, size).T))

    xx = xy_target[:,0].reshape(shape)
    yy = xy_target[:,1].reshape(shape)
    
    return xx, yy

def GetExtent(ds):
    """ Return list of corner coordinates from a gdal Dataset """
    xmin, xpixel,  ymax, ypixel = ds.GetGeoTransform()
    width, height = ds.RasterXSize, ds.RasterYSize
    xmax = xmin + width * xpixel
    ymin = ymax + height * ypixel
    return (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)

def ReprojectCoords(coords,src_srs,tgt_srs):
    """ Reproject a list of x,y coordinates. """
    trans_coords=[]
    transform = osr.CoordinateTransformation( src_srs, tgt_srs)
    for x,y in coords:
        x,y,z = transform.TransformPoint(x,y)
        trans_coords.append([x,y])
    return trans_coords


def main():

   # Parse Arguments
   argv = prepare_parser()

   instr = argv.instrument
   overpass = argv.overpass

   print('overpass: {}'.format(overpass))

   L1DIR = '/galaxy/data/arctic/' + instr + '/'

   #Get available overpass times
   if os.path.isdir(L1DIR):
       sat_ov_times=os.listdir(L1DIR)
   else:
       print('Directory {} does not exist. Exiting'.format(L1DIR))
       exit()

   #Sirt overpasses starting from the most recent one
   sat_ov_times.sort(reverse=True)

   if not overpass:
       sat_pass_date = sat_ov_times[0]
       print('Last dir with L1 files is {}.'.format(sat_pass_date))
   else:
       sat_pass_date = overpass
       print('Selected DB dir with L0 files is {}.'.format(sat_pass_date))

   print('sat_pass_date: {}'.format(sat_pass_date))

   #filestr = L1DIR + sat_pass_date + '/npp_viirs_?.tif'
   #print('filestr: {}'.format(filestr))

   if instr == "cris":
      tiffile = ''.join([os.path.basename(x) for x in glob.glob(L1DIR + sat_pass_date + '/fre*.tif')]) 
   else:
      print("Option nont implemented yet")
      exit

   print('tiffile: {}'.format(tiffile))

   resfile = L1DIR + overpass + '/mirto/results.nc'
   fovfile = L1DIR + overpass + '/mirto/fov.nc'
   imgfile = L1DIR + overpass + '/'  + tiffile


   #Paolo basicfile = '20200918_000000_BASIC_FULL_XX.nc'
   #Paolo trfile    = basicfile.replace('BASIC','TR')
   
   res   = Dataset(resfile,'r')
   fov   = Dataset(fovfile,'r')
   #Paolo tr    = Dataset(trfile,'r')
   #Paolo basic = Dataset(basicfile,'r')
    
   # Get lat lon extention from the geoTIFF
   ds=gdal.Open(imgfile)
   gt = ds.GetGeoTransform()
   proj = ds.GetProjection()

   print('Geotif projection: {}'.format(proj))
   print('Geotif geotransform: {}'.format(gt))

   xres = gt[1]
   yres = gt[5]

   # get the edge coordinates and add half the resolution 
   # to go to center coordinates
   xmin = gt[0] + xres * 0.5
   xmax = gt[0] + (xres * ds.RasterXSize) - xres * 0.5
   ymin = gt[3] + (yres * ds.RasterYSize) + yres * 0.5
   ymax = gt[3] - yres * 0.5
   
   # create a grid of xy coordinates in the original projection
   xy_source = np.mgrid[ymax+yres:ymin:yres, xmin:xmax+xres:xres]

   #ext=GetExtent(ds)
   #sr.ImportFromEPSG(9810)
   #ds.SetProjection(sr.ExportToWkt())

   inproj=osr.SpatialReference()
   inproj.ImportFromWkt(ds.GetProjection())
   #tgt_srs = src_srs.CloneGeogCS()
   
   #geo_ext=ReprojectCoords(ext, src_srs, tgt_srs)

   #print(geo_ext)
   print('inproj: {}'.format(inproj))

   #minx, miny = geo_ext[-1]
   #maxx, maxy = geo_ext[1]
   ext = (xmin-xmin,xmax-xmin,ymin-ymin,ymax-ymin)


   #print('min_latlong: {}, max_latlong: {}'.format(min_latlong,max_latlong))
   print('minx {}, miny {}, maxx {}, maxy {}'.format(xmin, ymin, xmax, ymax))
   
   # Read geoTIFF
   image = georaster.MultiBandRaster(imgfile,bands='all',load_data=True)
   print('Extent: {}'.format(image.extent))

   #pylab.close('all')
   #pylab.ion()

   # set Basemap with slightly larger extents
   # set resolution at intermediate level "i"
   #m = Basemap(width=5400000., height=5400000., projection='stere',boundinglat=75,lat_0=90,lon_0=0,resolution='h')
   #m = Basemap(projection='npstere',boundinglat=75,lat_0=90,lon_0=0,resolution='h')
   #m = Basemap(width=5400000.,height=5400000.,projection='stere',lat_0=90,lat_ts=71,lon_0=-180,resolution='h',epsg=3995)
   m = Basemap(width=5400000.,height=5400000.,projection='stere',lat_0=90,lat_ts=71,lon_0=0,resolution='h')

   outproj=osr.SpatialReference()
   #outproj.ImportFromEPSG(3995)
   outproj.ImportFromProj4(m.proj4string)

   print('outproj: {}'.format(outproj))

   # Convert from source projection to basemap projection
   #xx, yy = convertXY(xy_source, inproj, outproj)


   m.drawcoastlines(color="black",linewidth=.1)
   #m.fillcontinents(color='brown')
   #m.drawmapboundary(fill_color='aqua')
   
   # Plot RGB data (normalized to 0...1)
   #plt.imshow(image.r/255,extent=ext,alpha=1.0,interpolation='nearest')
   m.imshow(image.r/255,extent=image.extent,alpha=1.0,interpolation='nearest',origin='upper')
   
   #Paolo tr_lon = tr['longitude'][:]
   #Paolo tr_lat = tr['latitude'][:]
   #Paolo basic_lon = basic['longitude'][:]
   #Paolo basic_lat = basic['latitude'][:]
   lon = fov['Longitude'][:] 
   lat = fov['Latitude'][:]
   
   # Read rh
   from atmos.mirto_atmos_tools import mr2rh
   rh = mr2rh(res['p'][:],res['t'][:],res['q'][:],273.15)[0]
   
   rh_max = np.amax(rh,axis=1)
   good_indx = np.where( (res['d2'][:] < 5) & (rh_max < 100))[0]
   print('Number of good retrievals: {}'.format(len(good_indx)))
       
   
   x,y =lon,lat
   good_x,good_y = m(x[good_indx],y[good_indx])
   all_x, all_y  = m(x,y)
   
   # Plot all FOVS
   #plt.scatter(all_x, all_y, color='red',marker='x',s=8,label='BAD')
   
   # Plot good FOVS
   #plt.scatter(good_x, good_y, color='lightgreen',marker='x',s=.1,label='GOOD')
   
   #Paolo tr_x, tr_y = tr_lon,tr_lat
   #Paolo tr_x, tr_y = m(tr_x,tr_y)
   
   # Plot TR
   #plt.scatter(tr_x, tr_y, color='orange',marker='o',s=8,label='TR')
   #Paolo basic_x, basic_y = basic_lon,basic_lat
   #Paolo basic_x, basic_y = m(basic_x,basic_y)
   
   # Plot BASIC
   #plt.scatter(tr_x, tr_y, color='tab:orange',marker='o',s=8)
   
   
   # Draw FOVs ellipses #################################################
   
   sat_ape=[0.963, 0.897] # degrees
   sat_alt = 824; # km
   earth_radius = 6371 # km
   
   from obspy.geodetics.base import kilometers2degrees
   
   tand = lambda x: np.tan(np.radians(x))
   
   for fov_index, xlon in enumerate(lon):
   
       # Compute axis and angle
       L1=( tand(sat_ape[0]+fov['FOV_angle'][fov_index].data) - tand(fov['FOV_angle'][fov_index].data))*sat_alt
       L2= tand(sat_ape[1])*sat_alt
       
       # Convert 2 degrees
       L1 = 0.5*kilometers2degrees(L1,radius=earth_radius)
       L2 = 0.5*kilometers2degrees(L2,radius=earth_radius)
       
       angle = -90 - fov['Satellite_azimuth_angle'][fov_index]

       #ax = pylab.gca()
       ax = plt.gca()
   
       # Draw the patch
       if fov_index in good_indx:
           poly = m.ellipse(xlon, lat[fov_index], L1, L2, 300, edgecolor='green', facecolor='green', zorder=10,
            alpha=0.5, linewidth=.1)
           #e = Ellipse((xlon,lat[fov_index]),L1,L2,angle,facecolor='red',edgecolor='red',lw=.3)
           #e.set_fill(True)
           #print('xlon {},lat[fov_index] {},L1 {},L2 {}'.format(xlon,lat[fov_index],L1,L2))
           #plt.gca().add_patch(e)
       else:
           poly = m.ellipse(xlon, lat[fov_index], L1, L2, 300,  edgecolor='red', facecolor='red', zorder=10, alpha=0.5, linewidth=.1)
           #e = Ellipse((xlon,lat[fov_index]),L1,L2,angle,edgecolor='red',facecolor='none',lw=.3)
           #e.set_fill(False)
           #plt.gca().add_artist(e)
   
   ##########################################################################
   
   # Plot legend
   #plt.legend(loc='lower right')

   #plt.show()
   plt.savefig(L1DIR + overpass + '/' + overpass + '_mapplot.png',dpi=300)
  



if __name__ == '__main__':
    #main(sys.argv[1:])
    #main(sys.argv)
    main()

