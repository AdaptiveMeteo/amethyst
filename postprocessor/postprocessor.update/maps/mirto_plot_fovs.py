import georaster
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from mpl_toolkits.basemap import Basemap
import numpy as np
from osgeo import gdal, osr
import os, sys
from netCDF4 import Dataset
import glob

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


def GetCoordinates(FilePath):
    # get the existing coordinate system
    ds = gdal.Open(FilePath)
    old_cs= osr.SpatialReference()
    old_cs.ImportFromWkt(ds.GetProjectionRef())
    
    # create the new coordinate system
    wgs84_wkt = """
    GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563,
                AUTHORITY["EPSG","7030"]],
            AUTHORITY["EPSG","6326"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.01745329251994328,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4326"]]"""
    new_cs = osr.SpatialReference()
    new_cs.ImportFromWkt(wgs84_wkt)
    
    # create a transform object to convert between coordinate systems
    transform = osr.CoordinateTransformation(old_cs,new_cs) 
    
    #get the point to transform, pixel (0,0) in this case
    width = ds.RasterXSize
    height = ds.RasterYSize
    gt = ds.GetGeoTransform()
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*gt[5] 
    maxx = gt[0] + (gt[1] * width) - gt[1] * 0.5
    maxy = gt[3] - gt[5] * 0.5
    
    #get the coordinates in lat long
    min_latlong = transform.TransformPoint(minx,miny) 
    max_latlong = transform.TransformPoint(maxx,maxy)
    
    return min_latlong, max_latlong, new_cs

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

   filestr = L1DIR + sat_pass_date + '/npp_viirs_?.tif'
   #print('filestr: {}'.format(filestr))

   if instr == "cris":
      tiffile = ''.join([os.path.basename(x) for x in glob.glob(L1DIR + sat_pass_date + '/*.tif')]) 
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
   min_latlong, max_latlong, proj = GetCoordinates(imgfile)
   minx, miny = min_latlong[:2]
   maxx, maxy = max_latlong[:2]
   ext = (minx,maxx,miny,maxy)

   print('minx {}, miny {}, maxx {}, maxy {}, proj {}'.format(minx, miny, maxx, maxy, proj))
   
   # Read geoTIFF
   image = georaster.MultiBandRaster(imgfile,bands='all',load_data=True)
   
   # set Basemap with slightly larger extents
   # set resolution at intermediate level "i"
   #PaoloA m = Basemap( projection='cyl',
   #            llcrnrlon=minx,
   #            llcrnrlat=miny,
   #            urcrnrlon=maxx,
   #            urcrnrlat=maxy,
   #            resolution='h')
   m = Basemap( projection='cyl',
               llcrnrlon=minx,
               llcrnrlat=miny-90.0,
               urcrnrlon=maxx,
               urcrnrlat=maxy-90.0,
               resolution='h')
   
   #m.drawcoastlines(color="red")
   #m.fillcontinents(color='brown')
   #m.drawmapboundary(fill_color='aqua')
   
   # Plot RGB data (normalized to 0...1)
   plt.imshow(image.r/255,extent=ext,alpha=1.0,interpolation='nearest')
   
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
   # plt.scatter(good_x, good_y, color='lightgreen',marker='o',s=8,label='GOOD')
   
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
   
       # Draw the patch
       if fov_index in good_indx:
           e = Ellipse((xlon,lat[fov_index]),L1,L2,angle,facecolor='red',edgecolor='red',lw=.3)
           e.set_fill(True)
           plt.gca().add_patch(e)
       else:
           e = Ellipse((xlon,lat[fov_index]),L1,L2,angle,edgecolor='red',facecolor='none',lw=.3)
           e.set_fill(False)
           plt.gca().add_artist(e)
   
   ##########################################################################
   
   # Plot legend
   #plt.legend(loc='lower right')

   #plt.show()
   plt.savefig(L1DIR + overpass + '/' + overpass + '_mapplot.png',dpi=300)
  



if __name__ == '__main__':
    #main(sys.argv[1:])
    #main(sys.argv)
    main()

