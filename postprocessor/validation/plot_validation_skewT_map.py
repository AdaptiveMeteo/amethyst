import sys
import matplotlib.pyplot as plt
from matplotlib import cm
import glob
import logging
import os 
import numpy as np
from datetime import datetime
import matplotlib.patches as mpatches
from osgeo import gdal
import cartopy.crs as ccrs
from metpy.plots import SkewT
from metpy.units import pandas_dataframe_to_unit_arrays, units
from siphon.simplewebservice.wyoming import WyomingUpperAir
import metpy.calc as mpcalc
import matplotlib.gridspec as gridspec
from time import sleep

sys.path.insert(0,'/home/amethyst/amethyst/postprocessor/validation')
sys.path.insert(0,'/home/amethyst/amethyst/postprocessor/maps')

from utilities.atmos.mirto_atmos_tools import *
from utilities.data_reader.read_netcdf import netCDFReader

__author__    = "Paolo Antonelli and Paolo Scaccia"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__   = ["Paolo Antonelli","Paolo Scaccia"]
__license__   = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

if __name__ == '__main__':
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)

def plot_geotiff_and_fovs(argv):
    from utilities.data_reader.mirto_wrapper import MIRTO_joined_Wrapper

    l1dir    = argv['l1dir']
    overpass = argv['overpass']
    tiffile  = argv['tiff']
    plt = argv['axis']
    if 'special_fov' in argv:
        special_fov = argv['special_fov']
    else:
        special_fov = None
    if 'extra_points' in argv:
        extra_points = argv['extra_points']
    else:
        extra_points = None
    if 'ext_outfile' in argv:
        ext_outfile = argv['ext_outfile']
    else:
        ext_outfile = None
    if 'plot_var' in argv:
        plot_var = argv['plot_var']
    else:
        plot_var = None

    print('TIFF File: {}'.format(tiffile))

    #print('special_fovs index: {}'.format(special_fovs))
    if tiffile != None:
        try:
            ds=gdal.Open(tiffile) # Open GeoTiff File
            if ds == None:
                raise IOError
        except Exception as e:
            sys.exit("Cannot open geoTiff file: {}".format(e))
    
        log.debug("Geo Description:")
        log.debug(ds.GetDescription())
    
        data = ds.ReadAsArray()
        gt   = ds.GetGeoTransform()
        proj = ds.GetProjection()
        log.debug('Geotif projection:')
        log.debug('{}'.format(proj))
        log.debug("")
    else:
        import cartopy.feature

    # Define Projection
    if argv['world_area'] == 'arctic':
        log.debug('AREA: ARCTIC')
        myccrs=  ccrs.Stereographic(central_latitude=90.0, central_longitude=-150.0,
                         false_easting=0.0, false_northing=0.0,
                         true_scale_latitude=60)
    elif argv['world_area'] == 'italy':
        log.debug('AREA: Italy')
        globe = ccrs.Globe(semimajor_axis=6378137, flattening=1/298.257223563)
        myccrs=ccrs.PlateCarree(central_longitude=12.0, globe=globe)
    elif argv['world_area'] == 'pacific':
        log.debug('AREA: PACIFIC')
        globe = ccrs.Globe(semimajor_axis=6378137, flattening=1/298.257223563)
        myccrs=ccrs.PlateCarree(central_longitude=-160.0, globe=globe)

    
    # Define Map extention
    if tiffile != None:
        extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
                  gt[3] + ds.RasterYSize * gt[5], gt[3])
        log.debug('Geografic extent: {}'.format(extent))
        print('Geografic extent: {}'.format(extent))
    #else:
    #    #Redefine Extent using sonde lat lon
    #    if extra_points != None:
    #        extent = [ extra_points[0][1]-3, extra_points[0][1]+3, extra_points[0][0]-3, extra_points[0][0]+3 ]

    if argv['world_area'] == 'pacific' and  tiffile == None:
        extent = [-167, -147, 10, 30]
    elif argv['world_area'] == 'italy' and  tiffile == None:
        extent = [6.0, 18.0, 36.0, 48.0]
    elif argv['world_area'] == 'arctic' and  tiffile == None:
        extent = [-180, 180, 75, 90]

    if argv['world_area'] == 'italy':
        a , b = myccrs.transform_point(extent[0], extent[2], ccrs.Geodetic())
        c , d = myccrs.transform_point(extent[1], extent[3], ccrs.Geodetic())
        extent = (a, c, b, d)
        log.debug('Projected extent: {}'.format(extent))

    # Draw Map wit coast lines
    #plt.figure(figsize=(15, 15))
    ax = plt.subplot(1, 2, 2, projection=myccrs)
    #ax = plt.axes(projection=myccrs)
    gc=ax.coastlines(resolution='10m',linewidth=.75,color='orange')
    gl=ax.gridlines(draw_labels=True, linewidth=.75, color='black')
    gl.bottom_labels = False
    gl.left_labels = False

    if tiffile != None:
            
        log.debug('Data size in pixels: {}'.format(data.shape))
        if data.ndim == 3:
            img = ax.imshow(data[:3, :, :].transpose((1, 2, 0)),extent=extent,origin='upper')
            #img = ax.imshow(data[:3, :, :].transpose((1, 2, 0)),origin='upper')
        else:
            img = ax.imshow(data,extent=extent,cmap='gray',origin='upper')
    else:

        ax.add_feature(cartopy.feature.LAND)
        ax.add_feature(cartopy.feature.OCEAN)
        ax.add_feature(cartopy.feature.LAKES)
        ax.add_feature(cartopy.feature.RIVERS)

    # Read Retrieval data
    resfile = "/".join( [ l1dir , overpass , 'amethyst/amethyst_results.nc'])

    log.debug('Results file: {}'.format(resfile))

    retr = MIRTO_joined_Wrapper(resfile)
    lon = retr.longitude[:]
    lat = retr.latitude[:]

    # Filter FOVs
    rh_max = np.amax(retr.air_relative_humidity,axis=1)
    rh_min = np.min(retr.air_relative_humidity,axis=1)

    good_indx = np.where( (retr.d2[:] < 5) & (rh_max <= 100))[0]
    log.info('Number of good retrievals: {}'.format(len(good_indx)))

    #Get FOV parameters using satellite altitude and instrument apreture
    sat_ape=[0.963, 0.897] # degrees
    sat_alt = 824; # km
    earth_radius = 6371 # km

    tand = lambda x: np.tan(np.radians(x))


    # Define color map       
    if plot_var != None:
        from matplotlib import colors
        from matplotlib import cm
        if plot_var == 'air_relative_humidity':
           my_cmap = cm.get_cmap('BuPu')
        else:
           #my_cmap = cm.get_cmap('magma')
           my_cmap = cm.get_cmap('YlOrRd')

        if retr.__dict__[plot_var].ndim > 1:
           field = retr.__dict__[plot_var].max(axis=1)
        else:
           field = retr.__dict__[plot_var]
        var_min  = field.min()
        var_max  = field.max()
        print('var_min: {}., var_max: {}'.format(var_min,var_max))

        #scale_field = lambda x :  ( x - var_min) / (var_max - var_min)
        normalize = colors.Normalize(vmin=var_min, vmax=var_max, clip=True)
        mapper = cm.ScalarMappable(norm=normalize, cmap=my_cmap)


        sm = cm.ScalarMappable(cmap=my_cmap, norm=normalize)
        sm.set_array([])
        ax.colorbar(sm, label=plot_var, orientation="horizontal",fraction=0.046, pad=0.04)

    else:
        cmap = None
        normalize = None

    
    # Plot FOVs ellipses
    for fov_index, xlon in enumerate(lon):
        
        # Compute axis and angle
        # Double check if axis are semi-axis or full-axis
        L1=( tand(sat_ape[0]+retr.fov_angle[fov_index]) - tand(retr.fov_angle[fov_index]))*sat_alt
        L2= tand(sat_ape[1])*sat_alt
        try:
           angle = -90 - retr.azimuth_angle[fov_index]
        except:
           angle = -90

        #Transform lat and lon in cartesian coordinates
        projx1, projy1 = myccrs.transform_point(xlon, lat[fov_index], ccrs.Geodetic())
        #print(xlon,lat[fov_index],L1,L2,angle)

        # Define Ellipse Face and Edge Colour
        #color = cmap( scale_field(field[ fov_index ]) ) if plot_var != None else 'none'       
        color = mapper.to_rgba(field[fov_index]) if plot_var != None and fov_index in good_indx else 'none'
        edgecolor = 'lightgreen' if fov_index in good_indx else 'silver'
        lw=0.5
        
        if special_fov != None:
            if fov_index == special_fov:
                edgecolor = 'red'
                print('Got Special FOV: {}'.format(fov_index))
                lw=1.0
        
        # Draw the patch
        ax.add_patch(      mpatches.Ellipse(   xy        = [projx1, projy1],
                                               width     = L1*1000,
                                               height    = L2*1000,
                                               angle     = angle,
                                               edgecolor = edgecolor,
                                               alpha     = 1,
                                               facecolor = color,
                                               transform = myccrs,
                                               linewidth = 1 if tiffile !=None else 2.5 )
                     )
        #if fov_index not in good_indx and plot_var != None:
        #    ax.scatter(projx1, projy1,marker='x',c='black',s=2,zorder=10)  

        if extra_points != None:
            # Draw ROWINSONDE SITE
            
            for point in extra_points:
                # Define Ellipse Face and Edge Colour
                #Transform lat and lon in cartesian coordinates
                projx1, projy1 = myccrs.transform_point(point[1], point[0], ccrs.Geodetic())

                L1=5.0
                L2=5.0

                # Draw the patch
                ax.add_patch(      mpatches.Ellipse(   xy=[projx1, projy1],
                                                       width=L1*1000,
                                                       height=L2*1000,
                                                       angle=angle,
                                                       edgecolor='yellow',
                                                       alpha=1,
                                                       facecolor='none',
                                                       transform=myccrs,
                                                       linewidth=0.5 if tiffile != None else 2)
                            )

                L1=150.0
                L2=150.0

                # Draw the patch
                ax.add_patch(      mpatches.Ellipse(   xy=[projx1, projy1],
                                                       width=L1*1000,
                                                       height=L2*1000,

                                                       angle=angle,
                                                       edgecolor='yellow',
                                                       alpha=1,
                                                       facecolor='none',
                                                       transform=myccrs,
                                                       linewidth=0.5 if tiffile != None else 2)
                            )

    try:
        if point[0]>80:
            ax.set_extent([point[1]-6, point[1]+6, point[0]-3, point[0]+3])
        else:
            ax.set_extent([point[1]-3, point[1]+3, point[0]-3, point[0]+3])
    except:
        print("An error occurred while zooming on the rowinsonde site!")
        pass
    
    '''
    if ext_outfile != None:
        outfile = "/".join(  [ outdir , ext_outfile])
    elif plot_var != None:
        outfile = "/".join(  [ outdir , argv['world_area'] + '_' + overpass + '_mapplot_' + plot_var + '.png'])
    else:
        outfile = "/".join(  [ outdir , argv['world_area'] + '_' + overpass + '_mapplot.png'])
    #Note that plt.show opens a new figure. If image has to be saved, savefig command has to be issued before plt.show.
    
    try:
        plt.savefig(outfile,dpi=300)
        log.info("Map saved: {}".format(outfile))
    except Exception as e:
        log.info("Error in saving map plot: {}".format(e))
    '''
    #plt.show()    
    #plt.close()

def sonde_time_rounder(t):
    if t.hour >= 18:
        return t.replace(second=0, microsecond=0, minute=0, hour=0, day=t.day+1)
    elif t.hour <= 6:
        return t.replace(second=0, microsecond=0, minute=0, hour=0)
    else:
        return t.replace(second=0, microsecond=0, minute=0, hour=12)
    

def main():
   import argparse


   parser=argparse.ArgumentParser()
   parser.add_argument('-w',
                        '--workdir',
                        type=str,
                        default='/galaxy/data/hawaii/validation/test0006/',
                        help='Workdir')
   parser.add_argument('-o',
                        '--outdir',
                        type=str,
                        default='/galaxy/data/hawaii/validation/test0006/',
                        help='Output dir')
   parser.add_argument('--l1dir',
                        type=str,
                        default='/galaxy/data/hawaii/test/',
                        help='Outputdir')   
   parser.add_argument('-wa','--world_area',
                        type=str,
                        default='pacific',
                        help='World Area')   

   argv       = parser.parse_args()
   wdir       = argv.workdir
   l1dir      = argv.l1dir
   ncfiles    = glob.glob(wdir + '/*.nc')
   world_area = argv.world_area
   odir       = argv.outdir
   
   for file in ncfiles:
      try:
          nc_val_file = netCDFReader(file)
           
          snd_fmt = '%Y-%m-%dT%H:%M'
          print(file)
          sonde_date = os.path.basename(file).split('_')[4].replace('.nc','')
          print('sonde_date: {}'.format(sonde_date))  
          dt = datetime.strptime(sonde_date, snd_fmt)  
          
          station = os.path.basename(file).split('_')[3]
          print('station: {}'.format(station))
    
          # Read remote sounding data based on time (dt) and station
          df = WyomingUpperAir.request_data(dt, station)
          print('Downloaded Sonde for: {}'.format(dt.strftime("%m/%d/%Y, %H:%M:%S")))
                
          # Create dictionary of united arrays
          data = pandas_dataframe_to_unit_arrays(df)
    
          # Isolate united arrays from dictionary to individual variables
          p  = data['pressure']
          T  = data['temperature']
          Td = data['dewpoint']
          u  = data['u_wind']
          v  = data['v_wind']
            
          fov_lat = np.round(np.ma.getdata(nc_val_file.fov_lat).astype(float),2)
          fov_lon = np.round(np.ma.getdata(nc_val_file.fov_lon).astype(float),2)
    
          for index, val in enumerate(fov_lat):  
    
              overpass = nc_val_file.retr_name[index].replace('MIRTO_','')
              print('overpass: {}'.format(overpass))
              sat_fmt = '%Y%m%d_%H%M%S'
              t = datetime.strptime(overpass, sat_fmt)
    
              title_str = "MIRTO FOV Lat: {}, Lon: {}".format(fov_lat[index],                          
                                                              fov_lon[index])
              print(title_str)
    
              # Change default to be better for skew-T
              fig = plt.figure(figsize=(18, 11))
              #fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(18, 11))
    
    
              # Grid for plots
              gs = gridspec.GridSpec(1, 2)
              skew = SkewT(fig, rotation=45, subplot=gs[0, 0])
              
              # Initiate the skew-T plot type from MetPy class loaded earlier
              #skew = SkewT(fig, rotation=45,subplot(1,2,1))
    
              #PaoloS 16/07/2021
              retp = nc_val_file.fovs_levels[index]
              #retp = nc_val_file.fovs_levels[index]
              #p = np.ma.getdata(nc_val_file.pressure).astype(float)
              #print('pressure: {}'.format(retp))
    
              fgT = nc_val_file.prior_temp[index]
              retT = nc_val_file.post_temp[index]
              #sT = nc_val_file.sonde_temp
    
              fgMR = nc_val_file.prior_water_vapour[index]
              retMR = nc_val_file.post_water_vapour[index]
              #sMR = nc_val_file.sonde_water_vapour
    
              T0 = 273.15
    
              fgDP = mr2dp(retp,fgT,fgMR)-T0
              retDP = mr2dp(retp,retT,retMR)-T0
              #sDP = mr2dp(p,sT,sMR)
    
              # Plot the data using normal plotting functions, in this case using
              # log scaling in Y, as dictated by the typical meteorological plot
              skew.plot(p, T, 'k')
              skew.plot(retp * units.hPa, fgT * units.K, 'b')
              skew.plot(retp * units.hPa, retT * units.K, 'r')
    
              skew.plot(p, Td, 'grey')
              skew.plot(retp * units.hPa, fgDP * units.K, 'cornflowerblue')
              skew.plot(retp * units.hPa, retDP * units.K, 'orangered')
    
              #plot winds
              skew.plot_barbs(p[::3], u[::3], v[::3], y_clip_radius=0.03)
    
              # Set some appropriate axes limits for x and y
              skew.ax.set_xlim(-30, 40)
              skew.ax.set_ylim(1020, 200)
    
              # Add the relevant special lines to plot throughout the figure
              skew.plot_dry_adiabats(t0=np.arange(233, 533, 10) * units.K,
                               alpha=0.25, color='orangered')
              skew.plot_moist_adiabats(t0=np.arange(233, 400, 5) * units.K,
                                 alpha=0.25, color='tab:green')
              skew.plot_mixing_lines()
    
              # Calculate LCL height and plot as black dot. Because `p`'s first value is
              # ~1000 mb and its last value is ~250 mb, the `0` index is selected for
              # `p`, `T`, and `Td` to lift the parcel from the surface. If `p` was inverted,
              # i.e. start from low value, 250 mb, to a high value, 1000 mb, the `-1` index
              # should be selected.
              # lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
              # skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')
    
              # Calculate full parcel profile and add to plot as black line
              #prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
              #skew.plot(p, prof, 'k', linewidth=2)
    
              # Shade areas of CAPE and CIN
              #skew.shade_cin(p, T, prof, Td)
              #skew.shade_cape(p, T, prof)
    
              # Add some descriptive titles
              plt.title(title_str, loc='left')
              plt.title('Valid Time: {}'.format(t.strftime("%m/%d/%Y, %H:%M:%S")),loc='right')
    
              #argv.outdir = odir 
              tifffile = ''.join([os.path.basename(x) for x in glob.glob( "/".join( [l1dir , overpass , 'n*.tif']) )  ])
              if tifffile == '':
                  geotiff = None
              else:
                  geotiff = "/".join( [l1dir , overpass , tifffile] )

              print('GEOTIff: {}'.format(geotiff))
              
              #results_file = netCDFReader("/".join([l1dir,overpass,"amethyst/amethyst_results.nc"]))
              fov_index = np.ma.getdata(nc_val_file.fov_index[index]).astype(int)
              print('Fov index: {}'.format(fov_index))
    
              plot_var = None 
    
              # Rowinsonde Point
              extra_points = [ (np.ma.getdata(nc_val_file.sonde_lat),np.ma.getdata(nc_val_file.sonde_lon) ) ]
    
              outfile = 'MAP_' + os.path.basename(file).replace('nc','png')
              #outfile = outfile.replace(':','-') 
              #print('outfile: {}'.format(outfile))
    
              #Create argv dictionary for amethyst_sounder_mapper_on_imager_background
              argv = {'l1dir':l1dir, 'overpass':overpass, 'axis':plt,
                      'tiff':geotiff, 'world_area':world_area, 'verbose':'info', 'plot_var':plot_var, 
                      'ext_outfile':outfile, 'special_fov':fov_index, 'extra_points':extra_points
                      }
    
              plot_geotiff_and_fovs(argv)
              ofile = "/".join( [  odir, 
                                  'SKEWT_MAP_'+overpass+'_FOV_'+str(fov_index)+'.png' ] )
    
              #ofile = ofile.replace(':','-') 
              print("Saving " + ofile)
              plt.savefig(ofile, transparent = True, dpi=300)
              plt.close()
              
              sleep(4)
      except Exception as e:
         print("Skipped file {}: {}".format(file,e))
         pass
if __name__ == "__main__":
    log = logging.getLogger()
    main()
