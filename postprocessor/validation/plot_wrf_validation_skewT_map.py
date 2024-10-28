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
from utilities.atmos.mirto_atmos_tools import *
from utilities.data_reader.read_netcdf import netCDFReader
#from postprocessing.maps.mirto_sounder_mapper_on_imager_background import plot_geotiff_and_fovs

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

    wrf_data = argv['wrf']
    forecast = argv['forecast']
    #tiffile  = argv['tiff']
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

    #print('special_fovs index: {}'.format(special_fovs))

    lon = wrf_data.wrf_lon[:]
    lat = wrf_data.wrf_lat[:]

    # Define Projection
    if argv['world_area'] == 'arctic':
        log.debug('AREA: ARCTIC')
        myccrs=  ccrs.Stereographic(central_latitude=90.0, central_longitude=-150.0,
                         false_easting=0.0, false_northing=0.0,
                         true_scale_latitude=60)
    elif argv['world_area'] == 'pacific':
        log.debug('AREA: PACIFIC')
        globe = ccrs.Globe(semimajor_axis=6378137, flattening=1/298.257223563)
        myccrs=ccrs.PlateCarree(central_longitude=-160.0, globe=globe)

    #Redefine Extent using sonde lat lon
    #if extra_points != None:
    #    extent = [ extra_points[0][1]-3, extra_points[0][1]+3, extra_points[0][0]-3, extra_points[0][0]+3 ]

    if argv['world_area'] == 'pacific':
        a , b = lon.max(), lat.max()
        #myccrs.transform_point(extent[0], extent[2], ccrs.Geodetic())
        c , d = lon.min(), lat.min()
        # myccrs.transform_point(extent[1], extent[3], ccrs.Geodetic())
        extent = (a, c, b, d)
        log.debug('Projected extent: {}'.format(extent))

    # Draw Map wit coast lines
    #plt.figure(figsize=(15, 15))
    ax = plt.subplot(1, 2, 2, projection=myccrs)
    #ax = plt.axes(projection=myccrs)
    gc=ax.coastlines(resolution='10m',linewidth=1.5,color='orange')
    gl=ax.gridlines(draw_labels=True, linewidth=.75, color='black')
    #ax.drawmapboundary(fill_color='lightblue')
    gl.bottom_labels = False
    gl.left_labels = False
            
    from matplotlib import colors
    from matplotlib import cm

    # Define color map       
    cmap = None
    normalize = None

    
    # Plot FOVs ellipses
    for index, xlon,ylat in zip(range(lon.size),lon,lat):
        
        #Transform lat and lon in cartesian coordinates
        projx1, projy1 = myccrs.transform_point(xlon, ylat, ccrs.Geodetic())
        #print(xlon,lat[fov_index],L1,L2,angle)

        # Define Ellipse Face and Edge Colour
        #color = cmap( scale_field(field[ fov_index ]) ) if plot_var != None else 'none'       
        color     = 'none'
        edgecolor = 'silver'
        lw=0.5
        # Draw the patch
        ax.scatter( projx1, projy1,
                    edgecolor = edgecolor,
                    transform = myccrs,
                    s = 35,
                    color = 'tab:blue',
                    marker = 'X',
                    zorder=4,
                    label = 'WRF Grid center')
                     
        #if fov_index not in good_indx and plot_var != None:
        #    ax.scatter(projx1, projy1,marker='x',c='black',s=2,zorder=10)  

        if extra_points != None:
            # Draw ROWINSONDE SITE
            
            for point in extra_points:
                # Define Ellipse Face and Edge Colour
                #Transform lat and lon in cartesian coordinates
                projx1, projy1 = myccrs.transform_point(point[1], point[0], ccrs.Geodetic())

                angle = 90
                # Draw the patch
                ax.scatter( projx1, projy1,
                                edgecolor = edgecolor,
                                transform = myccrs,
                                s = 35,
                                color = 'red',
                                marker = 'o',
                                zorder=3,
                                label = 'Sonde Site')
                """
                ax.add_patch(      mpatches.Ellipse(   xy=[projx1, projy1],
                                                       width=L1*1000,
                                                       height=L2*1000,
                                                       angle=angle,
                                                       edgecolor='red',
                                                       alpha=1,
                                                       facecolor='none',
                                                       transform=myccrs,
                                                       linewidth=0.8,
                                                       label = '100 km')
                            )
                """
                L1=200.0
                L2=200.0
                #ax.text(projx1, projy1 + L1*1000,"150 km",transform=myccrs,c='red')
                # Draw the patch
                ax.add_patch(      mpatches.Ellipse(   xy=[projx1, projy1],
                                                       width=L1*1000,
                                                       height=L2*1000,

                                                       angle=angle,
                                                       edgecolor='red',
                                                       alpha=1,
                                                       facecolor='none',
                                                       transform=myccrs,
                                                       linewidth=0.8,
                                                       label = '100 km'))
    try:
        if point[0]>80:
            ax.set_extent([point[1]-6, point[1]+6, point[0]-3, point[0]+3])
        else:
            ax.set_extent([point[1]-3, point[1]+3, point[0]-3, point[0]+3])
    except:
        print("An error occurred while zooming on the rowinsonde site!")
        pass
    ax.legend(prop={'size':14})
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
   parser.add_argument('-wf','--wrf_file',
                        type=str,
                        default='/galaxy/data/hawaii/test/',
                        help='Outputdir')   
   parser.add_argument('-wa','--world_area',
                        type=str,
                        default='pacific',
                        help='World Area')   

   argv       = parser.parse_args()
   wdir       = argv.workdir
   wrf_file  = argv.wrf_file
   world_area = argv.world_area
   odir       = argv.outdir
   
   ncfiles    = glob.glob(wdir + '/*.nc')

   for file in ncfiles:
          #try:
          nc_val_file = netCDFReader(file)
           
          snd_fmt = '%Y%m%d%H%M'
          print("Plotting file ",file)
          sonde_date = os.path.basename(file).split('_')[6]
          print('sonde_date: {}'.format(sonde_date))  
          dt = datetime.strptime(sonde_date, snd_fmt)  
          
          station = os.path.basename(file).split('_')[4]
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
            
          wrf_lat = np.round(np.ma.getdata(nc_val_file.wrf_lat).astype(float),2)
          wrf_lon = np.round(np.ma.getdata(nc_val_file.wrf_lon).astype(float),2)
    
          for index, lon, lat in zip(range(wrf_lat.size),wrf_lon,wrf_lat):  
    
              forecast = nc_val_file.forecast[index].replace('wrfout_','')
              if 'analysis' in file:
                   print('Analysis: {}'.format(forecast))
              else:
                   print('Forecast: {} + {} '.format(forecast,os.path.basename(file).split('+')[-1].replace('.nc','')))
                   
              
              sat_fmt = '%Y-%m-%d_%H:%M:%S'
              t = datetime.strptime(forecast[4:], sat_fmt)
              title_str = "Lat: {}, Lon: {}".format(lon,lat)
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
              wrf_p = nc_val_file.wrf_levels[index]
              #retp = nc_val_file.fovs_levels[index]
              #p = np.ma.getdata(nc_val_file.pressure).astype(float)
              #print('pressure: {}'.format(retp))
              T0 = 273.15

              wrf_T = nc_val_file.temp[index] - T0
              #sT = nc_val_file.sonde_temp
    
              wrf_MR = nc_val_file.water_vapour[index]
    
    
              wrf_DP = mr2dp(wrf_p* 1e-2,wrf_T,wrf_MR)-T0
    
              # Plot the data using normal plotting functions, in this case using
              # log scaling in Y, as dictated by the typical meteorological plot
              skew.plot(p, T, 'k',label = 'Rowinsonde T')
              skew.plot(p, Td, 'grey',label='Rowinsonde DP')

              skew.plot(wrf_p *1e-2, wrf_T , 'b',label = 'WRF T')
              #skew.plot(wrf_p , wrf_T , 'b',label = 'WRF T')
              skew.plot(wrf_p * 1e-2, wrf_DP , 'cornflowerblue',label='WRF DP')
    
              #plot winds
              skew.plot_barbs(p[::3], u[::3], v[::3], y_clip_radius=0.03)
    
              # Set some appropriate axes limits for x and y
              skew.ax.set_xlim(-30, 40)
              skew.ax.set_ylim(1020, 100)
    
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
                  
              #results_file = netCDFReader("/".join([l1dir,overpass,"mirto/results.nc"]))
              fov_index = np.ma.getdata(nc_val_file.wrf_index[index]).astype(int)
              print('Grid Cell index: {}'.format(fov_index))
    
              plot_var = None 
              # Rowinsonde Point
              extra_points = [ (np.ma.getdata(nc_val_file.sonde_lat),np.ma.getdata(nc_val_file.sonde_lon) ) ]
    
              outfile = 'MAP_' + os.path.basename(file).replace('nc','png')
              #outfile = outfile.replace(':','-') 
              #print('outfile: {}'.format(outfile))
    
              #Create argv dictionary for mirto_sounder_mapper_on_imager_background
              argv = {'wrf':nc_val_file, 'forecast':forecast, 'axis':plt,
                      'tiff':None, 'world_area':world_area, 'verbose':'info', 'plot_var':plot_var, 
                      'ext_outfile':outfile, 'special_fov':nc_val_file.wrf_index, 'extra_points':extra_points
                      }
              fig.legend(prop={'size':11},bbox_to_anchor=(0.3, 0., 0.15, 0.83))
              
              plot_geotiff_and_fovs(argv)
              ofile = "/".join( [  odir, 
                                  'SKEWT_MAP_'+forecast+'_'+str(nc_val_file.wrf_index[index])+'.png' ] )
    
              #ofile = ofile.replace(':','-') 
              print("Saving " + ofile)
              #plt.show()
              plt.savefig(ofile, transparent = True, dpi=300)
              plt.close()
              
              sleep(4)
              #except Exception as e:
              #print("Skipped file {}: {}".format(file,e))
              #pass
if __name__ == "__main__":
    log = logging.getLogger()
    main()

