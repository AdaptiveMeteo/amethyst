    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, os
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from mpl_toolkits.basemap import Basemap
import numpy as np
import sys, os
from datetime import datetime
from matplotlib import ticker
import matplotlib.patches as mpatches

__author__    = "Paolo Antonelli and Paolo Scaccia"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__   = ["Paolo Antonelli","Paolo Scaccia"]
__license__   = "--"
__version__ = "0.0.2"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"


from data_reader.standard_data_format import STANDARD_FMT
mirto_color_maps = { STANDARD_FMT[x] : 'magma' for x in STANDARD_FMT.keys() }
mirto_color_maps['air_relative_humidity'] = 'BuPu'
mirto_color_maps['sfc_temperature']       = 'YlOrRd'


#### Grafical Tools for plotting FOVs ######################################################
############################################################################################
####                                                                                   #####
####     Used only with python version equal or less than 3.6                          #####
####                                                                                   ##### 
############################################################################################

if int(sys.version[2]) <= 6 : # check Python version
    import georaster
    from osgeo import gdal, osr
else:
    sys.exit("Python version > 3.6")
    
def GetTiffExtent(FilePath):

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
    
    return [ min_latlong[0], max_latlong[0], min_latlong[1], max_latlong[1] ]


class Tiff_Image(object):
    
    def __init__(self,overpass_path):
                
        TIF_N_HOURS = 4
        
        self.overpass_date = extract_date(overpass_path)
        ref_tif_date = datetime(self.overpass_date.year,self.overpass_date.month,self.overpass_date.day,12,0,0)
        
        if abs((self.overpass_date - ref_tif_date).total_seconds()//60) < 60*TIF_N_HOURS:
            tif_file     = [ file for file in os.listdir(overpass_path) if 'tif' in file and 'dnb' in file ][0]
        else:
            tif_file     = [ file for file in os.listdir(overpass_path) if 'tif' in file and 'color' in file ][0]
    
        self.tif_path = '/'.join([overpass_path,tif_file])
        
        # Get geoTIFF extention
        self.extent = GetTiffExtent(self.tif_path)

        # Read geoTIFF
        image  = georaster.MultiBandRaster(self.tif_path,bands='all',load_data=True)
        self.image = image.r
        
        
def plot_fovs(tif,fov,results):
    """
    Plot FOVs on 

    Parameters
    ----------
    tif_image : TIFF_Image class
        geoTIFF for the background
    tif_extent : geoTIFF extention
        corner of the GeoTIFF image
    fov : MIRTO_fov_Wrapper, MIRTO_tr_Wrapper or MIRTO_BASIC_Wrapper class
        Data structure with lat and lon (with standard format)
    results : MIRTO_fov_Wrapper, MIRTO_tr_Wrapper or MIRTO_BASIC_Wrapper class
        Data structure with info about relative humidity and d2 (with standard format)

    Returns
    -------
    None.

    """
    from obspy.geodetics.base import kilometers2degrees


    x,y = fov.longitude, fov.latitude    

    # Plot Map
    m = Basemap( projection='cyl',
            llcrnrlon=tif.extent[0], # min x
            llcrnrlat=tif.extent[2], # min y
            urcrnrlon=tif.extent[1], # max x
            urcrnrlat=tif.extent[3], # max y
            resolution='h')
    
    
    # Plot geoTIFF
    if tif.image.shape[2] > 1:
        plt.imshow(tif.image/255,extent=tif.extent,alpha=1.0,interpolation='nearest')
    else:
        plt.imshow(tif.image,extent=tif.extent,alpha=1.0,cmap='hot')
        
    # Tan in degrees
    tand = lambda x: np.tan(np.radians(x))

    # Satellite parameters
    sat_ape=[0.963, 0.897] # degrees
    sat_alt = 824; # km
    earth_radius = 6371 # km
    
    for fov_index, xlon in enumerate(x):
    
        # Compute axis and angle
        L1=( tand(sat_ape[0]+fov.fov_angle[fov_index]) - tand(fov.fov_angle[fov_index]))*sat_alt
        L2= tand(sat_ape[1])*sat_alt
        
        # Convert 2 degrees
        L1 = 0.5*kilometers2degrees(L1,radius=earth_radius)
        L2 = 0.5*kilometers2degrees(L2,radius=earth_radius)
        
        angle = -90 - fov.azimuth_angle[fov_index]
    
        # Draw the patch
        e = Ellipse((xlon,y[fov_index]),L1,L2,angle,facecolor="none",edgecolor='red')
        plt.gca().add_patch(e)
    
    plt.show()
    
    return

def read_overpass(overpass_path):
    from data_reader.mirto_wrapper import MIRTO_fov_Wrapper, \
                                          MIRTO_results_Wrapper
    
    res_file = overpass_path + '/mirto/results.nc'
    fov_file = res_file.replace('results','fov')  
    
    # Import data
    res = MIRTO_results_Wrapper(res_file)
    fov = MIRTO_fov_Wrapper(fov_file)
    tif = Tiff_Image(overpass_path)
    
    # Convert to standar data format
    res.to_standard()
    fov.to_standard()
    
    # Compute relative humidity
    res.compute_rh()
    
    return res, fov, tif
    
#######################################################################################    
####  End of Grafical Tools                                                     #######    
#######################################################################################    
def mirto_plot_native(input_path,outdir):
    import sys
    from datetime import datetime
    from netCDF4 import Dataset
    from data_reader.gc_wrapper import GCWrapper
    from piasi_reader.iasi_l1c_native_file import IasiL1cNativeFile
    from geometry.utilities.array_reshapers import array_1d

    MIN_LAT = 65
    
    overpasses = [ path for path in os.listdir(input_path) ] 
    nativepaths = overpasses
    nativepaths = [ path for path in overpasses if '27_0' in path]  # Date Filter - shitly hardcoded
    nativepaths.sort()
    
    ax = plt.gca()
    """
    # Color different sat
    
    sat_types  = np.array(['M01','M02','M03'])
    sat_colors = np.array(['tab:red','tab:blue','tab:green'])
    bool_label = np.array([False,False,False])
    """
    
    native_file_selector = { 'iasi'  : lambda x: '.nat' in x,
                             'cris'  : lambda x: '.h5' in x and 'GCRSO' in x,
                             'combined' : lambda x: ('.nat' in x) or ('.h5' in x and 'GCRSO' in x)
                             }
    fov_colors = { 'cris' : 'salmon', 'iasi': 'tab:blue'}
    if 'iasi' in input_path:
        mode = 'iasi'
        #sat_types = np.array(['M01','M02','M03'])

    elif 'cris' in input_path:
        mode = 'cris'
        # sat_types = np.array(['j01','npp','j02'])

    else:
        mode = 'combined'
    file_filter = native_file_selector[mode]
     
    tot_FOVs = 0
    nativepaths.sort()
    
    for native in nativepaths:
        filename = [ file for file in os.listdir('/'.join([input_path,native])) if file_filter(file)][0]
        if len(filename) == 0:
            print("No native file found in {}. Skipping dir...".format(native))
        else:
            if '.nat' in filename:
                #sat_type = filename.split('_')[3]

                native_data = IasiL1cNativeFile('/'.join([input_path,native,filename]))
                lat  = native_data.get_latitudes()
                lon  = native_data.get_longitudes() 
                instr = 'iasi'
                # Compute cloud fraction
                # avhrr_cloud_fraction = native_data.get_avhrr_cloud_fractions()
                # fov_indices = np.where( (lat > MIN_LAT ) & ( avhrr_cloud_fraction > -95))[0]
            else:
                #sat_type = filename.split('_')[1]
                native_data = GCWrapper('/'.join([input_path, native ,filename]))
                cloudmask_file = '/'.join([input_path, native ,'cloudmask.nc'])
                cloudmask      = Dataset(cloudmask_file,'r')
                lat = np.ma.masked_where(cloudmask['stats'][:,-1] > 2 ,array_1d(native_data.lats))
                lon = np.ma.masked_where(cloudmask['stats'][:,-1] > 2, array_1d(native_data.longs))
                instr = 'cris'
                
            #color    = sat_colors[ np.where(sat_types == sat_type)[0]  ]    
            fov_indices = np.where( lat > MIN_LAT )[0]
            tot_FOVs += fov_indices.size
            
            print('Plotting {} file {}...'.format(instr,native))
            """
            if ~bool_label[ np.where(sat_types == sat_type)[0] ]:
                label='IASI '+sat_type
                bool_label[ np.where(sat_types == sat_type)[0] ] = True
            else:
                label = ''
            """
            mirto_plot_fovs_in_polar_coordinates(lon[fov_indices], lat[fov_indices], 180, 90,label=instr.upper(),color=fov_colors[instr],size=0.1)

    if mode is 'combined':
        plt.legend(framealpha=1.0,shadow=True,loc='upper left').set_zorder(102)          
        
    all_dates = np.array([ datetime.strptime(date,'%Y%m%d_%H%M%S')  for date in nativepaths  ])
    plt.title('{} - {} tp {} - {} FOVs '.format(mode.upper(),all_dates.min().isoformat()[:-3],
                                                all_dates.max().isoformat()[:-3],
                                                tot_FOVs),size=10)
    
    filepath = outdir + '/MIRTO_native_FOVs_{}_{}_{}.png'.format(mode,nativepaths[0].replace('_',''),nativepaths[-1].replace('_',''))
    try:
        plt.savefig(filepath,dpi=600)
        print('Saved plot {}'.format(filepath))
    except Exception as error:
        sys.exit(error)
    return


def mirto_plot_fovs_in_polar_coordinates(lon, lat, lon_0, lat_0,size=1.0, marker='o',color='red',label=''):
    from geometry.earth_geometry import haversine
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt

    xmin = lon.min() 
    xmax = lon.max() 
    ymin = lat.min() 
    ymax = lat.max() 
    
    # width  = haversine(xmin,ymin,xmax,ymin)*1000
    # height = haversine(xmin,ymin,xmin,ymax)*1000
    
    width  = 4000000
    height = 4000000
    
    proj = 'npstere' if lat_0 > 0 else 'spstere'
    prj = Basemap(projection=proj, lon_0=lon_0, lat_0=lat_0,
                  boundinglat=0, resolution='c')
    
    m = Basemap(width=width,height=height,
                resolution='l',projection='stere',\
                lat_ts=80,lat_0=lat_0,lon_0=lon_0)    
    m.fillcontinents(color='palegoldenrod',lake_color='aqua',zorder=1)
    m.drawcoastlines(color='black')
    
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    
    m.drawmapboundary(fill_color='aqua')
    

    ax = plt.gca()
    x,y = m(lon,lat)

    plt.scatter(x, y, s=size,
                marker=marker, 
                c=color,
                label=label,
                zorder=11,
                alpha=1)

    return 

    
def mirto_plot_validation_profiles(x1,x2,x3,y,y_sonde,var,sup_title_str=None ,
                                   best_profile_index = None , outdir = None, file_string = None):
   """
        Input: 
                x1       - first guess profiles   (nfovs x nlev)
                x2       - posterior profiles     (nfovs x nlev)
                x1       - sonde profiles         (nfovs x nlev)
                y        - pressure grids         (nfovs x nlev)
                y_sonde  - sonde pressure grid    (nlev)

                var - var type               str
   """

   fig, ax = plt.subplots(1,2)
   fig.set_size_inches(16.5, 10.5)

   if var in 'T':
      xlabel_str_0 = 'Air Temperature [K]'
      xlabel_str_1 = 'Abs Delta Air Temperature [K]'
      ylabel_str = 'p [hPa]'
      title_str =  'Temperature'
   elif var in 'RH':
      xlabel_str_0 = 'Air Relative Humidity [%]'
      xlabel_str_1 = 'Abs Delta Air Relative Humidity [%]'
      ylabel_str = 'p [hPa]'
      title_str =  'Relative Humidity'
   elif var in 'WV':
      xlabel_str_0 = 'Air Water Vapor MR [g/kg]'
      xlabel_str_1 = 'Abs Delta Air Water Vapor MR [g/kg]'
      ylabel_str = 'p [hPa]'
      title_str =  'Water Vapor MR'
    
   sonde = x3
   sonde_good_levels = np.where(~sonde.mask)[0]
   ylim  = [ int(y[:,sonde_good_levels[0]].min() ) , int(y[:,sonde_good_levels[-1]].max())   ]

   for prof_index, profile in enumerate(zip(x1,x2,y)):
        if prof_index == best_profile_index:
            alpha = 1
            marker = 'o'
            order = 2

        else:
            alpha = 0.12
            marker = ''
            order = 2
            
        first_guess   = profile[0]
        posterior     = profile[1]
        pressure_grid = profile[2]       
        
        # Left plot
        ax[0].plot(first_guess, pressure_grid, 'b+', 
                   posterior, pressure_grid, 'ro',
                   sonde, y_sonde,'gx',
                   linestyle = 'solid',
                   marker    = marker,
                   alpha     = alpha,
                   zorder    = order)
        
        # Right plot
        d1 = np.abs( first_guess - sonde )
        d2 = np.abs( posterior   - sonde )
        ax[1].plot(d1, pressure_grid, 'b+', 
                   d2, pressure_grid, 'ro', 
                   linestyle='solid',
                   alpha=alpha,
                   marker=marker,
                   zorder=order)        


   # Set plot parameters
   ax[0].set_title(title_str, pad = 12, fontsize = 14)
   if var in 'WV':            # For WV fix the upper bound for the y axis  
       ylim[1] = 200          # to 200 hPa
   yticks  = np.flip( np.concatenate(([ylim[1]],np.arange(100,1100,100 ) ) ))
   ax[0].set_ylim(ylim)
   ax[0].set(xlabel = xlabel_str_0, ylabel = ylabel_str)
   ax[1].set(xlabel = xlabel_str_1, ylabel = ylabel_str)
   ax[1].set_title(title_str, pad = 12, fontsize = 14)
   #lgnd = ax[0].legend(['wrf fg', 'mirto ret', 'sonde'], loc='upper right', shadow=True)
   ax[0].grid(True)
   ax[0].set(yscale = 'log')
   
   labels = ['wrf fg', 'mirto ret', 'sonde']
   colors = ['blue','red','green']
   handles = [ mpatches.Patch(color=colors[i],label=l) for i,l in enumerate(labels)]
   lgnd = ax[0].legend(handles=handles, loc='upper right', shadow=True,framealpha = 1.0)

   labels = ['wrf fg - sonde', 'mirto ret - sonde']
   colors = ['blue','red']
   handles = [ mpatches.Patch(color=colors[i],label=l) for i,l in enumerate(labels)]
   ax[1].legend(handles=handles, 
                       loc='upper right', framealpha = 1.0,
                       shadow=True)

   ax[1].set_ylim(ylim)
   ax[1].grid(True)
   ax[1].set(yscale = 'log')
   fig.suptitle(sup_title_str, fontsize = 15)
   # plt.yticks(yticks,labels=yticks)

   ax[0].set_yticks( yticks )
   ax[1].set_yticks( yticks )

   ax[0].get_yaxis( ).set_major_formatter(ticker.LogFormatter())
   ax[1].get_yaxis( ).set_major_formatter(ticker.LogFormatter())
   
   ax[0].set_yticklabels([ repr(x) for x in yticks] )
   ax[1].set_yticklabels([ repr(x) for x in yticks] )

   # plt.yticks(yticks,labels=yticks)

   if outdir is None:
       plt.show()
   else:
       #plt.clf()
       plt.savefig('{}/MIRTO_validation_{}_{}.png'.format(outdir,file_string,var))
       print('Saved plot {}/MIRTO_validation_{}_{}.png'.format(outdir,file_string,var))

   return fig
    

def extract_date(string):
    
    string = string.split('/')[-1]
    
    try:
        return datetime.strptime(string[:15],"%Y%m%d_%H%M%S")
    except:
        datestring = "_".join(string.split('_')[4:6])
        return datetime.strptime(datestring,"%Y%m%d_%H%M%S")
        


###############################################################################
###                                                                         ###
###                 MAIN FUNCTION (test plot GeoTIFF)                       ###
###                                                                         ###
###############################################################################
def prepare_parser():
    from argparse import ArgumentParser

    v_levels = ['debug', 'info', 'warning']
    
    parser = ArgumentParser()
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='The level of verbosity of the software')
    parser.add_argument('--outdir','-o',default='/data/mirto/db/aggregated',type=str,
                        help='Output directory for TR and BASIC file')
    parser.add_argument('--mode','-m',default='plot_overpass',type=str,choices=['plot_overpass','plot_tr','plot_native'],
                        help='Plot mode')
    parser.add_argument('--input','-i',type=str,
                        help='Input file path')

    return parser.parse_args()


def main():
    
      argv = prepare_parser()

      if argv.mode == 'plot_overpass':
          res, fov, tif = read_overpass(argv.input) # Read geoTIFF, fov.nc and results.nc
          plot_fovs(tif, fov, res) 
      elif argv.mode == 'plot_native':
          mirto_plot_native(argv.input,argv.outdir)
          
      # Plot FOVs
      return
  
if __name__ == '__main__':
    sys.exit(main())
