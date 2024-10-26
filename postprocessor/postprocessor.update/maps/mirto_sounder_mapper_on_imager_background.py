import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os, sys
import glob
from osgeo import gdal
import cartopy.crs as ccrs
import sys
import logging
from data_reader.standard_data_format import mirto_colormaps, get_datatype

if __name__ == '__main__':
    log = logging.getLogger()
else:
    log = logging.getLogger(__name__)


def local_parser():
    """
        Terminal line parser   

    Returns
    -------
        argv : ArgumentParser
    """
    from argparse         import ArgumentParser
    
    v_levels = ['debug', 'info', 'warning']
    parser = ArgumentParser()
    parser.add_argument('--overpass','-ov',type=str,required=True,
                        help='String with overpass dirname %Y%m%d_%H%M%S')
    parser.add_argument('--basedir','-b',type=str,required=False,default='/home/oper/cris/common',
                        help='Base directory when running this script')
    parser.add_argument('--verbose','-v',type=str,choices=v_levels,default='info',
                        help='Logger verbosity')
    parser.add_argument('--l1dir',type=str,required=False,default= '/galaxy/data/arctic/cris/',
                        help='Logger verbosity')
    parser.add_argument('--tiff','-t',type=str,required=False,default= None ,
                        help='GeoTiff file path')
    parser.add_argument('--plot_var',type=str, default=None,
                        help='Retrieval field to be use for colouring FOVs')
    parser.add_argument('--outdir','-o',type=str, default=None,
                        help='Output directory')
    parser.add_argument('-wa','--world_area',type=str, required=True,choices=['arctic','pacific'],
                       help='World Area')
    parser.add_argument('--plot_var_mode',type=str, required=False,default=None ,
                       help='How to operate on retrieval field')
    parser.add_argument('--eps',type=bool, required=False,default=False ,
                   help='EPS Mode')

    parser.add_argument('--mode','-m',type=str,default='tif_and_fovs',choices=['tif_and_fovs','cloud_properties'],
                       help='Plot mode')
    argv = parser.parse_args()

    sys.path.insert(0,argv.basedir)
    if argv.tiff == None:
        argv.tiff = ''.join([os.path.basename(x) for x in glob.glob(argv.l1dir + argv.overpass + '/n*.tif')])
        argv.tiff = argv.l1dir + argv.overpass + '/' +argv.tiff
        print('tiff: {}'.format(argv.tiff))

    if argv.outdir == None:
        argv.outdir = "/".join([ argv.l1dir , argv.overpass ])
        
    return argv
  
def plot_geotiff_and_fovs(argv, plot_var = None):
    from data_reader.mirto_wrapper import MIRTO_joined_Wrapper

    l1dir    = argv.l1dir
    overpass = argv.overpass
    tiffile  = argv.tiff
    outdir   = argv.outdir

   
    try:
        ds=gdal.Open(tiffile) # Open GeoTiff File
        if ds == None:
            raise IOError
        tif_background = True
    except Exception as e:
        log.info("Cannot open geoTiff file: {}".format(e))                 
        tif_background = False
        
    if tif_background:
        log.debug("Geo Description:")
        log.debug(ds.GetDescription())

        data = ds.ReadAsArray()
        gt   = ds.GetGeoTransform()
        proj = ds.GetProjection()
        log.debug('Geotif projection:')
        log.debug('{}'.format(proj))
        log.debug("")

    
    #Proj4 string used to generate VIIRS truecolor image
    #proj4str='+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'   
 
    #Cartopy CRS definition based on proj4 parameters used to generate VIIRS truecolor image
    #defined in: /home/oper/cris/cspp/polar2grid_v_2_3/lib/python3.7/site-packages/polar2grid/grids/grids.conf
    
    #For the Arctic:
    #polar_alaska, proj4, +proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lat_ts=60.0 +lon_0=-150 +units=m, 400,-400
    #crs = ccrs.Stereographic(central_latitude=90.0, central_longitude=-150.0,
    #                 false_easting=0.0, false_northing=0.0,
    #                 true_scale_latitude=60, globe=None)
    
        
    # Define Projection
    if argv.world_area == 'arctic':
        log.debug('AREA: ARCTIC')
        myccrs=  ccrs.Stereographic(central_latitude=90.0, central_longitude=-150.0,
                                    false_easting=0.0,     false_northing=0.0,
                                    true_scale_latitude=60)
    elif argv.world_area == 'pacific':
        log.debug('AREA: PACIFIC')
        globe = ccrs.Globe(semimajor_axis=6378137, flattening=1/298.257223563)
        myccrs=ccrs.PlateCarree(central_longitude=-160.0, globe=globe)           

    # Draw Map wit coast lines
    plt.figure(figsize=(15, 15))
    ax = plt.axes(projection=myccrs)
    gc=ax.coastlines(resolution='10m',linewidth=.75,color='orange')
    #gl=ax.gridlines(linewidth=.5,color='black')
    gl=ax.gridlines(draw_labels=True, linewidth=.75, color='black')
    gl.bottom_labels = False
    gl.left_labels = False

    # Read Retrieval data
    resfile = "/".join( [ l1dir , overpass , 'mirto/results.nc'])
    log.debug('Results file: {}'.format(resfile))
    retr = MIRTO_joined_Wrapper(resfile)
    lon = retr.longitude[:] 
    lat = retr.latitude[:]

    if tif_background:
        # Define Map extention
        extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
                  gt[3] + ds.RasterYSize * gt[5], gt[3])
        log.debug('Geografic extent: {}'.format(extent))
    
        if argv.world_area == 'pacific':
            a , b = myccrs.transform_point(extent[0], extent[2], ccrs.Geodetic())
            c , d = myccrs.transform_point(extent[1], extent[3], ccrs.Geodetic())
            extent = (a, c, b, d)
            log.debug('Projected extent: {}'.format(extent))
    
        # Plot Tif Image    
        log.debug('Data size in pixels: {}'.format(data.shape))
        if data.ndim == 3:
            img = ax.imshow(data[:3, :, :].transpose((1, 2, 0)),extent=extent,origin='upper')
        else:
            img = ax.imshow(data,extent=extent,cmap='gray',origin='upper')
    else:
        xmin, xmax = lon.min(), lon.max()
        ymin, ymax = lat.min(), lat.max()

        ax.set_extent([ xmin + 12.5, xmax-12.5, ymin+12.5,ymax-12.5],crs=ccrs.PlateCarree())
        ax.background_img()

    # Filter FOVs
    rh_max = np.amax(retr.air_relative_humidity,axis=1)
    rh_min = np.min(retr.air_relative_humidity,axis=1)

    file_object = open(outdir + 'stats.txt', 'a+')

    good_indx = np.where( (retr.d2[:] < 5) & (rh_max < 100))[0]
    log.info('Number of good retrievals: {}'.format(len(good_indx)))
    file_object.write('Overpass: {}\n'.format(overpass))
    file_object.write('Number of good retrievals: {}\n'.format(len(good_indx)))
    good_indx = np.where( (retr.d2[:] < 5) & (rh_max <= 100))[0]
    log.info('Number of good retrievals including RH=100%: {}'.format(len(good_indx)))
    file_object.write('Number of good retrievals including RH=100%: {}\n'.format(len(good_indx)))
     
    file_object.close()

    
    #Get FOV parameters using satellite altitude and instrument apreture
    sat_ape=[0.963, 0.897] # degrees
    sat_alt = 824; # km
    earth_radius = 6371 # km

    tand = lambda x: np.tan(np.radians(x))

    # Define color map       
    if plot_var != None:
        from matplotlib import colors
        from matplotlib import cm

        # Read Proper colormap
        cmap = mirto_colormaps[plot_var]            

        # Import field as an atmospheric, surface or spectral datatype
        field = get_datatype(plot_var,retr.__dict__[plot_var]) 
        
        # If atmospheric datatype choose how to compute the scalar field
        if 'atmos' in str(type(field)):
            if argv.plot_var_mode != None:
                field.compute_scalar_field(argv.plot_var_mode)
            else:
                sys.exit("Specify how to compute the scalar field from atmospheric profiles")

        #field = np.copy(field.scalar_field)
        field = np.ma.masked_where(field.scalar_field >=  1e5, field.scalar_field)
        var_min  = field.min()
        var_max  = field.max()

        normalize = colors.Normalize(vmin=var_min, vmax=var_max, clip=True)
        mapper = cm.ScalarMappable(norm=normalize, cmap=cmap)

        sm = cm.ScalarMappable(cmap=cmap, norm=normalize)
        sm.set_array([])
        plt.colorbar(sm, label=plot_var, orientation="horizontal",fraction=0.046, pad=0.04)    

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
        except Exception as e:
            angle = -90
            log.info("Azimuth angle value not available, setting value to: {}".format(angle))

        #Transform lat and lon in cartesian coordinates
        projx1, projy1 = myccrs.transform_point(xlon, lat[fov_index], ccrs.Geodetic())
        
        # Define Ellipse Face and Edge Colour
        color = mapper.to_rgba(field[fov_index]) if plot_var != None else 'none'
        if fov_index in good_indx and rh_max[fov_index] < 100:
            edgecolor = 'lightgreen' if tif_background else 'darkslategray'
            talpha = 1
        elif fov_index in good_indx and rh_max[fov_index] == 100:
            edgecolor = 'red'
            talpha = 1
        else:
            edgecolor = 'yellow'
            talpha = 1 if argv.eps else 0.1
        zorder    = 1 if fov_index in good_indx else 2
        # Draw the patch
        ax.add_patch(      mpatches.Ellipse(   xy=[projx1, projy1], 
                                               width=L1*1000, 
                                               height=L2*1000, 
                                               angle=angle, 
                                               edgecolor= edgecolor,
                                               alpha=talpha,
                                               facecolor= color, 
                                               transform=myccrs,
                                               linewidth=0.4,
                                               zorder=zorder)
                     )  
        
        """
        if fov_index not in good_indx and plot_var != None:
            ax.scatter(projx1, projy1,marker='X',c='red',s=2,zorder=11)
        """
        
    #<----- cycle end
            
    if plot_var != None:
        outfile = "/".join(  [ outdir , argv.world_area + '_' + overpass + '_mapplot_' + plot_var + '.png'])
        if argv.plot_var_mode != None:
            outfile = outfile.replace('.png','_{}.png'.format(argv.plot_var_mode))
    else:
        outfile = "/".join(  [ outdir , argv.world_area + '_' + overpass + '_mapplot.png'])
    #Note that plt.show opens a new figure. If image has to be saved, savefig command has to be issued before plt.show.
        
    if argv.eps:
        outfile = outfile.replace('.png','.eps')
    try:
        plt.savefig(outfile,dpi=300)
        log.info("Map saved: {}".format(outfile))
    except Exception as e:
        log.info("Error in saving map plot: {}".format(e))

    #plt.show()    

def plot_cloud_properties(argv, plot_var):
    from data_reader.mirto_wrapper import MIRTO_fov_Wrapper
    from netCDF4 import Dataset
    from matplotlib import colors
    from matplotlib import cm

    l1dir    = argv.l1dir
    overpass = argv.overpass
    tiffile  = argv.tiff
    outdir   = argv.outdir
   
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
        
    # Define Projection
    if argv.world_area == 'arctic':
        log.debug('AREA: ARCTIC')
        myccrs=  ccrs.Stereographic(central_latitude=90.0, central_longitude=-150.0,
                                    false_easting=0.0,     false_northing=0.0,
                                    true_scale_latitude=60)
    elif argv.world_area == 'pacific':
        log.debug('AREA: PACIFIC')
        globe = ccrs.Globe(semimajor_axis=6378137, flattening=1/298.257223563)
        myccrs=ccrs.PlateCarree(central_longitude=-160.0, globe=globe)           

    # Draw Map wit coast lines
    plt.figure(figsize=(15, 15))
    ax = plt.axes(projection=myccrs)
    gc=ax.coastlines(resolution='10m',linewidth=.75,color='orange')
    #gl=ax.gridlines(linewidth=.5,color='black')
    gl=ax.gridlines(draw_labels=True, linewidth=.75, color='black')
    gl.bottom_labels = False
    gl.left_labels = False
    #gl.tick_params(labelsize=15) 
    
    # Define Map extention
    extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
              gt[3] + ds.RasterYSize * gt[5], gt[3])
    log.debug('Geografic extent: {}'.format(extent))

    if argv.world_area == 'pacific':
        a , b = myccrs.transform_point(extent[0], extent[2], ccrs.Geodetic())
        c , d = myccrs.transform_point(extent[1], extent[3], ccrs.Geodetic())
        extent = (a, c, b, d)
        log.debug('Projected extent: {}'.format(extent))

    log.debug('Data size in pixels: {}'.format(data.shape))
    if data.ndim == 3:
        img = ax.imshow(data[:3, :, :].transpose((1, 2, 0)),extent=extent,origin='upper')
    else:
        img = ax.imshow(data,extent=extent,cmap='gray',origin='upper')
        
    # Read Retrieval data
    fovfile = "/".join( [ l1dir , overpass , 'mirto/fov.nc'])
    cloudmask_file  = "/".join( [ l1dir , overpass , 'cloudmask.nc'])
    cloud_prop_file = "/".join( [ l1dir , overpass , 'cloud_properties.nc'])
    log.debug('FOV file: {}'.format(fovfile))   
    retr = MIRTO_fov_Wrapper(fovfile)
    retr.to_standard()
    lon = retr.longitude[:] 
    lat = retr.latitude[:]    
    cloudmask  = Dataset(cloudmask_file)
    cloud_prop = Dataset(cloud_prop_file) 
  
    #Get FOV parameters using satellite altitude and instrument apreture
    sat_ape=[0.963, 0.897] # degrees
    sat_alt = 824; # km
    earth_radius = 6371 # km
    tand = lambda x: np.tan(np.radians(x))

    # Define color map       
    # Read Proper colormap
    var_map = { 'cloud_top_temperature':{ 'color':'bwr',
                                          'label':'Cloud Top Temperature [K]',
                                          'vmin':200,
                                          'vmax':300},
                'cloud_top_pressure':{ 'color':'cool_r',
                                       'label':'Cloud Top Pressure [hPa]',                                                                          
                                       'vmin':100,
                                       'vmax':1015},
                'cloud_top_emissivity':{ 'color':'jet_r',
                                          'label':'Cloud Top Emissivity',
                                          'vmin':0,
                                          'vmax':1},
                'none':{'color':'none',
                        'label':'none'}
                }
    
    cmap = var_map[plot_var]['color']
    
    if plot_var != 'none':
        field = cloud_prop[plot_var][:]
        var_min  = field.min()
        var_max  = field.max()

        normalize = colors.Normalize(vmin=var_map[plot_var]['vmin'], vmax=var_map[plot_var]['vmax'], clip=True)
        mapper = cm.ScalarMappable(norm=normalize, cmap=cmap)
        sm = cm.ScalarMappable(cmap=cmap, norm=normalize)
        sm.set_array([])
        cbar = plt.colorbar(sm, label=var_map[plot_var]['label'], orientation="horizontal",fraction=0.046, pad=0.04)    
        cbar.ax.tick_params(labelsize=13) 
        cbar.set_label(var_map[plot_var]['label'], labelpad=10,size=20)

    # Plot FOVs ellipses
    for fov_index, fov in enumerate(zip(lon,lat)):
        xlon = fov[0]
        xlat = fov[1]
        cld_index = np.where(cloudmask['cris_lons'][:] == xlon  )[0][0]      
        # Compute axis and angle
        # Double check if axis are semi-axis or full-axis
        L1=( tand(sat_ape[0]+retr.FOV_angle[fov_index]) - tand(retr.FOV_angle[fov_index]))*sat_alt
        L2= tand(sat_ape[1])*sat_alt                    
        try:
            angle = -90 - retr.Satellite_azimuth_angle[fov_index]       
        except:
            angle = -90
            log.info("Azimuth angle value not available, setting value to: {}".format(angle))

        #Transform lat and lon in cartesian coordinates
        projx1, projy1 = myccrs.transform_point(xlon, xlat, ccrs.Geodetic())       
        # Define Ellipse Face and Edge Colour
        color = mapper.to_rgba(field[cld_index]) if plot_var != 'none' else 'none'
        edgecolor = 'white' 
        talpha = 1
        zorder = 1
        # Draw the patch
        ax.add_patch(      mpatches.Ellipse(   xy=[projx1, projy1], 
                                               width=L1*1000, 
                                               height=L2*1000, 
                                               angle=angle, 
                                               edgecolor= edgecolor,
                                               alpha=talpha,
                                               facecolor= color, 
                                               transform=myccrs,
                                               linewidth=0.4,
                                               zorder=zorder)
                     )          
        """
        if fov_index not in good_indx and plot_var != None:
            ax.scatter(projx1, projy1,marker='X',c='red',s=2,zorder=11)
        """
        
    #<----- cycle end
    plt.title('CLOUD TOP PROPERTIES',size=20,weight='bold',pad=10)
    outfile="/".join([outdir,argv.world_area+'_'+overpass+'_mapplot_'+plot_var+'.png'])
    #Note that plt.show opens a new figure. If image has to be saved, savefig   command has to be issued before plt.show.
    try:
        plt.savefig(outfile,dpi=300,bbox_inches='tight')
        log.info("Map saved: {}".format(outfile))
    except Exception as e:
        log.info("Error in saving map plot: {}".format(e))

    plt.clf()    
    return
    
def main():

    # Parse inline arguments    
    argv = local_parser()   
    # Prepare the log class
    verbosity = getattr(logging, argv.verbose.upper())
    log.setLevel(verbosity)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - '
                                  '%(funcName)s: %(message)s',
                                  datefmt='%m/%d/%Y %H:%M:%S')
    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(verbosity)
    streamhandler.setFormatter(formatter)
    log.addHandler(streamhandler)
    
    # Debug printouts
    log.debug("Overpass:    {}".format(argv.overpass))        
    log.debug("L1dir:       {}".format(argv.l1dir))               
    log.debug("Tiff file:   {}".format(argv.tiff))        
    log.debug("FOV var:     {}".format(argv.plot_var))        
    log.debug("Area:        {}".format(argv.world_area))        
    
    if argv.mode == 'tif_and_fovs':
        plot_geotiff_and_fovs(argv,plot_var = argv.plot_var)
    elif argv.mode == 'cloud_properties':
        plot_cloud_properties(argv,'cloud_top_temperature')
        plot_cloud_properties(argv,'cloud_top_pressure')
        plot_cloud_properties(argv,'cloud_top_emissivity')
        plot_cloud_properties(argv,'none')

    else:
        log.info("Wrong plot mode!")
    return

if __name__ == '__main__':
    sys.exit(main())
