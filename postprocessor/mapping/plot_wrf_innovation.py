import matplotlib.cm
from matplotlib import colors
from matplotlib import colorbar
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
import numpy as np
from datetime import datetime, timedelta
import georaster
from osgeo import gdal, osr
import wrf as wrf
import argparse
import os
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
from wrf import getvar
import matplotlib.patches as mpatches
import cartopy.feature

def contourn_fmt(x):
    s = f"{x:.2f}"
    #if s.endswith("0"):
    #    s = f"{x:.0f}"
    return rf"{s}" if plt.rcParams["text.usetex"] else f"{s}"

def plot_colormesh(lons,lats,field,title,outfile,world_area,
                   clabel="Temperature Difference [K]",
                   vmin=-1,vmax=1,colormap='bwr', tr_file = None):
        
    plt.clf()

    # Create a figure
    fig = plt.figure(figsize=(12,6))
    # Set the GeoAxes to the projection used by WRF
    if world_area == 'arctic':
        #proj = ccrs.Stereographic(central_latitude=90.0, central_longitude=-150.0,
        #                            false_easting=0.0,     false_northing=0.0,
        #                            true_scale_latitude=60)
        myccrs=  ccrs.Stereographic(central_latitude=90.0, central_longitude=-165.0,
                         false_easting=0.0, false_northing=0.0,
                         true_scale_latitude=60)
        #myccrs=  ccrs.NorthPolarStereo(central_longitude=-165.0, standard_longitude=15.0)
        #myccrs = get_cartopy(field)

        print('MYCCRS: {}'.format(myccrs))

    elif world_area == 'pacific':
        globe = ccrs.Globe(semimajor_axis=6378137, flattening=1/298.257223563)
        myccrs=ccrs.PlateCarree(central_longitude=-160.0, globe=globe)

    #PaoloA 08Feb2023
    if world_area == 'pacific': 
        extent = [-167, -147, 10, 30]
    elif world_area == 'arctic': 
         extent = [-180, 180, 75, 90]
         extent=None
    ax = plt.axes(projection=myccrs)
    ax.set_extent([-180, 180, 75, 90], ccrs.PlateCarree())

    gc=ax.coastlines(resolution='10m',linewidth=.75,color='orange')
    gl=ax.gridlines(draw_labels=True, linewidth=.75, color='black')
    gl.bottom_labels = False
    gl.left_labels = False


    #PaoloA 08Feb2023
    #ax.add_feature(cartopy.feature.LAND)
    #ax.add_feature(cartopy.feature.OCEAN)
    #ax.add_feature(cartopy.feature.LAKES)
    #ax.add_feature(cartopy.feature.RIVERS)

    
    # Download and add the states and coastlines
    #states = NaturalEarthFeature(category="cultural", scale="50m",
    #                             facecolor="none",
    #                             name="admin_1_states_provinces_shp")
    #ax.add_feature(states, linewidth=.5, edgecolor="black")
    #ax.coastlines('50m', linewidth=0.8,color='darkgoldenrod')
    
    # Make the contour outlines and filled contours for the smoothed sea level
    # pressure.
    #CS = plt.contour(to_np(lons), to_np(lats), to_np(field), 8, colors="black",linewidths=np.arange(0,1.125,0.125)[::-1],
    #            transform=proj,extent=extent)
    #PaoloA 08Feb2023            transform=proj)
    #ax.clabel(CS, CS.levels, inline=True,  fontsize=5,fmt=contourn_fmt)

    normalize = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    mapper = cm.ScalarMappable(norm=normalize, cmap=colormap)

    print('Field max position: {}'.format(np.argmax(to_np(field), axis=None)))
    print('Lat : {}'.format(np.min(to_np(lats))))
    print('Lon : {}'.format(np.min(to_np(lons))))
    #print('Lat max: {}'.format(to_np(lats(np.argmax(to_np(field), axis=None)))))
    #print('Lon max: {}'.format(to_np(lons(np.argmax(to_np(field), axis=None)))))

    mm = ax.pcolormesh(to_np(lons),\
                       to_np(lats),\
                       to_np(field),\
                       vmin=vmin,\
                       vmax=vmax,\
                       cmap=colormap)
    #                   transform=myccrs,cmap=colormap)
    #mm = plt.contourf(to_np(lons), to_np(lats), to_np(field), 8, colors="black",linewidths=np.arange(0,1.125,0.125)[::-1],
                #transform=myccrs,extent=extent,latlon=True)
    plt.colorbar(mm,label=clabel)
    #ax.set_xticks(np.arange(round(lons.min(),0),round(lons.max(),0),5.0))
    #ax.set_yticks(np.arange(round(lats.min(),0),round(lats.max(),0),5.0))
    #ax.gridlines(color='black',linestyle='dotted')
        
    # Set the map bounds
    #ax.set_xlim(lons.min(),lons.max())
    #ax.set_ylim(lats.min(),lats.max())
    sm = cm.ScalarMappable(cmap=colormap, norm=normalize)
    sm.set_array([])
 
    # Add the gridlines
    #ax.gridlines(color="black", linestyle="dotted")
    plt.title(title)
    ax.grid(True,which='minor',c='black',linestyle='-',linewidth=2,alpha=1.0)   
    #PaoloA 08Feb2023
    #ax.set_background()
    
    # Plot TR file if present
    if tr_file != None:

        tr_nc_file = Dataset(tr_file,'r')
        
        
        sat_ape      =  [0.963, 0.897] # degrees
        sat_alt      =  824; # km
        earth_radius = 6371 # km
        tand         = lambda x: np.tan(np.radians(x))
        n_fovs = tr_nc_file['longitude'][:].size
        for fov_index,lon,lat in zip(np.arange(n_fovs),
                                     tr_nc_file['longitude'][:],
                                     tr_nc_file['latitude'][:]):
                """            
                # Compute axis and angle
                # Double check if axis are semi-axis or full-axis
                L1=( tand(sat_ape[0]+tr_nc_file['field_of_view'][fov_index]) - tand(tr_nc_file['field_of_view'][fov_index]))*sat_alt
                L2= tand(sat_ape[1])*sat_alt
                             
                if not tr_nc_file['azimuth_angle'][fov_index].mask:
                   angle = -90 - tr_nc_file['azimuth_angle'][fov_index]
                else:
                   angle = -90 
                
                #Transform lat and lon in cartesian coordinates
                projx1, projy1 = myccrs.transform_point(lon,lat,
                                                      ccrs.Geodetic())
                
                # Define Ellipse Face and Edge Colour
                color = None
                edgecolor = 'grey'
                facecolor = edgecolor
                zorder    = 1 
                
                # Draw the patch
                ax.add_patch(      mpatches.Ellipse(   xy=[projx1, projy1], 
                                                       width=L1*1000, 
                                                       height=L2*1000, 
                                                       angle=angle, 
                                                       edgecolor= edgecolor,
                                                       alpha=1,
                                                       facecolor= facecolor, 
                                                       transform=proj,
                                                       linewidth=0.4,
                                                       zorder=zorder)
                             )  
                """
                L1=( tand(sat_ape[0]+tr_nc_file['field_of_view'][fov_index]) - tand(tr_nc_file['field_of_view'][fov_index]))*sat_alt
                L2= tand(sat_ape[1])*sat_alt
                             
                if not tr_nc_file['azimuth_angle'][fov_index].mask:
                   angle = -90 - tr_nc_file['azimuth_angle'][fov_index]
                else:
                   angle = -90 
               
                #PaoloA 08Feb2023
                plt.scatter(lon,lat,marker='+',c='red')
                #ax.add_patch(      mpatches.Ellipse(   xy=[lon, lat], 
                #                                       width=L1*0.01, 
                #                                       height=L2*0.01, 
                #                                       angle=angle, 
                #                                       edgecolor= 'darkolivegreen',
                #                                       alpha=1,
                #                                       facecolor= 'grey', 
                #                                       linewidth=0.4,
                #                                       zorder=1)
                #)
                plt.text(lon,lat,repr(fov_index+1),size=4,c='darkolivegreen')
    plt.savefig(outfile,dpi=300,bbox_inches='tight')
    print("Saved plot {}".format(outfile))   
    
    return

parser=argparse.ArgumentParser()
parser.add_argument('-o',
                    '--output',
                    type=str,
                    required=False,
                    default=None,
                    help='Output Path')
parser.add_argument('-aw', 
                    '--assimilationwindow', 
                    type    =str, 
                    required=True,
                    default ='20201126_000000', 
                    help    = 'Assimilation Window')
parser.add_argument('-w', 
                    '--wrf_path', 
                    type    =str, 
                    required=True,
                    default ='/galaxy/data/hawaii/cris/20201126_235326/wrfda/', 
                    help    = 'WRF output path')
parser.add_argument('--tr', 
                    type    =str, 
                    required=False,
                    default=None,
                    help    = 'TR file path')
parser.add_argument('--world_area', 
                    type    =str, 
                    required=False,
                    default='pacific',
                    choices = ['pacific','arctic'],
                    help    = 'World Area')

args = parser.parse_args()   

fg_file    = args.wrf_path + "/fg"
post_file  = args.wrf_path + "/wrfvar_output"

# Read All Dataset
fg                   = Dataset(fg_file)
post                 = Dataset(post_file)

# Save Overpass Name
assimilationwindow_title = os.path.basename(args.assimilationwindow)

figsize = (12,10)

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)


dif = getvar(post,'temp',units='K') - getvar(fg,'temp',units='K')
lats, lons = fg['XLAT'][0][:], fg['XLONG'][0][:]
#vmin = (dif.mean() - 3*dif.std())
#vmax = (dif.mean() + 3*dif.std())
vmin = -2.0
vmax = 2.0
# Plot Temp
for level in [1,10,20]:
    title = "T Innovations (WRFVAR-FG) Level {}".format(level)
    outfile = '{}/T_innovation_L{:02d}_{}.png'.format(args.output,level,assimilationwindow_title)
    plot_colormesh(lons,lats,dif[level-1,:,:],title,
                   outfile, args.world_area,clabel = 'Temperature Difference [K]',vmin=vmin,vmax=vmax, tr_file = args.tr)

# Plot WV
#dif = getvar(post,'temp',units='K') - getvar(fg,'temp',units='K')

dif = post['QVAPOR'][0][:,:,:] - fg['QVAPOR'][0][:,:,:]
dif *= 1e3
#vmin = (dif.mean() - 3*dif.std())
#vmax = (dif.mean() + 3*dif.std())
vmin = -2.0
vmax = 2.0
for level in [1,10,20]:        
    title = "QVAPOR Innovations (WRFVAR-FG) Level {}".format(level)
    outfile = '{}/Q_innovation_L{:02d}_{}.png'.format(args.output,level,assimilationwindow_title)
    plot_colormesh(lons,lats,dif[level-1,:,:],
                   title,outfile,args.world_area, clabel = 'Water Vapor Difference [g/Kg]',vmin=vmin,vmax=vmax, tr_file = args.tr)

# Plot Relative Humidity
    
dif = getvar(post,'rh') - getvar(fg,'rh')
#vmin = (dif.mean() - 3*dif.std())
#vmax = (dif.mean() + 3*dif.std())
vmin = -15
vmax = 15
for level in [1,10,20]:        
    title = "RH Innovations (WRFVAR-FG) Level {}".format(level)
    outfile = '{}/RH_innovation_L{:02d}_{}.png'.format(args.output,level,assimilationwindow_title)
    plot_colormesh(lons,lats,dif[level-1,:,:],title,outfile, args.world_area,
                   clabel = 'Relative Humidity Difference [%]',vmin=vmin,vmax=vmax, tr_file = args.tr)
