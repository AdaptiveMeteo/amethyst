from netCDF4 import Dataset
import wrf 
import numpy as np
import xarray as xr

import cartopy.crs as ccrs

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.cm import get_cmap
from matplotlib import colors
from matplotlib import colorbar

from datetime import datetime, timedelta
import georaster
from osgeo import gdal, osr
import argparse
import os
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
from wrf import getvar
import matplotlib.patches as mpatches
import cartopy.feature as cfe
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

from wrf import (to_np, getvar, smooth2d, get_cartopy, latlon_coords, cartopy_xlim, cartopy_ylim)
import numpy.ma as ma
import matplotlib.path as mpath

def contourn_fmt(x):
    s = f"{x:.2f}"
    #if s.endswith("0"):
    #    s = f"{x:.0f}"
    return rf"{s}" if plt.rcParams["text.usetex"] else f"{s}"

def plot_colormesh(lons,lats,field,projection,title,outfile,world_area,
                   clabel="Temperature Difference [K]",
                   vmin=-1,vmax=1, tr_file = None):

    colormap = 'RdBu_r'    
    if 'RH' in clabel:
       colormap = 'Greens'

    plt.clf()

    # Create a figure
    fig = plt.figure(figsize=(12,6))
    # Set the GeoAxes to the projection used by WRF
    if world_area == 'arctic':
        # Get the cartopy mapping object (use original data, rather than any processed data)
        myccrs = projection 

        # do masked-array on the lon_2d
        lon2d_greater = ma.masked_greater(wrf.to_np(lons), -0.01)
        lon2d_lesser = ma.masked_less(wrf.to_np(lons), 0)

        # apply masks to other associate arrays: lat_2d
        lat2d_greater = ma.MaskedArray(wrf.to_np(lats), mask=lon2d_greater.mask)
        lat2d_lesser = ma.MaskedArray(wrf.to_np(lats), mask=lon2d_lesser.mask)
        # apply masks to other associate arrays: delta
        field_2d_greater = ma.MaskedArray(wrf.to_np(field), mask=lon2d_greater.mask)
        field_2d_lesser = ma.MaskedArray(wrf.to_np(field), mask=lon2d_lesser.mask)


        print('MYCCRS: {}'.format(myccrs))

    elif world_area == 'pacific':
        globe = ccrs.Globe(semimajor_axis=6378137, flattening=1/298.257223563)
        myccrs=ccrs.PlateCarree(central_longitude=-160.0, globe=globe)
    #PaoloA 08Feb2023
    if world_area == 'pacific': 
        extent = [-167, -147, 10, 30]
    elif world_area == 'arctic': 
         extent = [-180, 180, 50, 90]
         #extent=None
   
    data_crs = crs.PlateCarree()

    ax = plt.axes(projection=myccrs)
    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)
    ax.set_extent([-180, 180, 50, 90], ccrs.PlateCarree())

    gc=ax.coastlines(resolution='10m',linewidth=.75,color='orange')
    #gl=ax.gridlines(draw_labels=True, linewidth=.75, color='black')
    #gl.bottom_labels = False
    #gl.left_labels = False
    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), lw=1, color="gray",
         y_inline=True, xlocs=range(-180,180,30), ylocs=range(50,90,5))
 
    gl.xlabel_style = {'size': 8, 'color': 'gray','rotation': 0, 'rotation_mode': 'anchor'}
    gl.ylabel_style = {'size': 8, 'color': 'gray','rotation': 0, 'rotation_mode': 'anchor'}
    #gl.xlabel_style = {'rotation': 0}
    #gl.ylabel_style = {'rotation': 0}

    gl.xlabel_style = {'size': 8, 'color': 'gray'}
   
    r_extent = 3331997.462
    r_extent *= 1.10      #increase a bit for better result

    # Prep circular boundary
    circle_path = mpath.Path.unit_circle()
    circle_path = mpath.Path(circle_path.vertices.copy() * r_extent,
                           circle_path.codes.copy())

    ax.set_boundary(circle_path)
    ax.set_frame_on(False)  #hide the boundary frame


    #gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    #              linewidth=0.5, color='black', alpha=1, linestyle='--')
    #gl.top_labels = False
    #gl.right_labels = False
    #gl.left_labels = True 
    #gl.xlines = True 
    #lon_range = 360
    #step = lon_range / 12  # Divide the range into 12 intervals
    #lon_pos = [-180 + step*interval for interval in range(13)]
    #gl.xlocator = mticker.FixedLocator(lon_pos)
    #lat_range = 30 
    #step = lat_range / 6  # Divide the range into 24 intervals
    #lat_pos = [60 + step*interval for interval in range(6)]
    #gl.ylocator = mticker.FixedLocator(lat_pos) 
    #gl.xformatter = LongitudeFormatter()
    #gl.yformatter = LatitudeFormatter()
    #gl.xlabel_style = {'size': 8, 'color': 'black','rotation': 0, 'rotation_mode': 'anchor'}
    #gl.ylabel_style = {'size': 8, 'color': 'black','rotation': 0, 'rotation_mode': 'anchor'}
    

    if vmin < 0.0:
         #norm = colors.TwoSlopeNorm(vmin=-extr, vcenter=0, vmax=extr)
         norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    else:  
         norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
         #norm = colors.Normalize(vmin=-extr, vmax=extr, clip=True)

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
    CS1 = plt.contour(lon2d_greater, lat2d_greater, field_2d_greater, np.linspace(vmin,vmax,5), colors="black", linewidths=0.6,transform=data_crs,extent=extent)
    ax.clabel(CS1, CS1.levels, inline=True,  fontsize=5,fmt=contourn_fmt)
    CS2 = plt.contour(lon2d_lesser, lat2d_lesser, field_2d_lesser, np.linspace(vmin,vmax,5), colors="black", linewidths=0.6, transform=data_crs,extent=extent)
    ax.clabel(CS2, CS2.levels, inline=True,  fontsize=5,fmt=contourn_fmt)

    normalize = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    #mapper = cm.ScalarMappable(norm=normalize, cmap=colormap)
    mapper = cm.ScalarMappable(norm=norm, cmap=colormap)

    if world_area == 'pacific':
        mm = ax.pcolormesh(to_np(lons),\
                       to_np(lats),\
                       to_np(field),\
                       vmin=vmin,\
                       vmax=vmax,\
                       transform=myccrs,cmap=colormap)
    elif world_area == 'arctic':
        mm = plt.contourf(lon2d_lesser, lat2d_lesser, field_2d_lesser, 100, vmin = vmin, vmax = vmax,
                transform=data_crs, cmap=colormap, norm=norm)
        mm = plt.contourf(lon2d_greater, lat2d_greater, field_2d_greater, 100, vmin = vmin, vmax = vmax,
                transform=data_crs, cmap=colormap, norm=norm)
        
    cbar = plt.colorbar(mm,label=clabel,format='%.2f')
#

    ticklabs = cbar.ax.get_yticklabels()
    cbar.ax.set_yticklabels(ticklabs, fontsize=8)

    #ax.set_xticks(np.arange(round(lons.min(),0),round(lons.max(),0),5.0))
    #ax.set_yticks(np.arange(round(lats.min(),0),round(lats.max(),0),5.0))
    #ax.gridlines(color='black',linestyle='dotted')
        
    # Set the map bounds
    #ax.set_xlim(lons.min(),lons.max())
    #ax.set_ylim(lats.min(),lats.max())
    #sm = cm.ScalarMappable(cmap=colormap, norm=normalize)
    sm = cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])
 
    # Add the gridlines
    #ax.gridlines(color="black", linestyle="dotted")
    plt.title(title,fontsize=10,loc='left')
    #ax.grid(True,which='minor',c='black',linestyle='-',linewidth=2,alpha=1.0)   
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
                #L1=( tand(sat_ape[0]+tr_nc_file['field_of_view'][fov_index]) - tand(tr_nc_file['field_of_view'][fov_index]))*sat_alt
                #L2= tand(sat_ape[1])*sat_alt
                             
                if not tr_nc_file['azimuth_angle'][fov_index].mask:
                   angle = -90 - tr_nc_file['azimuth_angle'][fov_index]
                else:
                   angle = -90 
               
                #PaoloA 08Feb2023
                try: 
                    plt.scatter(lon,lat,marker='.', s=.5, c='turquoise',alpha=0.75,transform=crs.PlateCarree())
                except Exception:
                    traceback.print_exc()

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
                #plt.text(lon,lat,repr(fov_index+1),size=4,c='darkolivegreen')
    try:
        plt.savefig(outfile,dpi=300,bbox_inches='tight')
        print("Saved plot {}".format(outfile))   
    except Exception:
        traceback.print_exc()

    plt.close()
    
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


lats, lons = fg['XLAT'][0][:], fg['XLONG'][0][:]

projection = get_cartopy(wrf.getvar(post, 'T', timeidx=0))

# Plot Temp Innovations
dif = getvar(post,'temp',units='K') - getvar(fg,'temp',units='K')
for level in [1,10,20]:
    title = "T Innovations (WRFVAR-FG) Level {}".format(level)
    outfile = '{}/T_innovation_L{:02d}_{}.png'.format(args.output,level,assimilationwindow_title)
    vmax, vmin = dif[level-1,:,:].max(), dif[level-1,:,:].min()
    plot_colormesh(lons,lats,dif[level-1,:,:],projection,title,
                   outfile, args.world_area,clabel = 'Temperature Difference [K]', vmin=vmin, vmax=vmax, tr_file = args.tr)

# Plot Potential Temp Innovations
dif = getvar(post,'th',units='K') - getvar(fg,'th',units='K')
for level in [1,10,20]:
    title = "$\Theta$ Innovations (WRFVAR-FG) Level {}".format(level)
    outfile = '{}/THETA_innovation_L{:02d}_{}.png'.format(args.output,level,assimilationwindow_title)
    vmax, vmin = dif[level-1,:,:].max(), dif[level-1,:,:].min()
    plot_colormesh(lons,lats,dif[level-1,:,:],projection,title,
                   outfile, args.world_area,clabel = 'Potential Temperature Difference [K]',vmin=vmin,vmax=vmax, tr_file = args.tr)


# Plot WV Innovations
dif = post['QVAPOR'][0][:,:,:] - fg['QVAPOR'][0][:,:,:]
dif *= 1e3
for level in [1,10,20]:        
    title = "QVAPOR Innovations (WRFVAR-FG) Level {}".format(level)
    outfile = '{}/Q_innovation_L{:02d}_{}.png'.format(args.output,level,assimilationwindow_title)
    vmax, vmin = dif[level-1,:,:].max(), dif[level-1,:,:].min()
    plot_colormesh(lons,lats,dif[level-1,:,:],projection,
                   title,outfile,args.world_area, clabel = 'Water Vapor Difference [g/Kg]',vmin=vmin,vmax=vmax, tr_file = args.tr)

# Plot Relative Humidity Innovations
dif = getvar(post,'rh') - getvar(fg,'rh')
for level in [1,10,20]:        
    title = "RH Innovations (WRFVAR-FG) Level {}".format(level)
    outfile = '{}/RH_innovation_L{:02d}_{}.png'.format(args.output,level,assimilationwindow_title)
    vmax, vmin = dif[level-1,:,:].max(), dif[level-1,:,:].min()
    plot_colormesh(lons,lats,dif[level-1,:,:],projection, title,outfile, args.world_area,
                   clabel = 'Relative Humidity Difference [%]',vmin=vmin,vmax=vmax, tr_file = args.tr)

# Plot FG Relative Humidity 
rh = getvar(fg,'rh')
for level in [1,10,20]:
    title = "RH(FG) Level {}".format(level)
    outfile = '{}/RH_FG_L{:02d}_{}.png'.format(args.output,level,assimilationwindow_title)
    vmax, vmin = rh[level-1,:,:].max(), rh[level-1,:,:].min()
    #vmax = 100.0
    #vmin = 0.0
    plot_colormesh(lons,lats,rh[level-1,:,:],projection, title,outfile, args.world_area,
                   clabel = 'Relative Humidity [%]',vmin=vmin,vmax=vmax, tr_file = args.tr)

# Plot FG Potential Vorticity
pv=wrf.g_vorticity.get_pvo(fg, timeidx=0)
for level in [1,10,20]:
    title = "PV(FG) Level {}".format(level)
    outfile = '{}/PV_FG_L{:02d}_{}.png'.format(args.output,level,assimilationwindow_title)
    vmax, vmin = pv[level-1,:,:].max(), pv[level-1,:,:].min()
    plot_colormesh(lons,lats,pv[level-1,:,:],projection, title,outfile, args.world_area,
                   clabel = 'Potential Vorticity [$m^2 s^{-1}K kg^{-1}$]',vmin=vmin,vmax=vmax, tr_file = args.tr)

# Plot WRFVAR_OUT  Potential Vorticity
pv=wrf.g_vorticity.get_pvo(post, timeidx=0)
for level in [1,10,20]:
    title = "PV(FG) Level {}".format(level)
    outfile = '{}/PV_POST_L{:02d}_{}.png'.format(args.output,level,assimilationwindow_title)
    vmax, vmin = pv[level-1,:,:].max(), pv[level-1,:,:].min()
    plot_colormesh(lons,lats,pv[level-1,:,:],projection, title,outfile, args.world_area,
                   clabel = 'Potential Vorticity [$m^2 s^{-1}K kg^{-1}$]',vmin=vmin,vmax=vmax, tr_file = args.tr)

# Plot Potential Vorticity Innovations
dif=wrf.g_vorticity.get_pvo(post, timeidx=0)-wrf.g_vorticity.get_pvo(fg, timeidx=0)
for level in [1,10,20]:
    title = "PV Innovations (WRFVAR-FG) Level {}".format(level)
    outfile = '{}/PV_innovation_L{:02d}_{}.png'.format(args.output,level,assimilationwindow_title)
    vmax, vmin = dif[level-1,:,:].max(), dif[level-1,:,:].min()
    plot_colormesh(lons,lats,dif[level-1,:,:],projection, title,outfile, args.world_area,
                   clabel = 'Potential Vorticity Difference[$m^2 s^{-1}K kg^{-1}$]',vmin=vmin,vmax=vmax, tr_file = args.tr)


