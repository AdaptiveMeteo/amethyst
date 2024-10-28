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

parser=argparse.ArgumentParser()
parser.add_argument('-vcf',
                    '--viirs_cloud_fraction',
                    type=str,
                    required=True,
                    help='VIIRS cloud fraction file')
parser.add_argument('-o',
                    '--output',
                    type=str,
                    required=False,
                    default=None,
                    help='Output Path')
parser.add_argument('-ov', 
                    '--overpass', 
                    type    =str, 
                    required=True,
                    default ='/galaxy/data/hawaii/cris/20201126_235326/', 
                    help    = 'Overpass path')
parser.add_argument('-w', 
                    '--wrf_path', 
                    type    =str, 
                    required=True,
                    default ='/galaxy/data/hawaii/cris/20201126_235326/wrfda/', 
                    help    = 'WRF output path')
parser.add_argument('-wa', 
                    '--world_area', 
                    type    =str, 
                    required=False,
                    choices=['pacific','arctic'],
                    default ='pacific', 
                    help    = 'World Area')
args = parser.parse_args()   

viirs_file = args.viirs_cloud_fraction
tif        = [ args.overpass + "/" + file for file in os.listdir(args.overpass) if '.tif' in file ][0]
fg_file    = args.wrf_path + "/fg"
post_file  = args.wrf_path + "/wrfvar_output"

# Read All Dataset
viirs_cloud_fraction = Dataset(viirs_file)
fg                   = Dataset(fg_file)
post                 = Dataset(post_file)

fg_cloud_fraction        = wrf.g_cloudfrac.get_cloudfrac(fg, meta = False)
#full_fg_cloud_fraction   = np.mean(fg_cloud_fraction,axis=0)
full_fg_cloud_fraction   = np.sum(fg_cloud_fraction,axis=0)
full_fg_cloud_fraction   = np.where(full_fg_cloud_fraction < 1, full_fg_cloud_fraction, 1)
post_cloud_fraction      = wrf.g_cloudfrac.get_cloudfrac(post, meta = False)
#full_post_cloud_fraction = np.mean(post_cloud_fraction,axis=0)
full_post_cloud_fraction = np.sum(post_cloud_fraction,axis=0)
full_post_cloud_fraction = np.where(full_post_cloud_fraction < 1, full_post_cloud_fraction, 1)

# Save Overpass Name
if args.overpass[-1] == "/":
    args.overpass = args.overpass[:-1]
overpass_title = os.path.basename(args.overpass)

figsize = (12,10)
ds = gdal.Open(tif)        

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
cs = osr.SpatialReference()
cs.ImportFromWkt(wgs84_wkt)
width = ds.RasterXSize
height = ds.RasterYSize
gt = ds.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5]
maxx = gt[0] + (gt[1]*width) - gt[1]*0.5
maxy = gt[3] - gt[5] * 0.5
#extent = [(minx, maxy), (maxx, maxy), (maxx, miny), (minx, miny)]
extent = (minx,maxx,miny,maxy)

parallels = np.arange(miny//1, maxy//1, 5)
meridians = np.arange(minx//1, maxx//1, 5)

image = georaster.MultiBandRaster(tif, 
#                                load_data=(minx, maxx,
#                                          miny, maxy), 
                                bands='all',
                                latlon=True)

if args.world_area == 'pacific':
    proj = ccrs.PlateCarree()
elif args.world_area == 'arctic':
    proj = ccrs.Stereographic(central_latitude=90.0, central_longitude=-150.0,
                              false_easting=0.0,     false_northing=0.0,
                              true_scale_latitude=60)

# Plot WRF Cloud Fraction before DA

#plt.imshow(image.r/255, extent=(minx, maxx,
#                                miny, maxy), alpha=0.98, interpolation = 'nearest')
plt.figure(figsize=figsize)
ax = plt.axes(projection=proj)
#ax.set_yticks(parallels)
#ax.set_yticklabels(parallels)
#ax.set_xticks(meridians)
#ax.set_xticklabels(meridians)
plt.imshow(image.r/255, alpha=1.0,extent=extent,interpolation='nearest')
ax.set_ylabel('Longitude (deg)',size=14,fontweight="bold")
ax.set_xlabel('Latitude (deg)',size=14,fontweight="bold")
ax.coastlines(color='yellow',resolution='10m')

projx = []
projy = []
for point in zip(fg['XLONG'][0][:,:].flatten(), fg['XLAT'][0][:,:].flatten()):
    tmp_x, tmp_y = proj.transform_point(point[0],point[1],ccrs.Geodetic())
    projx.append(tmp_x)
    projy.append(tmp_y)

fg_cloud_fraction = np.ma.masked_where( full_fg_cloud_fraction <= 0.5, full_fg_cloud_fraction)
projx = np.array(projx).reshape(fg_cloud_fraction.shape)
projy = np.array(projy).reshape(fg_cloud_fraction.shape)
mm = ax.pcolormesh( projx,
                    projy,
                   #fg['XLONG'][0][:,:],\
                    #fg['XLAT'][0][:,:],\
                   fg_cloud_fraction,\
                   vmin=0,\
                   vmax=1.0,\
                   #transform=proj,
                   cmap=cm.get_cmap('BuGn'))
plt.colorbar(mm,label="Cloud Fraction")
plt.xlim(minx,maxx)
plt.ylim(miny,maxy)

ofile = args.output + "/cloud_fraction_fg_{}.png".format( overpass_title )


plt.title("WRF Cloud Fraction Before DA",size=13,fontweight="bold")
if args.output != None:
    ofile = args.output + "/cloud_fraction_fg_{}.png".format( overpass_title )
    plt.savefig(ofile,dpi=300,bbox_inches='tight')
    print("Saved fg cloud fraction plot in {}".format(ofile))
else:
    plt.show()
plt.clf()

"""

# Plot WRF Cloud Fraction after DA
plt.figure(figsize=figsize)
ax = plt.axes(projection=proj)
ax.set_yticks(parallels)
ax.set_yticklabels(parallels)
ax.set_xticks(meridians)
ax.set_xticklabels(meridians)
ax.set_ylabel('Longitude (deg)',size=14,fontweight="bold")
ax.set_xlabel('Latitude (deg)',size=14,fontweight="bold")

plt.xlim(minx,maxx)
plt.ylim(miny,maxy)

plt.imshow(image.r/255, extent=(minx, maxx,
                                miny, maxy), alpha=0.98, interpolation = 'nearest')
if image.shape[2] > 1:
    plt.imshow(image/255,extent=tif.extent,alpha=1.0,interpolation='nearest')
else:
    plt.imshow(image,extent=tif.extent,alpha=1.0,cmap='hot')



ax.coastlines(color='yellow',resolution='10m')
post_cloud_fraction = np.ma.masked_where( full_post_cloud_fraction <= 0.5,full_post_cloud_fraction)
mm = ax.pcolormesh(post['XLONG'][0][:,:],\
                   post['XLAT'][0][:,:],\
                   post_cloud_fraction,\
                   vmin=0,\
                   vmax=1.0,\
                   transform=proj,cmap=cm.get_cmap('GnBu'),
 )
plt.colorbar(mm,label="Cloud Fraction")
plt.title("WRF Cloud Fraction After DA",size=13,fontweight="bold")
if args.output != None:
    ofile = args.output + "/cloud_fraction_post_{}.png".format( overpass_title )
    plt.savefig(ofile)
    print("Saved post da cloud fraction plot in {}".format(ofile))
else:
    plt.show()
plt.clf()
"""

viirs_cf = np.ma.masked_invalid(viirs_cloud_fraction["cloud_fraction"][:])
viirs_cloud_fraction_field = np.ma.masked_where( viirs_cf <= 0.5, viirs_cf)

#from geometry.earth_geometry import haversine,\
#                                    point_in_polygon,\
#                                    define_spatial_domain
#domain = define_spatial_domain(retrieval.latitude,retrieval.longitude)

# check if station is within domain  
#if not point_in_polygon([sonde.lon,sonde.lat],domain):

# Plot VIIRS Cloud Fraction
plt.figure(figsize=figsize)
ax = plt.axes(projection=proj)
#ax.set_yticks(parallels)
#ax.set_yticklabels(parallels)
#ax.set_xticks(meridians)
#ax.set_xticklabels(meridians)
plt.imshow(image.r/255, alpha=1.0,extent=extent,interpolation='nearest')
ax.set_ylabel('Longitude (deg)',size=14,fontweight="bold")
ax.set_xlabel('Latitude (deg)',size=14,fontweight="bold")
ax.coastlines(color='yellow',resolution='10m')



to_plot = ~ viirs_cloud_fraction_field.mask
plt.imshow(image.r/255, extent=extent, alpha=0.98, interpolation = 'nearest')
plt.xlim(minx,maxx)
plt.ylim(miny,maxy)
ax.coastlines(color='yellow',resolution='10m')
viirs_cloud_fraction_field[viirs_cloud_fraction_field.mask] = 0.0
mm = ax.pcolormesh( projx,
                    projy,
                    viirs_cloud_fraction_field,\
                    vmin=0,\
                    vmax=1.0,\
                    #transform=proj,
                    cmap=cm.get_cmap('Reds'))
plt.colorbar(mm,label="Cloud Fraction")
plt.title("VIIRS Cloud Fraction",size=13,fontweight="bold")
if args.output != None:
    ofile = args.output + "/cloud_fraction_viirs_{}.png".format(overpass_title)
    plt.savefig(ofile,dpi=300,bbox_inches='tight')
    print("Saved viirs cloud fraction plot in {}".format(ofile))
else:
    plt.show()
plt.clf()

"""
post_cloud_fraction = np.ma.masked_invalid(full_post_cloud_fraction)
fg_cloud_fraction = np.ma.masked_invalid(full_fg_cloud_fraction)


# Difference before DA
plt.figure(figsize=figsize)
ax = plt.axes(projection=proj)
ax.coastlines(color='yellow',resolution='10m')
ax.set_yticks(parallels)
ax.set_yticklabels(parallels)
ax.set_xticks(meridians)
ax.set_xticklabels(meridians)
ax.set_ylabel('Longitude (deg)',size=14,fontweight="bold")
ax.set_xlabel('Latitude (deg)',size=14,fontweight="bold")
plt.xlim(minx,maxx)
plt.ylim(miny,maxy)
plt.imshow(image.r/255, extent=(minx, maxx,
                                miny, maxy), alpha=0.98, interpolation = 'nearest')
mm = ax.pcolormesh(viirs_cloud_fraction['lons'][:,:],\
                   viirs_cloud_fraction['lats'][:,:],\
                   viirs_cf - fg_cloud_fraction,\
                   vmin=-1.0,\
                   vmax=1.0,\
                   transform=proj,cmap=cm.get_cmap('seismic'),
 )
plt.colorbar(mm,label="Cloud Fraction Difference")
plt.title("VIIRS - WRF Cloud Fraction (before DA)",size=13,fontweight="bold")
if args.output != None:
    ofile = args.output + "/cloud_fraction_difference_before_da_{}.png".format( overpass_title )
    plt.savefig(ofile)
    print("Saved viirs fg cloud fraction difference plot in {}".format(ofile))
else:
    plt.show()
plt.clf()

# Difference after DA
plt.figure(figsize=figsize)
ax = plt.axes(projection=proj)
ax.set_yticks(parallels)
ax.set_yticklabels(parallels)
ax.set_xticks(meridians)
ax.set_xticklabels(meridians)
ax.set_ylabel('Longitude (deg)',size=14,fontweight="bold")
ax.set_xlabel('Latitude (deg)',size=14,fontweight="bold")
ax.coastlines(color='yellow',resolution='10m')
plt.xlim(minx,maxx)
plt.ylim(miny,maxy)
plt.imshow(image.r/255, extent=(minx, maxx,
                                miny, maxy), alpha=0.98, interpolation = 'nearest')
mm = ax.pcolormesh(viirs_cloud_fraction['lons'][:,:],\
                   viirs_cloud_fraction['lats'][:,:],\
                   viirs_cf - post_cloud_fraction,\
                   vmin=-1.0,\
                   vmax=1.0,\
                   transform=proj,cmap=cm.get_cmap('seismic'),
 )
plt.colorbar(mm,label="Cloud Fraction Difference")
plt.title("VIIRS - WRF Cloud Fraction (after DA)",size=13,fontweight="bold")
if args.output != None:
    ofile = args.output + "/cloud_fraction_difference_after_da_{}.png".format( overpass_title )
    plt.savefig(ofile)
    print("Saved viirs post_da cloud fraction difference plot in {}".format(ofile))
else:
    plt.show()
plt.clf()

# Plot diff Histogram
dif_before =  (viirs_cf - fg_cloud_fraction ).reshape(fg_cloud_fraction.size)
dif_after = ( viirs_cf - post_cloud_fraction ).reshape(post_cloud_fraction.size)
plt.hist(dif_before,edgecolor='blue',label='before DA',align='mid',bins=np.arange(-1.0,1.1,0.1))
plt.hist(dif_after,color='wheat',edgecolor='orange',label='after DA',alpha=0.65,align='mid',bins=np.arange(-1.0,1.1,0.1))
plt.xlabel("Cloud Fraction Absolute Difference",size=14)
plt.ylabel("Histogram",size=14)
plt.xticks(np.arange(-1.0,1.1,0.1))
plt.legend()
plt.title("VIIRS - WRF Cloud Fraction",size=13,fontweight="bold")
if args.output != None:
    ofile = args.output + "/cloud_fraction_histogram_{}.png".format( overpass_title  )
    plt.savefig(ofile,dpi=300,bbox_inches='tight')
    print("Saved histogram plot in {}".format(ofile))
else:
    plt.show()
plt.clf()


# Plot Domain with GeoTIF 
plt.figure(figsize=figsize)
ax = plt.axes(projection=proj)
plt.xlim(minx,maxx)
plt.ylim(miny,maxy)
plt.imshow(image.r/255, extent=(minx, maxx,
                                miny, maxy), alpha=0.98, interpolation = 'nearest')
ax.coastlines(color='yellow',resolution='10m')
ax.set_yticks(parallels[1:])
ax.set_yticklabels(parallels[1:])
ax.set_xticks(meridians[1:])
ax.set_xticklabels(meridians[1:])
ax.set_ylabel('Longitude (deg)',size=14,fontweight="bold")
ax.set_xlabel('Latitude (deg)',size=14,fontweight="bold")

plt.title("DOMAIN",size=14,weight='bold')
if args.output != None:
    ofile = args.output + "/cloud_fraction_domain_{}.png".format(overpass_title)
    plt.savefig(ofile,dpi=300,bbox_inches='tight')
    print("Saved domain plot in {}".format(ofile))
else:
    plt.show()    
plt.clf()


n1_before = (dif_before <= 0.333).sum()
n1_after = (dif_after <= 0.333).sum()

n2_before = ((dif_before > 0.333) & (dif_before <= 0.666)).sum()
n2_after = ((dif_after > 0.333) & (dif_after <= 0.666)).sum()

n3_before = (dif_before > 0.666).sum()
n3_after = (dif_after > 0.666).sum()

print()
print(" Cloud Fraction Difference Occurrence (3 equidistant cloud fraction range )  ")
print("      BEFORE     AFTER      DIF")
print("N1 :    {}       {}         {}".format(n1_before,n1_after,n1_before-n1_after))
print("N2 :    {}       {}         {}".format(n2_before,n2_after,n2_before-n2_after))
print("N3 :    {}       {}         {}".format(n3_before,n3_after,n3_before-n3_after))



"""
