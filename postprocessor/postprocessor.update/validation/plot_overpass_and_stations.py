from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from argparse import ArgumentParser
from validation.read_stations_list import stations_list
from data_reader.mirto_wrapper import MIRTO_multi_overpass_Wrapper
from matplotlib.patches import Ellipse
from obspy.geodetics.base import kilometers2degrees
import os, sys
from geometry.earth_geometry import haversine
from postprocessing.profiles.mirto_plot_tools import mirto_plot_validation_profiles
from datetime import datetime

EARTH_RADIUS = 6371 # km
SAVE_PLOT    = False

def plot_points(lon, lat,marker, color='red',label='',size=1):
    #plt.title("Land Surface Temperature - Copernicus")
    x,y = m(lon,lat)

    plt.scatter(x, y, s=size,
                marker=marker, 
                c=color,
                label=label,
                zorder=1)
    return

# Prepare logger
v_levels = ['debug', 'info', 'warning']
parser = ArgumentParser()
parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                    help='The level of verbosity of the software')    
parser.add_argument('--input_dir', '-i',
                    help='Input dir with overpasses')    
parser.add_argument('--validation_file', '-vf', 
                    help='netCDF with validation results')
parser.add_argument('--outdir', '-o',default='',
                    help='Save plot')
argv = parser.parse_args()

if argv.outdir != '':
    SAVE_PLOT = True
    outdir = argv.outdir
    
# Read stations
stations = stations_list()

# Read validation file
valfile = Dataset(argv.validation_file,'r')
station_synop = os.path.basename(argv.validation_file).split('_')[3]
station = stations.get(station_synop)

# Read overpasses
overpass_list = os.listdir(argv.input_dir)
tmp1    = [ argv.input_dir+x+'/amethyst/amethyst_results.nc' for x in overpass_list if os.path.isdir(argv.input_dir+x+'/amethyst')]
tmp2    = [ argv.input_dir+x+'/amethyst_results.nc' for x in overpass_list if os.path.isfile(argv.input_dir+x+'/amethyst_results.nc')]
overpass_list = tmp1 + tmp2
print(f"Debug Overpass list: {overpass_list}")
overpass_list = [ path for path in overpass_list if 'MIRTO_'+path.split('/')[-2] in valfile['retr_name'][:] ]

if overpass_list == []:
    sys.exit("No overpass found")
retrievals = MIRTO_multi_overpass_Wrapper(overpass_list)

# Compute width and heigth

"""
lon_0 = retrievals.longitude.mean()
lat_0 = retrievals.latitude.mean()
"""

lon_0 = 180
lat_0 = 90

xmin = retrievals.longitude.min() 
xmax = retrievals.longitude.max() 
ymin = retrievals.latitude.min() 
ymax = retrievals.latitude.max() 

width  = haversine(xmin,ymin,xmax,ymin)*1000
height = haversine(xmin,ymin,xmin,ymax)*1000

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

# Plot all stations
lat = []
lon = []
for i,xlat in enumerate(stations.LAT):
    lat.append(xlat)
    lon.append(stations.LON[i])
#x, y = prj(lon, lat)
ax = plt.gca()
plot_points(lon,lat,color='tab:blue',marker='s',label='stations',size=4)

# Plot FOVs
plot_points(retrievals.longitude, retrievals.latitude, color='orange',marker='o',label='FOVs',size=0.8)

# Plot selected station
plot_points(station.lon, station.lat, color='red',marker='s',label='selected station',size=4)
R = kilometers2degrees(100,radius=EARTH_RADIUS)
m.tissot(station.lon,station.lat,R,500,zorder=3,facecolor='none',edgecolor='red') # plot circle

# Plot good FOVs
plot_points(valfile['fov_lon'][:], valfile['fov_lat'][:], color='red',marker='o',label='selected FOVS',size=0.8)

# Save or show coverage plot   
plt.title('STATION: {} DATE: {}'.format(argv.validation_file.split('_')[-2],
                                        argv.validation_file.split('_')[-1].replace('.nc','')),
          size=8)
plt.legend(framealpha = 1.0)
if SAVE_PLOT:
    plt.savefig('{}/MIRTO_validation_{}_{}_coverage.png'.format(outdir,argv.validation_file.split('_')[-2],argv.validation_file.split('_')[-1].replace('.nc','')),
                dpi=500)
    plt.clf()
    print('Saved plot {}/MIRTO_validation_{}_{}_coverage.png'.format(outdir,argv.validation_file.split('_')[-2],argv.validation_file.split('_')[-1].replace('.nc','')))
else:
    plt.show()

# Plot profiles
best_profile_index = valfile['post_rms_rh'][:].argmin()
plot_vars = {'temp':'T', 'rh':'RH', 'water_vapour':'WV'}
var_list = [ var for var in list(valfile.variables) if 'prior' in var and var.replace('prior_','') in plot_vars.keys() ]
for var in var_list:
    if SAVE_PLOT:
        outdir = outdir
    else:
        outdir = None
        
    # Compute temporal difference
    print(f"Debug plot_overpass_and_stations.py: {retrievals.fov_original_index(best_profile_index)[1]}") 
    retr_date = datetime.strptime(retrievals.fov_original_index(best_profile_index)[1].replace('MIRTO_',''),'%Y%m%d_%H%M%S')
    sonde_date = datetime.strptime(argv.validation_file.split('_')[-1].replace('.nc',''),'%Y-%m-%dT%H:%M')
    deltat = int((retr_date - sonde_date).total_seconds()//60)
    deltat_str = repr(abs(deltat))+' min AS' if deltat > 0 else repr(abs(deltat))+' min BS'
    mirto_plot_validation_profiles(valfile[var][:],
                                   valfile[var.replace('prior','post')][:],
                                   valfile[var.replace('prior','sonde')][:],
                                   valfile['fovs_levels'][:],
                                   valfile['sonde_levels'][:],
                                   plot_vars[var.replace('prior_','')],
                                   best_profile_index = best_profile_index,
                                   sup_title_str = 'STATION: {} ( {} km, {} ) - DATE: {}'.format(argv.validation_file.split('_')[-2],
                                                                                              int(valfile['distance'][best_profile_index]),
                                                                                              deltat_str,
                                                                                              argv.validation_file.split('_')[-1].replace('.nc','')),
                                   outdir = outdir,
                                   file_string = "{}_{}".format(argv.validation_file.split('_')[-2], argv.validation_file.split('_')[-1].replace('.nc',''))
                                   )
