from netCDF4 import Dataset
import wrf
import xarray as xr
import numpy as np

import cartopy.crs as crs

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.cm import get_cmap
from matplotlib import colors
from matplotlib import colorbar
import cartopy.feature as cfe

from wrf import (to_np, getvar, smooth2d, get_cartopy, latlon_coords, cartopy_xlim, cartopy_ylim)
import numpy.ma as ma


root_dir = '/home/adaptive/wrkdir/arctic/innovations'
fg = Dataset(root_dir+'/fg')
wout = Dataset(root_dir+'/wrfvar_output')
# t2 = wrf.getvar(nc, 'T', timeidx=wrf.ALL_TIMES)
fg_t2 = wrf.getvar(fg, 'T', timeidx=0) # extract 3rd time instance (t=2) - slow....
wo_t2 = wrf.getvar(wout, 'T', timeidx=0) 

print(fg_t2)
print(wo_t2)

delta=wo_t2[10,:]-fg_t2[10,:]
dd = wrf.to_np(delta)
vmax, vmin = dd.max(), dd.min()
#colormap='bwr'
#colormap='coolwarm'
colormap='RdBu'

# Get the latitude and longitude points (use original data, rather than any processed data)
lats, lons = wrf.latlon_coords(fg_t2)

# do masked-array on the lon_2d
lon2d_greater = ma.masked_greater(wrf.to_np(lons), -0.01)
lon2d_lesser = ma.masked_less(wrf.to_np(lons), 0)

# apply masks to other associate arrays: lat_2d
lat2d_greater = ma.MaskedArray(wrf.to_np(lats), mask=lon2d_greater.mask)
lat2d_lesser = ma.MaskedArray(wrf.to_np(lats), mask=lon2d_lesser.mask)
# apply masks to other associate arrays: delta 
delta_2d_greater = ma.MaskedArray(wrf.to_np(delta), mask=lon2d_greater.mask)
delta_2d_lesser = ma.MaskedArray(wrf.to_np(delta), mask=lon2d_lesser.mask)


delta.plot()
plt.savefig('TEST_DELTA.png',dpi=300,bbox_inches='tight')
plt.close()



# Get the cartopy mapping object (use original data, rather than any processed data)
cart_proj = wrf.get_cartopy(wo_t2)

data_crs = crs.PlateCarree()

# Create a figure
fig = plt.figure(figsize=(12,9))
# Set the GeoAxes to the projection used by WRF
ax = plt.axes(projection=cart_proj)

# Add coastlines
ax.coastlines('50m', linewidth=0.8)
ax.add_feature(cfe.NaturalEarthFeature('physical', 'antarctic_ice_shelves_lines', 
                                       '50m', linewidth=1.0, edgecolor='k', facecolor='none') )
#normalize = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
#mapper = cm.ScalarMappable(norm=normalize, cmap=colormap)
norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

# Plot contours
#plt.contour(lon2d_greater, lat2d_greater, delta_2d_greater, 10, vmin = -2, vmax = 2,
#                transform=data_crs, cmap=get_cmap("bwr"))
#plt.contour(lon2d_lesser, lat2d_lesser, delta_2d_lesser, 10, vmin = -2, vmax = 2,
#                transform=data_crs, cmap=get_cmap("bwr"))
plt.contourf(lon2d_greater, lat2d_greater, delta_2d_greater, 100, vmin = vmin, vmax = vmax, 
                transform=data_crs, cmap=colormap, norm=norm)
plt.contourf(lon2d_lesser, lat2d_lesser, delta_2d_lesser, 100, vmin = vmin, vmax = vmax,
                transform=data_crs, cmap=colormap, norm=norm)
#ax.contour(lons, lats, wo_t2, 10,
#                transform=crs.NorthPolarStereo(), cmap=get_cmap("Reds"))

# Add a color bar
cbar = plt.colorbar(ax=ax, shrink=.62)
cbar.set_label(fg_t2.units)

# Set the map limits.  Not really necessary, but used for demonstration.
ax.set_xlim(wrf.cartopy_xlim(fg_t2))
ax.set_ylim(wrf.cartopy_ylim(fg_t2))

# Add the gridlines
ax.gridlines(color="black", linestyle="dotted")

#plt.title('Δ'+fg_t2.description+'\n'+str(fg_t2.Time.values))
plt.title('Innovations for T @ Level 10\n'+str(fg_t2.Time.values))

print('DONE')

plt.savefig('TEST.png',dpi=300,bbox_inches='tight')


