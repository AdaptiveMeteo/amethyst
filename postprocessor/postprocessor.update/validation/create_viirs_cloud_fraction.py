#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 15:34:43 2021

@author: adrian
"""

import numpy as np
import os, sys
from   common.geometry.utilities.array_reshapers import array_1d, transform_index
import netCDF4 as nc
import wrf as wrf
from   mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm
from   matplotlib import colors
from multiprocessing import Process, Manager
from multiprocessing import Pool
from time import time
from datetime import datetime


def PARSER():
    
    import argparse

    parser=argparse.ArgumentParser()
    parser.add_argument('-ov', 
                        '--overpass', 
                        type    =str, 
                        default ='/galaxy/data/hawaii/cris/20201126_235326/', 
                        help    = 'Overpass path')
    parser.add_argument('-o', 
                    '--output', 
                    type    =str, 
                    default ='/galaxy/data/hawaii/validation/', 
                    help    = 'Overpass path')
    parser.add_argument('-w', 
                        '--wrf_path', 
                        type    =str, 
                        default ='/galaxy/data/hawaii/cris/20201126_235326/wrfda/', 
                        help    = 'WRF output path')
    parser.add_argument('-p',
                        '--picture_path',
                        type=str,
                        default='/galaxy/data/hawaii/validation/images/',
                        help='Path for output pictures')
    parser.add_argument('-n',
                        '--n_core',
                        type=int,
                        default=1,
                        help='Number of cores')
    
    args = parser.parse_args()     
    
    return args   

def get_wrf_cells_size(ncdataset):
    
    v2 = ncdataset['XLAT_U'][0]
    u2 = ncdataset['XLONG_U'][0]

    v3 = ncdataset['XLAT_V'][0]
    u3 = ncdataset['XLONG_V'][0]
    
    delta_lat =  np.abs(v3[1:,:] - v3[0:-1,:])
    delta_lon =  np.abs(u2[:,1:] - u2[:,0:-1])
        
    return delta_lon, delta_lat
    

class viirs_dataset_class(object):
    def __init__(self,filename):
        
        data = nc.Dataset(filename)
        self.lats      =  np.copy(data["Latitude"][:])
        self.longs     =  np.copy(data["Longitude"][:])
        self.cloudmask =  np.copy(data["CloudMask"][:])
        del(data)
        
def preprocess_geoloc_data(args):
    wrf_longs, wrf_lats, viirs_dataset_lon, viirs_dataset_lat, viirs_dataset_cm, wrf_delta_longs,wrf_delta_lats = args
            
    wrf_cell_position = np.array([ vector for vector in zip( wrf_longs,
                                                             wrf_lats )])
    wrf_cell_size     = np.array([ vector for vector in zip( wrf_delta_longs,
                                                             wrf_delta_lats )])
    
    viirs_lons = array_1d(viirs_dataset_lon)
    viirs_lats = array_1d(viirs_dataset_lat)
    
    viirs_shape = np.shape(viirs_dataset_lat)
    indices = [ [ ] for i in range(wrf_lats.size)] 

    for i, wrf_cell in enumerate(wrf_cell_position):
        
                # Check for cut at -180° -> 180°

                MAX_LON_DIFF = wrf_cell_size[i][0]
                MAX_LAT_DIFF = wrf_cell_size[i][1]

                lon_min = wrf_cell[0] - MAX_LON_DIFF/2.0
                lon_max = wrf_cell[0] + MAX_LON_DIFF/2.0
                lat_min = wrf_cell[1] - MAX_LAT_DIFF/2.0
                lat_max = wrf_cell[1] + MAX_LAT_DIFF/2.0
                            
                lon_ind = np.where(  (viirs_lons  >= lon_min)  & (viirs_lons <= lon_max)  )
                lat_ind = np.where(  (viirs_lats >= lat_min)  & (viirs_lats  <= lat_max) )
                
                indices[i] = [ transform_index(x,viirs_shape) for x in np.intersect1d(lon_ind,lat_ind)  
                                   if viirs_dataset_cm.mask[transform_index(x,viirs_shape)] == False ]  
    return indices


def merge_viirs_dataset(viirs_dataset,new_dataset ):
    viirs_dataset.longs = np.concatenate( (viirs_dataset.longs,
                                           new_dataset.longs ),
                                          axis = 0)

    viirs_dataset.lats  = np.concatenate( (viirs_dataset.lats,
                                           new_dataset.lats ),
                                          axis = 0)

    viirs_dataset.cloudmask = np.ma.concatenate( (viirs_dataset.cloudmask,
                                                  new_dataset.cloudmask ),
                                                  axis = 0 )
    if type(viirs_dataset.cloudmask.mask) == np.bool_:
                  viirs_dataset.cloudmask.mask = np.zeros_like(viirs_dataset.cloudmask,dtype=np.bool)

    return


def compute_viirs_cloudfraction_on_wrfgrid(index_array, viirs_dataset_cm, max_viirs_pixels):  

    THRESHOLD = 0.1 # Fraction of the maximum number of VIIRS pixels within 
                    # a wrf grid used as threshold for cloud fraction computation
    
    cloud_fraction    = np.empty(index_array.size)
    cloud_fraction[:] = np.nan
        
    for i,pixels_in_cell in enumerate(index_array):
        n_pixels = len(pixels_in_cell)
        if  n_pixels > THRESHOLD*max_viirs_pixels:
            counter = 0                   
            for viirs_fov in pixels_in_cell:             
                if viirs_dataset_cm[viirs_fov[0]][viirs_fov[1]] >= 2:
                    counter += 1
            cloud_fraction[i] = float(counter/n_pixels)

    return cloud_fraction

def new_compute_viirs_cloudfraction_on_wrfgrid(wrf_dataset, viirs_dataset, obs_filter):  

    THRESHOLD = 0.1 # Fraction of the maximum number of VIIRS pixels within 
                    # a wrf grid used as threshold for cloud fraction computation
    if obs_filter.shape != viirs_dataset.longs.shape:
        obs_filter = obs_filter.reshape(viirs_dataset.longs.shape)
        
    cloud_fraction    = np.zeros(wrf_dataset['XLONG'][0].shape,dtype=float)
    n_pixels          = np.zeros(wrf_dataset['XLONG'][0].shape,dtype=float)
    
    lat_indx, lon_indx = wrf.ll_to_xy(wrf_dataset,
                       viirs_dataset.lats,
                       viirs_dataset.longs,
                       meta=False)
    cm = viirs_dataset.cloudmask.flatten()
    
    for i,lon_i, lat_i in zip(range(lon_indx.size),lon_indx, lat_indx):
          n_pixels[lon_i,lat_i] += 1
          if cm[i] >= 2:
              cloud_fraction[lon_i,lat_i] += 1
              
    max_pixels = n_pixels.max()
    cloud_fraction = np.ma.masked_where( n_pixels < THRESHOLD*max_pixels, cloud_fraction)
    
    return np.ma.masked_invalid(cloud_fraction / n_pixels)

###################################### MAIN ####
if __name__ == "__main__":
        
    start = time()    
    args = PARSER()
        
    if not args.n_core and (202500 % args.n_core ) == 0:
               sys.exit("Invalid core number: it must be a divisor of 450x450!")

    fg_wrf_data   = nc.Dataset(args.wrf_path + "/fg") 
    post_wrf_data = nc.Dataset(args.wrf_path + "/wrfvar_output")
    
    minx = fg_wrf_data['XLONG'][0][:].min()
    miny = fg_wrf_data['XLAT'][0][:].min()
    maxx = fg_wrf_data['XLONG'][0][:].max()
    maxy = fg_wrf_data['XLAT'][0][:].max()
    
    maxx = max(minx,maxx)
    minx = min(minx,maxx)
    maxy = max(miny,maxy)
    miny = min(miny,maxy)
        
    wrf_delta_lon, wrf_delta_lat = get_wrf_cells_size(fg_wrf_data)

    viirs_list = [ "/".join([args.overpass,file]) for file in os.listdir(args.overpass) if "JRR-CloudMask_v2r0_" in file and '.nc' in file  ]    
    viirs_list.sort()
    
    # Merge CloudMask Files
    viirs_dataset = viirs_dataset_class(viirs_list[0])    
    if len(viirs_list) >  1:
        for viirs_file in viirs_list[1:]:
            new_viirs_dataset = viirs_dataset_class(viirs_file)
            viirs_miny = new_viirs_dataset.lats.min()
            viirs_maxy = new_viirs_dataset.lats.max()
            if (viirs_maxy >= maxy and viirs_miny  >= maxy) or \
               (viirs_miny <= miny and  viirs_maxy <= miny) :
                continue             # Filter Out Viirs Files outside WRF Domain
            else:                        
                merge_viirs_dataset( viirs_dataset, new_viirs_dataset)
            del(new_viirs_dataset)
    
    # Prepare data for indices selection
    n_wrf_points = int(fg_wrf_data['XLONG'][0].size // args.n_core)
    wrf_lons = [ array_1d(fg_wrf_data['XLONG'][0])[i*n_wrf_points:(i+1)*n_wrf_points]  for i in range(args.n_core) ]
    wrf_lats = [ array_1d(fg_wrf_data['XLAT'][0])[i*n_wrf_points:(i+1)*n_wrf_points]  for i in range(args.n_core) ]
    wrf_delta_lons = [ array_1d(wrf_delta_lon)[i*n_wrf_points:(i+1)*n_wrf_points] for i in range(args.n_core)]
    wrf_delta_lats = [ array_1d(wrf_delta_lat)[i*n_wrf_points:(i+1)*n_wrf_points] for i in range(args.n_core)]
    viirs_lons=[viirs_dataset.longs for i in range(args.n_core) ]
    viirs_lats=[viirs_dataset.lats for i in range(args.n_core) ]
    viirs_cm=[viirs_dataset.cloudmask for i in range(args.n_core) ]

    arg=zip(wrf_lons, wrf_lats,viirs_lons,viirs_lats, viirs_cm, wrf_delta_lons, wrf_delta_lats)
    
    print("Grouping viirs FOVs using {} cores...".format(args.n_core))
    pool = Pool(processes=args.n_core)
    multi_indices=pool.map(preprocess_geoloc_data,arg)
    indices = []
    for list in multi_indices:
        for tmp in list:
            indices.append(tmp)
    indices = np.array(indices)
    pool.close()

    print(indices.size)

    with nc.Dataset('./debug_viirs_cloud_fraction_indices.nc','w') as ncfile:
        ncfile.createDimension('n_indices',indices.size)
        outindx = ncfile.createVariable('indx',np.float64,('n_indices'))
        outindx[:] = indices
        print("Saved debug file ./debug_viirs_cloud_fraction_indices.nc")

    max_len = 0
    for l in indices:
        if len(l)>max_len:
            max_len = len(l)
    print("MAX VIIRS FOVS IN A WRF GRID CELL: {}".format(max_len))
    print("TIME FOR VIIRS FOVS SELECTION: {:.2f} s".format(time() - start))
    print("Computing viirs cloud fraction on wrf grid...")
    viirs_cloud_fraction = compute_viirs_cloudfraction_on_wrfgrid(indices,viirs_dataset.cloudmask,max_len)

    
    # Save NETCDF with viirs cloudfraction
    if args.overpass[-1] == '/':
        args.overpass = args.overpass[:-1]
    outnetcdf = args.output + "/viirs_cloud_fraction_{}.nc".format(os.path.basename(args.overpass))
    with nc.Dataset( outnetcdf, "w",format="NETCDF4") as ncfile:
        ncfile.createDimension('lon_size',fg_wrf_data['XLONG'][0].shape[0])
        ncfile.createDimension('lat_size',fg_wrf_data['XLONG'][0].shape[1])
        
        cf    = ncfile.createVariable('cloud_fraction',np.float64,('lon_size','lat_size'))
        nlats = ncfile.createVariable('lats',np.float64,('lon_size','lat_size'))
        nlons = ncfile.createVariable('lons',np.float64,('lon_size','lat_size'))
        
        nlons[:,:] = fg_wrf_data['XLONG'][0][:,:]
        nlats[:,:] = fg_wrf_data['XLAT'][0][:,:]
        try:
            #cf[:,:]    = viirs_cloud_fraction[:,:]
            cf[:,:]    = viirs_cloud_fraction.reshape(fg_wrf_data['XLAT'][0].shape)[:,:]

        except:
            cf[:,:]    = viirs_cloud_fraction.reshape(fg_wrf_data['XLAT'][0].shape)[:,:]
            
    print("Saved netcdf file {}".format(outnetcdf))
