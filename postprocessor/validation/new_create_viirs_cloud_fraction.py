#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:10:35 2022

@author: paolo
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 15:34:43 2021

@author: adrian
"""

import numpy as np
import os, sys
import netCDF4 as nc
import wrf as wrf
from   matplotlib import colors
from multiprocessing import Process, Manager
from multiprocessing import Pool
from time import time
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')


def max_available_number_of_cores(max_n,n_pixels):
    return max_n if n_pixels%max_n==0 else max(1, max_available_number_of_cores(max_n-1,n_pixels))

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
    parser.add_argument('--max_cores',
                        type=int,
                        default=1,
                        help='Max Number of cores')
    
    args = parser.parse_args()     
    
    return args   

class viirs_dataset_class(object):
    def __init__(self,filename):
        
        data = nc.Dataset(filename)
        self.lats      =  np.copy(data["Latitude"][:])
        self.longs     =  np.copy(data["Longitude"][:])
        self.cloudmask =  np.copy(data["CloudMask"][:])
        del(data)
        
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
                  viirs_dataset.cloudmask.mask = np.zeros_like(viirs_dataset.cloudmask,dtype=bool)

    return

"""
def compute_viirs_cloudfraction_on_wrfgrid(wrf_dataset, viirs_dataset):  

    THRESHOLD = 0.05 # Fraction of the maximum number of VIIRS pixels within 
                    # a wrf grid used as threshold for cloud fraction computation
        
    cloud_fraction    = np.zeros(wrf_dataset['XLONG'][0].shape,dtype=float)
    n_pixels          = np.zeros(wrf_dataset['XLONG'][0].shape,dtype=float)
    
    lat_indx, lon_indx = wrf.ll_to_xy(wrf_dataset,
                                      viirs_dataset.lats,
                                      viirs_dataset.longs,
                                      meta=False)
    cm = viirs_dataset.cloudmask.flatten()
    
    for i, lon_i, lat_i in zip(range(lon_indx.size),lon_indx, lat_indx):
          n_pixels[lon_i,lat_i] += 1

          if cm[i] >= 2:
               cloud_fraction[lon_i, lat_i] += 1
               
    max_pixels = n_pixels.max()
    cloud_fraction = np.ma.masked_where( n_pixels < 15, cloud_fraction)
    
    return np.ma.masked_invalid(cloud_fraction / n_pixels)
"""

def read_viirs_JRR_files(overpass):
    
    # Read VIIRS files
    viirs_list = [ "/".join([overpass,file]) for file in os.listdir(overpass) if "JRR-CloudMask_v2r0_" in file and '.nc' in file  ]    
    viirs_list.sort()
    
    # Merge CloudMask Files
    viirs_dataset = viirs_dataset_class(viirs_list[0])    
    if len(viirs_list) >  1:
        for viirs_file in viirs_list[1:]:
            merge_viirs_dataset( viirs_dataset, viirs_dataset_class(viirs_file))

    return viirs_dataset

def compute_viirs_cloudfraction_on_wrfgrid(args):  
    wrf_lat_indx, wrf_lon_indx, wrf_grid_shape, viirs_flatten_cm = args

    cloud_fraction    = np.zeros(wrf_grid_shape, dtype=float)
    n_pixels          = np.zeros(wrf_grid_shape, dtype=float)
    
    print('BEFORE',(n_pixels == 0).sum())
    np.add.at( n_pixels, [wrf_lon_indx, wrf_lat_indx], 1 )
    print('AFTER',(n_pixels == 0).sum())

    wrf_lon_indx_cloudy = wrf_lon_indx[viirs_flatten_cm >= 2]    
    wrf_lat_indx_cloudy = wrf_lat_indx[viirs_flatten_cm >= 2]    
    np.add.at( cloud_fraction, [wrf_lon_indx_cloudy, wrf_lat_indx_cloudy], 1 )
                   
    return cloud_fraction, n_pixels

###################################### MAIN ####
if __name__ == "__main__":
        
    start = time()    
    args  = PARSER()
        
    fg_wrf_data   = nc.Dataset(args.wrf_path + "/fg") 
    post_wrf_data = nc.Dataset(args.wrf_path + "/wrfvar_output")
    
    # Read VIIRS CloudMask
    viirs_dataset = read_viirs_JRR_files(args.overpass)
    
    # Map VIIRS pixels to WRF grid cells    
    wrf_lat_indx, wrf_lon_indx = wrf.ll_to_xy(fg_wrf_data,
                                              viirs_dataset.lats,
                                              viirs_dataset.longs,
                                              meta=False)
    
    # Find maximum divisor for VIIRS pixels given the max number of available cores
    n_cores = max_available_number_of_cores(args.max_cores,viirs_dataset.lats.size)
    n_wrf_points = int(fg_wrf_data['XLONG'][0].size // n_cores)
    
    # Prepare the input for multiprocessing
    print("Preparing the input for multiprocessing...")
    cm = viirs_dataset.cloudmask.flatten()
    input_wrf_lon_indx = [ wrf_lon_indx[i*n_wrf_points:(i+1)*n_wrf_points]  for i in range(n_cores) ]
    input_wrf_lat_indx = [ wrf_lat_indx[i*n_wrf_points:(i+1)*n_wrf_points]  for i in range(n_cores) ]
    input_viirs_cm     = [ cm[i*n_wrf_points:(i+1)*n_wrf_points]  for i in range(n_cores) ]
    input_wrf_shape    = [ fg_wrf_data['XLAT'][0].shape for i in range(n_cores)]
    multiproc_args     = zip(input_wrf_lon_indx, input_wrf_lat_indx, input_wrf_shape, input_viirs_cm)

    # Cloud Fraction MultiProcessing
    print("Computing viirs cloud fraction on wrf grid with {} cores...".format(n_cores))
    pool = Pool(processes = n_cores)
    output = pool.map( compute_viirs_cloudfraction_on_wrfgrid,\
                                                multiproc_args)
    print('OUTSIDE')
    print( (output[0][0] == 0).sum(), (output[1][0] == 0).sum())
    print( (output[0][1] == 0).sum(), (output[1][1] == 0).sum())

    cloud_frac_batches, n_pixels_batches = [ o[0] for o in output], [ o[1] for o in output]
    print( (n_pixels_batches[0] == 0).sum(), (n_pixels_batches[1] == 0).sum())

    #viirs_cloud_fraction    = np.zeros(fg_wrf_data['XLAT'][0].shape, dtype=float)
    #n_pixels          = np.zeros(fg_wrf_data['XLAT'][0].shape, dtype=float)
    n_pixels  = np.sum(n_pixels_batches, axis=0)
    viirs_cloud_fraction  = np.sum(cloud_frac_batches, axis=0)

    # Merge the output
    #for cf, npix in zip(cloud_frac_batches, n_pixels_batches):
    #    viirs_cloud_fraction += cf
    #    n_pixels += npix
        
    pool.close()
    del(cloud_frac_batches)    
    del(n_pixels_batches)
    print((n_pixels == 0).sum())
    viirs_cloud_fraction = np.ma.masked_invalid( viirs_cloud_fraction / n_pixels )
    #viirs_cloud_fraction.mask[ n_pixels < 15 ] = True
    
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
