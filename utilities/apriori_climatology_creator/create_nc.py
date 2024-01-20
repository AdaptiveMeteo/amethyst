#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 16:57:40 2023

@author: avalletti
"""
from netCDF4 import Dataset  
import numpy as np
import os
import pandas as pd
from netCDF4 import date2num
from datetime import datetime

def create_nc_file(GESDISC_Wrapper,  workdir = "/."):
    
    title = f"{GESDISC_Wrapper.time[0].year}{GESDISC_Wrapper.time[0].month}_{GESDISC_Wrapper.time[-1].year}{GESDISC_Wrapper.time[-1].month}_{GESDISC_Wrapper.variable}_GESDISC_DATA.nc" 
      
    var_bounds = {  
                  'T'   : {'min': 0, 'max': 1000},  # K
                  'O3'  : {'min': -50, 'max' : 50}, # log(kg/kg)
                  'H2O' : {'min': -50, 'max' : 50}  # log(kg/kg)
                 }
    ncfilepath = workdir+'/' + title
    if os.path.isfile(ncfilepath):
        os.system('rm {}'.format(ncfilepath))
    
    ncfile = Dataset(ncfilepath,mode="w",format="NETCDF4")
    
    coords_idxs = ncfile.createDimension('coords_idxs', len(GESDISC_Wrapper.latlons)) 
    time        = ncfile.createDimension('time',len(GESDISC_Wrapper.time))
    pres        = ncfile.createDimension('pres', len(GESDISC_Wrapper.pressure_levels))
    scalar      = ncfile.createDimension('scalar', 1)
    dim_lat     = ncfile.createDimension('lats', len(GESDISC_Wrapper.unique_latitude))
    dim_lon     = ncfile.createDimension('lons', len(GESDISC_Wrapper.unique_longitude))
    months      = ncfile.createDimension('months', 12)
    twodim      = ncfile.createDimension('tuple', 2)
    
    
    if GESDISC_Wrapper.variable == "Temperature":
        tag = "T"
    else:
        tag = GESDISC_Wrapper.variable
    # Define vars
    creation_date                 = ncfile.createVariable("creation_date",str,('scalar'))
    
    time                          = ncfile.createVariable("times",np.float64,('time'))
    time.units                    = 'seconds since 1970-01-01'
    time.long_name                = 'Dates in seconds'
    time.standard_name            = 'time'
    
    latlons                       = ncfile.createVariable("latlons", np.float32,('coords_idxs','tuple'))
    latlons.long_name             = "Couple of Coordinates: [latitude,longitude]"
    latlons.standard_name         = "Coordinates"
    
    latitudes                     = ncfile.createVariable("latitude", np.float32,('lats'))
    latitudes.long_name           = "Latitudes"
    latitudes.standard_name       = "Latitudes"
   
    longitudes                    = ncfile.createVariable("longitude", np.float32,('lons'))
    longitudes.long_name          = "Longitudes"
    longitudes.standard_name      = "Longitudes"
    
    pressure_levels               = ncfile.createVariable("pressure", np.float32, ('pres'))
    pressure_levels.units         = 'hPa'
    pressure_levels.long_name     = 'Pressure levels'
    pressure_levels.standard_name = 'Pressure levels'
    
    values                        = ncfile.createVariable(f"{tag}_values", np.float32, ('coords_idxs','time','pres'))
    values.units                  = GESDISC_Wrapper.unit
    values.long_name              = GESDISC_Wrapper.variable+" values"
    values.standard_name          = tag+"_values"
    
    precision                     = ncfile.createVariable(f"{tag}_precision", np.float32, ('coords_idxs','time','pres'))
    precision.units               = GESDISC_Wrapper.unit
    precision.long_name           = GESDISC_Wrapper.variable+" precision"
    precision.standard_name       = tag+"_precision"
    
    mm_values                     = ncfile.createVariable(f"mm_{tag}_values", np.float32, ('coords_idxs','months','pres'))
    mm_values.units               = GESDISC_Wrapper.unit
    mm_values.long_name           = GESDISC_Wrapper.variable+" values monthly means"
    mm_values.standard_name       = tag+"_values_monthly_means"
    mm_values.valid_min           = var_bounds[tag]['min']
    mm_values.valid_max           = var_bounds[tag]['max']

    mm_apriori_values                     = ncfile.createVariable(f"mm_apriori_{tag}_values", np.float32, ('coords_idxs','months','pres'))
    mm_apriori_values.units               = GESDISC_Wrapper.unit
    mm_apriori_values.long_name           = GESDISC_Wrapper.variable+" values monthly means"
    mm_apriori_values.standard_name       = tag+"_values_monthly_means"
    mm_apriori_values.valid_min           = var_bounds[tag]['min']
    mm_apriori_values.valid_max           = var_bounds[tag]['max']
    
    mm_precision                  = ncfile.createVariable(f"mm_{tag}_prec", np.float32, ('coords_idxs','months','pres'))
    mm_precision.units            = GESDISC_Wrapper.unit
    mm_precision.long_name        = GESDISC_Wrapper.variable+" precision monthly means"
    mm_precision.standard_name    = tag+"_prec_monthly_means"
    mm_precision.valid_min        = var_bounds[tag]['min']
    mm_precision.valid_max        = var_bounds[tag]['max']

    mm_apriori_precision                  = ncfile.createVariable(f"mm_apriori_{tag}_prec", np.float32, ('coords_idxs','months','pres'))
    mm_apriori_precision.units            = GESDISC_Wrapper.unit
    mm_apriori_precision.long_name        = GESDISC_Wrapper.variable+" precision monthly means"
    mm_apriori_precision.standard_name    = tag+"_prec_monthly_means"
    mm_apriori_precision.valid_min        = var_bounds[tag]['min']
    mm_apriori_precision.valid_max        = var_bounds[tag]['max']
    
    date_now              = datetime.now()
    date_now              = datetime.strftime(date_now,"%Y-%m-%d %H:%M:%S")
    creation_date[0]      = date_now
    
    times = date2num(GESDISC_Wrapper.time, time.units)
    time[:]           = times
    
    latlons[:]         = GESDISC_Wrapper.latlons
    latitudes[:]       = GESDISC_Wrapper.unique_latitude
    longitudes[:]      = GESDISC_Wrapper.unique_longitude
    pressure_levels[:] = GESDISC_Wrapper.pressure_levels
    values[:]          = GESDISC_Wrapper.values.values
    precision[:]         = GESDISC_Wrapper.precision.values
    mm_values[:]               = GESDISC_Wrapper.mm_values.values
    mm_apriori_values[:]       = GESDISC_Wrapper.mm_apriori_values.values
    mm_precision[:]    = GESDISC_Wrapper.mm_precision.values
    mm_apriori_precision[:]    = GESDISC_Wrapper.mm_ap_precision.values
       
    ncfile.close()
    
    return 
