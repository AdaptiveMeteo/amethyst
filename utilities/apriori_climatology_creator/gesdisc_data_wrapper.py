#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 13:51:45 2023

@author: avalletti
"""
import sys
from glob import glob
import h5py
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import xarray as xr
from create_nc import create_nc_file

MONTHS = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]

def centroid(vertexes):
     _x_list = [vertex[0] for vertex in vertexes]
     _y_list = [vertex[1] for vertex in vertexes]
     _len = len(vertexes)
     _x = sum(_x_list) / _len
     _y = sum(_y_list) / _len
     return(_x, _y)

def centroids(LATLONS):
    
    _centroids = [[],[]]
    for i in range(0,len(LATLONS[0][0])-1):
        for j in range(0,len(LATLONS[0])-1):        
            X,Y = centroid([[LATLONS[0][j][i],LATLONS[1][j][i]],
                            [LATLONS[0][j][i+1],LATLONS[1][j][i+1]],
                            [LATLONS[0][j+1][i],LATLONS[1][j+1][i]],
                            [LATLONS[0][j+1][i+1],LATLONS[1][j+1][i+1]]])
            _centroids[0].append(X)
            _centroids[1].append(Y)
    return np.array(_centroids)

def new_grid_coordinates(latmin,latmax,lonmin,lonmax,gap):
    LAT = np.arange(latmin,latmax,gap)
    LON = np.arange(lonmin,lonmax,gap)
    LATLONS = np.meshgrid(LAT,LON)

    CENTROIDS = centroids(LATLONS)
    return CENTROIDS

LATLONS = new_grid_coordinates(-82.5, 85, -180, 182.5, 2.5)
# TIME = np.unique([datetime(t.year, t.month, t.day) for t in all_GESDISC.time])

def found_neighbors(lat,lon,newgrid,k = 1,mode = "haversine"):
    from sklearn.neighbors import BallTree, KDTree

    olds = pd.DataFrame({"lat" : lat,"lon" : lon})
    news = pd.DataFrame({"lat" : newgrid[0],"lon" : newgrid[1]})
    
    if mode == "euclidean":
        kd = KDTree(news[["lat", "lon"]].values, metric='euclidean')
        distances, indices = kd.query(olds[["lat","lon"]], k = k)
    elif mode == "haversine":
        for column in news[["lat", "lon"]]:
            rad = np.deg2rad(news[column].values)
            news['{}_rad'.format(column)] = rad
        for column in olds[["lat", "lon"]]:
            rad = np.deg2rad(olds[column].values)
            olds['{}_rad'.format(column)] = rad
        ball = BallTree(news[["lat_rad", "lon_rad"]].values, metric='haversine')
        distances, indices = ball.query(olds[["lat_rad", "lon_rad"]].values, k = k)
    indices = [i[0] for i in indices]
    new_lat = news["lat"][indices]
    new_lon = news["lon"][indices]
    return np.array(indices), new_lat.values, new_lon.values

def get_means(masked_values,masked_precision,filt):
    v = masked_values[filt]
    p = masked_precision[filt]
    weights = 1/p**2
    v_mean, sum_of_weights  = np.ma.average(v,weights = weights, axis = 0, returned = True)

    # Precision as max daily precision
    #p_mean   = np.ma.max(p**2, axis = 0).data

    # Precision as sqrt of the summed weights
    p_mean   = 1/sum_of_weights

    return v_mean, np.ma.sqrt(p_mean)

def get_means_xr(xr_values,xr_precision):
    monthly_mean_values = xr_values.groupby("time.month").mean(dim="time")
    monthly_mean_precision = xr_precision.groupby("time.month").mean(dim="time")
    monthly_std_values  = xr_values.groupby("time.month").std(dim="time")
    p_mean  = monthly_mean_precision + monthly_std_values
    return monthly_mean_values, p_mean

def create_DataArray(data, coords,time,pressure):
    values = xr.DataArray(
        data=data,
        dims=["coords_indices", "time", "pressure"],
        coords={
            "coords_indices": range(len(coords)),
            "time": time,
            "pressure": pressure
        }
    )
    return values

class GESDISC_Wrapper(object):
    
    def __init__(self, filepath):
        
        tag         = filepath.split("/")[-1].split("-")[2].split("_")[0] 
        apriori_tag = tag +"-APriori"
        year        = int(filepath.split("/")[-1].split("-")[-1][4:8])
        nday        = int(filepath.split("/")[-1].split("-")[-1][9:12])
        x = h5py.File(filepath, "r")
        
        self.variable  = tag
        
        unique_coordinates = np.unique((LATLONS[0].flatten().tolist(), LATLONS[1].flatten().tolist()), axis=1).T
        pressure_levels = x["HDFEOS"]["SWATHS"][tag]["Geolocation Fields"]["Pressure"][:]
        pressure_levels[0] = 1100 # Fix first level (hPa)

        date = datetime(year = year, month = 1, day = 1) + timedelta(days = nday-1)
        
        TIME = np.array([date])
            
        num_time_steps = len(TIME)
        
        organized_values            = np.zeros((len(unique_coordinates), num_time_steps, len(pressure_levels)))*np.nan
        organized_precision         = organized_values.copy()
        organized_apriori_values    = organized_values.copy()
        organized_apriori_precision = organized_values.copy()

        time_idx = np.where(TIME == date)[0][0]
        lat = x["HDFEOS"]["SWATHS"][tag]["Geolocation Fields"]["Latitude"][:]
        lon = x["HDFEOS"]["SWATHS"][tag]["Geolocation Fields"]["Longitude"][:]
        indices, latitude, longitude = found_neighbors(lat,lon,LATLONS,k = 1,mode = "haversine")
        
        values            = x["HDFEOS"]["SWATHS"][tag]["Data Fields"]["L2gpValue"][:]
        precision         = x["HDFEOS"]["SWATHS"][tag]["Data Fields"]["L2gpPrecision"][:]
        apriori_values    = x["HDFEOS"]["SWATHS"][apriori_tag]["Data Fields"]["L2gpValue"][:]
        apriori_precision = x["HDFEOS"]["SWATHS"][apriori_tag]["Data Fields"]["L2gpPrecision"][:]
        precision[ precision<0 ] = apriori_precision[ precision < 0]

        quality           = x["HDFEOS"]["SWATHS"][tag]["Data Fields"]["Quality"][:]
        values[quality<1]            = 9999
        precision[quality<1]         = 9999
        apriori_values[quality<1]    = 9999
        apriori_precision[quality<1] = 9999
        
        values            = np.ma.masked_where(values == 9999, values)
        precision         = np.ma.masked_where(precision == 9999, precision)
        apriori_values    = np.ma.masked_where(apriori_values == 9999, apriori_values)
        apriori_precision = np.ma.masked_where(apriori_precision == 9999, apriori_precision)
        
        done_yet = []
        for i in indices:    
            if i in done_yet:
                continue
            filt = indices == i            
            v_mean,p_mean     = get_means(values,precision,filt)
            apv_mean,app_mean = get_means(apriori_values,apriori_precision,filt)
            
            organized_values[i, time_idx, :]           = v_mean
            organized_precision[i, time_idx, :]        = p_mean
            organized_apriori_values[i, time_idx, :]   = apv_mean
            organized_apriori_precision[i,time_idx, :] = app_mean
            
            done_yet.append(i)
            
        del(values, precision, done_yet, apriori_values, apriori_precision)
        
        self.latlons = unique_coordinates
        
        self.unique_latitude = np.unique(LATLONS[0])
        
        self.unique_longitude = np.unique(LATLONS[1])
        
        self.time = TIME
        
        self.pressure_levels = pressure_levels
        
        self.values = create_DataArray(organized_values, unique_coordinates, TIME, pressure_levels)

        self.precision = create_DataArray(organized_precision, unique_coordinates, TIME, pressure_levels)
        
        self.apriori_values = create_DataArray(organized_apriori_values, unique_coordinates, TIME, pressure_levels)
        
        self.apriori_precision = create_DataArray(organized_apriori_precision, unique_coordinates, TIME, pressure_levels)
        
        self.unit              = x['HDFEOS']['SWATHS'][tag]['Data Fields']['L2gpValue'].attrs['Units']
        
        # if tag == "H2O" or tag == "O3":
        #     logq  = np.log(organized_values)
        #     dlogq = organized_precision/organized_values
            
        #     self.logq              = create_DataArray(logq, unique_coordinates, TIME, pressure_levels)
        #     self.dlogq             = create_DataArray(dlogq, unique_coordinates, TIME, pressure_levels)
            
        del(organized_values, organized_precision, organized_apriori_values, organized_apriori_precision)
        
       
    def add_data(self, new_GESDISC_data):
        # we have to verify lat lon and pressures between files
        if self.variable == new_GESDISC_data.variable:
            self.time = np.concatenate((self.time,new_GESDISC_data.time))
            self.values            = xr.concat([self.values,new_GESDISC_data.values], "time")
            self.precision         = xr.concat([self.precision,new_GESDISC_data.precision], "time")
            self.apriori_values    = xr.concat([self.apriori_values,new_GESDISC_data.apriori_values], "time")
            self.apriori_precision = xr.concat([self.apriori_precision,new_GESDISC_data.apriori_precision], "time")
            # if self.variable == "H2O" or self.variable == "O3":
            #     self.logq          = xr.concat([self.logq,new_GESDISC_data.logq], "time")
            #     self.dlogq         = xr.concat([self.dlogq,new_GESDISC_data.dlogq], "time")
    
    def get_means_bymonth(self):

        self.mm_values,self.mm_precision  = get_means_xr(self.values,self.precision)
        self.mm_ap_values,self.mm_ap_precision  = get_means_xr(self.apriori_values,self.apriori_precision)
        # if self.variable == "H2O" or self.variable == "O3":
        #     self.mm_logq,_  = get_means_xr(self.logq,self.logq)
        #     self.mm_dlogq,_ = get_means_xr(self.dlogq,self.dlogq)
                        

def read_all_GESDISC_files(path,variable):
    # import matplotlib as plt
    filelist = glob(path+f"*{variable}*.he5")
    year       = [int(f.split("/")[-1].split("-")[-1][4:8]) for f in filelist]
    day        = [int(f.split("/")[-1].split("-")[-1][9:12]) for f in filelist]
    
    date = [datetime(year = year[i], month = 1, day = 1) + timedelta(days = day[i]-1) for i in range(0,len(day))]
    
    df = pd.DataFrame({'date':date,'filelist':filelist})
    df.sort_values("date",inplace = True,ignore_index=True)
    df.drop_duplicates("date", keep = "last",ignore_index=True,inplace = True)
    
    if len(filelist) != 0: 
        output = None
        # plt.figure(figsize = (20,20))
        for file in df["filelist"]: 
            print(datetime.now(),"Reading file: "+file)
            GESDISC = GESDISC_Wrapper(file)
            # plt.scatter(GESDISC.latitude,GESDISC.longitude)
            if output == None:
                output = GESDISC 
            else:
                output.add_data(GESDISC)
        
        # output.sort_by_date()
    else:
        sys.exit("Files not found!")

    return output 


def plot_monthly_profile(data,lat1,lon1):
    import matplotlib.pyplot as plt
    idx = np.where(np.logical_and(data.latlons[:,0]==lat1,data.latlons[:,1]==lon1))[0][0]
    
    plt.figure(figsize=(8,10))
    for m in range(0,12):
        month = MONTHS[m]
        plt.plot(data.values[idx,m,:],data.values["pressure"],label =month)
    plt.legend()
    plt.show()
    
    
# for i in range(0,len(LATLONS[0])):
    
#     lat1 = LATLONS[0][i]
#     lon1 = LATLONS[1][i]
#     plot_monthly_profile(Temperature, lat1, lon1)
def main():
    #path     = "/home/avalletti/MIRTO/GESDISC/data/"
    path = "/home/mirto/wrkdir/GESDISC_data/original_data/"
    outpath = "/home/mirto/wrkdir/GESDISC_data/"

    Temperature = read_all_GESDISC_files(path, "Temperature")
    H2O         = read_all_GESDISC_files(path, "H2O")
    O3          = read_all_GESDISC_files(path, "O3")
    
    Temperature.get_means_bymonth()
    H2O.get_means_bymonth()
    O3.get_means_bymonth()
    
    #print(Temperature)
    #print(O3)
    #print(H2O)
        
    create_nc_file(Temperature,outpath)
    create_nc_file(H2O,outpath)
    create_nc_file(O3,outpath)
    # medie mensili valori (logq H2O e O3) 
    # incertezza = sigma mese + valore incertezza medio

if __name__ ==  "__main__":
    
    sys.exit(main())
