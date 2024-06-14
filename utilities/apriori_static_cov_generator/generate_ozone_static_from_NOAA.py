#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 15:28:00 2024

@author: paoloscaccia
"""
import numpy as np
import pandas as pd
import os, sys
from progress.bar import Bar
from argparse import ArgumentParser
from netCDF4  import Dataset
from scipy.interpolate import interp1d

URL = { 
       "arctic"  : [ "https://gml.noaa.gov/aftp/data/ozwv/Ozonesonde/2_Field%20Projects/Barrow,%20AK/100%20m%20Average%20File/",
                     "https://gml.noaa.gov/aftp/data/ozwv/Ozonesonde/2_Field%20Projects/Fairbanks,%20AK/100%20m%20Average%20Files/" ],
       
       "pacific" : ["https://gml.noaa.gov/aftp/data/ozwv/Ozonesonde/4_Satellite%20Comparison/Hilo,%20Hawaii/Homogenized/"]
       }

parser = ArgumentParser()
parser.add_argument('--region','-r', type=str, required = True,choices=['pacific','arctic'], 
                    help='World region')
parser.add_argument('--outdir','-o', type=str, required = False,default='./',
                    help='Output folder')
parser.add_argument('--keep_csv',action='store_true', default = False, required = False,
                    help='Keep downloaded data.')
parser.add_argument('--pressure_file',type=str, default = None, required = False,
                    help='NetCDF file containing pressure grid.')
argv = parser.parse_args()
if argv.pressure_file is None:
       sys.exit("Case wihout the external pressure grid not yet implemented!")    
else:
    with Dataset(argv.pressure_file,'r') as external_source:
        pressure_grid = np.array(external_source["atmospheric_components"]["p"][:])
        if any(pressure_grid[1:] - pressure_grid[:-1] > 0): sys.exit("Pressure grid should be monotonically increasing")

out_df = pd.DataFrame( {'pressure' : pressure_grid})

print("CREATION OF OZONE STATIC COVARIANCE\n")
print("Region:         ",argv.region)
print("Outdir:         ",argv.outdir)
print("URL:            ","\n                 ".join(URL[argv.region]))
print("Pressure from:  ",argv.pressure_file)
print("Keep CSV:       ",argv.keep_csv)
print()

OUTDIR = argv.outdir

if not os.path.isdir(OUTDIR) and argv.keep_csv:
    os.system(f"mkdir -p {OUTDIR}")
    print("Created folder ",OUTDIR)
    
# Parse URL list for the given region
for isite, site in enumerate(URL[argv.region]):
    print(f"Scanning site {isite}...")
    df = pd.read_html(site)[0]
    
    # Init Progress bar
    bar = Bar("Downloading profiles: ",max=df.index.size)
    
    # Iterate on URL 
    for i,row in df.iterrows():
    
        try:
            # Parse URL for date info
            info = row["Name"].split('_')
            year,month,day,hour = info[1:5]
            hour = hour.split('.')[0]    
        except:
            continue
            
        # Download Data
        try:
            data = np.loadtxt(site + row["Name"],skiprows=29)
            O3_ppmv      = data[:,8]
            pressure_hPa = data[:,1]
        except:
            continue
        
        # Filter but data
        O3_ppmv[O3_ppmv > 90] = np.nan
        pressure_hPa[O3_ppmv > 90] = np.nan
        O3_ppmv = np.ma.masked_invalid(O3_ppmv)
        O3_kgkg      = O3_ppmv * 1.6571e-6
        pressure_hPa = np.ma.masked_invalid(pressure_hPa)
        
        # Drop Nan
        common_mask = np.logical_and( ~ O3_kgkg.mask, ~ pressure_hPa.mask)
        O3_kgkg = np.array(O3_kgkg[common_mask ]    )
        pressure_hPa = np.array(pressure_hPa[common_mask])
    
        # Keep single profiles if specified
        if argv.keep_csv:
            local_df = pd.DataFrame( {'O3 [kgkg]':O3_kgkg,'P [hPa]':pressure_hPa})
            local_df = local_df.dropna()
            outname = f"{OUTDIR}/O3_{year}_{month}_{day}_{hour}_site{isite}.csv"
            local_df.to_csv(outname)
            del(local_df)
    
        # Interpolate profiles
        f = interp1d(
                        np.log(pressure_hPa),
                        O3_kgkg,
                        kind='linear',
                        fill_value=np.nan,
                        bounds_error = False
                    )
        new_ozone = f(np.log(pressure_grid))
    
        # Update 
        label = f"{year}_{month}_{day}_{hour}_site{isite}"
        while label not in out_df:
            try:
                out_df.insert(0,label,new_ozone)
            except ValueError:
                label += "_1" 
        
        # Update Progress Bar
        bar.next()
    print()
print("Tot profiles: ",len(out_df.columns)-1)
print()

if argv.keep_csv:
    print("\nSingle profiles saved in ",OUTDIR)
    outfile =f"{OUTDIR}/O3_all_profiles.csv"
    out_df.to_csv(outfile,index=None)
    print("All profiles saved in ",outfile)

profiles_matrix = out_df.drop(columns=['pressure']).to_numpy()
profiles_matrix = np.ma.masked_invalid(profiles_matrix)
covariance = np.ma.cov( np.log(profiles_matrix))

results_nc = f"{OUTDIR}/ozone_static_apriori.nc"
with Dataset(results_nc,"w") as outfile:
    outfile.createDimension("nlev",covariance.shape[0])
    nc_oz = outfile.createVariable("ozone",np.float32,('nlev','nlev'),
                                   fill_value=1e9)
    nc_oz[:,:] = covariance
    nc_p = outfile.createVariable("pressure",np.float32,('nlev',),fill_value=1e9)
    nc_p[:] = pressure_grid
    
print("Results saved in ", results_nc)
print("Done!")
    
if os.path.isdir("./gml.noaa.gov"):
    os.system("rm -r gml.noaa.gov")
    print("Removed cache dir ./gml.noaa.gov")
    
    