import sys
import numpy as np
import logging
from datetime import datetime, timedelta
from utilities.atmos.atmos_tools import interpolate_pressure_grid
from utilities.data_reader.input_reader import InputData
from read_stations_list import stations_list
from download_wyoming_sondes import download_sondes
import getopt,os
from utilities.geometry.earth_geometry import haversine, point_in_polygon
from argparse import ArgumentParser
from data_reader.read_netcdf import netCDFReader

__author__ = "Paolo Scaccia, Paolo Antonelli and Stefano Piani"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__ = ["Paolo Scaccia", "Paolo Antonelli"]
__license__ = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

class DomainError(Exception):
    pass

def search_stations(stations, wrf_data):
    """
    INPUT
        - station:   stations_list class (defined in reads_station_list.py)
        - wrf_data:  InputData class (defined in input_reader.py)

    OUTPUT:
        - filtered_stations
    """
                                        
    # <----- IMPROVE SPATIAL DOMAIN USING GRID CORNER, NOT CENTER
                                        
    #define spatial domain    
    lon_min,lon_max = wrf_data.longitude.min(),wrf_data.longitude.max()
    lat_min,lat_max = wrf_data.latitude.min(),wrf_data.latitude.max()
    domain = [ (lon_min , lat_min ),    # bottom left
		        (lon_min , lat_max ),   # upper left
                (lon_max , lat_max ),   # upper right
		        (lon_max , lat_min),    # bottom right
	            (lon_min , lat_min )  ] # starting point
    
    out_indices = []
    for station_index, lon, lat in zip(range(len(stations.LON)),
                                       stations.LON,
                                       stations.LAT):
        # Check if station is within WRF domain  
        if point_in_polygon([lon,lat],domain):
            out_indices.append(station_index)
        
    return out_indices

class ValidationDataset(object):

    def __init__(self,temp,rh,wvmr,*arg,lat = None,lon = None):
        self.temp           = temp
        self.rh             = rh
        self.water_vapour   = wvmr
        self.pressure       = [] if len(arg) == 0 else arg[0]
        self.lat            = lat
        self.lon            = lon

def check_dates(wrf_times, MAX_MIN = 30):
    # wrf_time : np.array following according to wrapper 
    #            mirto.common.data_reader.input_reader
    # MAX_MIN : (optional) max temporal distance allowed [hours]
    wrf_times = np.array([ x.astype(datetime) for x in wrf_times])

    first_date = datetime(wrf_times[0].year,
                          wrf_times[0].month,
                          wrf_times[0].day,
                          0,0)
    last_date = datetime(wrf_times[-1].year,
                          wrf_times[-1].month,
                          wrf_times[-1].day,
                          12,0)
    n_steps = int((last_date - first_date).total_seconds()//43200) + 1
    ref_dates = [first_date + timedelta(hours=12*i) for i in range(n_steps)]

    out_dates   = []
    out_indices = []

    for date in ref_dates:
        dt = np.abs(wrf_times - date )
        if dt.min() <= timedelta(minutes=MAX_MIN):
            out_dates.append(date)
            out_indices.append(np.argmin(dt))
        del(dt)
    return np.array(out_dates), np.array(out_indices)

def WRF_SaveResults_nc(out_rms, wrf_data, sonde,
                       title, dists, workdir):
    
    from netCDF4 import Dataset
    
    ncfilepath = workdir+'/' + title
    if os.path.isfile(ncfilepath):
        os.system('rm {}'.format(ncfilepath))
    
    station = title.split('_')[4]
    forecast_hour = '0H' if 'analysis' in title else title.split('+')[-1].replace('.nc','')
    
    ncfile = Dataset(ncfilepath,mode="w",format="NETCDF4")

    n_fovs        = ncfile.createDimension('n_points',len(out_rms))
    nc_lev_dim    = ncfile.createDimension('n_lev',wrf_data.profiles[0].pressure.shape[1])
    nc_scalar_dim = ncfile.createDimension('scalar',1)

    # Define vars
    nc_sonde_lat = ncfile.createVariable('sonde_lat',np.float32,('scalar'))
    nc_sonde_lon = ncfile.createVariable('sonde_lon',np.float32,('scalar'))
    
    nc_pressure       = ncfile.createVariable('wrf_levels',np.float32,('n_points','n_lev'))
    nc_pressure.units = 'hPa'

    nc_sonde_pressure       = ncfile.createVariable('sonde_levels',np.float32,('n_lev',))
    nc_sonde_pressure.units = 'hPa'

    nc_dist           = ncfile.createVariable('distance',np.float32,('n_points',))
    nc_dist.units     = 'km'

    nc_fov_lat = ncfile.createVariable('wrf_lat',np.float32,('n_points',))
    nc_fov_lat.units     = 'north-degree'

    nc_fov_lon = ncfile.createVariable('wrf_lon',np.float32,('n_points',))
    nc_fov_lon.units     = 'east-degree'
    
    nc_temp  = ncfile.createVariable('temp',np.float32,('n_points','n_lev'))
    nc_temp.units = 'K'

    nc_sonde_temp       = ncfile.createVariable('sonde_temp',np.float32,('n_lev'))
    nc_sonde_temp.units = 'K'

    nc_rh         = ncfile.createVariable('rh',np.float32,('n_points','n_lev'))
    nc_rh.units = '%'

    nc_sonde_rh      = ncfile.createVariable('sonde_rh',np.float32,('n_lev',))
    nc_sonde_rh.units = '%'

    nc_wv      = ncfile.createVariable('water_vapour',np.float32,('n_points','n_lev'))
    nc_wv.units = 'g/kg'

    nc_sonde_wv      = ncfile.createVariable('sonde_water_vapour',np.float32,('n_lev',))
    nc_sonde_wv.units = 'g/kg'

    nc_rms_wv = ncfile.createVariable('rms_water_vapour',np.float32,('n_points','scalar'))
    nc_rms_rh = ncfile.createVariable('rms_rh',np.float32,('n_points','scalar'))
    nc_rms_temp = ncfile.createVariable('rms_temp',np.float32,('n_points','scalar'))

    nc_rms_profile_wv = ncfile.createVariable('rms_profile_water_vapour',np.float32,('n_points','n_lev'))
    nc_rms_profile_rh = ncfile.createVariable('rms_profile_rh',np.float32,('n_points','n_lev'))
    nc_rms_profile_temp = ncfile.createVariable('rms_profile_temp',np.float32,('n_points','n_lev'))

    nc_quality_control = ncfile.createVariable('quality_control','i1',('n_points',))
    nc_retrieval_name  = ncfile.createVariable('forecast',str,('n_points',))
    nc_id_string       = ncfile.createVariable('ID',str,('n_points',))
    nc_fov_index       = ncfile.createVariable('wrf_index',np.int64,('n_points',))

    nc_sonde_temp[:]     = sonde.temp
    nc_sonde_rh[:]       = sonde.rh
    nc_sonde_wv[:]       = sonde.water_vapour
    nc_sonde_pressure[:] = sonde.pressure
    nc_sonde_lat[:]      = sonde.lat
    nc_sonde_lon[:]      = sonde.lon

    wrf_profile  = wrf_data.profiles[0]
   
    for index, rms in enumerate(out_rms):

        fov_index = rms[-1]
        
        # FOVs vars
        nc_pressure[index,:]      = wrf_profile.pressure[fov_index,:]
        nc_dist[index]            = dists[fov_index]
        nc_fov_lat[index]         = wrf_data.latitude[fov_index]
        nc_fov_lon[index]         = wrf_data.longitude[fov_index]

        nc_temp[index,:]    = wrf_profile.temperature[fov_index,:]
        nc_rh[index,:]      = wrf_profile.rh[fov_index,:]     
        nc_wv[index,:]      = wrf_profile.water_vapour[fov_index,:]     
                  
        # Write Prior RMS
        nc_rms_temp[index,:] = rms[0].temp
        nc_rms_rh[index,:]   = rms[0].rh
        nc_rms_wv[index,:]   = rms[0].water_vapour
    
        # Write Prior RMS Profile
        nc_rms_profile_temp[index,:] = rms[1].temp
        nc_rms_profile_rh[index,:]   = rms[1].rh
        nc_rms_profile_wv[index,:]   = rms[1].water_vapour

        tmp = [ fov_index, wrf_data.name ]
        nc_fov_index[index]       = tmp[0]
        nc_retrieval_name[index]  = tmp[1]
        nc_id_string[index]       = '_'.join([station,tmp[1],forecast_hour])
        nc_quality_control[index] = 1
        if  np.isnan(rms[0].temp) or np.isinf(rms[0].rh) \
            or np.isnan(rms[0].water_vapour) or np.isinf(rms[0].water_vapour)   :
                nc_quality_control[index] = 0
    
    ncfile.close()

    return 
    


def WRF_Validation(wrf_data, sonde, station, 
                   dists, workdir, start_date, forecast_hour, timestep,
                   datatype = 'forecast'):

    out_rms = []
    index  = np.argmin(dists)
    dist = dists[index]
    
    wrf_profile = ValidationDataset(wrf_data.profiles[timestep].temperature[index,:],
                  wrf_data.profiles[timestep].rh[index,:],
                  wrf_data.profiles[timestep].water_vapour[index,:],
                  wrf_data.profiles[timestep].pressure[index,:]*1e-2)
    
    # merge pressure maps in a masked array   
    new_pressure_grid = np.ma.masked_invalid([ x if sonde.levels[-1]<x<sonde.levels[0]\
                                   else np.nan for x in wrf_profile.pressure ])
    
    # interpolated values of t,rh and mr
    int_sonde = ValidationDataset(
                        interpolate_pressure_grid(new_pressure_grid,
                                                  sonde.levels,sonde.temp),
                        interpolate_pressure_grid(new_pressure_grid,
                                                  sonde.levels,sonde.relh),
                        interpolate_pressure_grid(new_pressure_grid,
                                                  sonde.levels,sonde.mixr),
                        new_pressure_grid,
                        lat = sonde.lat,
                        lon = sonde.lon
                                )
    wrf_rms = ValidationDataset( 
                    np.sqrt(np.mean((int_sonde.temp - wrf_profile.temp)**2)),
                    np.sqrt(np.mean((int_sonde.rh   - wrf_profile.rh)**2)),
                    np.sqrt(np.mean((int_sonde.water_vapour - wrf_profile.water_vapour)**2))
                           )
    wrf_rms_profile = ValidationDataset( 
        np.sqrt(((int_sonde.temp - wrf_profile.temp)**2)),
        np.sqrt(((int_sonde.rh   - wrf_profile.rh)**2)),
        np.sqrt(((int_sonde.water_vapour - wrf_profile.water_vapour)**2)) 
                          )
    out_rms.append( [wrf_rms,  
                     wrf_rms_profile,
                     sonde.date,
                     index])
    
    #int( ( wrf_data.times[timestep[0]].astype(datetime) -\
    #                          start_date).total_seconds()//3600 )
           
    forecast_date_str = start_date.strftime('%Y%m%d%H%M')
            
    sonde_date_str    = sonde.date.replace(':','').replace('-','').replace('T','')
    title = 'WRF_Validation_{}_station_{}_r_{}_f_{}.nc'.format(datatype,
                                                       station.synop,
                                                       sonde_date_str,
                                                       forecast_date_str)
    if datatype != 'analysis':
        title = title.replace('.nc','+{}H.nc'.format(forecast_hour))
    WRF_SaveResults_nc(out_rms,
                        wrf_data,
                        int_sonde,
                        title,
                        dists,
                        workdir)
    print('Results saved in {}/{}'.format(workdir,title))
    
    return

def validate_single_station(station, dates,start_date, timestep_indx,
                            wrf_data, workdir, MAX_DIST = 100):
    """
    Parameters
    ----------
    station : station_class class (defined in read_stations_list.py)
    dates :   list containing rowinsonde dates
    timestep_indx : list containing timestep indices for each dates
    wrf_data : InputData class, WRF wrapped data (def. in input_reader.py)
    
    Returns
    -------
    None.

    """

    dist_f = lambda x : haversine(station.lon, station.lat,x[0],x[1] ) 
    dists  = np.ma.array(list(map(dist_f,
                               zip(wrf_data.longitude,wrf_data.latitude))
                              )
                         )
    #print('START DATE: ',start_date)
    # AGGIUNGERE CONTROLLO SU 
    
    for timestep, date in zip(timestep_indx, dates):
        sonde = download_sondes(station, date) # Download Rowinsonde
        if sonde == None:
            print("Warning: skipped station {}. "
                  "Rowinsonde data not found!".format(station.synop))
            continue


        forecast_hour = int( ( wrf_data.times[timestep].astype(datetime) -\
                               start_date).total_seconds()//3600 )
        print(wrf_data.times[timestep].astype(datetime), start_date,forecast_hour)
        WRF_Validation(wrf_data, sonde, station, 
                   dists, workdir, start_date, forecast_hour, timestep,
                   datatype = 'forecast' if date != start_date else 'analysis')
        
    return


def generate_statistics_file( outdir, timestep ):
    from netCDF4 import Dataset
    import re
    
    VALFILE_MASK = r'^WRF_Validation_(analysis|forecast)_station_\d+_r*'
    
    if timestep == 0:
        file_ext = 'analysis'
    elif timestep > 0 :
        file_ext = '+'+repr(timestep)+'H'
    else:
        file_ext = 'forecast'
    
    
    valfiles = [ outdir + '/' + file for file in os.listdir(outdir) if file_ext in file \
                                                                    and re.match(VALFILE_MASK,file) \
                                                                    and 'aggregated' not in file]
    
    if len(valfiles) == 0:
        print("No file found for the aggregation!")
        return
    
    for valfile in valfiles:
          if valfile == valfiles[0]:
              agg_nc_file = netCDFReader(valfile)
          else:
              tmp_val_file = netCDFReader(valfile)

              for var in tmp_val_file.__dict__.keys():        
                  if var not in ['n_points','n_lev','scalar'] \
                         and 'sonde' not in var:
                      agg_nc_file.__dict__[var] = np.ma.concatenate((agg_nc_file.__dict__[var],
                                                                  tmp_val_file.__dict__[var]))
                  elif 'sonde' in var:
                      agg_nc_file.__dict__[var] = np.ma.vstack((agg_nc_file.__dict__[var],
                                                                  tmp_val_file.__dict__[var]))
                      
              del(tmp_val_file)

    rms_temp_profile=np.sqrt(np.ma.mean(( agg_nc_file.sonde_temp - agg_nc_file.temp)**2,axis=0))
    rms_rh_profile=np.sqrt(np.ma.mean(( agg_nc_file.sonde_rh - agg_nc_file.rh)**2,axis=0))
    rms_wv_profile=np.sqrt(np.ma.mean((agg_nc_file.sonde_water_vapour-agg_nc_file.water_vapour)**2,axis=0))
        
    temp_rms = np.ma.mean( rms_temp_profile) 
    rh_rms = np.ma.mean( rms_rh_profile )
    wv_rms = np.ma.mean( rms_wv_profile)
    
    # Write NetCDF File
    
    outfile = outdir+'/WRF_Validation_{}_aggregated.nc'.format('analysis' if timestep == 0 else 'forecast')
    
    if 'forecast' in outfile and timestep > 0:
        outfile = outfile.replace('forecast','forecast+{}H'.format(timestep))
    elif 'forecast' in outfile and timestep < 0:
        outfile = outfile.replace('forecast','all_forecasts')

    if os.path.isfile(outfile):
        os.system('rm {}'.format(outfile))
    ncfile = Dataset(outfile,mode="w",format="NETCDF4")

    n_fovs        = ncfile.createDimension('n_cicles',len(valfiles))
    nc_lev_dim    = ncfile.createDimension('n_lev',agg_nc_file.sonde_levels.shape[1])
    nc_scalar_dim = ncfile.createDimension('scalar',1)

    # Define vars
    nc_sonde_lat = ncfile.createVariable('sonde_lat',np.float32,('n_cicles'))
    nc_sonde_lon = ncfile.createVariable('sonde_lon',np.float32,('n_cicles'))
    
    nc_pressure       = ncfile.createVariable('wrf_levels',np.float32,('n_cicles','n_lev'))
    nc_pressure.units = 'hPa'

    nc_sonde_pressure       = ncfile.createVariable('sonde_levels',np.float32,('n_cicles','n_lev',))
    nc_sonde_pressure.units = 'hPa'

    nc_fov_lat = ncfile.createVariable('wrf_lat',np.float32,('n_cicles',))
    nc_fov_lat.units     = 'north-degree'

    nc_fov_lon = ncfile.createVariable('wrf_lon',np.float32,('n_cicles',))
    nc_fov_lon.units     = 'east-degree'
    
    nc_temp  = ncfile.createVariable('temp',np.float32,('n_cicles','n_lev'))
    nc_temp.units = 'K'

    nc_sonde_temp       = ncfile.createVariable('sonde_temp',np.float32,('n_cicles','n_lev'))
    nc_sonde_temp.units = 'K'

    nc_rh         = ncfile.createVariable('rh',np.float32,('n_cicles','n_lev'))
    nc_rh.units = '%'

    nc_sonde_rh      = ncfile.createVariable('sonde_rh',np.float32,('n_cicles','n_lev',))
    nc_sonde_rh.units = '%'

    nc_wv      = ncfile.createVariable('water_vapour',np.float32,('n_cicles','n_lev'))
    nc_wv.units = 'g/kg'

    nc_sonde_wv      = ncfile.createVariable('sonde_water_vapour',np.float32,('n_cicles','n_lev',))
    nc_sonde_wv.units = 'g/kg'

    nc_rms_wv = ncfile.createVariable('rms_water_vapour',np.float32,('n_cicles','scalar'))
    nc_rms_rh = ncfile.createVariable('rms_rh',np.float32,('n_cicles','scalar'))
    nc_rms_temp = ncfile.createVariable('rms_temp',np.float32,('n_cicles','scalar'))

    nc_rms_profile_wv = ncfile.createVariable('rms_profile_water_vapour',np.float32,('n_cicles','n_lev'))
    nc_rms_profile_rh = ncfile.createVariable('rms_profile_rh',np.float32,('n_cicles','n_lev'))
    nc_rms_profile_temp = ncfile.createVariable('rms_profile_temp',np.float32,('n_cicles','n_lev'))

    nc_quality_control = ncfile.createVariable('quality_control','i1',('n_cicles',))
    nc_retrieval_name  = ncfile.createVariable('forecast',str,('n_cicles',))
    nc_id_string       = ncfile.createVariable('ID',str,('n_cicles',))

    nc_fov_index       = ncfile.createVariable('wrf_index',np.int64,('n_cicles',))

    nc_wv[:,:]   = agg_nc_file.water_vapour
    nc_rh[:,:]   = agg_nc_file.rh
    nc_temp[:,:] = agg_nc_file.temp
    
    nc_sonde_temp[:,:]     = agg_nc_file.sonde_temp
    nc_sonde_rh[:,:]       = agg_nc_file.sonde_rh
    nc_sonde_wv[:,:]       = agg_nc_file.sonde_water_vapour
    nc_sonde_pressure[:,:] = agg_nc_file.sonde_levels
    nc_sonde_lat[:]        = agg_nc_file.sonde_lat
    nc_sonde_lon[:]        = agg_nc_file.sonde_lon
    
    nc_pressure[:,:]  = agg_nc_file.wrf_levels
    nc_fov_lat[:]     = agg_nc_file.wrf_lat
    nc_fov_lon[:]     = agg_nc_file.wrf_lon
                      
    # Write RMS
    nc_rms_temp[:] =  temp_rms
    nc_rms_rh[:]   =  rh_rms
    nc_rms_wv[:]   =  wv_rms

    # Write RMS Profile
    nc_rms_profile_temp[:,:] = rms_temp_profile
    nc_rms_profile_rh[:,:]   = rms_rh_profile
    nc_rms_profile_wv[:,:]   = rms_wv_profile

    nc_fov_index[:]       = agg_nc_file.wrf_index
    nc_retrieval_name[:]  = agg_nc_file.forecast
    nc_id_string[:]       = agg_nc_file.ID
    
    nc_quality_control[:] = agg_nc_file.quality_control
    
    ncfile.close()
    
    print("Saved aggregated validation file {}".format(outfile))
    return

#####################################   MAIN  ##################  
def main(argv):
    parser = ArgumentParser()
    parser.add_argument('--verbose', '-v', default='info',
                    help='Verbose level')
    parser.add_argument('--outdir', '-o',default='./',
                    help='Output plot directory')
    parser.add_argument('--wrf_dir',required=True,
                    help='Output plot directory')
    parser.add_argument('--statistics',required=False,
                        default=False, type=bool,
                    help='Generate nc file with forecast statistics')
    parser.add_argument('--forecast_h','-fh',required=False,
                        default=-1, type=int,
                    help='Forecast hour for the statistics')
    parser.add_argument('--only_agg',required=False,
                    default=False, type=bool,
                help='Only aggregation Mode')


    print()
    print("Generating nc files for each station...")
    print()    
    argv = parser.parse_args()
    workdir   = argv.outdir
    verbose   = argv.verbose

    try:
        LOG.setLevel(verbose)
        streamhandler.setLevel(verbose)
    except:
        LOG.setLevel('INFO')
        streamhandler.setLevel('INFO')
    
    if not argv.only_agg:
        # Read all files in wrf_dir with 'wrfout' in the filename
        wrf_paths = [ argv.wrf_dir + '/' + x for x in os.listdir(argv.wrf_dir) if 'wrfout_d01' in x ]
        wrf_paths.sort()
        for wrf_path in wrf_paths:
            print('Reading ',wrf_path)
            # Wrap WRF Data
            wrf_data = InputData(wrf_path)

            if wrf_path == wrf_paths[0]:
                # Read all stations
                stations = stations_list()

                # Filter out those stations outside WRF domain
                stations_indx = search_stations(stations, wrf_data)
                print("Found {} stations within WRF domain".format(len(stations_indx)))

                # Save forecast start
                start_date = wrf_data.times[0].astype(datetime)
                #print('READ START DATE: ',start_date)
                #print('FILE: ',wrf_path)
                
            # Look for timesteps near rowinsonde dates
            sonde_dates, timestep_date_indx = check_dates(wrf_data.times)
            #print('WRF DATES ',wrf_data.times)
            #print('SONDE DATES ',sonde_dates)
        
            if len(sonde_dates) == 0:
                print("No Timestep available for WRF forecast {}".format(wrf_data.name))
                continue
            
            if stations_indx == []:
                return  
            
            valid_stations = [stations.get(stations.SYNOP[indx]) for indx in stations_indx]
            print(len(valid_stations))
            for station in valid_stations:
                validate_single_station(station, 
                                        sonde_dates, 
                                        start_date,
                                        timestep_date_indx,
                                        wrf_data,
                                        workdir)
          
    if argv.statistics or argv.only_agg:
            # AGGREGATE VALIDATION FILES
                    
            print()
            if argv.forecast_h != -1:
                print("Applyng statistics for forecast time +{}H...".format(argv.forecast_h ))
            else:
                print("Applyng statistics over all validation files...".format(argv.forecast_h ))
            print()    
            
            generate_statistics_file( argv.outdir,
                                      argv.forecast_h )
            
    return

# Global vars ###########
LOG = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - '
                              '%(funcName)s: %(message)s',
                               datefmt='%m/%d/%Y %H:%M:%S')
streamhandler = logging.StreamHandler()
streamhandler.setFormatter(formatter)

#########################
if __name__=="__main__":
    main(sys.argv[1:])
else:
    LOG = logging.getLogger(__name__)
    LOG.setLevel('INFO')
    LOG.addHandler(streamhandler)

