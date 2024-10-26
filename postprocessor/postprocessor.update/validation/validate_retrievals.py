import sys
import numpy as np
import logging
from datetime import datetime, timedelta

sys.path.append('/home/amethyst/amethyst/postprocessor')
from atmos.atmos_tools import interpolate_pressure_grid

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

class ValidationDataset(object):

    def __init__(self,temp,rh,wvmr,*arg,lat = None,lon = None):
        self.temp           = temp
        self.rh             = rh
        self.water_vapour   = wvmr
        self.pressure       = [] if len(arg) == 0 else arg[0]
        self.lat            = lat
        self.lon            = lon

def validate_retrieval_single_station(sonde,retrieval,station,workdir,MAX_DIST = 75):
    """
    INPUT
        - sonde:      raob_data class (defined in download_wyoming_sondes.py)
        - retrieval:  InputData class (defined in tools/input_reader.py)

    OUTPUT:
        - rms

    """
    from geometry.earth_geometry import haversine,\
                                        point_in_polygon,\
                                        define_spatial_domain

    #define spatial domain
    domain = define_spatial_domain(retrieval.latitude,retrieval.longitude)

    # check if station is within domain  
    if not point_in_polygon([sonde.lon,sonde.lat],domain):
        raise DomainError('Station outside retrieval\'s region:\n'
                         'ret domain (lat lon): \n'
                         '                      {} {}\n'
                         '                      {} {}\n'
                         '                      {} {}\n'
                         '                      {} {}\n'
                         'sonde lat lon: \n'
                         '                      {} {}'.format(
                                                       domain[0][1],domain[0][0],
                                                       domain[1][1],domain[1][0],
                                                       domain[2][1],domain[2][0],
                                                       domain[3][1],domain[3][0],
                                                       sonde.lat,sonde.lon)
                                                             )
        
    type_of_file = repr(type(retrieval))

    if 'SPS' in type_of_file:
        # find closest good FOV (quality_and_completeness == 1)
        sys.exit()
    elif 'MIRTO' in type_of_file:
        # find closest good FOV
        #dists = np.ma.masked_invalid( [ haversine(sonde.lon,sonde.lat,retrieval.longitude[i],retrieval.latitude[i]) \
        #                                 if retrieval.d2[i]<5 else np.nan      \
        #                               for i in range(0,len(retrieval.longitude))])
        
        dist_func = lambda x : haversine(sonde.lon, sonde.lat,x[0],x[1] ) 
        dists   = np.array(list(map(dist_func, zip(retrieval.longitude,retrieval.latitude))))
        dists   = np.ma.masked_where( (dists > MAX_DIST) & (retrieval.d2<5) , dists)
        indices = np.ma.masked_array(np.arange(retrieval.longitude.size))
        indices.mask = np.copy(dists.mask)
        
        # old version
        # dists = np.array( [haversine(sonde.lon,sonde.lat,rextrieval.longitude[i],retrieval.latitude[i]) \
        #            for i in range(0,len(retrieval.longitude))  ])
            
    if dists.mask.sum() == dists.size  :
        raise ValueError('No available retrieval')    

    # select what you are  validating
    if 'SPS' in type_of_file:  
        sys.exit('File list contains not only MIRTO files')
    elif 'MIRTO' in type_of_file:
        MIRTO_Validation(retrieval,sonde,station,dists,indices,workdir)
    elif 'WRF' in type_of_file:
        sys.exit('File list contains not only MIRTO files')

    return

def validate_retrieval_multi_station(retrievals, stations, dists, sonde_dates, workdir):
    """
    INPUT
        - retrieval:  InputData class (defined in tools/input_reader.py)
        - stations:   stations_list class (defined in read_stations_list.py)
        
    OUTPUT:
        - rms

    """
    # Convert sonde temperature from C to K
    # sonde.temp = sonde.temp + T_POINT_TEMP
    
    type_of_file = repr(type(retrievals))


    # select what you are  validating
    if 'SPS' in type_of_file:   
        sys.exit('File list contains not only MIRTO files')
    elif 'MIRTO' in type_of_file:
        MIRTO_MultiStation_Validation(retrievals, stations, dists, sonde_dates, workdir)
    elif 'WRF' in type_of_file:
        #WRF_Validation(retrieval,sonde,station,dists,min_index,workdir)
        sys.exit('File list contains not only MIRTO files')
    return

def SPS_Validation(retrieval,sonde,station,dists,min_index,workdir):
    from atmos.atmos_tools import interpolate_pressure_grid
    from atmos.mirto_atmos_tools import mr2rh
    
    prior = ValidationDataset(retrieval.profiles[0].temperature[min_index],
                              mr2rh(retrieval.profiles[0].pressure[min_index],
                                    retrieval.profiles[0].temperature[min_index],
                                    retrieval.profiles[0].water_vapour[min_index])[0],
                              retrieval.profiles[0].water_vapour[min_index],
                              retrieval.profiles[0].pressure[min_index])

    post  = ValidationDataset(retrieval.profiles[1].temperature[min_index],
                              mr2rh(retrieval.profiles[1].pressure[min_index],
                                    retrieval.profiles[1].temperature[min_index],
                                    retrieval.profiles[1].water_vapour[min_index])[0],
                              retrieval.profiles[1].water_vapour[min_index],
                              retrieval.profiles[1].pressure[min_index])


    # merge pressure maps in a masked array   
    new_pressure_grid = np.ma.masked_invalid([ x if sonde.levels[-1]<x<sonde.levels[0]\
                                       else np.nan for x in prior.pressure ])
    
    # interpolated values of t,rh and mr
    int_sonde = ValidationDataset(
                        interpolate_pressure_grid(new_pressure_grid,sonde.levels,sonde.temp),
                        interpolate_pressure_grid(new_pressure_grid,sonde.levels,sonde.relh),
                        interpolate_pressure_grid(new_pressure_grid,sonde.levels,sonde.mixr),
                        lat = sonde.lat,
                        lon = sonde.lon
                )
    
    if SPS_SaveResults(prior,
                   post,
                   int_sonde,
                   station.icao+'_'+sonde.date,
                   retrieval.name,
                   dists[min_index],
                   workdir):
        LOG.info('Results saved in {}/{}_{}_MIRTO_validation.nc'.format(workdir,station.icao,sonde.date))
    else:
        LOG.fatal('An error occurred in saving results!')

    return

def SPS_SaveResults(prior,post,sonde,string,name,dist,workdir):
    
    outfile = open(workdir+'/'+string+'_validation.txt','a')
#    LOG.info(prior.temp)
    print('===================================================',file=outfile)   
    print('      {}     '.format(name),file=outfile)
    print('===================================================',file=outfile)   
    
    print('\n===== RET INFO',file=outfile)
    print('dist: {0:.2f} km'.format(dist),file=outfile)    
    
    print('\n===== TEMPERATURE',file=outfile)    
    print('pressure prior post sonde ',file=outfile)
    for p,r,s,t,control in zip(prior.pressure,prior.temp,post.temp,sonde.temp,sonde.temp.mask):
        if control:
            print('{:.4}  {:.4} {:.4} --'.format(p,r,s),file=outfile)    
        else:
            print('{:.4}  {:.4} {:.4} {:.4}'.format(p,r,s,t),file=outfile)                
    print('\n===== REL. HUMIDITY',file=outfile)    
    print('pressure prior  post  sonde ',file=outfile)
    for p,r,s,t,control in zip(prior.pressure,prior.rh,post.rh,sonde.rh,sonde.rh.mask):
        if control:
            print('{:.4}  {:.4} {:.4}  --'.format(p,r,s),file=outfile)    
        else:
            print('{:.4}  {:.4} {:.4}  {:.4}'.format(p,r,s,t),file=outfile)                

    print('\n===== WATER VAPOR',file=outfile)    
    print('pressure prior post sonde ',file=outfile)
    for p,r,s,t,control in zip(prior.pressure,prior.water_vapour,post.water_vapour,sonde.water_vapour,sonde.water_vapour.mask):
        if control:
            print('  {:.4} {:.4} {:.4} --'.format(p,r,s),file=outfile)    
        else:
            print('  {:.4} {:.4} {:.4} {:.4}'.format(p,r,s,t),file=outfile)
            
    prior_rms = ValidationDataset( 
                        np.sqrt(np.mean((sonde.temp - prior.temp)**2)),
                        np.sqrt(np.mean((sonde.rh   - prior.rh)**2)),
                        np.sqrt(np.mean((sonde.water_vapour - prior.water_vapour)**2))
                )

    post_rms = ValidationDataset( 
                        np.sqrt(np.mean((sonde.temp  - post.temp)**2)),
                        np.sqrt(np.mean((sonde.rh   - post.rh)**2)),
                        np.sqrt(np.mean((sonde.water_vapour - post.water_vapour)**2))
                )


    # print prior and posterior rms
    print('\n===== PRIOR RMS:',file=outfile)
    print('% temp:        {0:.3f}'.format(prior_rms.temp),file=outfile)
    print('% rh:          {0:.3f}'.format(prior_rms.rh),file=outfile)
    print('% water vapor: {0:.3f}'.format(prior_rms.water_vapour),file=outfile)
   
    
    print('\n===== POSTERIOR RMS:',file=outfile)
    print('% temp:        {0:.3f}'.format(post_rms.temp),file=outfile)
    print('% rh:          {0:.3f}'.format(post_rms.rh),file=outfile)
    print('% water vapor: {0:.3f}'.format(post_rms.water_vapour),file=outfile)

    outfile.close()

    return 1

def MIRTO_Validation(retrieval,sonde,station,dists,indices,workdir):

    out_rms = []
    for fov_index in indices:
        
        if not np.ma.is_masked(fov_index):
            post = ValidationDataset(retrieval.profiles[1].temperature[fov_index],
                                     retrieval.profiles[1].rh[fov_index],
                                     retrieval.profiles[1].water_vapour[fov_index],
                                     retrieval.profiles[1].pressure[fov_index])

            prior = ValidationDataset(retrieval.profiles[0].temperature[fov_index],
                                      retrieval.profiles[0].rh[fov_index],
                                      retrieval.profiles[0].water_vapour[fov_index],
                                      retrieval.profiles[0].pressure[fov_index])
            
            # merge pressure maps in a masked array   
            new_pressure_grid = np.ma.masked_invalid([ x if sonde.levels[-1]<x<sonde.levels[0]\
                                               else np.nan for x in post.pressure ])
            
            # interpolated values of t,rh and mr
            int_sonde = ValidationDataset(
                                interpolate_pressure_grid(new_pressure_grid,sonde.levels,sonde.temp),
                                interpolate_pressure_grid(new_pressure_grid,sonde.levels,sonde.relh),
                                interpolate_pressure_grid(new_pressure_grid,sonde.levels,sonde.mixr)
                                        )
            
            prior_rms = ValidationDataset( 
                            np.sqrt(np.mean((int_sonde.temp - prior.temp)**2)),
                            np.sqrt(np.mean((int_sonde.rh   - prior.rh)**2)),
                            np.sqrt(np.mean((int_sonde.water_vapour - prior.water_vapour)**2))
                                   )
        
            post_rms = ValidationDataset( 
                            np.sqrt(np.mean((int_sonde.temp - post.temp)**2)),
                            np.sqrt(np.mean((int_sonde.rh   - post.rh)**2)),
                            np.sqrt(np.mean((int_sonde.water_vapour - post.water_vapour)**2))
                                   )

            prior_rms_profile = ValidationDataset( 
                                   (int_sonde.temp - retrieval.profiles[0].temperature[fov_index,:]),
                                   (int_sonde.rh   - retrieval.profiles[0].rh[fov_index,:]),
                                   (int_sonde.water_vapour - retrieval.profiles[0].water_vapour[fov_index,:])
                                   )
            
            post_rms_profile = ValidationDataset( 
                                   abs(int_sonde.temp         - retrieval.profiles[1].temperature[fov_index,:]),
                                   abs(int_sonde.rh           - retrieval.profiles[1].rh[fov_index,:]),
                                   abs(int_sonde.water_vapour - retrieval.profiles[1].water_vapour[fov_index,:])
                                   )
                    
            out_rms.append( [prior_rms, 
                             post_rms, 
                             prior_rms_profile,  
                             post_rms_profile,
                             sonde.date,
                             fov_index])
            
    #try:
    if  '--' not in station.icao:
        title = 'station_{}_{}_ov_{}'.format(station.icao, sonde.date, retrieval.name.replace('MIRTO_',''))
    else:
        title = 'station_{}_{}_ov_{}'.format(station.synop, sonde.date, retrieval.name.replace('MIRTO_',''))
    
    print(f"String: {title}; Workdir: {workdir}")    
    MIRTO_SaveMultiResults_nc(out_rms,
                                  retrieval,
                                  int_sonde,
                                  title,
                                  dists,
                                  workdir)
    LOG.info('Results saved in {}/MIRTO_validation_{}.nc'.format(workdir,title))
    #except Exception as e:
    #    LOG.info('An error occurred in saving results!\n{}'.format(e))

    return

def MIRTO_MultiStation_Validation(retrievals,stations,dists,sonde_dates,workdir):
    from atmos.atmos_tools import interpolate_pressure_grid
    from validation.download_wyoming_sondes import download_sondes
        
        
    # Save output stations used to create netcdf files
    out_stations      = []
    sonde_indx        = []
    out_sondes        = []
    for index, station in enumerate(stations):
        if station is not None:
            if station.synop not in out_stations and station.icao not in out_stations:
                out_stations.append(station.synop)
                """
                # Try to add the icao number of the station
                if '--' in station.icao:
                    out_stations.append(station.synop)
                else:
                    out_stations.append(station.icao)
                """
                sonde_indx.append(index)
    out_stations = np.array(out_stations)

    # Read rowinsondes for all the stations
    sondes       = [ ]
    for station_index, station in enumerate(stations):
        if station is None:
            sondes.append(None)
        else:
            sondes.append(download_sondes(station,sonde_dates[station_index]))
        if sondes[-1] is None and station is not None:
            stations[station_index] = None

    # Check if there is at least one station for the validation
    n_stations_found = retrievals.latitude.size - stations.count(None)
    if n_stations_found == 0:
        LOG.info("No valid station found for the validation")
        sys.exit()
    else:
        LOG.info("Number of FOVs used for the validation: {} / {}".format(n_stations_found, retrievals.latitude.size))                
       
    # Count how many FOVs share the same station
    nfovs_for_station = np.zeros( out_stations.size)
    for station in stations:
        if station is not None:
            nfovs_for_station[ out_stations == station ] += 1
    out_rms    =  { station : [] for station in out_stations }
    out_sondes =  { station : [] for station in out_stations }

    """
    for fov_index,station in enumerate(stations):
        if station is not None:
            LOG.debug('--> {} {} {} {}'.format(station.synop,
                                           sondes[fov_index].date,
                                           retrievals.fov_original_index(fov_index)[1],
                                           sondes[fov_index].temp[3]))
    """

    # Cycle on stations/FOVs
    for fov_index, station in enumerate(stations):
            if station is not None and retrievals.profiles[1].rh[fov_index,:].max() <= 100 : # Temporary

                prior = retrievals.profiles[0]
                post  = retrievals.profiles[1]
                sonde = sondes[fov_index]
                
                # merge pressure maps in a masked array   
                new_pressure_grid = np.ma.masked_invalid([ x if sonde.levels[-1]<x<sonde.levels[0]\
                                                           else np.nan for x in post.pressure[fov_index,:] ])
                average_indices = np.where(post.pressure[fov_index,:] > 200)[0]
                # interpolated values of t,rh and mr
                int_sonde = ValidationDataset(
                                        interpolate_pressure_grid(new_pressure_grid,sonde.levels,sonde.temp),
                                        interpolate_pressure_grid(new_pressure_grid,sonde.levels,sonde.relh),
                                        interpolate_pressure_grid(new_pressure_grid,sonde.levels,sonde.mixr),
                                        new_pressure_grid,
                                        lat = sonde.lat,
                                        lon = sonde.lon
                                            )
                prior_rms = ValidationDataset( 
                                np.sqrt(np.average((int_sonde.temp[average_indices] - prior.temperature[fov_index,average_indices])**2,weights=prior.pressure[fov_index,average_indices])),
                                np.sqrt(np.average((int_sonde.rh[average_indices]           - prior.rh[fov_index,average_indices])**2,weights=prior.pressure[fov_index,average_indices])),
                                np.sqrt(np.average((int_sonde.water_vapour[average_indices] - prior.water_vapour[fov_index,average_indices])**2,weights=prior.pressure[fov_index,average_indices]))
                                       )
            
                post_rms = ValidationDataset( 
                                np.sqrt(np.average((int_sonde.temp[average_indices]         - post.temperature[fov_index,average_indices])**2,weights=prior.pressure[fov_index,average_indices])),
                                np.sqrt(np.average((int_sonde.rh[average_indices]           - post.rh[fov_index,average_indices])**2,weights=prior.pressure[fov_index,average_indices])),
                                np.sqrt(np.average((int_sonde.water_vapour[average_indices] - post.water_vapour[fov_index,average_indices])**2,weights=prior.pressure[fov_index,average_indices]))
                                )
                 
                prior_rms_profile = ValidationDataset( 
                                       (int_sonde.temp - prior.temperature[fov_index,:]),
                                       (int_sonde.rh   - prior.rh[fov_index,:]),
                                       (int_sonde.water_vapour - prior.water_vapour[fov_index,:])
                                       )
                post_rms_profile = ValidationDataset( 
                                       abs(int_sonde.temp         - post.temperature[fov_index,:]),
                                       abs(int_sonde.rh           - post.rh[fov_index,:]),
                                       abs(int_sonde.water_vapour - post.water_vapour[fov_index,:])
                                        )
                out_station_index = np.where( (out_stations == station.synop) | (out_stations == station.icao) )[0]
                out_sondes[out_stations[out_station_index][0]].append( int_sonde )
                out_rms[out_stations[out_station_index][0]].append(   [ prior_rms, 
                                                                        post_rms, 
                                                                        prior_rms_profile, 
                                                                        post_rms_profile,
                                                                        sonde.date,
                                                                        fov_index
                                                                       ] )
    no_results_flag = True
    print('nstations: ',len(stations))
    for out_station_index,station in enumerate(out_stations):        # Cycle on output station
        # Select rowinsonde dates
        out_dates = []
        for data_structure in out_rms[station]:
            if data_structure[4] not in out_dates:
                out_dates.append(data_structure[4])
        """
        LOG.debug(' ')
        LOG.debug('station: {}'.format(station))
        LOG.debug('overpass:')
        for overpasses in out_rms[station]:
            LOG.debug('           {}'.format(retrievals.fov_original_index(overpasses[-1])[1]))
        LOG.debug('sonde date:')
        for overpasses in out_rms[station]:
            LOG.debug('           {}'.format(overpasses[-2]))
        LOG.debug('sonde temp:')
        for i,overpasses in enumerate(out_rms[station]):
            LOG.debug('           {}'.format( out_sondes[station][i].temp[:5] ))            
        LOG.debug(' ')
        """

        for date in out_dates:   # Cycle on rowinsonde date
            out_indices = [ ]
            for index, data_structure in enumerate(out_rms[station]):
                if data_structure[4] == date:
                    out_indices.append(index)
            filename = 'station_'+ station + '_' + date
            sonde = out_sondes[station][out_indices[0]]
            if len(out_indices)!=0:
                if no_results_flag:
                    no_results_flag = False 
                try:
                    MIRTO_SaveMultiResults_nc(out_rms[station][out_indices[0]:out_indices[-1]+1],
                                                 retrievals,
                                                 sonde,
                                                 filename,
                                                 dists,
                                                 workdir)
                    LOG.info('Results saved in {}/MIRTO_validation_{}.nc'.format(workdir,filename))
                    
                except Exception as error:
                        LOG.fatal('An error occurred in saving results for station {}!'.format(station))
                        LOG.debug('Error {}'.format(error))
      
    # Check if the code found at least one station
    if no_results_flag:
          LOG.info("No stations found for the validation!")

    return

def MIRTO_SaveResults(sonde,prior,post,string,name,dist,workdir):
    
    outfile = open(workdir+'/'+string+'_validation.txt','a')
    print('===================================================',file=outfile)   
    print('      {}     '.format(name),file=outfile)
    print('===================================================',file=outfile)   
    
    print('\n===== RET INFO',file=outfile)
    print('dist: {0:.2f} km'.format(dist),file=outfile)    
    
    print('\n===== TEMPERATURE',file=outfile)    
    print('pressure prior post sonde ',file=outfile)
    for p,t1,t2,s,control in zip(post.pressure,prior.temp,post.temp,sonde.temp,sonde.temp.mask):
        if control:
                print('{:.4} {:.4} {:.4}   --'.format(p,t1,t2),file=outfile)    
        else:
                print('{:.4} {:.4} {:.4}  {:.4}'.format(p,t1,t2,s),file=outfile)    

    print('\n===== REL. HUMIDITY',file=outfile)    
    print('pressure prior post  sonde ',file=outfile)
    for p,rh1,rh2,s,control in zip(post.pressure,prior.rh,post.rh,sonde.rh,sonde.rh.mask):
        if control:
            print('{:.4}  {:.4}  {:.4}  --'.format(p,rh1,rh2),file=outfile)    
        else:
            print('{:.4}  {:.4}  {:.4}  {:.4}'.format(p,rh1,rh2,s),file=outfile)    

    print('\n===== WATER VAPOR',file=outfile)    
    print('pressure prior post sonde',file=outfile)
    for p,q1,q2,s,control in zip(post.pressure,
                                 prior.water_vapour,
                                 post.water_vapour,
                                 sonde.water_vapour,
                                 sonde.water_vapour.mask):
        if control:
            print('{:.4} {:.4} {:.4}  --'.format(p,q1,q2),file=outfile)    
        else:
            print('{:.4} {:.4} {:.4} {:.4}'.format(p,q1,q2,s),file=outfile)    
            
            
    prior_rms = ValidationDataset( 
                    np.sqrt(np.mean((sonde.temp - prior.temp)**2)),
                    np.sqrt(np.mean((sonde.rh   - prior.rh)**2)),
                    np.sqrt(np.mean((sonde.water_vapour - prior.water_vapour)**2))
                           )

    post_rms = ValidationDataset( 
                    np.sqrt(np.mean((sonde.temp - post.temp)**2)),
                    np.sqrt(np.mean((sonde.rh   - post.rh)**2)),
                    np.sqrt(np.mean((sonde.water_vapour - post.water_vapour)**2))
                           )
     
    # print prior and posterior rms
    print('\n===== PRIOR RMS:',file=outfile)
    print('% temp:        {0:.2f}'.format(prior_rms.temp),file=outfile)
    print('% rh:          {0:.2f}'.format(prior_rms.rh),file=outfile)
    print('% water vapor: {0:.2f}'.format(prior_rms.water_vapour),file=outfile)

    print('\n===== POST RMS:',file=outfile)
    print('% temp:        {0:.2f}'.format(post_rms.temp),file=outfile)
    print('% rh:          {0:.2f}'.format(post_rms.rh),file=outfile)
    print('% water vapor: {0:.2f}'.format(post_rms.water_vapour),file=outfile)
        
    # DA CAMBIAREEEE!
    try:
        if  not post_rms.temp.count() or not post_rms.rh.count() or not post_rms.water_vapour.count():
            print("\nNOT GOOD RETRIEVAL!\n",file=outfile)
    except:
        pass
    
    outfile.close()

    return


def MIRTO_SaveResults_nc(sonde,prior,post,prior_rms,post_rms,string,name,dist,workdir):
    from netCDF4 import Dataset
    
    ncfilepath = workdir+'/'+string+'_MIRTO_validation.nc'
    
    ncfile = Dataset(ncfilepath,mode="w",format="NETCDF4")

    nc_lev_dim = ncfile.createDimension('n_lev',post.pressure.size)
    nc_scalar_dim = ncfile.createDimension('scalar',1)
    
    # Write profiles    
    nc_pressure       = ncfile.createVariable('pressure',np.float32,('n_lev',))
    nc_pressure.units = 'hPa'
    nc_pressure[:]    = post.pressure
    
    nc_dist           = ncfile.createVariable('distance',np.float32,('scalar',))
    nc_dist[:]        = dist
    nc_dist.units     = 'km'
    
    nc_prior_temp       = ncfile.createVariable('prior_temp',np.float32,('n_lev',))
    nc_prior_temp[:]    = prior.temp
    nc_prior_temp.units = 'K'

    nc_post_temp        = ncfile.createVariable('post_temp',np.float32,('n_lev',))
    nc_post_temp[:]     = post.temp
    nc_post_temp.units  = 'K'

    nc_sonde_temp       = ncfile.createVariable('sonde_temp',np.float32,('n_lev',))
    nc_sonde_temp[:]    = sonde.temp
    nc_sonde_temp.units = 'K'

    nc_prior_rh         = ncfile.createVariable('prior_rh',np.float32,('n_lev',))
    nc_prior_rh[:]    = prior.rh
    nc_prior_rh.units = '%'

    nc_post_rh       = ncfile.createVariable('post_rh',np.float32,('n_lev',))
    nc_post_rh[:]    = post.rh
    nc_post_rh.units = '%'

    nc_sonde_rh      = ncfile.createVariable('sonde_rh',np.float32,('n_lev',))
    nc_sonde_rh[:]    = sonde.rh
    nc_sonde_rh.units = '%'

    nc_prior_wv      = ncfile.createVariable('prior_water_vapour',np.float32,('n_lev',))
    nc_prior_wv[:]    = prior.water_vapour
    nc_prior_wv.units = 'g/kg'

    nc_post_wv      = ncfile.createVariable('post_water_vapour',np.float32,('n_lev',))
    nc_post_wv[:]    = post.water_vapour
    nc_post_wv.units = 'g/kg'

    nc_sonde_wv      = ncfile.createVariable('sonde_water_vapour',np.float32,('n_lev',))
    nc_sonde_wv[:]    = sonde.water_vapour
    nc_sonde_wv.units = 'g/kg'
                            
    # Write Prior RMS
    nc_prior_rms_temp = ncfile.createVariable('prior_rms_temp',np.float32,('scalar',))
    nc_prior_rms_temp[:] = prior_rms.temp
    nc_prior_rms_rh = ncfile.createVariable('prior_rms_rh',np.float32,('scalar',))
    nc_prior_rms_rh[:] = prior_rms.rh
    nc_prior_rms_wv = ncfile.createVariable('prior_rms_water_vapour',np.float32,('scalar',))
    nc_prior_rms_wv[:] = prior_rms.water_vapour

    # Write Post RMS
    nc_post_rms_temp = ncfile.createVariable('post_rms_temp',np.float32,('scalar',))
    nc_post_rms_temp[:] = post_rms.temp
    nc_post_rms_rh = ncfile.createVariable('post_rms_rh',np.float32,('scalar',))
    nc_post_rms_rh[:] = post_rms.rh
    nc_post_rms_wv = ncfile.createVariable('post_rms_water_vapour',np.float32,('scalar',))
    nc_post_rms_wv[:] = post_rms.water_vapour

    nc_quality_control = ncfile.createVariable('quality_control','i1',('scalar',))
    nc_quality_control[:] = 1    
    try:
        if  not post_rms.temp.count() or not post_rms.rh.count() or not post_rms.water_vapour.count():
            nc_quality_control[:] = 0
    except:
        pass
    
    ncfile.close()

    return 

def MIRTO_SaveMultiResults_nc(out_rms, retrievals, sonde, string ,dists,workdir):
    from netCDF4 import Dataset
    import os
    
    ncfilepath = workdir+'/MIRTO_validation_{}.nc'.format(string)
    if os.path.isfile(ncfilepath):
        os.system('rm {}'.format(ncfilepath))
    
    ncfile = Dataset(ncfilepath,mode="w",format="NETCDF4")

    n_fovs        = ncfile.createDimension('n_fovs',len(out_rms))
    nc_lev_dim    = ncfile.createDimension('n_lev',retrievals.profiles[0].pressure.shape[1])
    nc_scalar_dim = ncfile.createDimension('scalar',1)
    
    # Define vars
    nc_sonde_lat = ncfile.createVariable('sonde_lat',np.float32,('scalar'))
    nc_sonde_lon = ncfile.createVariable('sonde_lon',np.float32,('scalar'))
    
    nc_pressure       = ncfile.createVariable('fovs_levels',np.float32,('n_fovs','n_lev'))
    nc_pressure.units = 'hPa'

    nc_sonde_pressure       = ncfile.createVariable('sonde_levels',np.float32,('n_lev',))
    nc_sonde_pressure.units = 'hPa'

    nc_dist           = ncfile.createVariable('distance',np.float32,('n_fovs',))
    nc_dist.units     = 'km'

    nc_fov_lat           = ncfile.createVariable('fov_lat',np.float32,('n_fovs',))
    nc_fov_lat.units     = 'north-degree'

    nc_fov_lon           = ncfile.createVariable('fov_lon',np.float32,('n_fovs',))
    nc_fov_lon.units     = 'east-degree'
    
    nc_prior_temp       = ncfile.createVariable('prior_temp',np.float32,('n_fovs','n_lev'))
    nc_prior_temp.units = 'K'

    nc_post_temp        = ncfile.createVariable('post_temp',np.float32,('n_fovs','n_lev'))
    nc_post_temp.units  = 'K'

    nc_sonde_temp       = ncfile.createVariable('sonde_temp',np.float32,('n_lev'))
    nc_sonde_temp.units = 'K'

    nc_prior_rh         = ncfile.createVariable('prior_rh',np.float32,('n_fovs','n_lev'))
    nc_prior_rh.units = '%'

    nc_post_rh       = ncfile.createVariable('post_rh',np.float32,('n_fovs','n_lev'))
    nc_post_rh.units = '%'

    nc_sonde_rh      = ncfile.createVariable('sonde_rh',np.float32,('n_lev',))
    nc_sonde_rh.units = '%'

    nc_prior_wv      = ncfile.createVariable('prior_water_vapour',np.float32,('n_fovs','n_lev'))
    nc_prior_wv.units = 'g/kg'

    nc_post_wv      = ncfile.createVariable('post_water_vapour',np.float32,('n_fovs','n_lev'))
    nc_post_wv.units = 'g/kg'

    nc_sonde_wv      = ncfile.createVariable('sonde_water_vapour',np.float32,('n_lev',))
    nc_sonde_wv.units = 'g/kg'

    nc_prior_rms_wv = ncfile.createVariable('prior_rms_water_vapour',np.float32,('n_fovs','scalar'))
    nc_prior_rms_rh = ncfile.createVariable('prior_rms_rh',np.float32,('n_fovs','scalar'))
    nc_prior_rms_temp = ncfile.createVariable('prior_rms_temp',np.float32,('n_fovs','scalar'))

    nc_prior_rms_profile_wv = ncfile.createVariable('prior_rms_profile_water_vapour',np.float32,('n_fovs','n_lev'))
    nc_prior_rms_profile_rh = ncfile.createVariable('prior_rms_profile_rh',np.float32,('n_fovs','n_lev'))
    nc_prior_rms_profile_temp = ncfile.createVariable('prior_rms_profile_temp',np.float32,('n_fovs','n_lev'))

    nc_post_rms_wv = ncfile.createVariable('post_rms_water_vapour',np.float32,('n_fovs','scalar'))
    nc_post_rms_rh = ncfile.createVariable('post_rms_rh',np.float32,('n_fovs','scalar'))
    nc_post_rms_temp = ncfile.createVariable('post_rms_temp',np.float32,('n_fovs','scalar'))

    nc_post_rms_profile_wv   = ncfile.createVariable('post_rms_profile_water_vapour',np.float32,('n_fovs','n_lev'))
    nc_post_rms_profile_rh   = ncfile.createVariable('post_rms_profile_rh',np.float32,('n_fovs','n_lev'))
    nc_post_rms_profile_temp = ncfile.createVariable('post_rms_profile_temp',np.float32,('n_fovs','n_lev'))
    
    nc_quality_control = ncfile.createVariable('quality_control','i1',('n_fovs',))
    nc_retrieval_name  = ncfile.createVariable('retr_name',str,('n_fovs',))
    nc_fov_index       = ncfile.createVariable('fov_index',np.int64,('n_fovs',))


    # Write sonde vars
    try:
        nc_sonde_temp[:]     = sonde.temp
        nc_sonde_rh[:]       = sonde.rh
        nc_sonde_wv[:]       = sonde.water_vapour
        nc_sonde_pressure[:] = sonde.pressure
    except:  # Single Station Version - Improve this part !!
        nc_sonde_temp     = sonde.temp
        nc_sonde_rh       = sonde.rh
        nc_sonde_wv       = sonde.water_vapour
        nc_sonde_pressure = sonde.pressure
    nc_sonde_lat[:]        = sonde.lat
    nc_sonde_lon[:]        = sonde.lon

    # Write retrievals vars
    post  = retrievals.profiles[1]
    prior = retrievals.profiles[0]
    
    for index, rms in enumerate(out_rms):

        fov_index = rms[-1]
        
        # FOVs vars
        nc_pressure[index,:]      = post.pressure[fov_index,:]
        nc_dist[index]            = dists[fov_index]
        nc_fov_lat[index]          = retrievals.latitude[fov_index]
        nc_fov_lon[index]          = retrievals.longitude[fov_index]

        nc_prior_temp[index,:]    = prior.temperature[fov_index,:]
        nc_post_temp[index,:]     = post.temperature[fov_index,:]
        nc_prior_rh[index,:]      = prior.rh[fov_index,:]    
        nc_post_rh[index,:]       = post.rh[fov_index,:]     
        nc_prior_wv[index,:]      = prior.water_vapour[fov_index,:]     
        nc_post_wv[index,:]       = post.water_vapour[fov_index,:]
                  
        # Write Prior RMS
        nc_prior_rms_temp[index,:] = rms[0].temp
        nc_prior_rms_rh[index,:]   = rms[0].rh
        nc_prior_rms_wv[index,:]   = rms[0].water_vapour
    
        # Write Post RMS
        nc_post_rms_temp[index,:] = rms[1].temp
        nc_post_rms_rh[index,:]   = rms[1].rh
        nc_post_rms_wv[index,:]   = rms[1].water_vapour

        # Write Prior RMS Profile
        nc_prior_rms_profile_temp[index,:] = rms[2].temp
        nc_prior_rms_profile_rh[index,:]   = rms[2].rh
        nc_prior_rms_profile_wv[index,:]   = rms[2].water_vapour

        # Write Post RMS Profile
        nc_post_rms_profile_temp[index,:] = rms[3].temp
        nc_post_rms_profile_rh[index,:]   = rms[3].rh
        nc_post_rms_profile_wv[index,:]   = rms[3].water_vapour

        # Save fov index and original file name
        try:
            tmp = retrievals.fov_original_index(fov_index)
        except:  # Single station Version - Improve this part
            tmp = [ fov_index, retrievals.name ]
        nc_fov_index[index]       = tmp[0]
        nc_retrieval_name[index]  = tmp[1]
        nc_quality_control[index] = 1
        try:
            if np.isnan(rms[1].temp) or np.isinf(rms[1].rh) or np.isnan(rms[1].water_vapour) or np.isinf(rms[1].water_vapour):
                nc_quality_control[index] = 0
        except:
            pass
    
    ncfile.close()

    return 
    
def WRF_Validation(retrieval,sonde,station,dists,min_index,workdir):
    """
    To be added...
    """
    return


def check_date(name,sonde_date,exitmode=True):
    from os.path import basename
    
    # Max temporal difference (hours)
    MAX_TIME = 1
        
    # extract date from the name
    if 'L2VDP' in name:
        # SPS
        retr_date = basename(name).split('-')[2]
        retr_date = datetime.strptime(retr_date,"%Y%m%d%H%M%SZ")
    elif 'wrf' in name:
        """
        To be added...
        """
        raise ValueError('WRF control To be added...')
    else:
        """
        MIRTO
        """
        name      = name.split('/')[-1]
        retr_date = datetime.strptime(name,"%Y%m%d_%H%M%S")

    time_diff = int(np.abs((retr_date - sonde_date).total_seconds())/3600.)
    if exitmode:
        if time_diff >= MAX_TIME:
            raise DomainError('Time difference between retrieval and rowinsonde'
                             'is bigger than {} hours:\n'
                             'time diff: {} h\n'
                             'ret date: {}\n'
                             'sonde date: {}\n'.format(
                                     MAX_TIME,
                                     time_diff,
                                     retr_date,
                                     sonde_date))
    else:
        return False if time_diff >= MAX_TIME else True
    
    
def main(argv):
    import getopt,os
    from validation.read_stations_list import stations_list
    from data_reader.input_reader import InputData

    retrieval  = ''
    station    = ''
    sonde_date = ''
    singlemode = True
    
    # get options
    opts ,args = getopt.getopt(argv,"r:s:d:w:v:",["retrieval=","station=","date=","workdir=","verbose="])
    for opt,arg in opts:
        if opt in ('-r','--retrieval'):
            retr_name = arg
        elif opt in ('-s','--station'):
            station = arg
        elif opt in ('-d','--date'):
            sonde_date = datetime.strptime(arg,'%Y-%m-%dT%H:%M')
        elif opt in ('-w','--workdir'):
            workdir = arg
        elif opt in ('-v','--verbose'):
            verbose = arg.upper()
            
    try:
        LOG.setLevel(verbose)
        streamhandler.setLevel(verbose)
    except:
        LOG.setLevel('INFO')
        streamhandler.setLevel('INFO')

    # LOG.addHandler(streamhandler)
    date = datetime.now().isoformat()[:-10]
    hdlr = logging.FileHandler('/var/tmp/validate_retrievals_{}.log'.format(date))
    LOG.addHandler(hdlr)
    LOG.addHandler(streamhandler)
    
    # Check the input received
    if retr_name == '':
        sys.exit("Insert retrieval's path (with --retrieval= or -r)")

    # check workdir string
    if workdir == '':
        workdir = os.getcwd()
    elif workdir[-1] == '/':
        workdir = workdir[:-1]

    # Import station list
    stations = stations_list()
    
    if station == '':
        # Multi-station mode:
        # Version without a particular station.
        # The script search for the nearest stations to each FOV
        singlemode = False
    else:
        # station is a station_class    
        # class defined in read_stations_list.py
        station = stations_list().get(station)
        
    if singlemode:    
        # Single-station mode
        from validation.download_wyoming_sondes import download_sondes

        # Build retr_list 
        if os.path.isdir(retr_name):
            
            file_list = os.listdir(retr_name)
            file_list.sort()
            if retr_name[-1] != '/':
                retr_name = retr_name+'/'
            
            
            try:
                AMETHYST_file_list = [ retr_name+x for x in file_list if 'amethyst' in os.listdir(retr_name+x)]
            except:
                AMETHYST_file_list = []

            AMETHYST_file_list.sort()
                
            retr_list = AMETHYST_file_list
            
        else:
            retr_list = [retr_name]
            
        for retr_name in retr_list :
     
            LOG.info('Validating '+os.path.basename(retr_name)+'...')
    
            # retrieval is an InputData class 
            # defined in tools/input_reader.py
            try:
                # check date
                check_date(retr_name,sonde_date)
                                
                # import retrieval
                retrieval = InputData(retr_name) if 'wrf' in retr_name or       \
                                                    'L2VDP' in retr_name else   \
                            InputData(retr_name+'/amethyst/results.nc')
                            
                # download rawinsonde from the given station
                sonde = download_sondes(station,sonde_date)
                
                if sonde != None:
                    # validate retrieval with the given station
                    validate_retrieval_single_station(sonde,retrieval,station,workdir)
    
            # Save errors in the logfile
            except OSError as os_err:
                LOG.debug(os_err)
            except DomainError as v_err:
                LOG.debug(v_err)

    else:
        # Multi-station mode
        # Note (P.Scaccia): Some redundancies are still present in the code...to be improved
        from data_reader.mirto_wrapper import MIRTO_multi_overpass_Wrapper

        # Build retr_list 
        if os.path.isdir(retr_name):
            file_list = os.listdir(retr_name)
            file_list.sort()
            if retr_name[-1] != '/':
                retr_name = retr_name+'/'
            AMETHYST_file_list = [ retr_name+x+'/amethyst/results.nc' for x in file_list if os.path.isdir(retr_name+x+'/amethyst')]
            AMETHYST_file_list_2 = [ retr_name+x+'/results.nc' for x in file_list if os.path.isfile(retr_name+x+'/results.nc')]
            AMETHYST_file_list.sort()
            AMETHYST_file_list_2.sort()

            retr_list = AMETHYST_file_list + AMETHYST_file_list_2
        else:
            retr_list = [retr_name]
        if retr_list == []:
            sys.exit('No overpass found')

        try:
            # Store all overpasses in a single dataset 
            retrievals = MIRTO_multi_overpass_Wrapper(retr_list)
            dists    = np.zeros( retrievals.latitude.size )
            validating_stations = []    

            # Read stations nearest to each FOV
            for fov_index,coord in enumerate(zip(retrievals.longitude,retrievals.latitude)):
                lon = coord[0]
                lat = coord[1]

                station, dist =  stations.get_nearest_station(lat,lon , R = 75) 
                if station is None:
                    LOG.debug('No station found for FOV n. {} with coordinates: ({:.3f},{:.3f})'.format(fov_index,lat,lon))
                validating_stations.append(station)
                dists[fov_index] = dist
            dists = np.ma.masked_invalid(dists)

            l = list(dict.fromkeys(validating_stations))
            LOG.info("Found {} stations after the spatial and temporal filter".format( len(l) - l.count(None)))
    

            # Time filter
            retr_dates = [ retrievals.fov_original_index(fov_index)[1].replace('MIRTO_','') for fov_index in range(retrievals.latitude.size) ]            
            if sonde_date != '':  
                sonde_dates = [ sonde_date for x in range(retrievals.longitude.size)]
            else:   
                sonde_dates = []
                for retr_date in retr_dates:        
                        date = datetime.strptime(retr_date.replace('MIRTO_',''),'%Y%m%d_%H%M%S')
                        if (date.hour - 12) < -6:
                            sonde_dates.append(datetime(date.year, date.month, date.day,0,0,0 ))
                        elif abs(date.hour - 12) <= 6:
                            sonde_dates.append(datetime(date.year, date.month, date.day,12,0,0 ))
                        elif (date.hour - 12) > 6:
                            sonde_dates.append(datetime(date.year, date.month, date.day ,0,0,0 ) + timedelta(hours=24) )
                        else:
                            sys.exit("Error in date format!")
            for fov_index, sonde_date in enumerate(sonde_dates):
                    try:
                        if not check_date(retr_dates[fov_index],sonde_date, exitmode = False):
                            LOG.debug("FOV n.{} removed from the validation.".format(fov_index))
                            LOG.debug("Overpass date {} too distant from rowinsonde time".format(retrievals.fov_original_index(fov_index)[1]))
                            validating_stations[fov_index] = None   # Remove station for fov outside the temporal window
                    except Exception as error:
                        LOG.info(error)
                        validating_stations[fov_index] = None
                
   
            # Validate retrievals in the aggregated dataset with all the stations
            validate_retrieval_multi_station(retrievals, validating_stations, dists, sonde_dates, workdir)

        # Save errors in the logfile
        except OSError as os_err:
                sys.exit(os_err)
        except DomainError as v_err:
                sys.exit(v_err)


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

