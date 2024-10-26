import logging
import numpy as np

__author__    = "Paolo Scaccia and Paolo Antonelli"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__   = ["Paolo Scaccia", "Paolo Antonelli"]
__license__   = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

LOGGER = logging.getLogger(__name__)

import logging
#import datetime
from argparse import ArgumentParser
from sys import exit as sysexit
from data_reader.cris_wrapper import CRIS_Wrapper
from data_reader.viirs_wrapper import VIIRS_Wrapper
from datetime import datetime, timedelta
    
if __name__ == '__main__':
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)

def read_epoch(date,time):
    return datetime(        int(date[:4]),
                            int(date[4:6]),
                            int(date[6:8]),
                            int(time[:2]),
                            int(time[2:4]),
                            int(time[4:6])
                                             ).timestamp()

def search_files(gcrso_path,filetype):
    import os
   
    #filetype can be: GMODO, GITCO, GMCTO, IICMO
 
    # Check input
    if os.path.isfile(gcrso_path) is False:
        sysexit("File GCRSO not found!")
    elif 'GCRSO' not in gcrso_path:
        sysexit("Not a GCRSO file!")
        
    if '/' in gcrso_path:
        gcrso_name = gcrso_path.split('/')[-1]    
        gcrso_path = gcrso_path.strip(gcrso_name)
    else:
        gcrso_name = gcrso_path
        gcrso_path = './'
    
    
    start = read_epoch( gcrso_name[11:19],gcrso_name[21:27])
    end   = read_epoch( gcrso_name[11:19],gcrso_name[30:37])
    
    #filelist = [ name for name in os.listdir(gcrso_path) 
    #                         if filetype in name 
    #                         and ( start < read_epoch( name[11:19],name[21:27]) < end or
    #                               start < read_epoch( name[11:19],name[30:37]) < end )  
    #                         ]   
    
    filelist = [ name for name in os.listdir(gcrso_path) 
                             if filetype in name 
                             and ( abs(start - read_epoch( name[11:19],name[21:27])) < 15*60)  
                            ]   
    
    return gcrso_name , filelist
    

def merge_viirs_dataset(viirs_dataset,new_dataset ):

    viirs_dataset.longs = np.concatenate( (viirs_dataset.longs, 
                                           new_dataset.longs ),
                                          axis = 0)
    
    viirs_dataset.lats  = np.concatenate( (viirs_dataset.lats, 
                                           new_dataset.lats ),
                                          axis = 0)
    
    viirs_dataset.height = np.concatenate( (viirs_dataset.height, 
                                            new_dataset.height ),
                                            axis = 0)

    viirs_dataset.cloudmask = np.ma.concatenate( (viirs_dataset.cloudmask, 
                                                  new_dataset.cloudmask ),
                                                  axis = 0 )
    
    viirs_dataset.scan_position = np.concatenate(  (viirs_dataset.scan_position,
                                                    new_dataset.scan_position ),
                                                    axis = 0  )

def write_output(cris_dataset, viirs_dataset, indices, cloud_stats):    
    from netCDF4 import Dataset
    from geometry.utilities.array_reshapers import array_1d, transform_index
    
    outfile = 'cloudmask.nc'
    
    results = Dataset(outfile, 'w' , format = 'NETCDF4')

    # Find max VIIRS FOVs number inside a signle CrIS FOV (necessary to define viirs lat/lon field)
    max_viirs_ind = 0
    for x in indices:
        if len(x) > max_viirs_ind:
            max_viirs_ind = len(x)

    coor_dim      = results.createDimension('CrIS_Size', cris_dataset.longs.size )
    stat_dim      = results.createDimension('n_prob', 4)
    viirs_ind_dim = results.createDimension('viirs_index', max_viirs_ind)
    
    cris_lats = results.createVariable('cris_lats',np.float64 ,('CrIS_Size',) )
    cris_lons = results.createVariable('cris_lons',np.float64 ,('CrIS_Size',) )
    stats = results.createVariable('stats',float,('CrIS_Size','n_prob') )
    viirs_lats = results.createVariable('viirs_lats',np.float64,('CrIS_Size','viirs_index'))
    viirs_lons = results.createVariable('viirs_lons',np.float64,('CrIS_Size','viirs_index'))
    viirs_cloudmask = results.createVariable('viirs_cloudmask',np.int,('CrIS_Size','viirs_index'))
    
    results.description = 'Cloudmask statistics for CrIs file: {}'.format(cris_dataset.path.split('/')[-1])
    results.hystory = "Created {}".format(
                                          datetime.now().isoformat().replace('T',' ')[:-10]
                                          )
    cris_lats.units = 'degree_north'
    cris_lons.units = 'degree_east'
    
    viirs_lats.units = 'degree_north'
    viirs_lons.units = 'degree_east'
    stats.units = 'percent'    
    
    cris_lats[:] = array_1d(cris_dataset.lats)[:]
    cris_lons[:] = array_1d(cris_dataset.longs)[:]
    
    lat = array_1d(viirs_dataset.lats)
    lon = array_1d(viirs_dataset.longs)
    cmask = array_1d(viirs_dataset.cloudmask)
    
    # define index convertion function
    f = lambda x: transform_index(x, viirs_dataset.lats.shape)

    # for each CrIS FOV read currispondent VIIRS fields    
    for i_row,viirs_indices in enumerate(indices):
        viirs_indices = np.array([ f(ind) for ind in viirs_indices ])
        n = len(viirs_indices)
        #Paolo Antonelli  20.Nov.2019
        #print('i_row= %s viirs_indices= %s' % (i_row,viirs_indices))
        #viirs_lats[i_row,:n] = lat[viirs_indices]
        #viirs_lons[i_row,:n] = lon[viirs_indices]
        #viirs_cloudmask[i_row,:n] = cmask[viirs_indices]
        if n > 0:
           viirs_lats[i_row,:n] = lat[viirs_indices]
           viirs_lons[i_row,:n] = lon[viirs_indices]
           viirs_cloudmask[i_row,:n] = cmask[viirs_indices]
        else:
           n = 100
           viirs_lats[i_row,:n] = -999 
           viirs_lons[i_row,:n] = -999 
           viirs_cloudmask[i_row,:n] = 999
          

    stats[:] = cloud_stats.reshape( (cris_dataset.longs.size,4) )[:]
    
    results.close()
    return

def write_ASCII(cris_dataset, viirs_dataset, indices, cloud_stats):
    from geometry.utilities.array_reshapers import array_1d, transform_index
        
        # Output files
    cris_file  = open('cris_table.dat','w')
    stat_file  = open('stats_table.dat','w')
    
    # Comment first lines with file info
    cris_file.write("% file: {}\n".format(cris_dataset.path.split('/')[-1]))
    cris_file.write("% lat lon\n")
    
    cris_shape = np.shape(cris_dataset.longs)

    lons = array_1d(cris_dataset.longs)
    
    for i, lon in enumerate(lons):
        cris_indices = transform_index(i,cris_shape)
        cris_file.write(" {} {} \n".format(lon,cris_dataset.lats[cris_indices]) )
        stat_file.write(" {} {} {} {} \n".format(*cloud_stats[cris_indices]))


    cris_file.close()
    stat_file.close()
    
    return


def main():
    from atmos.atmos_tools import cris_viirs_cloudmask
    
    v_levels = ['debug', 'info', 'warning']

    # Read external options
    parser = ArgumentParser()
    parser.add_argument('GCRSO', type=str,
                        help='A path to GCRSO file to be used to geolocalize the data')
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='The level of verbosity of the software')
    argv = parser.parse_args()

    # Configure the LOGGER object
    verbosity = getattr(logging, argv.verbose.upper())
    LOGGER.setLevel(verbosity)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - '
                                  '%(funcName)s: %(message)s',
                                  datefmt='%m/%d/%Y %H:%M:%S')
    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(verbosity)
    streamhandler.setFormatter(formatter)
    LOGGER.addHandler(streamhandler)

    # Read gcrso name and search for gmodo or gmcto, or gitco files
    gcrso_name , viirs_geo_list = search_files(argv.GCRSO,'GMODO')
    viirs_iicmo_list = []
    if viirs_geo_list == []:
        gcrso_name , viirs_geo_list = search_files(argv.GCRSO,'GMTCO')
        gcrso_name , viirs_iicmo_list = search_files(argv.GCRSO,'IICMO')
    if viirs_geo_list == []:
        gcrso_name , viirs_geo_list = search_files(argv.GCRSO,'GITCO')
        gcrso_name , viirs_iicmo_list = search_files(argv.GCRSO,'IICMO')

    LOGGER.debug("File cris: {}".format(gcrso_name))
 
    # Exit if empty list is returned
    if viirs_geo_list == []:
        sysexit("No VIIRS file found for the given CrIS file")

    LOGGER.debug("Reading and merging VIIRS Datasets")        
    # Read and merge togheter all viirs dataset
    LOGGER.debug("Reading {}".format(viirs_geo_list[0]))   
    if viirs_iicmo_list == []:
        viirs_dataset = VIIRS_Wrapper(viirs_geo_list[0])
        if len(viirs_geo_list) > 1:
            for geo_filename in viirs_geo_list[1:]:
                LOGGER.debug("Reading {}".format(geo_filename))
                new_dataset = VIIRS_Wrapper(geo_filename)
                merge_viirs_dataset(viirs_dataset,new_dataset)
    else:
        viirs_dataset = VIIRS_Wrapper(viirs_geo_list[0],viirs_iicmo_list[0])
        if len(viirs_geo_list) > 1:
            for geo_filename, iic_filename in zip(viirs_geo_list[1:],viirs_iicmo_list[1,:]):
                LOGGER.debug("Reading {}".format(geo_filename))
                new_dataset = VIIRS_Wrapper(geo_filename,iic_filename)
                merge_viirs_dataset(viirs_dataset,new_dataset)

    # Release last dataset read to save memory 
    new_dataset = None

    LOGGER.debug("Reading CrIS Datasets...")            
    # Reda cris dataset
    cris_dataset = CRIS_Wrapper(gcrso_name)
    
    LOGGER.debug("Computing cloudmask statistics...")        
    # Get collocation and compute cloudmask statistics
    cloud_stats , indices = cris_viirs_cloudmask(cris_dataset, viirs_dataset)

    LOGGER.debug("Writing results")        
    # Write results in a NETCDF4 
    write_output(cris_dataset, viirs_dataset, indices, cloud_stats )   
    # write_ASCII(cris_dataset,viirs_dataset, indices, cloud_stats)
    
    return

    
    
if __name__ == '__main__':
    sysexit(main())
