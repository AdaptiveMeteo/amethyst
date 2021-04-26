#!/usr/bin/env python

# This file is part of Cris2observations.
#
# Cris2observations is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Cris2observations is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Cris2observation. If not, see <http://www.gnu.org/licenses/>.

"""
VIIRS/CRIS 2 Amethyst Cloud Mask

This script is an open source software written to convert

:copyright: 2021 by AdaptieMeteo S.r.l.
"""
import logging
from argparse import ArgumentParser
from sys import exit as sysexit
from traceback import format_exc
import numpy as np
from preprocessor.fov_generator.cris.fov_generator_cris_utilities.gclos_wrapper import GCLOS_Wrapper
from preprocessor.fov_generator.cris.fov_generator_cris_utilities.scris_wrapper import SCrISWrapper
from preprocessor.fov_generator.cris.fov_generator_cris_utilities.viirs_wrapper import VIIRS_Wrapper
from preprocessor.fov_generator.cris.fov_generator_cris_utilities.atmos_tools   import cris_viirs_cloudmask
from preprocessor.fov_generator.fov_generator_utilities.boxes                   import Rectangle
from datetime import datetime

__author__ = 'Paolo Antonelli<paolo.antonelli@adaptivementeo.com>'
__copyright__ = "Copyright 2021, AdaptiveMeteo S.r.l."
__credits__ = ["Paolo Antonelli", "Paolo Scaccia"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = [ "Paolo Antonelli", "Paolo Scaccia"]
__email__ = "paolo.antonelli@adaptivemeteo.com"


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

    #filetype can be: GMODO, GITCO, GMCTO, IICMO, JRR-CloudMask
    
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
    try:
        filelist = [ name for name in os.listdir(gcrso_path)
                                 if filetype in name
                                 and ( abs(start - read_epoch( name[11:19],name[21:27])) < 15*60)
                                ]
    except:
        filelist = [ name for name in os.listdir(gcrso_path)
                                 if filetype in name
                                 and ( abs(start - read_epoch(name.split('_')[3][1:9],
                                                              name.split('_')[3][9:-1])) < 15*60)
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
    if type(viirs_dataset.cloudmask.mask) == np.bool_:
                  viirs_dataset.cloudmask.mask = np.zeros_like(viirs_dataset.cloudmask,dtype=np.bool)
    viirs_dataset.scan_position = np.concatenate(  (viirs_dataset.scan_position,
                                                    new_dataset.scan_position ),
                                                    axis = 0  )
    return

def write_output(cris_dataset, viirs_dataset, indices, cloud_stats, outfile):
    from netCDF4 import Dataset
    from preprocessor.fov_generator.cris.fov_generator_cris_utilities.array_reshapers import array_1d, transform_index

    #outfile = 'cloudmask.nc'

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
    #PaoloA#viirs_lats = results.createVariable('viirs_lats',np.float64,('CrIS_Size','viirs_index'))
    #PaoloA#viirs_lons = results.createVariable('viirs_lons',np.float64,('CrIS_Size','viirs_index'))
    #PaoloA#viirs_cloudmask = results.createVariable('viirs_cloudmask',np.int,('CrIS_Size','viirs_index'))

    results.description = 'Cloudmask statistics for CrIs file: {}'.format(cris_dataset.path.split('/')[-1])
    results.hystory = "Created {}".format(
                                          datetime.now().isoformat().replace('T',' ')[:-10]
                                          )
    cris_lats.units = 'degree_north'
    cris_lons.units = 'degree_east'

    #PaoloA#viirs_lats.units = 'degree_north'
    #PaoloA#viirs_lons.units = 'degree_east'
    stats.units = 'percent'

    cris_lats[:] = array_1d(cris_dataset.lats)[:]
    cris_lons[:] = array_1d(cris_dataset.longs)[:]
    stats[:] = cloud_stats.reshape( (cris_dataset.longs.size,4) )[:]

    results.close()
    return

def write_ASCII(cris_dataset, viirs_dataset, indices, cloud_stats):
    from preprocessor.fov_generator.cris.fov_generator_cris_utilities.array_reshapers import array_1d, transform_index

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


if __name__ == '__main__':
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


def main():
    v_levels = ['debug', 'info', 'warning']

    parser = ArgumentParser()
    parser.add_argument('GCRSO', type=str,
                        help='A GCRSO file to be used to geolocalize the data')
    parser.add_argument('SCRIS', type=str,
                        help='A SCRIS file with the CrIS data')
    parser.add_argument('--outfile', '-o', type=str,
                        help='The netcdf file that must be created as output')
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='The level of verbosity of the software')
    parser.add_argument('--index-cm-path', '-m', type=str, default=None,
                        help='The path of the file where the index map will '
                             'be saved.')
    parser.add_argument('--observation-number', '-n', type=int, default=None,
                        help='The maximum number of observations that will be '
                             'saved. This is useful to test the software using '
                             'less observations')
    parser.add_argument('--latmin', type=float, default=None,
                        help='All the FOV with latitude lesser than this value '
                             'will be discarded')
    parser.add_argument('--latmax', type=float, default=None,
                        help='All the FOV with latitude greater than this '
                             'value will be discarded')
    parser.add_argument('--lonmin', type=float, default=None,
                        help='All the FOV with longitude lesser than this '
                             'value will be discarded')
    parser.add_argument('--lonmax', type=float, default=None,
                        help='All the FOV with longitude greater than this '
                             'value will be discarded')

    argv = parser.parse_args()

    # Configure the LOGGER object
    verbosity = getattr(logging, argv.verbose.upper())
    LOGGER.setLevel(verbosity)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - '
                                  '%(filename)s: %(message)s',
                                  datefmt='%m/%d/%Y %H:%M:%S')

    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(verbosity)
    streamhandler.setFormatter(formatter)
    LOGGER.addHandler(streamhandler)

    LOGGER.info('Reading the GCRSO file')
    try:
        gcrso = GCLOS_Wrapper(argv.GCRSO)
    except:
        LOGGER.error('Read of the GCRSO file failed!')
        LOGGER.debug(format_exc())
        return 1

    LOGGER.info('Reading the SCRIS file')
    try:
        scris = SCrISWrapper(argv.SCRIS)
    except:
        LOGGER.error('Read of the SCrIS file failed!')
        LOGGER.debug(format_exc())
        return 2

    # Create the boxes for the observations
    lats = gcrso.lats.flatten()
    lons = gcrso.longs.flatten()

    if argv.latmin is None:
        latmin = np.min(lats)
    else:
        latmin = argv.latmin
    if argv.latmax is None:
        latmax = np.max(lats)
    else:
        latmax = argv.latmax
    if argv.lonmin is None:
        lonmin = np.min(lons)
    else:
        lonmin = argv.lonmin
    if argv.lonmax is None:
        lonmax = np.max(lons)
    else:
        lonmax = argv.lonmax
    LOGGER.debug('Creating a box of longitude from {} to {} and latitude from '
                 '{} to {}'.format(lonmin, lonmax, latmin, latmax))
    box = Rectangle(lonmin, lonmax, latmin, latmax)

    LOGGER.debug('Checking which observations are inside the box')
    inside_box = box.contains(lons, lats)

    LOGGER.debug('Checking which observations are related to a turned off '
                 'decoder')
    turned_off_decoder = scris.valid_radiances().flatten()

    obs_filter = np.bool_(np.logical_and(inside_box, turned_off_decoder))
    LOGGER.debug('Keeping {} observations'.format(obs_filter))

    if argv.observation_number is not None:
        obs_max = argv.observation_number
        LOGGER.debug('Decreasing the number of the observations to {}'
                     ''.format(obs_max))
        true_entries = np.where(obs_filter)[0]
        if len(true_entries) > obs_max:
            clear_after = true_entries[obs_max]
            obs_filter[clear_after:] = False

    # We count the number of FOVs inside the box to get the number of
    # the FOVs we will save on the output file
    numobs = np.count_nonzero(obs_filter)
    LOGGER.debug('Keeping {} observations'.format(numobs))


    # Read gcrso name and search for gmodo or gmcto, or gitco files
    gcrso_name , viirs_geo_list = search_files(argv.GCRSO,'GMODO')

    overpass_path = argv.GCRSO.replace(gcrso_name,'')

    viirs_cm_list = []
    if viirs_geo_list == []:
        gcrso_name , viirs_geo_list = search_files(argv.GCRSO,'GMTCO')
        gcrso_name , viirs_cm_list = search_files(argv.GCRSO,'IICMO')
    if viirs_geo_list == []:
        gcrso_name , viirs_geo_list = search_files(argv.GCRSO,'GITCO')
        gcrso_name , viirs_cm_list = search_files(argv.GCRSO,'IICMO')
    if viirs_cm_list == []:
    	gcrso_name , viirs_cm_list = search_files(argv.GCRSO,'JRR-CloudMask')
    
    LOGGER.debug("File cris: {}".format(gcrso_name))

    # Exit if empty list is returned
    if viirs_geo_list == []:
        sysexit("No VIIRS file found for the given CrIS file")

    LOGGER.debug("Reading and merging VIIRS Datasets")
    # Read and merge togheter all viirs dataset

    #Add overpass path to list of viirs geo files
    viirs_geo_list = [overpass_path +'/' + string for string in viirs_geo_list]

    #Paolo added 24/03/2021
    viirs_geo_list.sort()
    viirs_cm_list.sort()

    LOGGER.debug("Reading {}".format(viirs_geo_list[0]))
    if viirs_cm_list == []:
        viirs_dataset = VIIRS_Wrapper(viirs_geo_list[0])
        if len(viirs_geo_list) > 1:
            for geo_filename in viirs_geo_list[1:]:
                LOGGER.debug("Reading {}".format(geo_filename))
                new_dataset = VIIRS_Wrapper(geo_filename)
                merge_viirs_dataset(viirs_dataset,new_dataset)
    else:
        viirs_cm_list = [overpass_path +'/' + string for string in viirs_cm_list]
        viirs_dataset = VIIRS_Wrapper(viirs_geo_list[0],viirs_cm_list[0])
        if len(viirs_geo_list) > 1:
            for geo_filename, cm_filename in zip(viirs_geo_list[1:],viirs_cm_list[1:]):
                LOGGER.debug("Reading {}".format(geo_filename))
                new_dataset = VIIRS_Wrapper(geo_filename,cm_filename)
                merge_viirs_dataset(viirs_dataset,new_dataset)

    # Release last dataset read to save memory
    new_dataset = None

    #sysexit("Temporary Stop")

    LOGGER.debug("Reading CrIS Datasets...")
    # Reda cris dataset
    # cris_dataset = CRIS_Wrapper(gcrso_name)

    cris_dataset = gcrso

    LOGGER.debug("Computing cloudmask statistics...")
    # Get collocation and compute cloudmask statistics
    cloud_stats , indices = cris_viirs_cloudmask(cris_dataset, viirs_dataset, obs_filter)

    LOGGER.debug("Writing results")
    # Write results in a NETCDF4 
    write_output(cris_dataset, viirs_dataset, indices, cloud_stats, argv.outfile)
    # write_ASCII(cris_dataset,viirs_dataset, indices, cloud_stats)
    

    LOGGER.info('Execution complete')
    return 0


if __name__ == '__main__':
    sysexit(main())
