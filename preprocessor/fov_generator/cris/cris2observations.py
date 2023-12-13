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
CRIS 2 Mirto Observation file converter

This script is an open source software written to convert

:copyright: 2016 by eXact-lab
"""

import logging
from argparse import ArgumentParser
from sys import exit as sysexit
from traceback import format_exc

import numpy as np
import netCDF4 as nc

from preprocessor.fov_generator.cris.fov_generator_cris_utilities.gc_wrapper    import GCWrapper
from preprocessor.fov_generator.cris.fov_generator_cris_utilities.scris_wrapper import SCrISWrapper
from preprocessor.fov_generator.cris.fov_generator_cris_utilities.obs_indices   import ObsMap
from preprocessor.fov_generator.fov_generator_utilities.boxes                   import Rectangle
from preprocessor.fov_generator.fov_generator_utilities.fov_file                import FovFile

__author__ = 'Stefano Piani'
__copyright__ = "Copyright 2021, Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Stefano Piani"
__email__      = "paolo.antonelli@adaptivemeteo.com"


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
    parser.add_argument('output', type=str,
                        help='The netcdf file that must be created as output')
    parser.add_argument('--CLOUDMASK','-cmf', type=str,
                        help='A file with the VIIRS collocated CM data')
    parser.add_argument('--cm_treshold','-cmt', type=float, default=None,
                        help='Clear sky fraction needed to process FOV')
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='The level of verbosity of the software')
    parser.add_argument('--index-map-file', '-m', type=str, default=None,
                        help='The path of the file where the index map will '
                             'be saved. If not specified, the file will not '
                             'be created')
    parser.add_argument('--observation_number', '-n', type=int, default=None,
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
        gcrso = GCWrapper(argv.GCRSO)
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

    if argv.cm_treshold is None:
        cm_treshold=.95
    else:
        cm_treshold = argv.cm_treshold

    if argv.CLOUDMASK is None:
        LOGGER.info('NO CLOUDMASK file available')
    else:
        LOGGER.info('Reading the CLOUDMASK file')
        src = nc.Dataset(argv.CLOUDMASK)
        cm = src.variables['stats'][:]
        cloud_free_index = cm[:,0] + cm[:,1] >= cm_treshold
        obs_filter = np.bool_(np.logical_and(obs_filter, cloud_free_index))
        LOGGER.info('Decreasing the number of the observations to {}'
                     ''.format(np.sum(obs_filter)))
       

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
    LOGGER.debug('Saving {} observations'.format(numobs))
    num_channels = scris.radiances.shape[2]
    LOGGER.info('Writing the output file')
    try:
        fov_file = FovFile(argv.output, numobs, num_channels,'cris')
    except:
        LOGGER.error('Writing of the output file failed!')
        LOGGER.debug(format_exc())
        return 3

    sol_azimuth_angle = gcrso.solar_azimuth_angle
    sol_zenith_angle = gcrso.solar_zenith_angle

    with fov_file:
        fov_file.save_latitude(lats, obs_filter)
        fov_file.save_longitude(lons, obs_filter)
        fov_file.save_time(gcrso.times, obs_filter)
        fov_file.save_radiance(scris.radiances, obs_filter)
        fov_file.save_wavenumbers(scris.wavenumbers)
        fov_file.save_fov_angle(gcrso.sat_zenith_angle, obs_filter)
        fov_file.save_sat_azimuth_angle(gcrso.sat_azimuth_angle, obs_filter)
        fov_file.save_sat_zenith_angle(gcrso.sat_zenith_angle, obs_filter)
        fov_file.save_solar_azimuth_angle(sol_azimuth_angle, obs_filter)
        fov_file.save_solar_zenith_angle(sol_zenith_angle, obs_filter)

    if argv.index_map_file is not None:
        LOGGER.debug('Saving index map on file {}'.format(argv.index_map_file))

        if argv.index_map_file == argv.output:
            mode = 'a'
        else:
            mode = 'w'

        with ObsMap(numobs, argv.index_map_file, mode) as obs_map:
            obs_map.save_index_map(scris.generate_index_map(), obs_filter)

    LOGGER.info('Execution complete')
    return 0


if __name__ == '__main__':
    sysexit(main())
