#!/usr/bin/env python

# This file is part of Iasi2observations.
#
# Iasi2observations is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Iasi2observations is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Iasi2observations. If not, see <http://www.gnu.org/licenses/>.

"""
IASI L1C 2 Mirto Observation file converter

This script is an open source software written to convert the IASI L1C files
to a format that Mirto can read.

:copyright: 2016 by eXact-lab and Paolo Antonelli
"""
from __future__ import print_function

import logging
from argparse import ArgumentParser
from sys import stderr, exit as sysexit
from traceback import format_exc

import numpy as np

from preprocessor.fov_generator.fov_generator_utilities.fov_file import FovFile
from preprocessor.fov_generator.fov_generator_utilities.boxes import Rectangle

# Check if the piasi_reader library is installed
try:
    import piasi_reader
except:
    print(format_exc(), file=stderr)
    print('\nERROR: piasi_reader library not found. You can download it from\n'
          '  https://github.com/spiani/piasi_reader\n'
          'EXECUTION ABORTED!\n', file=stderr)
    sysexit(100)

# Get the piasi_reader version and check that is at least 0.9.7
try:
    pr_version = piasi_reader.version
except:
    pr_version = (0, 0, 0)

if pr_version < (0, 9, 7):
    print('\nERROR: The found piasi_reader library is too old. This '
          'software requires at least the version 0.9.7. You can download the '
          'last version from:\n'
          '  https://github.com/spiani/piasi_reader\n'
          'EXECUTION ABORTED!\n', file=stderr)
    sysexit(101)

from piasi_reader.iasi_l1c_native_file import IasiL1cNativeFile

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

if __name__ == '__main__':
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


def prepare_parser():
    v_levels = ['debug', 'info', 'warning']

    parser = ArgumentParser()

    parser.add_argument('iasifile', type=str,
                        help='A IASI L1C file to be converted')
    parser.add_argument('output', type=str,
                        help='The netcdf file that must be created as output')
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='The level of verbosity of the software')
    parser.add_argument('--observation-number', '-n', type=int, default=None,
                        help='The maximum number of observations that will be '
                             'saved. This is useful to test the software using '
                             'less observations')
    parser.add_argument('--cm_treshold','-cmt', type=float, default=None,
                        help='Clear sky fraction needed to process FOV')
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
    return parser


def prepare_logger(verbosity_level):
    verbosity = getattr(logging, verbosity_level)
    LOGGER.setLevel(verbosity)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - '
                                  '%(funcName)s: %(message)s',
                                  datefmt='%m/%d/%Y %H:%M:%S')

    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(verbosity)
    streamhandler.setFormatter(formatter)
    LOGGER.addHandler(streamhandler)


def main():
    parser = prepare_parser()
    argv = parser.parse_args()
    verbosity_level = argv.verbose.upper()
    prepare_logger(verbosity_level)

    LOGGER.info('Reading the IASI file')
    try:
        iasi_file = IasiL1cNativeFile(argv.iasifile)
    except:
        LOGGER.error('Read of the IASI file failed!')
        LOGGER.debug(format_exc())
        return 1

    # Get the latitudes of the observations
    lats = iasi_file.get_latitudes()
    lons = iasi_file.get_longitudes()

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

    LOGGER.info('Total numberof FOVs: {}'.format(len(lats)))
    # obs_filter is a variable that is true when the observation
    # must be saved and is false otherwise
    obs_filter = inside_box[:]


    #Remove Cloudy FOV (according to AVHRR) 
    avhrr_cloud_fraction = iasi_file.get_avhrr_cloud_fractions()

    #PaoloA 26102020
    if argv.cm_treshold:
        cloud_free_mask = avhrr_cloud_fraction < argv.cm_treshold
        clear_entries = np.where(cloud_free_mask)[0]
        LOGGER.info('Number of clear FOVs accroding to AVHRR cloud mask (cf < {}%): {}'
                     ''.format(argv.cm_treshold,len(clear_entries)))
        obs_filter = np.bool_(np.logical_and(obs_filter, cloud_free_mask))
        
    true_entries = np.where(obs_filter)[0]
    LOGGER.info('Keeping {} clear FOVs inside AOI'
                  ''.format(len(true_entries)))

    obs_max = argv.observation_number
    if obs_max is not None and obs_max < lats.size:
        LOGGER.info('Decreasing the number of the observations to {}'
                     ''.format(obs_max))
        if len(true_entries) > obs_max:
            clear_after = true_entries[obs_max]
            obs_filter[clear_after:] = False

    # We count the number of FOVs inside the box to get the number of
    # the FOVs we will save on the output file
    numobs = np.count_nonzero(obs_filter)
    LOGGER.debug('Saving {} observations'.format(numobs))

    # Multiply radiances by 1e5 to be consistent with Mirto units
    radiances = iasi_file.get_radiances() * 1e5
    num_channels = radiances.shape[1]

    LOGGER.info('Writing the output file')
    try:
        fov_file = FovFile(argv.output, numobs, num_channels,'iasi')
    except:
        LOGGER.error('Write of the output file failed!')
        LOGGER.debug(format_exc())
        return 3

    # Convert the times in seconds from epoc
    times = iasi_file.get_obs_times().astype('datetime64[ms]').astype(np.int64)

    wavenumbers = iasi_file.get_channels()
    zenith_angle = iasi_file.get_zenith_angles()
    sol_azimuth_angle = iasi_file.get_solar_azimuth_angles()
    sol_zenith_angle = iasi_file.get_solar_zenith_angles()

    with fov_file:
        fov_file.save_latitude(lats, obs_filter)
        fov_file.save_longitude(lons, obs_filter)
        fov_file.save_time(times, obs_filter)
        fov_file.save_radiance(radiances, obs_filter)
        fov_file.save_wavenumbers(wavenumbers)
        fov_file.save_fov_angle(zenith_angle, obs_filter)
        fov_file.save_sat_zenith_angle(zenith_angle, obs_filter)
        fov_file.save_solar_azimuth_angle(sol_azimuth_angle, obs_filter)
        fov_file.save_solar_zenith_angle(sol_zenith_angle, obs_filter)
        fov_file.save_avhrr_cloud_mask(avhrr_cloud_fraction, obs_filter)

    LOGGER.info('Execution complete')
    return 0


if __name__ == '__main__':
    sysexit(main())
