#!/usr/bin/env python

# This file is part of Wrf2firstguess.
#
# Wrf2firstguess is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or1
# (at your option) any later version.

# Wrf2firstguess is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Wrf2firstguess. If not, see <http://www.gnu.org/licenses/>.

"""
WRF output model To Mirto FirstGuess converter

This script is an open source software written to use the data generated
from a wrf model as a first guess for Mirto, an open source software for
elaborating meteorological interferometer data

:copyright: 2016, eXact-lab and Paolo Antonelli
"""

import logging
from argparse import ArgumentParser
from sys import exit as sysexit
from traceback import format_exc

import numpy as np
from netCDF4 import Dataset

from sources.local_file import LocalFile
from amethyst_config import preprocessor_vars

__author__ = 'Stefano Piani'
__copyright__ = "Copyright 2016, Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Stefano Piani"

if __name__ == '__main__':
    log = logging.getLogger()
else:
    log = logging.getLogger(__name__)

# The name of the NetCDF tables that will be read by this script
LATITUDE = 'Latitude'
LONGITUDE = 'Longitude'
TIME = 'Time'

def main():
    v_levels = ['debug', 'info', 'warning']

    parser = ArgumentParser()
    parser.add_argument('observations', type=str,
                        help='A valid observation file for Mirto')
    parser.add_argument('output', type=str,
                        help='The first guess NetCDF file that will be'
                             ' generated')
    parser.add_argument('--source', '-s', type=str, default='local_file',
                        help='Where are the Wrf model data. At the moment, the'
                             ' only accepted value is "local_file"')
    parser.add_argument('--input', '-i', type=str, default=None,
                        help='If the source is "local_file", please specify the'
                             ' path of the input file')
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='the level of verbosity of the software')
    parser.add_argument('--levels', '-l', type=int, default=81,
                        help='The total number of levels')
    parser.add_argument('--top', '-t', type=float, default=0.005,
                        help='The pressure of the heighest level of the '\
                             'output first guess')

    argv = parser.parse_args()

    # Prepare the log class
    verbosity = getattr(logging, argv.verbose.upper())
    log.setLevel(verbosity)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - '
                                  '%(filename)s: %(message)s',
                                  datefmt='%m/%d/%Y %H:%M:%S')

    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(verbosity)
    streamhandler.setFormatter(formatter)
    log.addHandler(streamhandler)

    # Read the input
    if argv.source == 'local_file' and argv.input is None:
        log.error('No input file specified even if the source is "local_file"')
        return 11

    log.info('Opening observation file')
    try:
        with Dataset(argv.observations, 'r') as obs_file:
            log.debug('Reading longitudes')
            lons = np.array(obs_file.variables[LONGITUDE][:], dtype=np.float32)
            log.debug('Reading latitudes')
            lats = np.array(obs_file.variables[LATITUDE][:], dtype=np.float32)
    except:
        log.error('Read of the latitude or longitude failed!')
        log.debug(format_exc())
        return 2

    try:
        with Dataset(argv.observations, 'r') as position_file:
            log.debug('Reading time')
            obs_times = position_file.variables[TIME][:]
            # Convert to numpy datetime
            obs_times.dtype = 'datetime64[ms]'
    except:
        log.error('Read of observation time failed!')
        log.debug(format_exc())
        return 4

    max_date = str(np.max(obs_times))
    min_date = str(np.min(obs_times))
    log.debug('Observations from {} to {}'.format(min_date, max_date))

    if argv.source not in ('local_file',):
        log.error('Unknown source for the WRF files')
        return 3

    if argv.source == 'local_file':
        wrf_source = LocalFile(argv.input)

    n_levs = argv.levels
    top_lev = np.float32(argv.top)

    #try:
    wrf_source.read_and_save(obs_times, lons, lats,
                                 n_levs, top_lev, argv.output)
    #except:
    #    log.error('Error converting data!')
    #   log.debug(format_exc())
    #    return 100

    log.info('Execution complete')
    return 0


if __name__ == '__main__':
    sysexit(main())
