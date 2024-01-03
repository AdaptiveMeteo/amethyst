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
import pandas as pd
from netCDF4 import Dataset

from amethyst_config import common_vars, preprocessor_vars
from preprocessor.fg_generator.fg_generator_utilities.climatology import ClimatologyGrid
 
__author__ = [ 'Paolo Scaccia <paolo.scaccia@adaptivemeteo.com>']
__copyright__ = "Copyright 2023, Adaptive Meteo S.r.l."
__credits__ = ["Paolo Antonelli","Paolo Scaccia"]
__license__ = "GPL"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"

if __name__ == '__main__':
    log = logging.getLogger()
else:
    log = logging.getLogger(__name__)

# The name of the NetCDF tables that will be read by this script
LATITUDE  = 'Latitude'
LONGITUDE = 'Longitude'
TIME      = 'Time'

def main():
    v_levels = ['debug', 'info', 'warning']

    parser = ArgumentParser()
    parser.add_argument('observations', type=str,
                        help='A valid observation file for Mirto')
    parser.add_argument('output', type=str,
                        help='The first guess NetCDF file that will be'
                             ' generated')
    parser.add_argument('--source', '-s', type=str, required = False, default = None,
                        help='Where are the Wrf model data.')
    parser.add_argument('--input', '-i', type=str, default=None,
                        help='If the source is "local_file", please specify the'
                             ' path of the input file')
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='the level of verbosity of the software')
    parser.add_argument('--levels', '-l', type=int, default=common_vars['levels'],
                        help='The total number of levels')
    parser.add_argument('--top', '-t', type=float, default=preprocessor_vars['first_guess']['top_pressure'],
                        help='The pressure of the heighest level of the '\
                             'output first guess')
    parser.add_argument('--bottom', '-b', type=float, default=preprocessor_vars['first_guess']['bottom_pressure'],
                        help='The pressure of the lower level of the '\
                             'output first guess')
    parser.add_argument('--h2o', type=str, default=preprocessor_vars["climatology"]["h2o"],
                        help='Water Vapor climatology file')
    parser.add_argument('--temp', type=str, default=preprocessor_vars["climatology"]["temperature"],
                        help='Temperature climatology file')
    parser.add_argument('--o3', type=str, default=preprocessor_vars["climatology"]["o3"],
                        help='Ozone climatology file')

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
            day_of_year = pd.DatetimeIndex(obs_times).dayofyear[0] - 1
    except:
        log.error('Read of observation time failed!')
        log.debug(format_exc())
        return 4

    max_date = str(np.max(obs_times))
    min_date = str(np.min(obs_times))
    log.debug('Observations from {} to {}'.format(min_date, max_date))
    log.debug('Climatology Files: \n'
              '                   {} (Temperature) \n'
              '                   {} (Water Vapor) \n'
              '                   {} (Ozone)'.format(argv.temp,argv.h2o,argv.o3))

    # Define First Guess Pressure Grid
    if argv.source != None:
        # SOURCE READING: To be Implemented!
        log.error('Not yet implemented!')
        return 99
        # wrf_source = LocalFile(argv.input, argv.h2o, argv.temperature, argv.o3)
    else:
        # Otherwise, introduce new levels equispaced in the
        # log space
        pressure_grid = np.linspace(np.log(argv.top),
                                    np.log(argv.bottom),
                                    argv.levels + 1)
        pressure_grid = np.exp(pressure_grid)[::-1]
        
    # Read input Climatology grid: argv.temp for temperature, argv.h20 
    #                        for water vapor and argv.o3 for ozone.
    known_climatology = ClimatologyGrid(argv.temp, argv.h2o, argv.o3)
    
    # Check consistency between the given pressure grid and climatology profiles
    # in order to avoid extrapolations.
    known_climatology.check_pressure_bounds(argv.top, argv.bottom)
    
    # Extract Climatology profiles at each observation site (lat, lon) 
    # and save first guess file (argv.output).
    known_climatology.read_and_save(obs_times, 
                                    day_of_year, 
                                    lons, lats, 
                                    argv.output)

    log.info('Execution complete')
    return 0


if __name__ == '__main__':
    sysexit(main())
