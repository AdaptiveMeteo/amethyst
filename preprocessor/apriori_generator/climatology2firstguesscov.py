#!/usr/bin/env python

# This file is part of CovTable2FirstGuessCov.
#
# CovTable2FirstGuessCov is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# CovTable2FirstGuessCov is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with CovTable2FirstGuessCov. If not, see <http://www.gnu.org/licenses/>.

"""
Covariance table To Mirto FirstGuessCovariance converter

This script generates the error covariances of the first guess data that Mirto
uses. It starts from a covariance table: a netcdf file that associates to
each point of a grid a particular covariance matrix. The output of this script
is a netcdf file that is suitable to be used with Mirto.

:copyright: 2016 by eXact-lab and Paolo Antonelli
"""

import logging
from argparse import ArgumentParser
from sys import exit as sysexit
import sys
from traceback import format_exc
from os import path

from netCDF4 import Dataset
import numpy as np
import pandas as pd

from preprocessor.fg_generator.utilities.climatology import ClimatologyGrid
from utilities.geometry import min_distance_indx, dist_on_earth
from amethyst_config import common_vars, preprocessor_vars

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

if __name__ == '__main__':
    LOG = logging.getLogger()
    AMETHYST_PATH = [ x for x in sys.path if path.basename(x) == 'amethyst' ][0]

else:
    LOG = logging.getLogger(__name__)

LONGITUDE    = 'Longitude'
LATITUDE     = 'Latitude'
NATMLEVELS   = 'number_of_atmospheric_levels'
ATMGROUP     = 'atmospheric_components'
COVGROUP     = 'Covariances'
NFOVS        = 'number_of_FOVs'
ASSOCIATIONS = 'Associations'
PROFILES     = 'Profiles'
PRESSVAR     = 'p'
ZOVERHVAR    = 'z_over_H'
TIME         = 'Time'

def main():
    v_levels = ['debug', 'info', 'warning']

    parser = ArgumentParser()
    parser.add_argument('observations', type=str,
                        help='A valid observation file for Mirto')
    parser.add_argument('output', type=str,
                        help='The first guess NetCDF file that will be '
                             'generated')
    parser.add_argument('--firstguess', '-f', required = True,
                        help='Amethyst First Guess')
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='the level of verbosity of the software')
    parser.add_argument('--h2o', type=str, default=preprocessor_vars["climatology"]["h2o"],
                        help='Water Vapor climatology file')
    parser.add_argument('--temp', type=str, default=preprocessor_vars["climatology"]["temperature"],
                        help='Temperature climatology file')
    parser.add_argument('--o3', type=str, default=preprocessor_vars["climatology"]["o3"],
                        help='Ozone climatology file')
    parser.add_argument('--compression', '-c', type=int, default=4,
                        help='The level of compression of the output file. By '
                             'default is 4, 0 means "no compression", 9 is the '
                             'maximum')
    argv = parser.parse_args()

    # Prepare the log class
    verbosity = getattr(logging, argv.verbose.upper())
    LOG.setLevel(verbosity)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - '
                                  '%(filename)s: %(message)s',
                                  datefmt='%m/%d/%Y %H:%M:%S')

    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(verbosity)
    streamhandler.setFormatter(formatter)
    LOG.addHandler(streamhandler)

    # Read the compression level
    if argv.compression < 0:
        LOG.error('The compression level must be positive')
        return 10
    elif argv.compression > 9:
        LOG.error('The compression level must be less than 9')
        return 11
    else:
        cmp_level = argv.compression
        if cmp_level == 0:
            enable_cmp = False
            LOG.debug('Compression disabled')
        else:
            enable_cmp = True
            LOG.debug('Compression enabled with level {}'.format(cmp_level))

        
    # Read the input
    LOG.info('Opening observation file')
    try:
        with Dataset(argv.observations, 'r') as obs_file:
            LOG.debug('Reading longitudes')
            lons = obs_file.variables[LONGITUDE][:]
            LOG.debug('Reading latitudes')
            lats = obs_file.variables[LATITUDE][:]
            LOG.debug('Reading time')
            obs_times = obs_file.variables[TIME][:]
            # Convert to numpy datetime
            obs_times.dtype = 'datetime64[ms]'
            month = pd.DatetimeIndex(obs_times).month[0] - 1

    except:
        LOG.error('Read of the latitude or longitude failed!')
        LOG.debug(format_exc())
        return 2

    # Read the first guess
    LOG.info('Opening first guess file')
    try:
        with Dataset(argv.firstguess, 'r') as obs_file:
            LOG.debug('Reading pressure')
            pressure_grid = obs_file.groups[ATMGROUP].variables[PRESSVAR][:]

    except:
        LOG.error('Read of the first guess pressure failed!')
        LOG.debug(format_exc())
        return 2

    LOG.info('Saving output on file {}'.format(argv.output))
    mode = 'w'
    if path.exists(argv.output):
        LOG.info('File is already present, the content will be appended')
        mode = 'a'

    # Open the output and the covtable
    with Dataset(argv.output, mode) as output_f:
        # Get the number of the FOVs
        numobs = lats.size
        # If the dimension NFOVS is already present, check that it is
        # consistent. Otherwise, create it
        if NFOVS in output_f.dimensions:
                output_file_nfovs = len(output_f.dimensions[NFOVS])
                if output_file_nfovs != numobs:
                    raise ValueError('The number of FOVs in the {} file is '
                                     '{} while this software is preparing {} '
                                     'FOVs'.format(argv.output,
                                                   output_file_nfovs,
                                                   numobs))
                LOG.debug('{} already present. Using that one'.format(NFOVS))
        else:
            LOG.debug('Creating dimension {} of size {}'.format(NFOVS,
                                                                numobs))
            output_f.createDimension(NFOVS, numobs)

        # Do the same for the variable LATITUDE
        if LATITUDE in output_f.variables:
            LOG.debug('{} variable already present in file {}: '
                      'comparing values...'.format(LATITUDE, argv.output))
            lats_var = output_f.variables[LATITUDE][:]
            if not np.allclose(lats_var, lats):
                error_mess = 'The values of the {} variable on the file '\
                             '{} are different from the expected ones.'\
                             ''.format(LATITUDE, argv.output)
                LOG.error(error_mess)
                raise IOError(error_mess)
        else:
            LOG.debug('Creating variable {}'.format(LATITUDE))
            lats_var = output_f.createVariable(LATITUDE,
                                               'f4',
                                               (NFOVS,),
                                               zlib=True,
                                               fill_value=1e9
                                               )
            lats_var[:] = lats

        # The same for the longitude
        if LONGITUDE in output_f.variables:
            LOG.debug('{} variable already present in file {}: '
                      'comparing values...'.format(LONGITUDE, argv.output))
            lons_var = output_f.variables[LONGITUDE][:]
            if not np.allclose(lons_var, lons):
                error_mess = 'The values of the {} variable on the file '\
                             '{} are different from the expected ones.'\
                             ''.format(LONGITUDE, argv.output)
                LOG.error(error_mess)
                raise IOError(error_mess)
        else:
            LOG.debug('Creating variable {}'.format(LONGITUDE))
            lons_var = output_f.createVariable(LONGITUDE,
                                               'f4',
                                               (NFOVS,),
                                               zlib=True,
                                               fill_value=1e9
                                               )
            lons_var[:] = lons

        # If does not exist, create the atmospheric group
        if ATMGROUP in output_f.groups:
            LOG.debug('{} group already present in the file {}'
                      ''.format(ATMGROUP, argv.output))
            atm = output_f.groups[ATMGROUP]
        else:
            LOG.debug('Creating group {}'.format(ATMGROUP))
            atm = output_f.createGroup(ATMGROUP)

        # Create a COVGROUP inside the atmospheric one
        if COVGROUP in atm.groups:
            error_msg = '{} already present in the output file. Execution'\
                        ' aborted'.format(COVGROUP)
            LOG.error(error_msg)
            raise IOError(error_msg)
        else:
            LOG.debug('Creating group {}'.format(COVGROUP))
            cov_group = atm.createGroup(COVGROUP)

        nlevs = common_vars['levels']
        mols = ['T', 'q', 'O3']

        if NATMLEVELS in atm.dimensions:
            output_file_nlevs = len(atm.dimensions[NATMLEVELS])
            if output_file_nlevs != nlevs:
                err_msg = 'The number of FOVs in the {} file is {} while '\
                          'this software is preparing {} FOVs'\
                          ''.format(argv.output, output_file_nfovs, numobs)
                LOG.error(err_msg)
                raise IOError(err_msg)
            LOG.debug('{} already present. Using that one'
                      .format(NATMLEVELS))
        else:
            LOG.debug('Creating dimension {} of size {}'.format(NATMLEVELS, nlevs))
            output_f.createDimension(NATMLEVELS, nlevs)

        output_tables = {}
        for mol in mols:
            LOG.debug('Saving table {} in the output file'.format(mol))
            table_name = mol
            table_type = 'f4'  # floating point, 32 bytes
            table_dims = (NFOVS, NATMLEVELS, NATMLEVELS)
            first_dim_chunk_size = min(50, numobs)
            chunks = (first_dim_chunk_size, nlevs, nlevs)
            new_table = cov_group.createVariable(
                                                 table_name,
                                                 table_type,
                                                 table_dims,
                                                 zlib=enable_cmp,
                                                 chunksizes=chunks,
                                                 complevel=cmp_level
                                                 )
            output_tables[mol] = new_table


        output_tables['T_q'] = cov_group.createVariable('T_q',
                                                      'f4',
                                                      (NFOVS, NATMLEVELS, NATMLEVELS),
                                                      zlib=enable_cmp,
                                                      chunksizes=chunks,
                                                      complevel=cmp_level)
        output_tables['T_q'][:] = np.zeros_like(output_tables['T_q'][:])

        # Read Known Climatology
        climatology = ClimatologyGrid(argv.temp, argv.h2o, argv.o3, precision = True)
        
        # Now, for each FOV, the closest precision is extracted from
        # the climatology grid
        for i in range(lats.size):
            # For each molecule, look for the profile number that must
            # be used in the association tables
            
            # Read Climatology Precision
            temp_precision, \
            wv_precision, \
            ozone_precision = climatology.get_precision(month, 
                                                        lons[i],lats[i],
                                                        pressure_grid = np.array(pressure_grid[i]))
            
            for precision, mol in zip([temp_precision, wv_precision, ozone_precision],mols):
                # Save Covariance Matrix as diagonal matrix using climatology precision
                output_tables[mol][i, :] = np.diag(precision**2)

    LOG.info('Execution complete')
    return 0



if __name__ == '__main__':
    sysexit(main())

