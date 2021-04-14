#!/usr/bin/env python

# This file is part of Emiss2FirstGuess.
#
# Emiss2FirstGuess is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Emiss2FirstGuess is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Emiss2FirstGuess. If not, see <http://www.gnu.org/licenses/>.

"""
Emissivity climatology To Mirto FirstGuess converter

This script is an open source software written to use climatological data
for emissivity for generating a first guess for Mirto, an open source software
for elaborating metereological interferometer data

:copyright: 2016 by eXact-lab and Paolo Antonelli
"""

import logging
from argparse import ArgumentParser
from sys import exit, stderr
from traceback import format_exc

import numpy as np
from netCDF4 import Dataset

from constants import OBSERVATIONS
from igbp.igbp import Igbp, IGBP_CLASSES
from land.land_climatology import LandClimatology
from sea.M_climatology import MasudaClimatology
from output_formats.memory import MemoryOutput as Output

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

if __name__ == '__main__':
    log = logging.getLogger()
else:
    log = logging.getLogger(__name__)


def main():
    v_levels = ['debug', 'info', 'warning']

    parser = ArgumentParser()
    parser.add_argument('observations', type=str,
                        help='A valid observation file for Mirto')
    parser.add_argument('igbp', type=str,
                        help='The file with the IGBP map')
    parser.add_argument('land', type=str,
                        help='The file with the climatology for the land')
    parser.add_argument('sea', type=str,
                        help='The file with the climatology for the sea')
    parser.add_argument('output', type=str,
                        help='The first guess NetCDF file that will be'
                             ' generated')
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='the level of verbosity of the software')
    argv = parser.parse_args()
    
    # Prepare the log class
    verbosity = getattr(logging, argv.verbose.upper())
    log.setLevel(verbosity)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - '
                                  '%(funcName)s: %(message)s',
                                  datefmt='%m/%d/%Y %H:%M:%S')

    sh = logging.StreamHandler()
    sh.setLevel(verbosity)
    sh.setFormatter(formatter)
    log.addHandler(sh)

    log.info('Opening observation file')
    try:
        with Dataset(argv.observations, 'r') as obs_file:
            log.debug('Reading longitudes')
            lons = obs_file.variables[OBSERVATIONS.LONGITUDEFIELD][:]
            log.debug('Reading latitudes')
            lats = obs_file.variables[OBSERVATIONS.LATITUDEFIELD][:]
            log.debug('Reading fov angles')
            fovs = obs_file.variables[OBSERVATIONS.FOVANGLEFIELD][:]
    except:
        log.error('Read of the observations failed!')
        log.debug(format_exc())
        return 1

    log.info('Opening IGBP file')
    try:
        igbp = Igbp(argv.igbp)
    except:
        log.error('Read of the IGBP map failed!')
        log.debug(format_exc())
        return 2

    log.info('Opening the file for the land climatology')
    try:
        land_climatology = LandClimatology(argv.land)
    except:
        log.error('Read of the file for the land climatology failed!')
        log.debug(format_exc())
        return 3

    log.info('Opening the file for the sea climatology')
    try:
        masuda_climatology = MasudaClimatology(argv.sea)
    except:
        log.error('Read of the file for the sea climatology failed!')
        log.debug(format_exc())
        return 4

    # Check that the climatology for the sea and for the land use the same
    # wavenumbers to represent the functions
    try:
        same_wn = np.allclose(land_climatology.wavenumbers,
                              masuda_climatology.wavenumbers)
    except ValueError:
        log.debug(format_exc())
        same_wn = False
    if not same_wn:
        log.error('The land climatology emissivity and the sea one use '
                  'different wavenumbers. They are not compatible')
        return 10
    
    # Check that they use the same number of eigenvalues
    n_of_eigenvalues = land_climatology.n_of_eigenvalues
    log.debug('Land climatology uses {} eigenvalues'.format(n_of_eigenvalues))
    sea_n_of_eigenvalues = masuda_climatology.n_of_eigenvalues
    log.debug('Sea climatology usese {} eigenvalues'
              ''.format(sea_n_of_eigenvalues))
    if n_of_eigenvalues != sea_n_of_eigenvalues:
        log.error('Land climatology and sea climatology use a different '
                  'number of eigenvectors')
        return 11

    log.info('Preparing the output file')
    obsnum = lons.size
    try:
        output = Output(
                        argv.output,
                        obsnum,
                        land_climatology.wavenumbers,
                        n_of_eigenvalues,
                        lats,
                        lons
                        )
    except:
        log.error('Impossible to create the output file')
        log.debug(format_exc())
        return 5

    # Create an array with the igbp classes of each point of the observation
    log.debug('Creating IGBP associations')
    igbp_classes = igbp(lons, lats)

    # Open the output object where the data will be saved
    with output:
        # Save each observation
        for i in range(lons.size):
            fov = fovs[i]
            igbp_class = igbp_classes[i]
            if igbp_class == IGBP_CLASSES.OCEAN_WATER:
                log.debug('Saving observation {} as it was on water'.format(i))
                data = masuda_climatology.associate(fov)
                output.save(i, 0, *data)
            else:
                log.debug('Saving observation {} as it was on land'.format(i))
                emiss = land_climatology.associate(igbp_class)
                output.save(i, 1, *emiss.data())


if __name__ == '__main__':
    exit(main())
