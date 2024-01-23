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

import logging
from os import path
import numpy as np
from netCDF4 import Dataset

from output_formats.output import Output
from constants import FIRSTGUESS

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

log = logging.getLogger(__name__)


class MemoryOutput(Output):
    """
    Memory output is an output object that stores all the information in
    memory. When the object is closed, all the data are copyed on a netcdf
    file on the disk. This significantly improves the execution time because
    it allows to compress the data just once.
    """
    def __init__(self, output_file, obsnum, wavenumbers, n_of_eigenvalues,
                 lats, lons):
        logging.info('Writing {}'.format(output_file))

        self.file = output_file
        self.filepointer = None

        mode = 'w'
        if path.exists(output_file):
            log.debug('File is already present, the content will be appended')
            mode = 'a'

        with Dataset(output_file, mode) as f:
            # Check if the number of fovs is consistent
            if FIRSTGUESS.FOVNUM in f.dimensions:
                fovs = len(f.dimensions[FIRSTGUESS.FOVNUM])
                log.debug('{} already present. Using that one'
                          ''.format(FIRSTGUESS.FOVNUM))
                if fovs != obsnum:
                    raise ValueError('The number of FOVs in the {} file is '
                                     '{} while this software is preparing {} '
                                     'FOVs'.format(output_file, fovs, obsnum))
            else:
                log.debug('Creating dimension {} of size {}'
                          ''.format(FIRSTGUESS.FOVNUM, obsnum))
                f.createDimension(FIRSTGUESS.FOVNUM, obsnum)

            # Do the same for the variable LATITUDE
            if FIRSTGUESS.LATITUDE in f.variables:
                log.debug('{} variable already present in file {}: '
                          'comparing values...'.format(FIRSTGUESS.LATITUDE,
                                                       output_file))
                lats_var = f.variables[FIRSTGUESS.LATITUDE][:]
                if not np.allclose(lats_var, lats):
                    raise IOError('The values of the {} variable on the file '
                                  '{} are different from the expected ones.'
                                  ''.format(FIRSTGUESS.LATITUDE, output_file))
            else:
                log.debug('Creating variable {}'.format(FIRSTGUESS.LATITUDE))
                lats_var = f.createVariable(FIRSTGUESS.LATITUDE,
                                            'f4',
                                            (FIRSTGUESS.FOVNUM,),
                                            zlib=True,
                                            fill_value=1e9
                                            )
                lats_var[:] = lats

            # and for the longitude
            if FIRSTGUESS.LONGITUDE in f.variables:
                log.debug('{} variable already present in file {}: '
                          'comparing values...'.format(FIRSTGUESS.LONGITUDE,
                                                       output_file))
                lons_var = f.variables[FIRSTGUESS.LONGITUDE][:]
                if not np.allclose(lons_var, lons):
                    raise IOError('The values of the {} variable on the file '
                                  '{} are different from the expected ones.'
                                  ''.format(FIRSTGUESS.LONGITUDE,
                                            output_file))
            else:
                log.debug('Creating variable {}'.format(FIRSTGUESS.LONGITUDE))
                lons_var = f.createVariable(FIRSTGUESS.LONGITUDE,
                                            'f4',
                                            (FIRSTGUESS.FOVNUM,),
                                            zlib=True,
                                            fill_value=1e9
                                            )
                lons_var[:] = lons

            # If the SURFGROUP is already present, then the user is
            # requested to delete the file
            if FIRSTGUESS.SURFGROUP in f.groups:
                log.error('Group {} already present! I can not write on '
                          'file {}'.format(FIRSTGUESS.SURFGROUP, output_file))
                raise IOError('Group {} already present! I can not write on '
                              'file {}'
                              ''.format(FIRSTGUESS.SURFGROUP, output_file))

            # Do not remove this sync()! For some strange reason, without the
            # sync all the dimensions that should be FOVNUM become
            # n_of_eigenvectors in the netCDF file
            f.sync()

            log.debug('Creating group {}'.format(FIRSTGUESS.SURFGROUP))
            surf = f.createGroup(FIRSTGUESS.SURFGROUP)

            n_of_wavenumbers = wavenumbers.size
            log.debug('Creating dimension {} of size {}'
                      ''.format(FIRSTGUESS.N_OF_WN, n_of_wavenumbers))
            surf.createDimension(FIRSTGUESS.N_OF_WN, n_of_wavenumbers)

            log.debug('Creating dimension {} of size {}'
                      ''.format(FIRSTGUESS.N_OF_EIGENVALUES, n_of_eigenvalues))
            surf.createDimension(FIRSTGUESS.N_OF_EIGENVALUES, n_of_eigenvalues)

            log.debug('Saving wavenumbers')
            wn = surf.createVariable(
                                     FIRSTGUESS.WAVENUMBERS,
                                     'f4',
                                     (FIRSTGUESS.N_OF_WN,),
                                     zlib=True,
                                     complevel=9,
                                     )
            wn[:] = wavenumbers

            log.debug('Creating table {}'.format(FIRSTGUESS.FUNCTIONS))
            dims = (
                    FIRSTGUESS.FOVNUM,
                    FIRSTGUESS.N_OF_EIGENVALUES,
                    FIRSTGUESS.N_OF_WN
                    )
            surf.createVariable(
                                FIRSTGUESS.FUNCTIONS,
                                'f4',
                                dims,
                                zlib=True,
                                complevel=9,
                               )
            self.functions = np.empty((obsnum,
                                       n_of_eigenvalues,
                                       n_of_wavenumbers),
                                      dtype=np.float32)

            log.debug('Creating table {}'.format(FIRSTGUESS.BIAS))
            dims = (
                    FIRSTGUESS.FOVNUM,
                    FIRSTGUESS.N_OF_WN
                    )
            surf.createVariable(
                                FIRSTGUESS.BIAS,
                                'f4',
                                dims,
                                zlib=True,
                                complevel=9,
                                )
            self.bias = np.empty((obsnum,
                                  n_of_wavenumbers),
                                 dtype=np.float32)

            log.debug('Creating table {}'.format(FIRSTGUESS.COVARIANCE))
            dims = (
                    FIRSTGUESS.FOVNUM,
                    FIRSTGUESS.N_OF_EIGENVALUES,
                    )
            surf.createVariable(
                                FIRSTGUESS.COVARIANCE,
                                'f4',
                                dims,
                                zlib=True,
                                complevel=9,
                                )
            self.covariance = np.empty((obsnum,
                                        n_of_eigenvalues),
                                       dtype=np.float32)

            log.debug('Creating table {}'.format(FIRSTGUESS.LANDWATER))
            surf.createVariable(
                                FIRSTGUESS.LANDWATER,
                                'i4',
                                (FIRSTGUESS.FOVNUM,),
                                zlib=True,
                                complevel=9,
                                fill_value=-1,
                                )
            self.land_or_water = np.empty((obsnum,), dtype=np.bool_)
            f.sync()

    def __enter__(self):
        if self.filepointer is None:
            self.filepointer = Dataset(self.file, 'a')
            return self
        raise IOError('File already opened!')

    def __exit__(self, *args):
        if self.filepointer is not None:
            surf = self.filepointer.groups[FIRSTGUESS.SURFGROUP]

            functions_table = surf.variables[FIRSTGUESS.FUNCTIONS]
            bias_table = surf.variables[FIRSTGUESS.BIAS]
            covariance_table = surf.variables[FIRSTGUESS.COVARIANCE]
            land_or_water_table = surf.variables[FIRSTGUESS.LANDWATER]

            land_or_water_table[:] = self.land_or_water[:]
            functions_table[:] = self.functions
            bias_table[:] = self.bias
            covariance_table[:] = self.covariance

            self.filepointer.close()

        self.filepointer = None

    def save(self, pos, land_or_water, bias, functions, covariance):
        if self.filepointer is None:
            raise IOError('Can not read or write on the file'
                          ' while the file is closed')

        self.land_or_water[pos] = land_or_water
        self.functions[pos, :] = functions
        self.bias[pos, :] = bias
        self.covariance[pos, :] = covariance
