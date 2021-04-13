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
The observation indices module contains the class ObsMap, which is used to
save the original position of the observations in the CrIS file.
"""

from os import path
import logging
from traceback import format_exc

from netCDF4 import Dataset

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

LOGGER = logging.getLogger(__name__)

# The names of the dimensions, groups and variables
# in the netcdf file
FOVNUM = 'number_of_FOVs'
INDICES = 'indices'
INDEXTABLE = 'index_map'

class ObsMap(object):
    """
    An ObsMap is a wrapper around a NetCDF file that contains an
    index map as described in function generate_index_map.

    Args:
        - *numobs*: The number of observations that will be saved
        - *filepath*: The path of the file that will be written
    """

    def __init__(self, numobs, filepath, mode='w'):
        logging.debug('Writing {}'.format(filepath))

        self.file = filepath
        self.filepointer = None

        if mode != 'a' and path.exists(self.file):
            LOGGER.error('File {} is already present'.format(self.file))
            raise IOError('File already exists')

        with Dataset(self.file, mode) as obsf:

            # Save the dimensions
            try:
                obsf.createDimension(FOVNUM, numobs)
            except RuntimeError:
                error_trace = format_exc().replace('\n', ' | ')
                LOGGER.debug('Dimension {} not created for the following '
                             'reason: '.format(FOVNUM) + error_trace)
            obsf.createDimension(INDICES, 3)

            # Create a function to save the variables
            LOGGER.debug('Creating NetCDF table {}'.format(INDEXTABLE))
            obsf.createVariable(INDEXTABLE, 'i4', (FOVNUM, INDICES),
                                zlib=True, complevel=9, fill_value=-1)
            obsf.sync()

    def __enter__(self):
        if self.filepointer is None:
            self.filepointer = Dataset(self.file, 'a')
            return self
        raise IOError('File already opened!')

    def __exit__(self, *args):
        if self.filepointer is not None:
            self.filepointer.close()
        self.filepointer = None

    def save_index_map(self, index_map, filter=None):
        """
        Save the content of an index map on the netcdf file
        """
        if self.filepointer is None:
            raise IOError('Can not read or write the index map'
                          ' while the file is closed')
        LOGGER.debug('Saving index map on {}'.format(self.file))
        coords = index_map.shape[-1]
        numobs = index_map.size // coords
        index_map = index_map.reshape(numobs, coords)
        self.filepointer.variables[INDEXTABLE][:] = index_map[filter]
