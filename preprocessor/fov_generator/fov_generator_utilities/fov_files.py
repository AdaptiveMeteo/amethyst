# This file is part of Cris2observations and Iasi2observations.
#
# Cris2observations and IASI2observation are free softwares: you can redistribute 
# it and/or modify it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Both software are distributed in the hope that they will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Cris2observation and Iasi2observation. If not, see <http://www.gnu.org/licenses/>.

"""
Inside the fov_file module there is a class (FovFile) which is
an interface to write nicely a FOV file for Mirto
"""

import logging
from os import path

from netCDF4 import Dataset
import numpy as np

__author__ = 'Stefano Piani'
__copyright__ = "Copyright 2021, Adaptive Meteo S.r.l"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "paolo.scaccia@adaptivemeteo.com"

LOGGER = logging.getLogger(__name__)

# The names of the dimensions, groups and variables
# in the netcdf file
FOVNUM = 'number_of_FOVs'
NCHANNEL = 'number_of_channels'
LATTABLE = 'Latitude'
LONTABLE = 'Longitude'
TIMETABLE = 'Time'
RADTABLE = 'Radiance'
WAVENUMBERS = 'Wavenumbers'
FOVANGLETABLE = 'FOV_angle'
AZIMUTHTABLE = 'Satellite_azimuth_angle'
ZENITHTABLE = 'Satellite_zenith_angle'
SOLAZIMUTHTABLE = 'Solar_azimuth_angle'
SOLZENITHTABLE = 'Solar_zenith_angle'
AVHRRCLOUDFRACTIONTABLE = 'AVHRR_cloud_fraction'



class FovFile(object):
    """
    A FovFile is an object that stores all the data from a CrIS interferometer
    needed by Mirto.

    Before saving the data, the object must be open and when all the data are
    saved, the object must be closed.

    To open and close the object, it is also possible to use the "with"
    statement.

    The variables of a FovFile are array-like objects that store, for each fov,
    information about the position of the satellite and of the Sun, the
    latitude and longitudes of the FOV and the radiances measured by the
    CrIS instrument.

    :ivar pressure_levels: The pressure values (in hPa) that define the levels
    :ivar temperature: The temperature (in K) for each level
    :ivar water_vapour: The water vapour for each level
    :ivar ozone: The ozone for each level
    
    Args:
        - *netcdf_file*
        - *fov_num*: 
        - *n_of_channels*
    """

    def __init__(self, netcdf_file, fov_num, n_of_channels):
        logging.debug('Writing {}'.format(netcdf_file))

        self.file = netcdf_file
        self.filepointer = None

        if path.exists(self.file):
            LOGGER.error('File already present! Execution aborted!')
            raise IOError('File already present!')

        with Dataset(self.file, 'w') as fovf:

            # Save the dimensions
            fovf.createDimension(FOVNUM, fov_num)
            fovf.createDimension(NCHANNEL, n_of_channels)

            # Create a function to save the variables
            def create_var(name, type='f4', dim=(FOVNUM,), fill_val=1e9):
                LOGGER.debug('Creating NetCDF table {}'.format(name))
                table = fovf.createVariable(name, type, dim, zlib=True,
                                            complevel=9, fill_value=fill_val)
                return table

            create_var(LATTABLE)
            create_var(LONTABLE)
            create_var(TIMETABLE, type='i8', fill_val = -1)
            create_var(RADTABLE, dim=(FOVNUM, NCHANNEL))
            create_var(WAVENUMBERS, dim=(NCHANNEL,))
            create_var(FOVANGLETABLE)
            create_var(ZENITHTABLE)
            create_var(SOLAZIMUTHTABLE)
            create_var(SOLZENITHTABLE)
            
            if 'cris' in path.__file__:
                create_var(AZIMUTHTABLE)
            elif 'iasi' in path.__file__:
                # iasi
                create_var(AVHRRCLOUDFRACTIONTABLE, type='u1', fill_val = 0)
            else:
                IOError("Cannot choose between IASI and CRIS. Error with script path!")
            fovf.sync()

    def __enter__(self):
        if self.filepointer is None:
            self.filepointer = Dataset(self.file, 'a')
            return self
        raise IOError('File already opened!')

    def __exit__(self, *args):
        if self.filepointer is not None:
            self.filepointer.close()
        self.filepointer = None

    def save_latitude(self, lats, filter=None):
        if self.filepointer is None:
            raise IOError('Can not read or write the latitude'
                          ' while the file is closed')        
        LOGGER.debug('Saving latitude on {}'.format(self.file))
        self.filepointer.variables[LATTABLE][:] = lats.flatten()[filter]

    def save_longitude(self, lons, filter=None):
        if self.filepointer is None:
            raise IOError('Can not read or write the longitude'
                          ' while the file is closed')        
        LOGGER.debug('Saving longitude on {}'.format(self.file))
        self.filepointer.variables[LONTABLE][:] = lons.flatten()[filter]

    def save_time(self, time, filter=None):
        if self.filepointer is None:
            raise IOError('Can not read or write the time'
                          ' while the file is closed')        
        LOGGER.debug('Saving time on {}'.format(self.file))
        # Ensure that the times is saved as milliseconds
        time_msec = np.array(time ,dtype='datetime64[ms]').flatten()
        # Convert it in integer (before was saved as a date)
        time_int = time_msec.astype(np.int64)

        self.filepointer.variables[TIMETABLE][:] = time_int[filter]

    def save_radiance(self, radiances, filter=None):
        if self.filepointer is None:
            raise IOError('Can not read or write the radiance'
                          ' while the file is closed')        
        LOGGER.debug('Saving radiance on {}'.format(self.file))
        n_elements = radiances.size
        n_channels = radiances.shape[-1]
        n_fovs = n_elements // n_channels
        # Reshape the array so that we have just one dimension for the
        # FOVs and one dimension for the channels
        radiances = radiances.reshape(n_fovs, n_channels)

        self.filepointer.variables[RADTABLE][:] = radiances[filter]

    def save_wavenumbers(self, wavenumbers, filter=None):
        if self.filepointer is None:
            raise IOError('Can not read or write the wavenumbers'
                          ' while the file is closed')        
        LOGGER.debug('Saving wavenumbers on {}'.format(self.file))
        wavenumbers_table = self.filepointer.variables[WAVENUMBERS]
        wavenumbers_table[:] = wavenumbers.flatten()[filter]

    def save_fov_angle(self, angl, filter=None):
        if self.filepointer is None:
            raise IOError('Can not read or write the FOV angle'
                          ' while the file is closed')
        LOGGER.debug('Saving FOV angle on {}'.format(self.file))
        self.filepointer.variables[FOVANGLETABLE][:] = angl.flatten()[filter]

    def save_sat_azimuth_angle(self, angl, filter=None):
        if self.filepointer is None:
            raise IOError('Can not read or write the satellite azimuth angle'
                          ' while the file is closed')        
        LOGGER.debug('Saving satellite azimuth angle on {}'.format(self.file))
        self.filepointer.variables[AZIMUTHTABLE][:] = angl.flatten()[filter]

    def save_sat_zenith_angle(self, angl, filter=None):
        if self.filepointer is None:
            raise IOError('Can not read or write the satellite zenith angle'
                          ' while the file is closed')        
        LOGGER.debug('Saving satellite zenith angle on {}'.format(self.file))
        self.filepointer.variables[ZENITHTABLE][:] = angl.flatten()[filter]

    def save_solar_azimuth_angle(self, angl, filter=None):
        if self.filepointer is None:
            raise IOError('Can not read or write the solar azimuth angle'
                          ' while the file is closed')        
        LOGGER.debug('Saving solar azimuth angle on {}'.format(self.file))
        self.filepointer.variables[SOLAZIMUTHTABLE][:] = angl.flatten()[filter]

    def save_solar_zenith_angle(self, angl, filter=None):
        if self.filepointer is None:
            raise IOError('Can not read or write the solar zenith angle'
                          ' while the file is closed')        
        LOGGER.debug('Saving solar zenith angle on {}'.format(self.file))
        self.filepointer.variables[SOLZENITHTABLE][:] = angl.flatten()[filter]
        
    def save_avhrr_cloud_mask(self, cloud_fraction, filter=None):
        if self.filepointer is None:
            raise IOError('Can not read or write the AVHRR cloud mask'
                          ' while the file is closed')        
        LOGGER.debug('Saving AVHRR cloud mask on {}'.format(self.file))
        mask = cloud_fraction.flatten()[filter]
        self.filepointer.variables[AVHRRCLOUDFRACTIONTABLE][:] = mask