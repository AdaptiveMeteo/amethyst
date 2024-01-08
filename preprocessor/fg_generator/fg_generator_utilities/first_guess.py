# This file is part of Wrf2firstguess.
#
# Wrf2firstguess is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Wrf2firstguess is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Wrf2firstguess. If not, see <http://www.gnu.org/licenses/>.

from os import path
import logging

from netCDF4 import Dataset
import numpy as np

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

log = logging.getLogger(__name__)

# The names of the dimensions, groups and variables
# in the netcdf file
FOVNUM = 'number_of_FOVs'
LATITUDE = 'Latitude'
LONGITUDE = 'Longitude'
TIME = 'Time'
LEVNUM = 'number_of_atmospheric_levels'
ATMGROUP = 'atmospheric_components'
TEMPERATURE = 'T'
PRESSURELEVELS = 'p'
WATERVAPOUR = 'q'
OZONE = 'O3'
SKINTEMPERATURE = 'skT'
SURFACEPRESSURE = 'sp'

class AtmosphericFirstGuess(object):
    """
    An AtmosphericFirstGuess is an object that stores all the data required to
    generate the atmospheric part of a first guess. It could write the data
    on the disk or save them in memory.

    Before saving the data, the object must be open and when all the data are
    saved, the object must be closed.

    To open and close the object, it is also possible to use the "with"
    statement.

    The variables of an AtmosphericFirstGuess are array-like objects that
    store, for each observation of the input, the values of a particular
    physical quantity over the different levels defined by the pressure.

    :ivar pressure_levels: The pressure values (in hPa) that define the levels
    :ivar temperature: The temperature (in K) for each level
    :ivar water_vapour: The water vapour for each level
    :ivar ozone: The ozone for each level
    """
    
    def __init__(self):
        raise NotImplementedError
    
    def __enter__(self):
        self.open()
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        pass

    def open(self):
        pass

    @property
    def pressure_levels(self):
        raise NotImplementedError

    @property
    def temperature(self):
        raise NotImplementedError

    @property
    def water_vapour(self):
        raise NotImplementedError

    @property
    def ozone(self):
        raise NotImplementedError


class NetcdfAtmosphericFirstGuess(AtmosphericFirstGuess):
    """
    An NetcdfAtmosphericFirstGuess is an AtmosphericFistGuess that saves all
    the data on a NetCDF file.

    Args:
        - *lons*: the longitude the profiles refers to
        - *lats*: the latitude the profiles refers to
        - *times*: the date and time of the moment when the profiles where
          generated
        - *levsnum*: the number of levels of each profile
        - *output_file*: the path of the NetCDF file where all data will be
          saved
    """
    def __init__(self, lons, lats, times, levsnum, output_file):
        logging.info('Writing {}'.format(output_file))

        self.file = output_file
        self.filepointer = None

        mode = 'w'
        if path.exists(output_file):
            log.info('File is already present, the content will be appended')
            mode = 'a'

        with Dataset(output_file, mode) as f:
            # If the dimension FOVNUM is already present, check that it is
            # consistent. Otherwise, create it
            if FOVNUM in f.dimensions:
                fovs = len(f.dimensions[FOVNUM])
                if fovs != lons.size:
                    raise ValueError('The number of FOVs in the {} file is '
                                     '{} while this software is preparing {} '
                                     'FOVs'.format(output_file,fovs,lons.size))
                log.debug('{} already present. Using that one'.format(FOVNUM))
            else:
                log.debug('Creating dimension {} of size {}'.format(FOVNUM,
                                                                    lons.size))
                f.createDimension(FOVNUM, lons.size)

            # Do the same for the variable LATITUDE
            if LATITUDE in f.variables:
                log.debug('{} variable already present in file {}: '
                          'comparing values...'.format(LATITUDE, output_file))
                lats_var = f.variables[LATITUDE][:]
                if not np.allclose(lats_var, lats):
                    raise IOError('The values of the {} variable on the file '
                                  '{} are different from the expected ones.'
                                  ''.format(LATITUDE, output_file))
            else:
                log.debug('Creating variable {}'.format(LATITUDE))
                lats_var = f.createVariable(LATITUDE,
                                           'f4',
                                           (FOVNUM,),
                                           zlib=True,
                                           fill_value=1e9
                                           )
                lats_var[:] = lats

            # For the longitude
            if LONGITUDE in f.variables:
                log.debug('{} variable already present in file {}: '
                          'comparing values...'.format(LONGITUDE, output_file))
                lons_var = f.variables[LONGITUDE][:]
                if not np.allclose(lons_var, lons):
                    raise IOError('The values of the {} variable on the file '
                                  '{} are different from the expected ones.'
                                  ''.format(LONGITUDE, output_file))
            else:
                log.debug('Creating variable {}'.format(LONGITUDE))
                lons_var = f.createVariable(LONGITUDE,
                                           'f4',
                                           (FOVNUM,),
                                           zlib=True,
                                           fill_value=1e9
                                           )
                lons_var[:] = lons

            # and for the time (after it has been converted it in milliseconds)
            times_msec = np.array(times, dtype='datetime64[ms]')
            times_int = times_msec.astype(np.int64)
            if TIME in f.variables:
                log.debug('{} variable already present in file {}: '
                          'comparing values...'.format(TIME, output_file))
                times_var = f.variables[TIME][:]
                if not np.allclose(times_int, times_var):
                    raise IOError('The values of the {} variable on the file '
                                  '{} are different from the expected ones.'
                                  ''.format(TIME, output_file))
            else:
                log.debug('Creating variable {}'.format(TIME))
                times_var = f.createVariable(TIME,
                                            'i8',
                                            (FOVNUM,),
                                            zlib=True,
                                            fill_value=-1
                                            )
                times_var[:] = times_int

            if ATMGROUP not in f.groups:
                log.debug('Creating group {}'.format(ATMGROUP))
                atm = f.createGroup(ATMGROUP)
            else:
                log.debug('Group {} already present!'.format(ATMGROUP))
                atm = f.groups[ATMGROUP]

            # In the atm group, save the number of levels as a dimension, and
            # the variables of the profile as variables
            if LEVNUM in atm.dimensions:
                levsnum_file = len(atm.dimensions[LEVNUM])
                if levsnum_file != levsnum:
                    error_string = 'The number of levels of the file {} ({}) '\
                                   'is different than the one expected ({})'\
                                   ''.format(levsnum_file, levsnum)
                    log.error(error_string)
                    raise IOError(error_string)
            else:
                log.debug('Creating dimension {} of size {}'.format(LEVNUM,
                                                                    levsnum))
                atm.createDimension(LEVNUM, levsnum)

            # A function that create a variable inside the atmospheric group.
            # If a variable with that name already exists, it returns an error
            def create_var_in_atm(name, unit, type='f4', dims=(FOVNUM, LEVNUM),
                                  zlib=True, fill_value=-99999, min=None, max=None):
                if name in atm.variables:
                    error_string_raw = 'Variable {} already defined in group {}'
                    error_string = error_string_raw.format(name, ATMGROUP)
                    log.error(error_string)
                    raise IOError(error_string)

                log.debug('Creating variable {} in group {}'.format(name,
                                                                    ATMGROUP))
                netcdf_var = atm.createVariable(
                                                name,
                                                type,
                                                dims,
                                                zlib=zlib,
                                                fill_value=fill_value
                                               )
                netcdf_var.units=unit
                if min is not None:
                    if type=='f4':
                        netcdf_var.valid_min = np.float32(min)
                    else:
                        netcdf_var.valid_min = min
                if max is not None:
                    if type=='f4':
                        netcdf_var.valid_max = np.float32(max)
                    else:
                        netcdf_var.valid_max = max

                return netcdf_var

            # Save pressure levels
            create_var_in_atm(PRESSURELEVELS, 'hPa', max=1500.)

            # Save temperature
            create_var_in_atm(TEMPERATURE, 'K', min=170., max=360.)

            # Save water vapour
            create_var_in_atm(WATERVAPOUR, 'kg/kg', min = 0., max = 1.0)

            # Save ozone
            create_var_in_atm(OZONE, 'kg/kg', min=0., max=1)

            # Save skin temperature
            create_var_in_atm(SKINTEMPERATURE, 'K', dims=(FOVNUM,),
                              min=170, max=360)

            # Save surface pressure
            create_var_in_atm(SURFACEPRESSURE, 'hPa', dims=(FOVNUM,),
                              min=100, max=1500)
            f.sync()
            

    def __enter__(self):
        if self.filepointer is None:
            self.filepointer = Dataset(self.file, 'a')
            return self
        raise IOError('File already opened!')

    def __exit__(self, *args):
        if self.filepointer is not None:
            self.filepointer.close()
        self.filepointer = None

    @property
    def pressure_levels(self):
        if self.filepointer is None:
            raise IOError('Can not read or write the pressure levels'
                          ' while the file is closed')

        atm = self.filepointer.groups[ATMGROUP]
        return atm.variables[PRESSURELEVELS]

    @property
    def temperature(self):
        if self.filepointer is None:
            raise IOError('Can not read or write the temperature variable'
                          ' while the file is closed')

        atm = self.filepointer.groups[ATMGROUP]
        return atm.variables[TEMPERATURE]

    @property
    def water_vapour(self):
        if self.filepointer is None:
            raise IOError('Can not read or write the water vapour variable'
                          ' while the file is closed')

        atm = self.filepointer.groups[ATMGROUP]
        return atm.variables[WATERVAPOUR]


    @property
    def ozone(self):
        if self.filepointer is None:
            raise IOError('Can not read or write the ozone variable'
                          ' while the file is closed')

        atm = self.filepointer.groups[ATMGROUP]
        return atm.variables[OZONE]


    @property
    def skin_temperature(self):
        if self.filepointer is None:
            raise IOError('Can not read or write the skin temperature'
                          ' while the file is closed')

        atm = self.filepointer.groups[ATMGROUP]
        return atm.variables[SKINTEMPERATURE]

    @property
    def surface_pressure(self):
        if self.filepointer is None:
            raise IOError('Can not read or write the surface pressure'
                          ' while the file is closed')

        atm = self.filepointer.groups[ATMGROUP]
        return atm.variables[SURFACEPRESSURE]
