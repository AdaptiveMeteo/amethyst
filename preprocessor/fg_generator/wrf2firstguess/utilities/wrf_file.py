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

import logging
from traceback import format_exc

from netCDF4 import Dataset
from numpy import array, datetime64, rollaxis, power, float32

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

log = logging.getLogger(__name__)

class WrfProfile(object):
    """
    A WrfProfile is a profile of temperature, pressure, water vapour and
    ozone on a particular point located in space and time of a WtrMode

    Args:
        - *lon* : the longitude of the point where the profile is located
        - *lat* : the latitude of the point where the profile is located
        - *n_of_levels*: the number of levels of the WrfModel
        - *pressure_levels*: the pressure levels of the profiles
        - *temperature*: the profile of the temperature over the pressure
          levels
        - *water_vapour*: the profile of the water vapour over the pressure
          levels
        - *skin_temperature*: the value of the temperature on the ground
        - *surface_pressure*: the value of the pressure on the ground

    :ivar lon: The longitude of the point where the profile is taken
    :ivar lat: The latitude of the point where the profile is taken
    :ivar n_of_levels: The numbers of the levels of the profile
    :ivar pressure_levels: The pressure values that define the levels of the
                           profile
    :ivar temperature: An array with the values of the temperature in Kelvin
                       for each level of the profile
    :ivar water_vapour: An array with the values of the water vapour for each
                        level of the profile
    :ivar skin_temperature: The temperature of the ground
    :ivar surface_pressure: The value of the pressure at the ground level
    """

    def __init__(self, lon, lat, n_of_levels, pressure, temperature,
                 water_vapour, skin_temperature, surface_pressure):

        self.lon = lon
        self.lat = lat
        self.n_of_levels = n_of_levels
        self.pressure_levels = pressure
        self.temperature = temperature
        self.water_vapour = water_vapour
        self.skin_temperature = skin_temperature
        self.surface_pressure = surface_pressure


class WrfFile(object):
    """
    A wrapper around a file genertated by a WRF model.

    This object is able to return different parameters of the WRF model as
    numpy arrays.

    This class has been written trying to minimize the disk access. For this
    reason, each table is read just once (the first time somebody tries to
    access it) and then it is cached.

    Moreover, the class is able to tell if the underlying netCDF file handler
    is already opened or not. You are therefore free to ask for a particular
    table in any moment: if the file is opened (with the usual syntax "with
    WrfFile(file_name): ...") then that file handler will be used. Otherwise,
    a new file handler will be created. This is usefull because if more than
    one table is needed it is possible to open the file and to read all of
    them using a single file handler (and therefore avoiding to open and
    close the file multiple times).

    Args:
        - *file_name*: The file with the data of the WRF model

    :ivar times: An array of datetime64 objects that are the time step of the
                 model
    :ivar lons: A 2D array with the longitude of the points of the model grid
    :ivar lats: A 2D array with the latitude of the points of the model grid
    """

    def __init__(self, file_name):
        self.__filename = file_name
        self.__mode = 'r'
        self.__filepointer = None
        self.clean_cache()

    def clean_cache(self):
        """
        The WrfFile keeps in memory all the information that it reads to recall
        it faster later. This function wipes all the data allowing to save
        memory
        """

        self.__dict__['cache'] = {
                                  'times': None,
                                  'lats' : None,
                                  'lons' : None,
                                  'pressure' : None,
                                  'temperature' : None,
                                  'water_vapour' : None,
                                  'skin_temperature' : None,
                                  'surface_pressure' : None,
                                  }

    def __enter__(self):
        log.debug('Opening file {}'.format(self.__filename))
        self.__filepointer = Dataset(self.__filename, self.__mode)
        return self

    def __exit__(self, exc_type, exc_value, tb):
        if self.__filepointer is not None:
            try:
                log.debug('Closing file {}'.format(self.__filename))
                self.__filepointer.close()
            except:
                log.warning(format_exc())

        self.__filepointer = None

    @property
    def times(self):
        if self.cache['times'] is None:
            if self.__filepointer is not None:
                log.debug('Reading "Times" table')
                raw_time = self.__filepointer.variables['Times']
            else:
                log.debug('Opening file {}'.format(self.__filename))
                with Dataset(self.__filename, 'r') as f:
                    log.debug('Reading "Times" table')
                    raw_time = f.variables['Times']

            # The raw times are array of bytes; we have to convert them in
            # strings
            string_times = []
            for t in raw_time:
                t_str = ''.join(k.decode('ASCII') for k in t)
                string_times.append(t_str)
            # Replace the _ with a T to divide date from time
            string_times = [t.replace('_', 'T') + u'Z' for t in string_times]
            self.cache['times'] = array(string_times, dtype=datetime64)
        return self.cache['times']

    @property
    def lats(self):
        if self.cache['lats'] is None:
            if self.__filepointer is not None:
                log.debug('Reading "XLAT" table')
                fp = self.__filepointer
                self.cache['lats'] = array(
                                           fp.variables['XLAT'][:],
                                           dtype=float32
                                           )
            else:
                log.debug('Opening file {}'.format(self.__filename))
                with Dataset(self.__filename, 'r') as f:
                    log.debug('Reading "XLAT" table')
                    self.cache['lats'] = array(
                                               f.variables['XLAT'][:],
                                               dtype=float32,
                                               )
        return self.cache['lats']

    @property
    def lons(self):
        if self.cache['lons'] is None:
            if self.__filepointer is not None:
                log.debug('Reading "XLONG" table')
                fp = self.__filepointer
                self.cache['lons'] = array(
                                           fp.variables['XLONG'][:],
                                           dtype=float32,
                                           )
            else:
                log.debug('Opening file {}'.format(self.__filename))
                with Dataset(self.__filename, 'r') as f:
                    log.debug('Reading "XLONG" table')
                    self.cache['lons'] = array(
                                               fp.variables['XLONG'][:],
                                               dtype=float32,
                                               )
        return self.cache['lons']

    @property
    def water_vapour(self):
        if self.cache['water_vapour'] is None:
            if self.__filepointer is not None:
                log.debug('Reading "QVAPOR" table')
                fp = self.__filepointer
                vapour = array(
                               fp.variables['QVAPOR'][:],
                               dtype=float32
                               )
                # Multiply by 1000 to return values in g/Kg                
            else:
                log.debug('Opening file {}'.format(self.__filename))
                with Dataset(self.__filename, 'r') as f:
                    log.debug('Reading "QVAPOR" table')
                    vapour = array(
                                   fp.variables['QVAPOR'][:],
                                   dtype=float32
                                   )

            vapour *= 1E03

            # move the levels on the last position
            vapour = rollaxis(vapour, 2, 1)
            vapour = rollaxis(vapour, 3, 2)

            self.cache['water_vapour'] = vapour

        return self.cache['water_vapour']

    @property
    def pressure(self):
        if self.cache['pressure'] is None:
            if self.__filepointer is not None:
                log.debug('Reading "P" table')
                p = self.__filepointer.variables['P'][:]
                log.debug('Reading "PB" table')
                pb = self.__filepointer.variables['PB'][:]
            else:
                log.debug('Opening file {}'.format(self.__filename))
                with Dataset(self.__filename, 'r') as f:
                    log.debug('Reading "P" table')
                    p = f.variables['P'][:]
                    log.debug('Reading "PB" table')
                    pb = f.variables['PB'][:]

            pressure = array(p + pb, dtype=float32)

            # Convert to hPa
            pressure *= 0.01

            # move the levels on the last position
            pressure = rollaxis(pressure, 2, 1)
            pressure = rollaxis(pressure, 3, 2)

            self.cache['pressure'] = pressure

        return self.cache['pressure']

    @property
    def temperature(self):
        if self.cache['temperature'] is None:
            if self.__filepointer is not None:
                log.debug('Reading "T" table')
                temperature = self.__filepointer.variables['T'][:] + 300
            else:
                log.debug('Opening file {}'.format(self.__filename))
                with Dataset(self.__filename, 'r') as f:
                    log.debug('Reading "T" table')
                    temperature = f.variables['T'][:] + 300

            # move the levels on the last position
            temperature = rollaxis(temperature, 2, 1)
            temperature = rollaxis(temperature, 3, 2)

            temperature *= power(self.pressure * 1E-3, 2 / 7.)

            self.cache['temperature'] = array(temperature, dtype=float32)

        return self.cache['temperature']

    @property
    def skin_temperature(self):
        if self.cache['skin_temperature'] is None:
            if self.__filepointer is not None:
                log.debug('Reading "TSK" table')
                skin_temperature = self.__filepointer.variables['TSK'][:]
            else:
                log.debug('Opening file {}'.format(self.__filename))
                with Dataset(self.__filename, 'r') as f:
                    log.debug('Reading "TSK" table')
                    skin_temperature = f.variables['TSK'][:]

            self.cache['skin_temperature'] = array(
                                                   skin_temperature,
                                                   dtype=float32,
                                                   )

        return self.cache['skin_temperature']

    @property
    def surface_pressure(self):
        if self.cache['surface_pressure'] is None:
            if self.__filepointer is not None:
                log.debug('Reading "PSFC" table')
                surface_pressure = self.__filepointer.variables['PSFC'][:]
            else:
                log.debug('Opening file {}'.format(self.__filename))
                with Dataset(self.__filename, 'r') as f:
                    log.debug('Reading "PSFC" table')
                    surface_pressure = f.variables['PSFC'][:]

            # Convert to hPa
            surface_pressure *= 0.01

            self.cache['surface_pressure'] = array(
                                                   surface_pressure,
                                                   dtype=float32,
                                                   )

        return self.cache['surface_pressure']

    @property
    def lev_num(self):
        return self.water_vapour.shape[-1]

    def get_profile(self, t, i, j):
        """
        Return the profile for the point of coordinate t, i, j in the model.
        Be carefull that t, i and j are not physical quantities, but just the
        index of the particular point in the model

        Args:
            - *t*: the time-step of the point in the model
            - *i*: the first index of the spacial grid of the model
            - *j*: the second index of the spacial grid of the model

        Return:
            A WrfProfile over that point
        """

        lon = self.lons[t, i, j]
        lat = self.lats[t, i, j]
        n_of_levels = self.lev_num
        pressure = self.pressure[t, i, j]
        temperature = self.temperature[t, i, j]
        water_vapour = self.water_vapour[t, i, j]
        skin_temperature = self.skin_temperature[t, i, j]
        surface_pressure = self.surface_pressure[t, i, j]

        p = WrfProfile(lon, lat, n_of_levels, pressure, temperature,
                       water_vapour, skin_temperature, surface_pressure)

        return p
