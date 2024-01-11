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
import numpy as np
from sys import exit as sysexit

from preprocessor.fg_generator.fg_generator_utilities.wrf_file  import WrfFile
from preprocessor.fg_generator.fg_generator_utilities.geometry  import min_distance_indx, dist_on_earth

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

# Get the main dir of the software (we expect that the current file
# is one level inside the tree directory)
SCRIPT_DIR =path.realpath(__file__).split('amethyst/')[0]+'amethyst'
TEMP_EXTRAPOLATOR_COEFFICIENTS = path.join(
                                           SCRIPT_DIR,
                                           'ancillary/atmosphere/fg_temperature_coefficients.nc'
                                           )

log = logging.getLogger(__name__)


class SourceFile(object):
    """
    A LocalFile source retrieves information for the first guess from a single
    file on the disk (usually, the output of a WRF model).
    
    It divides the array of observations in strides. Each stride refers to the
    same time-step of the model (which is the closer). For each observation,
    the LocalFile looks for the closest point of the model in space inside
    that time-step.
    
    After that, the object has to interpolate the profile of the model on the
    levels that are defined in the level files. For the temperature, this
    is done by the
    :meth:`utilities.level_interpolations.interp_temperature_over_levels`
    function while for the water vapour this is done by the
    :meth:`utilities.level_interpolations.interp_water_vapour_over_levels`
    function.
    
    The coefficients for the temperature extrapolator required by the
    :meth:`utilities.level_interpolations.interp_temperature_over_levels`
    function are read from the file amethyst/ancillary/atmosphere/fg_temperature_coefficients.nc.
    
    Args:
        - *wrf_file*: the path of the file with the data from the WRF model
    """

    def __init__(self, wrf_file):
        self.path             = wrf_file
        if not path.exists(self.path):
            raise ValueError('File {} does not exist'.format(self.path))
        if not path.isfile(self.path):
            raise ValueError('{} is not a regular file'.format(self.path))

        with  WrfFile(self.path) as open_wrf_file:
            self.times    = open_wrf_file.times
            self.wrf_file = open_wrf_file
            self.skin_temperature = open_wrf_file.skin_temperature
            self.surface_pressure = open_wrf_file.surface_pressure

    def __exit__(self):
        del(self.wrf_file)
        return
    
    def get_obs_profile(self, time_step, lon, lat):
        """
        Given an observation and its position, read its profile
        
        Args:
            - *lon*: The longitude of the observation
            - *lat*: The latitude of the observation
            - *time_step*: The time step of the wrf_model
            - *wrf_file*: The output of a wrf model which contains the data
              that must be saved

        Returns:
            A WrfProfile over that point
        """
        wrf_file = self.wrf_file
        wrf_lons = wrf_file.lons[time_step, :]
        wrf_lats = wrf_file.lats[time_step, :]

        indx1, indx2 = min_distance_indx(
                                         lon,
                                         lat,
                                         wrf_lons,
                                         wrf_lats,
                                         )

        reference_lon = wrf_lons[indx1, indx2]
        reference_lat = wrf_lats[indx1, indx2]

        dist = dist_on_earth(
                             lon,
                             lat,
                             reference_lon,
                             reference_lat
                             )

        if dist < 100:
            log.debug('Using point ({:.2f},{:.2f}) of the model (indices '
                      '({}, {})) for the point ({:.2f},{:.2f}) which is '
                      '{:.2f} Km far.'.format(reference_lat, reference_lon,
                                              indx1, indx2, lat, lon, dist))

        if dist >= 100:
            log.warning('Using point ({:.2f},{:.2f}) of the model for the '
                        'point ({:.2f},{:.2f}) which is {:.2f} Km far.'
                        ''.format(reference_lat, reference_lon, lat, lon, dist))
        
        profile = wrf_file.get_profile(time_step, indx1, indx2)
        
        profile.lon = lon
        profile.lat = lat
        
        return profile
    
    def get_closest_timestep(self, time):
        if len(self.times) == 1:
            return 0
        else:
            sysexit("WRF Time selection to be implemented!")
