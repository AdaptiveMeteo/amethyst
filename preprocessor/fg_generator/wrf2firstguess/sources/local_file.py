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

from sources.source import Source
from preprocessor.fg_generator.wrf2firstguess.wrf2firstguess_utilities.wrf_file             import WrfFile
from preprocessor.fg_generator.wrf2firstguess.wrf2firstguess_utilities.split_by_time        import split_by_time
from utilities.geometry                                                                     import min_distance_indx, dist_on_earth
from preprocessor.fg_generator.wrf2firstguess.wrf2firstguess_utilities.temp_extrapolator    import TempExtrapolator
from preprocessor.fg_generator.wrf2firstguess.wrf2firstguess_utilities.level_interpolations import interp_temperature_over_levels, \
                                           interp_water_vapour_over_levels
from preprocessor.fg_generator.wrf2firstguess.wrf2firstguess_utilities.climatology          import generate_ozone_profile, TemperatureClimatology, WaterVaporClimatology
from preprocessor.fg_generator.wrf2firstguess.wrf2firstguess_utilities.first_guess          import NetcdfAtmosphericFirstGuess

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


class LocalFile(Source):
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

    @staticmethod
    def get_obs_profile(time_step, lon, lat, wrf_file, temperature_climatology, water_vapor_climatology):
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
        
        temp_indx1 = min_distance_indx(
                                         lon,
                                         lat,
                                         temperature_climatology.lons,
                                         temperature_climatology.lats
                                         )
        
                
        
        return profile

    def __init__(self, wrf_file, h20_climatology , temperature_climatology, o3_climatology ):
        self.path             = wrf_file
        self.h2o_file         = h20_climatology
        self.temperature_file = temperature_climatology
        self.o3_file          = o3_climatology
        
        for file in [self.path, self.h2o_file, self.temperature_file, self.o3_file]:
            if not path.exists(file):
                raise ValueError('File {} does not exist'.format(file))
            if not path.isfile(file):
                raise ValueError('{} is not a regular file'.format(file))


    def read_and_save(self, obs_time, lons, lats,
                      n_levs, top_lev, first_guess_file):
        log.info('Opening WRF file')

        with  WrfFile(self.path) as wrf_file, \
              TemperatureClimatology(self.temperature_file) as temp_climatology, \
              WaterVaporClimatology(self.h2o_file)          as wv_climatology:
                  
            wrf_times = wrf_file.times
            wrf_levels = wrf_file.lev_num

            log.debug('Found {} timesteps'.format(len(wrf_times)))
            split = split_by_time(obs_time, wrf_times)

            num_obs = lons.size

            # Prepare an extrapolator for the temperature
            t_extr = TempExtrapolator(TEMP_EXTRAPOLATOR_COEFFICIENTS)
   
            # Prepare the space where the data will be saved
            first_guess = NetcdfAtmosphericFirstGuess(lons, lats, obs_time,
                                                      n_levs, first_guess_file)

            # Open the first_guess object and prepare it for saving
            # the read data
            with first_guess:
                for time_step, start, end in split:
                    log.debug('Using step {} for observations from {} to '
                              '{}'.format(time_step, start, end))

                    wrf_lats = wrf_file.lats[time_step]
                    wrf_lons = wrf_file.lons[time_step]
                    wrf_month = wrf_file.times[time_step].astype(object).month

                    for obs in range(start, end):
                        log.debug('Looking for the position of the '
                                  'observation {}'.format(obs))
                        p = LocalFile.get_obs_profile(time_step, lons[obs], lats[obs],
                                                      wrf_file, temp_climatology, wv_climatology)

                        # Now we need to create the levels for the first guess
                        m_levs = p.pressure_levels
                        assert p.n_of_levels <= n_levs, 'The model contains '\
                                                        'more levels than the'\
                                                        ' requested ones'
                        assert np.min(m_levs) > top_lev, 'The top_level is '\
                                                         'too low'
                        if p.n_of_levels == n_levs:
                            # If the users want the same number of levels of the
                            # model, use its levels
                            levs = m_levs
                        else:
                            # Otherwise, introduce new levels equispaced in the
                            # log space
                            m_top = np.min(m_levs)
                            levs_to_add = np.linspace(np.log(top_lev),
                                                      np.log(m_top),
                                                      n_levs - p.n_of_levels+1)
                            # We need to revert them back from the log space
                            # and to reverse their order. I also remove the
                            # first element because it is already in the model
                            # (it is the m_top level)
                            levs_to_add = np.exp(levs_to_add)[::-1][1:]
                            levs = np.array(np.hstack([m_levs, levs_to_add]),
                                            dtype=np.float32)

                        t  = interp_temperature_over_levels(p, levs, t_extr)
                        wv = interp_water_vapour_over_levels(p, levs)
                        o3 = generate_ozone_profile(p, levs, wrf_month)

                        # Save the profiles on the first guess object
                        first_guess.pressure_levels[obs, :] = levs[:]
                        first_guess.temperature[obs, :] = t[:]
                        first_guess.water_vapour[obs, :] = wv[:]
                        first_guess.ozone[obs, :] = o3[:]
                        first_guess.skin_temperature[obs] = p.skin_temperature
                        first_guess.surface_pressure[obs] = p.surface_pressure

