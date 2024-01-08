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
from netCDF4 import Dataset
from preprocessor.fg_generator.fg_generator_utilities.first_guess import NetcdfAtmosphericFirstGuess
from preprocessor.fg_generator.wrf2firstguess.wrf2firstguess_utilities.geometry  import min_distance_indx, dist_on_earth
from scipy.interpolate import interp1d
import numpy as np

__author__     = [ 'Paolo Scaccia <paolo.scaccia@adaptivemeteo.com>']
__copyright__  = "Copyright 2023, Adaptive Meteo S.r.l."
__credits__    = ["Paolo Antonelli","Paolo Scaccia"]
__license__    = "GPL"
__maintainer__ = "Paolo Scaccia"
__email__      = "paolo.scaccia@adaptivemeteo.com"

log = logging.getLogger(__name__)

class ClimatologyReadingError(Exception):
    pass

class ClimatologyBoundsError(Exception):
    pass

class Profile(object):
    def __init__(self, temperature = None, 
                       pressure    = None, 
                       water_vapor = None,
                       ozone       = None):
        self.temperature = temperature
        self.pressure    = pressure
        self.water_vapor = water_vapor
        self.ozone       = ozone
        
        return
        

class ClimatologyGrid(Profile):
    def convert_units(values, var):
        
        return 
    def __init__(self, temp_climatology, h2o_climatology, o3_climatology, precision = False):
        """
            This class contains the climatology grid
            readining the three input files (temperature, water vapor
            and ozone) as input, togheter with the method to extract 
            a profile and the corresponding error near a given point (lat,lon).
            
            Parameters
            ----------
            temp_climatology : str
                Netcdf file with the temperature climatology.
            h2o_climatology : str
                Netcdf file with the water vapor climatology.
            o3_climatology : str
                Netcdf file with the ozone climatology.
            precision : bool, optional
                Flag for precision reading. The default is False.
    
        """
        super().__init__(self)
        try:
            with Dataset(temp_climatology,'r') as temp_file:
                self.temperature     = temp_file['mm_T_values'][:]
                self.n_lats          = temp_file.dimensions['lats'].size
                self.n_lons          = temp_file.dimensions['lons'].size
                self.latlons         = temp_file['latlons'][:]
                self.lats            = self.latlons[:,0].reshape(self.n_lats,self.n_lons)
                self.lons            = self.latlons[:,1].reshape(self.n_lats,self.n_lons)
                self.pressure        = temp_file['pressure'][:]
                if precision:
                    self.temp_precision  = np.abs(temp_file['mm_T_prec'][:])
                
        except:
            raise ClimatologyReadingError("Error in reading temperature climatology")
            
        try:
            with Dataset(h2o_climatology,'r') as wv_file:
                self.water_vapor     = wv_file['mm_H2O_values'][:]*1e4
                if precision:
                    self.wv_precision  = np.abs(wv_file['mm_H2O_prec'][:])*1e4
        except:
            raise ClimatologyReadingError("Error in reading water vapor climatology")

        try:  
            with Dataset(o3_climatology,'r') as ozone_file:
                self.ozone     = ozone_file['mm_O3_values'][:]*1e4
                if precision:
                    self.ozone_precision  = np.abs(ozone_file['mm_O3_prec'][:])*1e4
        except:
            raise ClimatologyReadingError("Error in reading temperature climatology")

        return
    
    def check_pressure_bounds(self, top, bottom):
        """
            This method checks wether the given pressure bounds (top and bottom [hPa])
            fall outside the climatology pressure grid.
        """
        if self.pressure.min() > top:
            raise ClimatologyBoundsError("The given pressure top is outside the climatology pressure grid.\n")
        elif self.pressure.max() < bottom:
            raise ClimatologyBoundsError('The given pressure bottom is outside the climatology pressure grid.\n')
        return
    
    def get_closest_cell(self,lon,lat):
        """
            This method returns the index of the cell closest to the given point.
            
            Args:
                - *lon*: The longitude of the observation
                - *lat*: The latitude of the observation
            Returns:
                A Climatology Grid cell index
        """
        
        indx = min_distance_indx(
                                        lon,
                                        lat,
                                        self.lons,
                                        self.lats,
                                )

        reference_lon = self.lons[indx]
        reference_lat = self.lats[indx]

        dist = dist_on_earth(
                             lon,
                             lat,
                             reference_lon,
                             reference_lat
                             )

        if dist < 100:
            log.debug('Using point ({:.2f},{:.2f}) of the climatology (indices '
                      '({})) for the point ({:.2f},{:.2f}) which is '
                      '{:.2f} Km far.'.format(reference_lat, reference_lon,
                                              indx, lat, lon, dist))

        return np.ravel_multi_index(indx, self.lons.shape)

    
    def get_obs_profile(self, month, lon, lat, 
                        pressure_grid = None, keep_top_climatology = False):
        """
        Given an observation and its position, read its profile
        
        Args:
            - *month*: Month of the year (0-11)
            - *lon*: The longitude of the observation
            - *lat*: The latitude of the observation
            - *pressure_grid*: (optional) Reference pressure grid above which 
                                          to interpolate the extracted profiles
            - *keep_top_climatology*: (optional) Boolean to keep the climatology 
                                      profile above the reference pressure grid
        Returns:
            A Climatology profile over that point
        """

        # Get closest grid cell
        indx = self.get_closest_cell(lon,lat)
        
        if pressure_grid is None:
            # Case with no pressure grid in input: just read the profile
            return Profile(temperature = self.temperature[indx,month,:],
                           water_vapor = self.water_vapor[indx,month,:],
                           ozone       = self.ozone[indx,month,:],
                           pressure    = self.pressure )

        elif not keep_top_climatology:
            # Interpolate above the given pressure grid
            # and cut the climatology profiles above the top pressure
            temp_interp = interp1d(
                                     np.log(self.pressure[::-1]),
                                     self.temperature[indx,month,:][::-1],
                                     kind='linear',
                                     copy=False,
                                     bounds_error=True,
                                     )
            wv_interp = interp1d(
                                     np.log(self.pressure[::-1]),
                                     self.water_vapor[indx,month,:][::-1],
                                     kind='linear',
                                     copy=False,
                                     bounds_error=True,
                                     )
            ozone_interp = interp1d(
                                     np.log(self.pressure[::-1]),
                                     self.ozone[indx,month,:][::-1],
                                     kind='linear',
                                     copy=False,
                                     bounds_error=True,
                                     )
            return Profile(temperature = temp_interp( np.log(pressure_grid)[::-1] )[::-1],
                           water_vapor = wv_interp(   np.log(pressure_grid)[::-1] )[::-1]/1000,  # kg/kg
                           ozone       = ozone_interp( np.log(pressure_grid)[::-1] )[::-1]/1000, # kg/kg
                           pressure    = pressure_grid )

        else:
            # Otherwise, interpolate above the given pressure grid,
            # mantain the climatology above the top and smooth
            # the values to ensure the continuity at the merging point

            raise ClimatologyBoundsError("Not implemented yet!")


    def get_precision(self, month, lon, lat, 
                        pressure_grid = None, keep_top_climatology = False):
        """
        Given an observation and its position, read its precision
        
        Args:
            - *month*: Month of the year (0-363)
            - *lon*: The longitude of the observation
            - *lat*: The latitude of the observation
            - *pressure_grid*: (optional) Reference pressure grid above which 
                                          to interpolate the extracted profiles
            - *keep_top_climatology*: (optional) Boolean to keep the climatology 
                                      profile above the reference pressure grid
        Returns:
            A Climatology precision over that point
        """
        
        # Check if the climatology precision has been loaded
        if not all( [ hasattr(self, 'temp_precision'), 
                      hasattr(self, 'wv_precision'), 
                      hasattr(self, 'ozone_precision')]):
            raise ClimatologyReadingError("Climatology precision not loaded!")
        
        # Get closest grid cell
        indx = self.get_closest_cell(lon,lat)
        
        if pressure_grid is None:
            # Case with no pressure grid in input: just read the precision

            return self.temp_precision[indx,month,:],\
                   self.wv_precision[indx,month,:],\
                   self.ozone_precision[indx,month,:]
        elif not keep_top_climatology:
            # Interpolate above the given pressure grid
            # and cut the climatology precision profiles above the top pressure
            temp_precision_interp = interp1d(
                                     np.log(self.pressure[::-1]),
                                     self.temp_precision[indx,month,:][::-1],
                                     kind='linear',
                                     copy=False,
                                     bounds_error=True,
                                     )
            # Interpolate Relative error: STDEV / Q 
            wv_precision_interp = interp1d(
                                     np.log(self.pressure[::-1]),
                                     (self.wv_precision/self.water_vapor)[indx,month,:][::-1],
                                     kind='linear',
                                     copy=False,
                                     bounds_error=True,
                                     )
            # Interpolate Relative error: STDEV / O3 
            ozone_precision_interp = interp1d(
                                     np.log(self.pressure[::-1]),
                                     (self.ozone_precision/self.ozone)[indx,month,:][::-1],
                                     kind='linear',
                                     copy=False,
                                     bounds_error=True,
                                     )

            return  temp_precision_interp( np.log(pressure_grid)[::-1] )[::-1],\
                    wv_precision_interp( np.log(pressure_grid)[::-1] )[::-1],\
                    ozone_precision_interp( np.log(pressure_grid)[::-1] )[::-1],

        else:
            # Otherwise, interpolate above the given pressure grid,
            # mantain the climatology above the top and smooth
            # the values to ensure the continuity at the merging point

            raise ClimatologyBoundsError("Not implemented yet!")

    def read_and_save_profiles(self, obs_time, month, lons, lats, first_guess_file,
                               pressure_grid = None, keep_top_climatology = False, 
                               surface_pressure = 1013):

        n_lev = self.pressure.size if pressure_grid is None else  pressure_grid.size

        # Prepare the space where the data will be saved
        first_guess = NetcdfAtmosphericFirstGuess(lons, lats, obs_time, n_lev, first_guess_file)

        # Open the first_guess object and prepare it for saving
        # the read data
        with first_guess:
            
                for obs, lat, lon in zip(range(lats.size), lats, lons):
                    log.debug('Looking for the position of the '
                              'observation {}'.format(obs))
                    p = self.get_obs_profile(month, lon, lat, pressure_grid = pressure_grid)

                    # Save the profiles on the first guess object
                    first_guess.pressure_levels[obs, :] = p.pressure[:]
                    first_guess.temperature[obs, :]     = p.temperature[:]
                    first_guess.water_vapour[obs, :]    = p.water_vapor[:]
                    first_guess.ozone[obs, :]           = p.ozone[:]
                    first_guess.skin_temperature[obs]   = p.temperature[0]
                    first_guess.surface_pressure[obs]   = surface_pressure
                    
