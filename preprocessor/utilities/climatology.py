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
from utilities.geometry  import min_distance_indx, dist_on_earth
from preprocessor.fg_generator.fg_generator_utilities.first_guess import NetcdfAtmosphericFirstGuess
from preprocessor.fg_generator.fg_generator_utilities.source  import SourceFile
from scipy.interpolate import interp1d
import numpy as np
from pandas import DatetimeIndex
from supersmoother import SuperSmoother

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

class AdditionalFileReadingError(Exception):
    pass

class Profile(object):
    def __init__(self, temperature = None, 
                       pressure    = None, 
                       water_vapor = None,
                       ozone       = None,
                       skin_temperature = None,
                       surface_pressure = None):
        self.temperature = temperature
        self.pressure    = pressure
        self.water_vapor = water_vapor
        self.ozone       = ozone
        self.skin_temperature = skin_temperature
        self.surface_pressure = surface_pressure
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
                self.water_vapor     = wv_file['mm_H2O_values'][:] # kg/kg
                if precision:
                    self.wv_precision  = np.abs(wv_file['mm_H2O_prec'][:])/self.water_vapor # log(kg/kg)
        except:
            raise ClimatologyReadingError("Error in reading water vapor climatology")

        try:  
            with Dataset(o3_climatology,'r') as ozone_file:
                self.ozone     = np.abs(ozone_file['mm_O3_values'][:]) # kg/kg
                if precision:
                    self.ozone_precision  = np.abs(ozone_file['mm_O3_prec'][:])/self.ozone # log( kg/kg )
        except Exception as e:
            raise ClimatologyReadingError("Error in reading temperature climatology: {}".format(e))

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

    
    def get_obs_profile(self, month, time, lon, lat, 
                        pressure_grid = None,  source_data = None):
        """
        Given an observation and its position, read its profile
        
        Args:
            - *month*: Month of the year (0-11)
            - *time*: The time (np.datetime64[ms]) of the observation
            - *lon*: The longitude of the observation
            - *lat*: The latitude of the observation
            - *pressure_grid*: (optional) Reference pressure grid above which 
                                          to interpolate the extracted profiles
            - *source_data*: (optional) Additional data from wich to read 
                                  atmospheric profiles for the first guess
        Returns:
            A Climatology profile over that point
        """

        # Get closest climatology grid cell
        indx = self.get_closest_cell(lon,lat)
        
        if pressure_grid is None and source_data is None:
            # Case with no pressure grid in input: just read the profile
            return Profile(temperature = self.temperature[indx,month,:],
                           water_vapor = self.water_vapor[indx,month,:], # kg/kg
                           ozone       = self.ozone[indx,month,:],       # kg/kg
                           pressure    = self.pressure )
        else:
            # Case in which a pressure_grid and/or source data are given
            
            # Create climatology interpolators
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
            
            if source_data is None:
                # If an external file is not given then use the
                # return the interpolated climatology profiles

                return Profile(temperature = temp_interp( np.log(pressure_grid)[::-1] )[::-1],
                               water_vapor = wv_interp(   np.log(pressure_grid)[::-1] )[::-1],  # kg/kg
                               ozone       = ozone_interp( np.log(pressure_grid)[::-1] )[::-1], # kg/kg
                               pressure    = pressure_grid )

            else:
                # Otherwise, read the source profile and merge it with the climatology,
                # and then smooth the values to ensure the continuity at the merging point
                
                # Get source timestep
                timestep = source_data.get_closest_timestep(time)
                
                # Read source profile
                source_profile = source_data.get_obs_profile(timestep, lon, lat)

                #print('DEBUG')
                #print('temp clima',self.temperature[indx,month,:])
                #print('ozone clima',self.ozone[indx,month,:])
                #print('pressure clima',self.pressure)

                # Define climatology pressure levels above the source
                n_climatology_levels  = pressure_grid.shape[-1] - source_profile.n_of_levels
                source_top            = np.min(source_profile.pressure_levels)
                climatology_top       = np.min(pressure_grid)
                log_top_pressure_grid = np.linspace( np.log(climatology_top),
                                                     np.log(source_top),
                                                     n_climatology_levels )[::-1]
                top_pressure_grid   = np.exp(log_top_pressure_grid)


                # Retrieve climatology profile for the levels avove the source
                climatology_profile =  Profile(temperature = temp_interp(  log_top_pressure_grid[::-1] )[::-1],
                                               water_vapor = wv_interp(    log_top_pressure_grid[::-1] )[::-1], # kg/kg
                                               ozone       = ozone_interp( log_top_pressure_grid[::-1] )[::-1], # kg/kg
                                               pressure    = top_pressure_grid )

                # Merge all profiles
                merged_pressure    = np.concatenate((source_profile.pressure_levels, top_pressure_grid))
                merged_temperature = np.concatenate((source_profile.temperature, 
                                                     climatology_profile.temperature))
                merged_water_vapor = np.concatenate((source_profile.water_vapour, 
                                                     climatology_profile.water_vapor))

                if merged_pressure.max() > self.pressure.max() or merged_pressure.min() < self.pressure.min() :
                       # If WRF pressure grid starts below the climatology pressure
                       # fill the gap with constant ozone levels
                       ozone_indx = np.logical_and(  merged_pressure <= self.pressure.max(), merged_pressure >= self.pressure.min()) 
                       merged_ozone = np.empty_like(merged_pressure)*np.nan
                       merged_ozone[ ~ ozone_indx  ] = self.ozone[indx,month,0]
                       merged_ozone[ozone_indx] = ozone_interp( np.log(merged_pressure[ozone_indx])[::-1])[::-1] # kg/kg
                       #print('interp pressure',merged_pressure[ozone_indx])

                else:
                       merged_ozone = ozone_interp( np.log(merged_pressure)[::-1])[::-1] # kg/kg
                #print('interp ozone', merged_ozone)

                # Smooth merged profiles
                sm = SuperSmoother() # Define Smoother
                def smooth(x, cut, n_points = 2):
                   sm.fit(merged_pressure, x)
                   
                   return np.concatenate(( x[:cut - n_points:] ,
                                           sm.predict(merged_pressure[cut - n_points:cut + n_points ]),
                                           x[cut + n_points:] ))
                smoothed_temperature = smooth(merged_temperature, cut = source_profile.n_of_levels)
                smoothed_water_vapor = smooth(merged_water_vapor, cut = source_profile.n_of_levels)
                return Profile(temperature = smoothed_temperature,
                               water_vapor = smoothed_water_vapor,
                               ozone       = merged_ozone,
                               pressure    = merged_pressure,
                               skin_temperature = source_profile.skin_temperature,
                               surface_pressure = source_profile.surface_pressure
                               )
            
    def get_precision(self, month, lon, lat):
        """
        Given an observation and its position, read its precision
        
        Args:
            - *month*: Month of the year (0-363)
            - *lon*: The longitude of the observation
            - *lat*: The latitude of the observation
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
        return  self.temp_precision[indx,month,:],  \
                self.wv_precision[indx,month,:],    \
                self.ozone_precision[indx,month,:]
        """
        # Deprecated Version: the interpolation is not 
        # correct neither necessary
        
        # if pressure_grid is None:
        #     # Case with no pressure grid in input: just read the precision

        #     return self.temp_precision[indx,month,:],\
        #            self.wv_precision[indx,month,:],\
        #            self.ozone_precision[indx,month,:]
                   
        # else:
        #     # Interpolate above the given pressure grid
        #     temp_precision_interp = interp1d(
        #                              np.log(self.pressure[::-1]),
        #                              self.temp_precision[indx,month,:][::-1],
        #                              kind='linear',
        #                              copy=False,
        #                              bounds_error=True,
        #                              )
            
        #     # Interpolate Log error: STDEV / Q 
        #     wv_precision_interp = interp1d(
        #                              np.log(self.pressure[::-1]),
        #                              (self.wv_precision/self.water_vapor)[indx,month,:][::-1],
        #                              kind='linear',
        #                              copy=False,
        #                              bounds_error=True,
        #                              )
        #     # Interpolate Log error: STDEV / O3 
        #     ozone_precision_interp = interp1d(
        #                              np.log(self.pressure[::-1]),
        #                              (self.ozone_precision/self.ozone)[indx,month,:][::-1],
        #                              kind='linear',
        #                              copy=False,
        #                              bounds_error=True,
        #                              )

        #     return  temp_precision_interp( np.log(pressure_grid)[::-1] )[::-1],\
        #             wv_precision_interp( np.log(pressure_grid)[::-1] )[::-1],\
        #             ozone_precision_interp( np.log(pressure_grid)[::-1] )[::-1],
        """

    def read_and_save_profiles(self, obs_times, lons, lats, first_guess_file,
                               pressure_grid = None,  keep_top_climatology = False, 
                               source_file = None,    default_surface_pressure = 1013 ):

        
        obs_month = DatetimeIndex(obs_times).month[0] - 1

        n_lev = self.pressure.size if pressure_grid is None else  pressure_grid.size

        # Prepare the space where the data will be saved
        first_guess = NetcdfAtmosphericFirstGuess(lons, lats, obs_times, n_lev, first_guess_file)

        # If the additional file is given, try to read it
        # using different wrappers (right now WrfFile is the only implemented).
        if source_file is None:
            source_data = None
        else:
            try:
                # using different wrappers (right now only WrfFile is implemented)
                # Try to read externel file with the WRF wrapper
                source_data = SourceFile(source_file)
            except:
                raise AdditionalFileReadingError("Source is not a WRF File. Specific wrapper not implemented yet!")
            
        # Open the first_guess object and prepare it for saving
        # the read data
        with first_guess:
            
                for obs, time, lat, lon in zip(range(lats.size), obs_times, lats, lons):
                    log.debug('Looking for the position of the '
                              'observation {}'.format(obs))

                    #print('preparing obs n',obs)
                    # Retrieve profile closest to the observation
                    p = self.get_obs_profile(obs_month, time, lon, lat, 
                                             pressure_grid = pressure_grid,
                                             source_data = source_data)
                    
                    # Read Superificial Values
                    skin_temperature = p.temperature[0] if p.skin_temperature is None else p.skin_temperature
                    surface_pressure = default_surface_pressure if p.surface_pressure is None else p.surface_pressure

                    # Save the profiles on the first guess object
                    first_guess.pressure_levels[obs, :] = p.pressure[:]
                    first_guess.temperature[obs, :]     = p.temperature[:]
                    first_guess.water_vapour[obs, :]    = p.water_vapor[:]
                    first_guess.ozone[obs, :]           = p.ozone[:]
                    first_guess.skin_temperature[obs]   = skin_temperature
                    first_guess.surface_pressure[obs]   = surface_pressure
