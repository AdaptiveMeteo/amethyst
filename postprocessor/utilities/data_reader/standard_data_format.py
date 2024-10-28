import numpy as np
import re

__author__    = "Paolo Antonelli and Paolo Scaccia"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__   = ["Paolo Antonelli","Paolo Scaccia"]
__license__   = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"




## STANDARD VARIABLES  #########################################################
"""
    Dictionary with the STANDARD format for MIRTO data
"""
###################  INTERNAL NAME              OUTUPUT NAME
STANDARD_FMT = {   'pressure'               : 'air_pressure' ,
                   'temperature'            : 'air_temperature',
                   'water vapor'            : 'air_water_vapor_mr',
                   'relative humidity'      : 'air_relative_humidity',
                   'ozone'                  : 'air_ozone',
                   'surface temperature'    : 'sfc_temperature',
                   'surface pressure'       : 'sfc_pressure',
                   'latitude'               : 'latitude',
                   'longitude'              : 'longitude',
                   'sat azimuth angle'      : 'azimuth_angle',
                   'solar zenith angle'     : 'solar_zenith_angle',
                   'solar azimuth angle'    : 'solar_azimuth_angle',
                   'fov angle'              : 'fov_angle',
                   'first guess residuals'  : 'fg_residuals',
                   'residuals'              : 'residuals',
                   'quality control'        : 'qc',
                   'obs minus back'         : 'omb',
                   'emissivity scores'      : 'ems_coeff',
                   'TR eigenvalues'         : 'DA_Lambda'

                   # Other fields to be added....
                   }
## DATATYPE  ######################################################################
def get_datatype(variable,profile):
    if variable in ['air_pressure' ,
                    'air_temperature',
                    'air_water_vapor_mr',
                    'air_relative_humidity',
                    'air_ozone']:
        # Get Atmospheric Data    
        return atmospheric_data(profile)
    
    elif variable in [ 'sfc_temperature', 'sfc_pressure','ems_scores' ]:
        # Get Surface Data
        return surface_data(profile)
    elif variable in ['fg_residuals','residuals']:
        # Get Spectral Data
        return spectral_data(profile)      
    else:
        raise IOError('Variable not belonging to any datatype')

class atmospheric_data(object):
    def __init__(self,profile):
        self.full_profile = np.copy(profile)
        self.scalar_field = None
        return
    
    def compute_scalar_field(self,mode,axis = 1,lev_axis=1,lev1=0,lev2=None):
        if mode == 'min':
            self.get_min(axis)
        elif mode == 'max':
            self.get_max(axis)
        elif bool(re.match(r'level_(?P<lev>\d+)',mode)):
            lev = int(re.match(r'level_(?P<lev>\d+)',mode).group('lev'))
            self.get_level(lev)
        elif mode == 'vertical_average':
            self.get_level_average(lev1, lev2,lev_axis)
        elif '_average_' in mode:
            match = re.match(r'vertical_average_(?P<lev1>\d+)_(?P<lev2>\d+)',mode)
            if not bool(match):
                raise IOError("Wrong input mode for custom vertical average: vertical_average_lev1_lev2 ")
            lev1 = int(match.group('lev1'))
            lev2 = int(match.group('lev2'))
            self.get_level_average(lev1,lev2,lev_axis)
        
    def get_min(self,axis=1):
        self.scalar_field = self.full_profile.min(axis=axis)
    def get_max(self,axis=1):
        self.scalar_field = self.full_profile.max(axis=axis)
    def get_level(self,level,lev_axis=1):
        if lev_axis == 0:
            self.scalar_field = self.full_profile[level]
        elif lev_axis == 1:
            self.scalar_field = self.full_profile[:,level]
        else:
            raise IOError('Axis not implemented')
    def get_level_average(self,lev1,lev2,lev_axis=1):
        if lev_axis == 1:
            if lev2 == None:
                self.scalar_field = np.mean( self.full_profile[:,lev1:],axis=lev_axis)
            else:
                self.scalar_field = np.mean( self.full_profile[:,lev1:lev2],axis=lev_axis)
                
        else:
            raise IOError('Mean over axis != 1 not yet implemented')
  
    
class surface_data(object):
    def __init__(self,profile):
        self.scalar_field = np.copy(profile)
        return
    
class spectral_data(object):
    def __init__(self,profile):
        self.full_profile = np.copy(profile)
    
    
## COLORMAPS  #####################################################################

mirto_colormaps = { STANDARD_FMT[x] : 'YlOrRd' for x in STANDARD_FMT.keys() }
mirto_colormaps['air_relative_humidity'] = 'BuPu'
