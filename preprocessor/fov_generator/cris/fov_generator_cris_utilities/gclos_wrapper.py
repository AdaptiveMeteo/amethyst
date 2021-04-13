# This tool is adapted from Cris2Observation input tools.

"""
cris_wrapper is a module that contains the class CRIS_Wrapper, which is a wrapper
to handle the GCRSO files containing the data regarding the geolocationing of
the satellite
"""

import h5py
import logging
import numpy as np

from preprocessor.fov_generator.cris.fov_generator_cris_utilities.array_reshapers import array_twist, array_spread

__author__    = "Paolo Scaccia and Paolo Antonelli"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__   = ["Paolo Scaccia", "Paolo Antonelli"]
__license__   = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

LOGGER = logging.getLogger(__name__)

class GCLOS_Wrapper(object):
    """
    A wrapper around the GCRSO files. A CRIS_Wrapper reads all the information
    inside a specified GCRSO file and loads it into memory. After that, it
    expose all the data through its interface.

    :ivar str path: The path of the original GCRSO file
    :ivar longs: The longitude of each observation
    :ivar lats: The latitude of each observation
    :ivar times: The time when each observation has been taken
    :ivar sat_azimuth_angle: The azimuth angle of the satellite
    :ivar sat_zenith_angle: The zenith angle of the satellite
    :ivar solar_azimuth_angle: The azimuth angle of the Sun
    :ivar solar_zenith_angle: The zenith angle of the Sun
    """

    def __init__(self, *gcrso_path):
        
        # Input control
        if len(gcrso_path) == 0:
            gcrso_path = './'           # default path
        else:
            gcrso_path = gcrso_path[0]
        
        
        self.path = gcrso_path
        with h5py.File(self.path, 'r') as gc_file:
            geo_all = gc_file['All_Data/CrIS-SDR-GEO_All']


            LOGGER.debug('Reading longitudes')
            self.longs = np.array(geo_all['Longitude'][:], dtype=np.float32)
            LOGGER.debug('Twisting longitudes')
            self.longs = array_twist(self.longs)
            # /-------------------------/
            LOGGER.debug('Reading latitudes')
            self.lats = np.array(geo_all['Latitude'][:], dtype=np.float32)
            LOGGER.debug('Twisting latitudes')
            self.lats = array_twist(self.lats)
            # /-------------------------/                        
            """
            # TIMES OMITTED
            LOGGER.debug('Reading times')
            # Times are saved as microseconds from 01-01-1958
            # Compute the shift in us from Unix Epoc
            unix_epoc = np.datetime64('1970-01-01T00:00:00Z', 'us')
            ietc_epoc = np.datetime64('1958-01-01T00:00:00Z', 'us')
            shift = unix_epoc - ietc_epoc
            # Read the times as pure integers
            time_us = np.array(geo_all['FORTime'], dtype='int64')
            # Convert the table in microseconds from unix epoc
            times = np.array(time_us, dtype='datetime64[us]')
            # Shift the time to make it start from 1958
            self.times = times - shift

            LOGGER.debug('Spreading times')
            self.times = array_spread(self.times)
            # /-------------------------/
            LOGGER.debug('Computing LOS in ECEF coordinates...')
            self.los  = compute_LOS_ECEF(self.sat_range,
                                         self.sat_zenith_angle,
                                         self.sat_azimuth_angle,
                                         self.lats,
                                         self.longs)
            # /-------------------------/
            """
            LOGGER.debug('Reading SatelliteRange')
            sat_range = np.array(
                                         geo_all['SatelliteRange'][:],
                                         dtype=np.float32
                                         )
            LOGGER.debug('Twisting SatelliteRange')
            self.sat_range = array_twist(sat_range)
            # /-------------------------/            
            LOGGER.debug('Reading SatelliteAzimuthAngle')
            sat_azimuth_angle = np.array(
                                         geo_all['SatelliteAzimuthAngle'][:],
                                         dtype=np.float32
                                         )
            LOGGER.debug('Twisting SatelliteAzimuthAngle')
            self.sat_azimuth_angle = array_twist(sat_azimuth_angle)
            # /-------------------------/
            LOGGER.debug('Reading SatelliteZenithAngle')
            sat_zenith_angle = np.array(
                                         geo_all['SatelliteZenithAngle'][:],
                                         dtype=np.float32
                                         )
            LOGGER.debug('Twisting SatelliteZenithAngle')
            self.sat_zenith_angle = array_twist(sat_zenith_angle)
            # /-------------------------/
            LOGGER.debug('Reading SolarAzimuthAngle')
            solar_azimuth_angle = np.array(
                                           geo_all['SolarAzimuthAngle'][:],
                                           dtype=np.float32
                                           )
            LOGGER.debug('Twisting SolarAzimuthAngle')
            self.solar_azimuth_angle = array_twist(solar_azimuth_angle)
            # /-------------------------/
            LOGGER.debug('Reading SolarZenithAngle')
            solar_zenith_angle = np.array(
                                          geo_all['SolarZenithAngle'][:],
                                          dtype=np.float32
                                          )
            LOGGER.debug('Twisting SolarZenithAngle')
            self.solar_zenith_angle = array_twist(solar_zenith_angle)
            # /-------------------------/            
            LOGGER.debug('Reading Scan Position')   # Doesn't require twisting
            self.scan_position      = np.array(
                                          geo_all['SCPosition'][:],
                                          dtype=np.float32
                                          )
            # /-------------------------/            
            LOGGER.debug('Reading Height')
            height      = np.array(
                                          geo_all['Height'][:],
                                          dtype=np.float32
                                          )
            self.height = array_twist(height)
            # /-------------------------/    
