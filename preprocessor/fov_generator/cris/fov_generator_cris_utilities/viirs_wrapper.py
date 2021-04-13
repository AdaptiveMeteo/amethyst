"""
viirs_wrapper is a module that contains the class VIIRS_Wrapper, which is a wrapper
to handle the GITCO/GMTCO files containing the data regarding the geolocationing of
the satellite
"""
import h5py
import numpy as np
import logging

__author__    = "Paolo Scaccia and Paolo Antonelli"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__   = ["Paolo Scaccia", "Paolo Antonelli"]
__license__   = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

LOGGER = logging.getLogger(__name__)

FILE_TYPE = {'GMTCO' : 'VIIRS-MOD-GEO-TC_All' ,
             'GITCO' : 'VIIRS-IMG-GEO-TC_All' ,
             'GMODO' : 'VIIRS-MOD-GEO_All' ,
             'IICMO' : 'VIIRS-MOD-GEO_All' }


class VIIRS_Wrapper(object):
    """
    A wrapper around the GITCO/GMTCO files. A VIIRS_Wrapper reads all the information
    inside the specified file and loads it into memory. After that, it
    expose all the data through its interface.

    :ivar str path: The path of the original GITCO/GMTCO file
    :ivar longs: The longitude of each observation
    :ivar lats: The latitude of each observation
    :ivar times: The time when each observation has been taken
    :ivar sat_azimuth_angle: The azimuth angle of the satellite
    :ivar sat_zenith_angle: The zenith angle of the satellite
    :ivar solar_azimuth_angle: The azimuth angle of the Sun
    :ivar solar_zenith_angle: The zenith angle of the Sun
    """

    def __init__(self,*viirs_files_path):

        # Input control
        if len(viirs_files_path) == 0:
            viirs_path = './'          # default path
        elif len(viirs_files_path) == 1:
            viirs_path = viirs_files_path[0]
            print('GMODO file selected: all VIIRS CM information collected from single file')
        else:
            viirs_path = viirs_files_path[0]
            cm_path = viirs_files_path[1]
            self.cm_path  = cm_path
            cm_filename  = cm_path.split('/')[-1]
            print('VIIRS CM information collected from {} file'.format(cm_filename))

            
        self.viirs_los = np.array([])
        self.path = viirs_path
        
        
        filename  = viirs_path.split('/')[-1]
        self.type = filename[:5]
                                   
        with h5py.File(viirs_path,'r') as file:

            viirs_file = file['All_Data'][FILE_TYPE[self.type]]

            LOGGER.debug('Reading longitudes')
            self.longs = np.array(viirs_file['Longitude'][:], dtype=np.float32)
            # /-------------------------/
            LOGGER.debug('Reading latitudes')
            self.lats = np.array(viirs_file['Latitude'][:], dtype=np.float32)
            # /-------------------------/
            LOGGER.debug('Reading SatelliteAzimuthAngle')
            self.sat_azimuth_angle = np.array(
                                         viirs_file['SatelliteAzimuthAngle'][:],
                                         dtype=np.float32
                                         )
            # /-------------------------/
            LOGGER.debug('Reading SatelliteZenithAngle')
            self.sat_zenith_angle = np.array(
                                         viirs_file['SatelliteZenithAngle'][:],
                                         dtype=np.float32
                                         )
            # /-------------------------/
            LOGGER.debug('Reading SolarAzimuthAngle')
            self.solar_azimuth_angle = np.array(
                                           viirs_file['SolarAzimuthAngle'][:],
                                           dtype=np.float32
                                           )
            # /-------------------------/
            LOGGER.debug('Reading SolarZenithAngle')
            self.solar_zenith_angle = np.array(
                                          viirs_file['SolarZenithAngle'][:],
                                          dtype=np.float32
                                          )
            # /-------------------------/
            LOGGER.debug('Reading Scan Position')
            self.scan_position      = np.array(
                                          viirs_file['SCPosition'][:],
                                          dtype=np.float32
                                          )
            # /-------------------------/
            LOGGER.debug('Reading Heights')
            self.height      = np.array(
                                          viirs_file['Height'][:],
                                          dtype=np.float32
                                          )
            # /-------------------------/
            """
            LOGGER.debug('Reading Times')
            st_times      = np.array(
                                          viirs_file['StartTime'][:],
                                          dtype=np.int64
                                          )
            mid_times      = np.array(
                                          viirs_file['MidTime'][:],
                                          dtype=np.int64
                                          )
            self.times = np.ndarray( shape = ( 192 , 200) )
            for n_line,start in enumerate(st_times):
                dt = mid_times[n_line] / 100       
                self.times[n_line] = [ start + n_obs*dt for n_obs in range(200)]
            
            # /-------------------------/
            self.los = compute_LOS_ECEF(         viirs_file['SatelliteRange'],
                                                 viirs_file['SatelliteZenithAngle'],
                                                 viirs_file['SatelliteAzimuthAngle'],
                                                 viirs_file['Latitude'],
                                                 viirs_file['Longitude'] )
            """
            # /-------------------------/
            if self.type == 'GMODO':
                self.read_gmodo_cloudmask(self.path)  # read cloudmask from GMODO file
            elif self.type == 'GITCO' or self.type == 'GMTCO':
                self.read_cm_cloudmask(self.cm_path)  # read cloudmask from IICMO file
            else:
                self.cloudmask = []                   # default empty mask
            
    def read_gmodo_cloudmask(self,path):
        from utilities.array_reshapers import array_1d, transform_index
        """
        Read cloudmask from GMODO file. Read Ref for bytes interpretation. 
        Ref: https://www.star.nesdis.noaa.gov/jpss/documents/ATBD/D0001-M01-S01-011_JPSS_ATBD_VIIRS-Cloud-Mask_E.pdf
        """
        self.cloudmask = np.zeros(shape = self.longs.shape ,  dtype = np.int8)
        self.CMQ       = np.zeros( shape = self.longs.shape , dtype = np.int8)
        self.daynight  = np.zeros( shape = self.longs.shape , dtype = np.int8)
        
        with h5py.File(path,'r') as file:
                      LOGGER.debug('Reading Cloud Mask')
                      tmp_mask = file['All_Data']['VIIRS-CM-EDR_All']['QF1_VIIRSCMEDR'][:]
                      
                      mask_shape = tmp_mask.shape
                      
                      # 1d mask
                      tmp_mask = array_1d(tmp_mask)

                      for index, value in enumerate(tmp_mask):

                          # convert from uint8 to bin
                          bin_string = np.binary_repr(value, width = 8)
                          
                          matrix_indices = transform_index(index,mask_shape) 
                          quality = bin_string[-2:]

                          # Assign to cloudmask's element an integer from 0 to 3 (four possible combinations
                          #  with second and third bit), with the following meaning:    0 - CONF. CLEAR
                          #                                                             1 - PROB. CLEAR
                          #                                                             2 - PROB. CLOUDY
                          #                                                             3 - CONF. CLOUDY
                          self.CMQ[ matrix_indices ] = 0 if quality == '00' else 1
                          self.cloudmask[ matrix_indices ] = int(bin_string[-3]) + 2*int(bin_string[-4]) \
                                                             if self.CMQ[ matrix_indices ] == 1          \
                                                             else -128
                                                             
        self.cloudmask = np.ma.masked_where(self.cloudmask == -128,self.cloudmask,copy=False).astype(int)
        if type(self.cloudmask.mask) == np.bool_:
                  self.cloudmask.mask = np.zeros_like(self.cloudmask,dtype=np.bool)

    def read_cm_cloudmask(self,path):
        """
        Read cloudmask from IICMO or JRR (netcdf) file. Read Ref for bytes interpretation. 
        Ref: https://www.star.nesdis.noaa.gov/jpss/documents/ATBD/D0001-M01-S01-011_JPSS_ATBD_VIIRS-Cloud-Mask_E.pdf
        """

        from utilities.array_reshapers import array_1d, transform_index

        self.cloudmask = np.zeros(shape = self.longs.shape , dtype  = np.int8)
        self.CMQ       = np.zeros( shape = self.longs.shape , dtype = np.int8)
        self.daynight  = np.zeros( shape = self.longs.shape , dtype = np.int8)

        LOGGER.debug('Reading Cloud Mask')
    
        if 'IICMO' in path:
                # Read Cloud Mask from IICMO file
                        
                with h5py.File(path,'r') as file:
     
                      tmp_mask = file['All_Data']['VIIRS-CM-IP_All']['QF1_VIIRSCMIP'][:]
    
                      mask_shape = tmp_mask.shape
    
                      # 1d mask
                      tmp_mask = array_1d(tmp_mask)
                      
                      for index, value in enumerate(tmp_mask):
    
                          # convert from uint8 to bin
                          bin_string = np.binary_repr(value, width = 8)
    
                          matrix_indices = transform_index(index,mask_shape)
                          quality = bin_string[-2:]
        
                          # Assign to cloudmask's element an integer from 0 to 3 (four possible combinations
                          #  with second and third bit), with the following meaning:    0 - CONF. CLEAR
                          #                                                             1 - PROB. CLEAR
                          #                                                             2 - PROB. CLOUDY
                          #                                                             3 - CONF. CLOUDY
                          self.CMQ[ matrix_indices ] = 0 if quality == '00' else 1
    
                          self.cloudmask[ matrix_indices ] = int(bin_string[-3]) + 2*int(bin_string[-4])  \
                                                             if self.CMQ[ matrix_indices ] == 1 \
                                                             else -128  # Invalid Value
                # Mask invalid values
                self.cloudmask = np.ma.masked_where(self.cloudmask == -128,self.cloudmask)
                if type(self.cloudmask.mask) == np.bool_:
                              self.cloudmask.mask = np.zeros_like(self.cloudmask,dtype=np.bool)


        elif 'JRR' in path:
                # Read Cloud Mask from JRR-CloudMask file
                from netCDF4 import Dataset

                self.cloudmask = np.zeros(shape = self.longs.shape , dtype=np.int8)
                self.CMQ       = np.zeros( shape = self.longs.shape , dtype=bool)

                with Dataset(path,'r') as file:
 
                    self.cloudmask = np.copy(file['CloudMask'][:])
                    self.CMQ       = np.copy(file['CloudMaskQualFlag'][:])

                self.cloudmask = np.ma.masked_where(self.cloudmask==-128,self.cloudmask)
                self.CMQ       = np.ma.masked_where(self.CMQ != 0 , self.CMQ)
                if type(self.cloudmask.mask) == np.bool_:
                         self.cloudmask.mask = np.zeros_like(self.cloudmask,dtype=np.bool)
                print('--->',self.cloudmask.mask.sum())
        else:
            raise ValueError("Unknown Cloud Mask file")
                            
            
        return
