"""
viirs_wrapper is a module that contains the class VIIRS_Wrapper, which is a wrapper
to handle the GITCO/GMTCO files containing the data regarding the geolocationing of
the satellite
"""
import h5py
import numpy as np
import logging
import sys

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
            LOGGER.debug('GMODO file selected: all VIIRS CM information collected from single file')
        else:
            viirs_path = viirs_files_path[0]
            cm_path = viirs_files_path[1]
            self.cm_path  = cm_path
            cm_filename  = cm_path.split('/')[-1]
            LOGGER.debug('Information on VIIRS CM collected from {} file'.format(cm_filename))
            
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
                try:
                        self.read_cm_cloudmask(self.cm_path)  # read cloudmask from IICMO file
                except:
                        LOGGER.info('CloudMask not found in {}'.format(self.path))
                        self.cloudmask = []
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


        elif 'JRR' in path:
                # Read Cloud Mask from JRR-CloudMask file
                from netCDF4 import Dataset

                self.cloudmask = np.zeros(shape = self.longs.shape , dtype=np.int8)
                self.CMQ       = np.zeros( shape = self.longs.shape , dtype=bool)

                with Dataset(path,'r') as file:
 
                    # self.cloudmask = np.ma.masked_where(tmp_mask<0,tmp_mask)
                    # self.CMQ       = np.ma.masked_where(tmp_QC!=0,tmp_QC)

                    self.cloudmask = np.copy(file['CloudMask'][:])
                    self.CMQ       = np.copy(file['CloudMaskQualFlag'][:])
                    
                    self.cloudmask = np.ma.masked_where(self.cloudmask==-128,self.cloudmask)
                    self.CMQ       = np.ma.masked_where(self.CMQ!=0,self.CMQ)
                    
                    """
                    mask_shape = tmp_mask.shape
                    
                    # 1d mask
                    tmp_mask = array_1d(tmp_mask)
                    tmp_QC   = array_1d(tmp_QC)
 
                    for index, value in enumerate(tmp_mask):
    
                              matrix_indices = transform_index(index,mask_shape)
                              
                              self.cloudmask[ matrix_indices ] = value
                              self.CMQ[ matrix_indices ]       = np.int(tmp_QC[index])
                    """
        else:
            raise ValueError("Unknown Cloud Mask file")
                            
            
        return



def search_VIIRS_CM_files(cris_geo_path):
    """
    Function for merging 

    Parameters
    ----------
    cris_geo_path : str
        Path to GMTCO file

    Returns
    -------
        List of files with VIIRS CloudMask corresponding to the given GMTCO

    """
    import os 
    from datetime import datetime
    
    cris_geo_file = os.path.basename(cris_geo_path)
    start_date = datetime.strptime(cris_geo_file.split('_')[2].replace('d','') + cris_geo_file.split('_')[3].replace('t','')[:-3],'%Y%m%d%H%M')
    end_date   = datetime.strptime(cris_geo_file.split('_')[2].replace('d','') + cris_geo_file.split('_')[4].replace('e','')[:-3],'%Y%m%d%H%M')
    out_files = []
    viirs_cm_files =[ name for name in os.listdir(cris_geo_path.replace(cris_geo_file,'')) if 'CloudMask' in name ]
    for cm_file in viirs_cm_files:
        if 'AGGREGATED' not in cm_file:
            if datetime.strptime(cm_file.split('_')[3].replace('s','')[:-3],'%Y%m%d%H%M') >= start_date and datetime.strptime(cm_file.split('_')[4].replace('e','')[:-3],'%Y%m%d%H%M') <= end_date: 
                out_files.append(cris_geo_path.replace(cris_geo_file,'') + '/' + cm_file)
    
    return out_files

def read_single_JRR_cloudmask(cm_file):
    from netCDF4 import Dataset
    
    with Dataset(cm_file,'r') as file:
        binary_data = file['CloudMask'][:]
        quality_data = file['CloudMaskQualFlag'][:]
        lat          = file['Latitude'][:]
        lon          = file['Longitude'][:]
    return binary_data, quality_data, lat , lon
    
def merge_VIIRS_CM_files(viirs_cm_files):
    from netCDF4 import Dataset
    from datetime import datetime
    
    # Sort list
    viirs_cm_files.sort()

    # Parse file name and file path
    first_filename = viirs_cm_files[0].split('/')[-1]
    files_path     = viirs_cm_files[0].replace('/{}'.format(first_filename),'')
    # Read starting and ending date
    start_date = datetime.strptime(first_filename.split('_')[3][:-3],'s%Y%m%d%H%M') 
    end_date   = datetime.strptime(viirs_cm_files[-1].split('/')[-1].split('_')[4].replace('e','')[:-3],'%Y%m%d%H%M')
    
    # Read and aggregate cloud masks
    binary_data, quality_data, lat_data, lon_data = read_single_JRR_cloudmask( viirs_cm_files[0] )
    for file in viirs_cm_files[1:]:
        binary_buffer, quality_buffer, lat, lon = read_single_JRR_cloudmask(file)
        
        binary_data = np.concatenate(( binary_data, binary_buffer),   axis=0)
        quality_data = np.concatenate(( quality_data, quality_buffer),axis=0)
        lat_data = np.concatenate(( lat_data, lat),axis=0)
        lon_data = np.concatenate(( lon_data, lon),axis=0)
       
    
    out_file_name = '{}/{}_{}_{}_s{}_e{}_AGGREGATED.nc'.format(files_path,   
                                                               first_filename.split('_')[0],
                                                               first_filename.split('_')[1],
                                                               first_filename.split('_')[2],
                                                               start_date.strftime('%Y%m%d%H%M%S'), 
                                                               end_date.strftime('%Y%m%d%H%M%S')).replace('///','/')
    
    try:
        # Try to write the aggregated netcdf
        with Dataset(out_file_name,'w',format='NETCDF4') as ncfile:
            
            ncfile.createDimension('rows_dim',binary_data.shape[0])
            ncfile.createDimension('columns_dim',binary_data.shape[1])
            
            nc_binary_matrix = ncfile.createVariable('CloudMask',np.int8,('rows_dim','columns_dim'))
            nc_lons = ncfile.createVariable('Longitude',np.float,('rows_dim','columns_dim'))
            nc_lats = ncfile.createVariable('Latitude',np.float,('rows_dim','columns_dim'))

            nc_quality_matrix = ncfile.createVariable('CloudMaskQualFlag',np.int8,('rows_dim','columns_dim'))
            
            nc_binary_matrix[:] = binary_data
            nc_quality_matrix[:] = quality_data
            nc_lats[:]           = lat_data
            nc_lons[:]           = lon_data
            
        LOGGER.info('Created aggregated Cloud Mask file {}'.format(out_file_name))
    except Exception as error:
        sys.exit('Cannot create aggregated Cloud Mask file: {}'.format(error))

    return  out_file_name