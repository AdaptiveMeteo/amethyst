import logging
from data_reader.read_netcdf import netCDFReader
import numpy as np
from data_reader.input_reader import MIRTOInputData

__author__    = "Paolo Antonelli and Paolo Scaccia"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__   = ["Paolo Antonelli","Paolo Scaccia"]
__license__   = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

LOGGER = logging.getLogger(__name__)


class MIRTO_results_Wrapper(netCDFReader):
    """
    A wrapper around the MIRTO "results.nc" files. A MIRTO_results_Wrapper reads all the information
    inside a specified MIRTO "results.nc" file and loads it into memory. After that, it
    expose all the data through its interface.

    """
    def __init__(self,file_path):
       netCDFReader.__init__(self,file_path)
       
    def to_standard(self):
        netCDFReader.to_standard(self,{'pressure':'p','temperature':'t','water vapor':'q'})
        
    def compute_rh(self):
        try:
            netCDFReader.compute_rh(self,{'pressure':'p','temperature':'t','water vapor':'q'})
        except:
            netCDFReader.compute_rh(self,{'pressure':'air_pressure','temperature':'air_temperature','water vapor':'air_water_vapor_mr'})
        netCDFReader.to_standard(self,{'relative humidity':'rh'})
        

class MIRTO_fov_Wrapper(netCDFReader):
    """
    A wrapper around the MIRTO "fov.nc" files. A MIRTO_results_Wrapper reads all the information
    inside a specified MIRTO "fov.nc" file and loads it into memory. After that, it
    expose all the data through its interface.

    """

    def __init__(self,file_path):
       netCDFReader.__init__(self,file_path)
       
    def to_standard(self):
        """
        convert to Standard format
        """
        netCDFReader.to_standard(self,{'latitude'              : 'Latitude',
                                       'longitude'             : 'Longitude',
                                       'temperature'           : 'T',
                                       'fov angle'             : 'FOV_angle',
                                       'sat azimuth angle'     : 'Satellite_azimuth_angle',
                                       'solar azimuth angle'   : 'Solar_azimuth_angle',
                                       'solar zenith angle'    : 'Solar_zenith_angle'
                                       }
                                 )
class MIRTO_fg_Wrapper(netCDFReader):
    """
    A wrapper around the MIRTO "fg.nc" files. A MIRTO_results_Wrapper reads all the information
    inside a specified MIRTO "fov.nc" file and loads it into memory. After that, it
    expose all the data through its interface.
    """

    def __init__(self,file_path):
       netCDFReader.__init__(self,file_path)

    def compute_rh(self):
        netCDFReader.compute_rh(self,{'pressure':'atmospheric_components_p','temperature':'atmospheric_components_T','water vapor':'atmospheric_components_q'})
        netCDFReader.to_standard(self,{'relative humidity':'rh'})
               
    def to_standard(self):
        """
        convert to Standard format
        """
        netCDFReader.to_standard(self,{'latitude' : 'Latitude',
                                       'longitude': 'Longitude',
                                       'pressure':'atmospheric_components_p',
                                       'temperature':'atmospheric_components_T',
                                       'water vapor':'atmospheric_components_q'
                                       }
                                       
                                 )

class MIRTO_joined_Wrapper(netCDFReader):
        """
        A wrapper joining the MIRTO "fov.nc" and "results.nc" files. 
        """
        def __init__(self,file_path):
           netCDFReader.__init__(self,file_path)
               
           # Add fov or results file
           if 'results' in file_path:
               self.add_fov(file_path)
           elif 'fov.nc' in file_path:
               self.add_results(file_path)   
               
           # Compute RH
           netCDFReader.compute_rh(self,{'pressure':'p','temperature':'t','water vapor':'q'})
           
           # Convert all data to standard format
           netCDFReader.to_standard(self,{
                                      'pressure':'p',
                                      'temperature':'t',
                                      'water vapor':'q',
                                      'latitude' : 'Latitude',
                                      'longitude': 'Longitude',
                                      'sat azimuth angle'     : 'Satellite_azimuth_angle',
                                      'solar azimuth angle'   : 'Solar_azimuth_angle',
                                      'solar zenith angle'    : 'Solar_zenith_angle',
                                      'fov angle'             : 'FOV_angle',
                                      'relative humidity'     : 'rh',
                                      'surface temperature'   : 'skt'
                                      }
                                    )
           
               

        def add_fov(self,file_path):
            file = file_path.replace('results','fov')
            netCDFReader.__init__(self,file)
            
            return
        
        def add_results(self,file):
            file = file.replace('fov','results')
            netCDFReader.__init__(self,file)
            return
        
        

            
class dimension(object):
    def __init__(self,size):
        self.__dict__ = {'size':size}
        
        
class WRFDA_fov_Wrapper(netCDFReader):
    """
    A wrapper around the WRDDA "diags_metop-2-iasi_r_yyymmddHHMM.nc" files.
    A WRFDA_fov_Wrapper reads all the information inside a specified WRFDA
    "diags_metop-2-iasi_r_yyymmddHHMM.nc" file and loads it into memory.
    After that, it expose all the data through its interface.

    """

    def __init__(self,file_path):
       netCDFReader.__init__(self,file_path)

    def to_standard(self):
        """
        convert to Standard format
        """
        netCDFReader.to_standard(self,{'latitude'          : 'lat',
                                       'longitude'         : 'lon',
                                       'sat azimuth angle' : 'satazi',
                                       'sat zenith angle'  : 'satzen',
                                       'quality control'   : 'tb_qc',
                                       'obs minus bac'     : 'tb_inv'}
                                 )

class MIRTO_multi_overpass_Wrapper(object):
    def __init__(self,filepath):

           starting_file = filepath if  type(filepath) is str else filepath[0]
           dataset = MIRTOInputData(starting_file)

           self.latitude   = np.copy(dataset.latitude)
           self.longitude  = np.copy(dataset.longitude)
           self.d2         = np.copy(dataset.d2)
           self.profiles   = dataset.profiles
           self.nfovs      =    [  self.latitude.size ]
           self.name       = dataset.name   
           if type(filepath) is list:
               self.name = [self.name.split('/')[-1]]
               self.nfovs = [  dataset.latitude.size ]
               LOGGER.info('Combining MIRTO datasets...')
               for path in filepath[1:]:
                           try:
                               self.combine_dataset(path)
                           except:
                               LOGGER.info('Skipped file {}'.format(path))
           return

    def combine_dataset(self,input_filename):
        dataset  = MIRTOInputData(input_filename)
        self.nfovs.append(dataset.latitude.size)

        self.latitude = np.concatenate((self.latitude, dataset.latitude[:]), axis=0)
        self.longitude = np.concatenate((self.longitude, dataset.longitude[:]), axis=0)
        self.d2 = np.concatenate((self.d2, dataset.d2[:]), axis=0)

        for i in range(0,2):

            self.profiles[i].pressure     = np.concatenate((self.profiles[i].pressure,dataset.profiles[i].pressure[:]),axis=0)
            self.profiles[i].temperature  = np.concatenate((self.profiles[i].temperature,dataset.profiles[i].temperature[:]),axis=0)
            self.profiles[i].water_vapour = np.concatenate((self.profiles[i].water_vapour,dataset.profiles[i].water_vapour[:]),axis=0)
            if dataset.profiles[0].rh.shape[0] != 81:
                self.profiles[i].rh           = np.concatenate((self.profiles[i].rh,dataset.profiles[i].rh[:]),axis=0)
            else:
                self.profiles[i].rh           = np.concatenate((self.profiles[i].rh,dataset.profiles[i].rh.T),axis=0)
          
        self.name.append(dataset.name.split('/')[-1])
        dataset = None
        return
    
    def to_standard(self):
        """
                TO BE ADDED...
        Returns
        -------
        None.

        """
        return
    
    def fov_original_index(self,index):
        counter = 0
        for idim,dim in enumerate(self.nfovs):
            counter += dim
            if index // counter == 0:
                break
    
        return int(index - np.sum(self.nfovs[:idim])), self.name[idim]
    
    def close(self):
       self.variables  = {}
       self.groups     = {}
       self.dimensions = {}
       return

