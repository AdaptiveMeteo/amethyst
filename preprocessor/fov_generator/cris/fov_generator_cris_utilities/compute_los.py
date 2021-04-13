"""
This module belongs to the 'collocation' tool of AdaptiveTools
"""
import numpy as np
import logging

__author__     = "Paolo Scaccia and Paolo Antonelli"
__copyright__  = "Copyright 2016, AdaptiveMeteo S.r.l."
__credits__    = ["Paolo Scaccia", "Paolo Antonelli"]
__license__    = "--"
__version__    = "1.0"
__maintainer__ = "Paolo Scaccia"
__email__      = "paolo.scaccia@adaptivemeteo.com"
__status__     = "Development"

LOGGER = logging.getLogger(__name__)


def compute_LOS_ECEF(R,Zenith,Azimuth,Lat,Lon):
        """ 
        Function to compute the line-of-sight (LOS) in Earth-Centered-Earth-Fixed 
        coordinates, using Local Spherical Coordinates and Geodetic Latitude-Longitude-Altitude Coordinates.
        """
        
        datashape  = np.shape(Zenith)

        out = np.zeros( shape = datashape + (3,) )
        
        # import as np.array
        Lat = np.array(Lat)
        Lon = np.array(Lon)
        Zenith = np.array(Zenith)
        Azimuth = np.array(Azimuth)
        
        cos_lon     = np.cos(Lon*np.pi/180)
        cos_azimuth = np.cos(Azimuth*np.pi/180)
        cos_lat     = np.cos(Lat*np.pi/180)
        cos_zenith  = np.cos(Zenith*np.pi/180)

        sin_lat     = np.sin(Lat*np.pi/180)
        sin_lon     = np.sin(Lon*np.pi/180)
        sin_zenith  = np.sin(Zenith*np.pi/180)
        sin_azimuth = np.sin(Azimuth*np.pi/180)
        
        E = R*sin_zenith*sin_azimuth
        N = R*sin_zenith*cos_azimuth
        U = R*cos_zenith
        
        X = -E*(sin_lon) - N*(cos_lon*sin_lat) + U*(cos_lon*cos_lat)
        Y =  E*cos_lon   - N*sin_lon*sin_lat   + U*sin_lon*cos_lat
        Z =                N*cos_lat           + U*sin_lat 
        
        out = np.array([ vector for vector in zip(  -X.reshape( Zenith.size ),
                                                    -Y.reshape( Zenith.size ),
                                                    -Z.reshape( Zenith.size )   
                                                  ) 
                      ])
        return out.reshape( datashape + (3,) )   # return it with dataset shape
    
def confront_LOS(cris_los,viirs_los):
    """
    LOS confrontation. Method described in Ref.
    """
    
    COS_THRESHOLD = 0.9999646886138214
    
    cris_los  = np.array(cris_los)
    viirs_los = np.array(viirs_los)
    
    norm = np.sqrt( cris_los.dot(cris_los) * viirs_los.dot(viirs_los) )
    
    cosin = cris_los.dot(viirs_los) / norm
    
    if cosin > COS_THRESHOLD:
        return True
    else:
        return False
    
    
def select_LOS(cris_dataset, viirs_dataset, indices):
    """
    This function select from cris FOV's neigneighborhood, VIIRS FOV satisfing the 
    LOS confront. For a detailed description of the method used see Ref.
    Ref DOI:  doi:10.3390/rs8010076

    Dev NOTE:  adapt function to different methods
    """
    from preprocessor.fov_generator.cris.fov_generator_cris_utilities.array_reshapers import transform_index
    from preprocessor.fov_generator.cris.fov_generator_cris_utilities.earth_geometry  import compute_FOV_position_ECEF

    cris_shape     = cris_dataset.longs.shape
    
    # Compute obs. position vectors
    G_viirs   = compute_FOV_position_ECEF(viirs_dataset)
    G_cris    = compute_FOV_position_ECEF(cris_dataset)    
    
    
    for cris_index,index_row in enumerate(indices):
        # Read satellite position (same scan position every 3 scan lines) and compute CRIS LOS
        P_cris     =  cris_dataset.scan_position[ transform_index(cris_index,cris_shape)[0] // 3 ]  
        cris_los   =  G_cris[transform_index(cris_index,cris_shape)] - P_cris

        trash_indices = [   index for index in index_row  
                            if confront_LOS(cris_los, G_viirs[index] - P_cris ) is False 
                        ]
        
        # Remove trash indices from indices line 
        indices[cris_index] = list( set(indices[cris_index]) - set(trash_indices) )
        #if len(indices[cris_index]) != 0:
        #        print('size indices[{}]: {}'.format(cris_index,len(indices[cris_index])))
                    
    return indices
