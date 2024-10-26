"""

"""

import logging
import numpy as np
import os


__author__    = "Paolo Scaccia and Paolo Antonelli"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__   = ["Paolo Scaccia", "Paolo Antonelli"]
__license__   = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

LOGGER = logging.getLogger(__name__)

def check_dir(dir_path,instr):
    """
    Check whether the directory contains a CRIS/VIIRS geolocalization file.
    """
    tag = {   'cris'    : 'GCRSO',
               'viirs I': 'GITCO',
               'viirs M': 'GMTCO'}

    if instr not in tag.keys():
        LOGGER.error('Wrong input instrument')
        raise IOError
        
    filelist = os.listdir(dir_path)
    for filename in filelist:
        if tag[instr] in filename:
            return True
        
    return False    
    
def preprocess_geoloc_data(geoloc_cris,geoloc_viirs):
    from geometry.utilities.array_reshapers import array_1d, transform_index
    
    MAX_LON_DIFF = 0.25  # max degree diff. 
    MAX_LAT_DIFF = 0.25  # max degree diff.     
            
    cris_fov_position = np.array([ vector for vector in zip( array_1d(geoloc_cris.longs),
                                                             array_1d(geoloc_cris.lats)  ) ] )
       
    viirs_lons = array_1d(geoloc_viirs.longs)
    viirs_lats = array_1d(geoloc_viirs.lats)

    cris_shape  = np.shape(geoloc_cris.lats)
    viirs_shape = np.shape(geoloc_viirs.lats)

    indices = np.ndarray( shape = cris_shape )
    indices = [ [ ] for i in range(geoloc_cris.longs.size)] 
    
    for i,cris_fov in enumerate(cris_fov_position):

        lon_min = cris_fov[0] - MAX_LON_DIFF
        lon_max = cris_fov[0] + MAX_LON_DIFF        
        lat_min = cris_fov[1] - MAX_LAT_DIFF
        lat_max = cris_fov[1] + MAX_LAT_DIFF
        
        lon_ind = np.where(  (viirs_lons>= lon_min)   & (viirs_lons <= lon_max)  )
        lat_ind = np.where(  (viirs_lats >= lat_min)  & (viirs_lats  <= lat_max)  )
        
        # Add indices satisfing VIIRS cloudmask's quality control
        indices[i] = [ transform_index(x,viirs_shape) for x in np.intersect1d(lon_ind,lat_ind)  
                       if geoloc_viirs.cloudmask.mask[transform_index(x,viirs_shape)] == False ]  
        
    return indices
    

def get_collocation(geoloc_cris, geoloc_viirs):
    from geometry.utilities.compute_los import select_LOS
    
    # IMPROVE THIS PART MERGING THE TWO FUNCTIONS (preprocc.. and select_LOS)
    LOGGER.debug("Preprocessing geolocation VIIRS data...")
    indices = preprocess_geoloc_data(geoloc_cris,geoloc_viirs)
    
    LOGGER.debug("Selecting VIIRS FOVs using LOS criterion... ")
    indices = select_LOS(geoloc_cris, geoloc_viirs, indices)
    
    return indices


def test(cris_file_path,viirs_file_path):
    from geometry.utilities.array_reshapers import transform_index
    from data_reader.cris_wrapper import CRIS_Wrapper
    from data_reader.viirs_wrapper import VIIRS_Wrapper

    # Read files
    geoloc_cris  = CRIS_Wrapper(cris_file_path)
    geoloc_viirs = VIIRS_Wrapper(viirs_file_path)
    
    # Sub-sample the data
    N_CRIS_LINES  = 30
    N_VIIRS_LINES = 2200

    geoloc_cris.longs  = geoloc_cris.longs[:N_CRIS_LINES,:]
    geoloc_cris.lats   = geoloc_cris.lats[:N_CRIS_LINES,:]
    geoloc_cris.height = geoloc_cris.height[:N_CRIS_LINES,:]
    
    geoloc_viirs.longs  = geoloc_viirs.longs[:N_VIIRS_LINES,:]
    geoloc_viirs.lats   = geoloc_viirs.lats[:N_VIIRS_LINES,:]
    geoloc_viirs.height = geoloc_viirs.height[:N_VIIRS_LINES,:]
    geoloc_viirs.cloudmask = geoloc_viirs.cloudmask[:N_VIIRS_LINES,:]
    geoloc_viirs.CMQ = geoloc_viirs.CMQ[:N_VIIRS_LINES,:]

    indices = get_collocation(geoloc_cris, geoloc_viirs)
    
    """
    RAPID TEST: COMMENT UPPER LINE AND UNCOMMENT THIS PART
    from data_reader.cris_wrapper import CRIS_Wrapper
    from data_reader.viirs_wrapper import VIIRS_Wrapper
    geoloc_cris = CRIS_Wrapper(cris_file_path)
    geoloc_viirs = VIIRS_Wrapper(viirs_file_path)    
    indices =  [ [] for i in range(geoloc_cris.longs.size)]
    indices[0] = [(1,1), (2,2),(3,3)]
    """
    
    # ------------   PRINTOUTS   ---------------
    
    # Output files
    cris_file  = open('cris_table.dat','w')
    lon_file   = open('viirs_lon_table.dat','w')
    lat_file   = open('viirs_lat_table.dat','w')

    # Comment first lines with file info
    cris_file.write("% file: {}\n".format(geoloc_cris.path.split('/')[-1]))
    cris_file.write("% lat lon zenith_angle azimuth_angle\n")
    lon_file.write( "% file: {}\n".format(geoloc_viirs.path.split('/')[-1]))
    lat_file.write( "% file: {}\n".format(geoloc_viirs.path.split('/')[-1]))
    
    cris_shape = np.shape(geoloc_cris.longs)
    
    for index,viirs_indices in enumerate(indices):
        
        # Convert index for 1dimensional cris_data to 2dimensional format
        cris_index = transform_index(index,cris_shape)
        
        # Print cris lat, lon, zenith_angle, azimuth_angle
        cris_file.write("{} {} {} {}\n".format( geoloc_cris.lats[cris_index],
                                                geoloc_cris.longs[cris_index],
                                                geoloc_cris.sat_zenith_angle[cris_index],
                                                geoloc_cris.sat_azimuth_angle[cris_index] 
                                              ) 
                        )
        for i in range(300):
            if i < len(viirs_indices):
                viirs_index = viirs_indices[i]
                lon_file.write(" {} ".format(geoloc_viirs.longs[viirs_index] ))
                lat_file.write(" {} ".format(geoloc_viirs.lats[viirs_index] ))  
            else:
                lon_file.write(" Nan ")
                lat_file.write(" Nan ")
        lon_file.write("\n")
        lat_file.write("\n")
        
    cris_file.close()
    lon_file.close()
    lat_file.close()
