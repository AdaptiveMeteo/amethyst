from netCDF4 import Dataset
import logging
from data_reader.standard_data_format import STANDARD_FMT

__author__    = "Paolo Antonelli and Paolo Scaccia"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__   = ["Paolo Antonelli","Paolo Scaccia"]
__license__   = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

LOGGER = logging.getLogger(__name__)

class netCDFReader(object):
    def __init__(self,filepath):
        """
        Parameters
        ----------
        filepath : string
            Path to netcdf file

        Returns
        -------
        None.

        """
        if type(filepath) is str:
            nc_fid = Dataset(filepath) 
        else:
            nc_fid = filepath                # read already built Dataset
        
        
        # Read groups        
        list_groups = nc_fid.groups.keys()
            
        # Read variables from each group
        for grp in list_groups:
           list_dim = nc_fid.groups[grp].dimensions.keys()
           try:
               for dim in list_dim:
                   LOGGER.debug('Reading {} group dimension {}'.format(grp,dim))
                   cmd = 'self.{}_{} = nc_fid.groups["{}"].dimensions["{}"].size'.format(grp,dim,grp,dim)
                   exec(cmd)
           except:
               pass
           
           list_var = nc_fid.groups[grp].variables.keys()
           for var in list_var:
               LOGGER.debug('Reading variable {} from group {}'.format(var,grp))
               cmd = grp + '_' + var + ' = nc_fid.groups["' + grp + '"].variables["' + var + '"][:]'
               cmd = 'self.{}_{} = nc_fid.groups["{}"].variables["{}"][:]'.format(grp,var,grp,var)
               exec(cmd)
           
        # Read var dimensions         
        list_dim = nc_fid.dimensions.keys()
        for dim in list_dim:
              LOGGER.debug('Reading dimension {}'.format(dim))
              cmd = 'self.{} = nc_fid.dimensions["{}"].size'.format(dim,dim)
              exec(cmd)

        list_var = nc_fid.variables.keys()
        for var in list_var:
              LOGGER.debug('Reading variable {}'.format(var))
              cmd = 'self.{} = nc_fid.variables["{}"][:]'.format(var,var)
              exec(cmd)
        
        nc_fid.close()
        return
    
    def compute_rh(self,profile_names):
        """
        Compute Relative Humidity from pressure, temperature and water vapor
        
        Parameters
        ----------
        profile_names : dict
                        Pressure, temperature and water vapor variable names

        Returns
        -------
            None. 
        """
        from atmos.mirto_atmos_tools import mr2rh
        
        FIX_TEMP = 273.15
        cmd = "self.rh = mr2rh(self.{},self.{},self.{},{})[0]".format(profile_names['pressure'],
                                                                      profile_names['temperature'],
                                                                      profile_names['water vapor'],
                                                                      FIX_TEMP)
        exec(cmd)
        
        return
    
    def to_standard(self,variables):
        """
        Convert profile names

        Parameters
        ----------
        variables : dict
                    Profile names to change with standard names

        Returns
        -------
        None.

        """
        for varname in variables.keys():
            cmd1 = "self.{} = self.{}".format(STANDARD_FMT[varname],variables[varname])  # Rename attribute
            cmd2 = "delattr(self,'{}')".format(variables[varname])
            try:
                exec(cmd1)
                exec(cmd2)
            except:
                pass   # Skip unknown variables

        return
