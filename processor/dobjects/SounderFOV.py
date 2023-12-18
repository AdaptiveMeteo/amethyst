"""
Class to retrieve data from a Measured Radiances File.
"""

from netCDF4 import Dataset
import numpy as np

def FOVCount(datafile):
    """
    This function is a wrapper around a netCDF data file
    It retrieves the number of observations in the file
    """
    with Dataset(datafile, mode='r') as df:
        obsnum = len(df.dimensions['number_of_FOVs'])

    return obsnum

class SounderFOV(object):

    def __init__(self, conf_vars):
        """
        Initialize the masured data file
        """
        self.df = Dataset(conf_vars["fov_file"], mode='r')
        self.single_fov = False

        self.obsnum = len(self.df.dimensions['number_of_FOVs'])
        self.channels = len(self.df.dimensions['number_of_channels'])

        self.fov_type = conf_vars["instrument"].upper()
        
        # Check if the right channel list file is given in configuration.xml
        if self.fov_type.lower() not in conf_vars["instr_chan_list"]:
            raise ValueError("Wrong channel list file for instrument {}.".format(self.fov_type))
            
        self.indx =  np.loadtxt(conf_vars["instr_chan_list"], dtype=np.int32) - 1  # Read Channel Indices
        self.selchannels = len(self.indx)
        self.wnR = self.df.variables['Wavenumbers'][:]
        self.wnR = self.wnR[self.indx]

    def __del__(self):
        self.df.close()

    def fov(self, obs):
        """
         Get one of the Observation spectra from file.
        """

        Latitude  = self.df.variables['Latitude'][obs]
        Longitude = self.df.variables['Longitude'][obs]
        FOVangle  = self.df.variables['FOV_angle'][obs]
        Solar_zenith_angle = self.df.variables['Solar_zenith_angle'][obs]
        Solar_azimuth_angle = self.df.variables['Solar_azimuth_angle'][obs]

        Rad = self.df.variables['Radiance'][obs, Ellipsis]
        Rad = Rad[self.indx]

        return [Latitude, Longitude, None, FOVangle, Solar_zenith_angle, Solar_azimuth_angle, np.copy(self.wnR), Rad]


class sounderfov(object):
    """
    This class is a wrapper around a SounderFOV object.
    It is used to read just an observation and keep it in memory
    until the object is deleted.
    """
    def __init__(self, fov, obs):
        """
        Initialize the object using a SounderFOV and the obs index
        """
        self.indx = np.copy(fov.indx)
        [self.Latitude,
         self.Longitude,
         self.TimeFracDay,
         self.FOVangle,
         self.Solar_zenith_angle,
         self.Solar_azimuth_angle,
         self.wnR,
         self.Rad] = fov.fov(obs)
        self.obs = obs

    def fov(self, obs):
        """
        Return the in-memory object
        """
        if obs != self.obs:
            raise IndexError('Object initilized with obs='+repr(self.obs)+
                             ' but obs='+repr(obs)+' requested !')
        return [self.Latitude,
                self.Longitude,
                self.TimeFracDay,
                self.FOVangle,
                self.Solar_zenith_angle,
                self.Solar_azimuth_angle,
                self.wnR,
                self.Rad]

if __name__ == '__main__':
    import amethyst_config

    if __package__ is None:
        raise ImportError('The file "SounderFOV.py" is embedded into '
                          'the dobjects package\nTo lauch it, use:\n'
                          '"python -m dobjects.SounderFOV" '
                          'from the main directory of this project.')

    measure = SounderFOV(amethyst_config)
    measure.test()
