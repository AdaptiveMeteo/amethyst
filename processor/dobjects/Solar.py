#!/usr/bin/env python

"""
This class export an object which contains data from a netCDF file.
File expected content is:

      float WVNM(nspectral) ;
        WVNM:units = "cm-1" ;
      float IRRADIANCE(nspectral) ;
        IRRADIANCE:units = "mW m-2 / cm-1" ;

The Spectral Irradiance given by this class are

  ### HERE IS MISSING A REFERENCE. ASTM ??? ###

The solar provides also a linear interpolator which
can interpolate irradiances on an ordered array of input wave numbers.
"""


from numpy import array, interp, linspace
from netCDF4 import Dataset

# Change the following to reflect the file naming
# If a variable named as the value of WVNM would
# not be found, the program will automatically
# look for a variable named as the value of
# WVNM_2TRY
WVNM = 'WVNM'
WVNM_2TRY = 'FREQ'
IRRADIANCE = 'IRRADIANCE'

class Solar(object):
    """Solar irradiance dataset"""

    def __init__(self, datafile):
        """Initialize object from a netCDF file"""
        solar_file = Dataset(datafile, mode='r')

        try:
            self.wvnm = array(solar_file.variables[WVNM][:])
        except KeyError:
            self.wvnm = array(solar_file.variables[WVNM_2TRY][:])
        self.irr = array(solar_file.variables[IRRADIANCE][:])
        solar_file.close()


    def get(self, wvnm=None):
        """Linear interpolation on selected frequencies"""
        if wvnm is None:
            return self.irr
        irr = interp(wvnm, self.wvnm, self.irr)
        return irr


