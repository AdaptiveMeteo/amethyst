# This file is part of Emiss2FirstGuess.
#
# Emiss2FirstGuess is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Emiss2FirstGuess is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Emiss2FirstGuess. If not, see <http://www.gnu.org/licenses/>.

import logging
import numpy as np
from netCDF4 import Dataset

from constants import SEAEMISSIVITY

log = logging.getLogger(__name__)

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"


class MasudaClimatology(object):
    """
    A wrapper around the file that contains the masuda climatology for the sea
    
    Args:
        - *f_path*: the path to the file
    """
    def __init__(self, f_path):
        with Dataset(f_path, 'r') as f:
            self.bias = f.variables[SEAEMISSIVITY.BIAS][:]
            self.functions = f.variables[SEAEMISSIVITY.FUNCTIONS][:]
            self.covariances = f.variables[SEAEMISSIVITY.COVARIANCE][:]
            self.wavenumbers = f.variables[SEAEMISSIVITY.WAVENUMBERS][:]
            self.z_angles = f.variables[SEAEMISSIVITY.ZENITHANGLEVECTOR][:]
            
            self.n_of_eigenvalues = self.covariances.shape[1]

    def associate(self, zenith_angle):
        """
        Return functions, bias, and covariances of a point with a specific
        zenith angle
         
        Args:
            - *zenith_angle*: The zenith angle of the observation
        
        Returns:
            - *bias*
            - *functions*
            - *covariance*
        """
        log.debug('Looking for a known angle near {}'
                  ''.format(zenith_angle))
        dist = np.abs(self.z_angles - np.abs(zenith_angle))
        indx = np.argmin(dist)
        
        log.debug('Using angle {} ({} degrees far)'
                  ''.format(self.z_angles[indx], dist[indx]))
        
        bias = self.bias[indx]
        functions = self.functions[indx]
        covariance = self.covariances[indx]

        return bias, functions, covariance
