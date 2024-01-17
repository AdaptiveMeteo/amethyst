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
from netCDF4 import Dataset

from constants import LANDEMISSIVITY

log = logging.getLogger(__name__)

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"


class LandClimatologyGroup(object):
    """
    A land climatology group contains all values that are
    associated with one or more IGBP classes.
    
    Args:
        - *netcdf_group*: A netcdf group: in the file that contains the
          emissivity there must be some netcdf groups. To create an object
          of the LandClimatologyGroup, plese use the abstraction of one
          of those groups provided by the netCDF4 library
    """
    def __init__(self, netcdf_group):
        self.functions = netcdf_group.variables[LANDEMISSIVITY.FUNCTIONS][:]
        self.bias = netcdf_group.variables[LANDEMISSIVITY.BIAS][:]
        self.covariance = netcdf_group.variables[LANDEMISSIVITY.COVARIANCE][:]
        
        self.n_of_eigenvalues = self.functions.shape[0]
    
    def data(self):
        return self.bias, self.functions, self.covariance


class LandClimatology(object):
    """
    A wrapper around the climatology emissivity file for the land
    """
    def __init__(self, f_path):
        # Save all the groups that are inside the netcdf file
        emiss_groups = {}
        with Dataset(f_path, 'r') as f:
            for grp in f.groups:
                if grp.startswith('group_'):
                    log.debug('Saving group {}'.format(grp))
                    grp_number = int(grp.split('_')[1])
                    grp_content = LandClimatologyGroup(f.groups[grp])
                    emiss_groups[grp_number] = grp_content
        
            # Check that in the file there is at least one group
            if len(emiss_groups.items()) == 0:
                raise ValueError('There are no valid groups inside the '
                                 '{} file'.format(f_path))
        
            log.debug('Reading the values of the vector that associates the '
                      'IGBP classes to the groups')
            self.mappingvector = f.variables[LANDEMISSIVITY.MAPPINGVECTOR][:]
            
            # The number of the allowed IGBP classes
            self.n_of_IGBP_classes = self.mappingvector.size
            log.debug('There are {} IGBP classes in the land file'
                      ''.format(self.n_of_IGBP_classes))
            
            # Create a group that associate to each IGBP classes its group
            self.groups = {}
            for i in range(self.n_of_IGBP_classes):
                self.groups[i+1] = emiss_groups[self.mappingvector[i]]
            
            # Save the wavenumbers
            self.wavenumbers = f.variables[LANDEMISSIVITY.WAVENUMBERS][:]

            # Check that each group has the same number of eigenvalues
            # Get one group to use as reference
            ref_group = list(self.groups.keys())[0]
            self.n_of_eigenvalues = self.groups[ref_group].n_of_eigenvalues
            # and compare all the others
            for grp_num, grp in self.groups.items():
                if self.n_of_eigenvalues != grp.n_of_eigenvalues:
                    raise ValueError('Different number of eigenvalues in '
                                     'among groups')

    def associate(self, i):
        return self.groups[i]
