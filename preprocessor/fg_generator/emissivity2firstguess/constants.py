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

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"


class SimpleNamespace(object):
    def __init__(self):
        raise Exception('{} is just a namespace'
                        ''.format(self.__class__.__name__))


class OBSERVATIONS(SimpleNamespace):
    LONGITUDEFIELD = 'Longitude'
    LATITUDEFIELD = 'Latitude'
    FOVANGLEFIELD = 'FOV_angle'


class LANDEMISSIVITY(SimpleNamespace):
    FUNCTIONS = 'LandClimaModelModelFunctions'
    BIAS = 'LandClimaModelModelFunctionsBias'
    COVARIANCE = 'LandClimaModelCovariance'
    WAVENUMBERS = 'LandClimaModelWnModelFunctions'

    MAPPINGVECTOR = 'MappingIGBP2groups'


class SEAEMISSIVITY(SimpleNamespace):
    FUNCTIONS = 'ModelFunctions'
    BIAS = 'ModelFunctionsBias'
    COVARIANCE = 'ModelCovariance'
    WAVENUMBERS = 'ModelWaveNumbers'

    ZENITHANGLEVECTOR = 'ZenithAngles'


class FIRSTGUESS(SimpleNamespace):
    FOVNUM = 'number_of_FOVs'
    LATITUDE = 'Latitude'
    LONGITUDE = 'Longitude'
    SURFGROUP = 'surface_components'
    N_OF_WN = 'n_of_wavenumbers'
    N_OF_EIGENVALUES = 'n_of_eigenvalues'

    FUNCTIONS = 'ModelFunctions'
    BIAS = 'ModelFunctionsBias'
    COVARIANCE = 'ModelCovariance'
    WAVENUMBERS = 'ModelWaveNumbers'
    
    LANDWATER = 'LandOrWater'
