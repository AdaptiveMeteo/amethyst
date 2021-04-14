# This file is part of Wrf2firstguess.
#
# Wrf2firstguess is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Wrf2firstguess is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Wrf2firstguess. If not, see <http://www.gnu.org/licenses/>.

import logging
from os import path

import numpy as np
from scipy.interpolate import interp1d
from netCDF4 import Dataset

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

# Get the main dir of the software (we expect that the current file
# is one level inside the tree directory)
SCRIPT_DIR =path.realpath(__file__).split('amethyst/')[0]+'amethyst'

# The path of a valid file for the ozone profiles
OZONE_PROFILE_FILES = path.join(SCRIPT_DIR, 'ancillary/atmosphere/fg_ozone_profiles.nc')

log = logging.getLogger(__name__)


class OzoneClimatologyMatrix(np.ndarray):
    """
    A Ozone ClimatologyMatrix is a numpy matrix whose columns
    are climatological profiles of the ozone on different belts
    of the Earth. In particular:
    - the first column is the ozone profile for the tropical belt
    - the second one is the ozone profile for the mid-latitude belt during the
      summer
    - the third one is the ozone profile for the subartic belt during the
      summer
    - the fourth one is the ozone profile for the mid-latitude belt during the
      winter
    - the fifth one is the ozone profile for the subartic belt during the
      winter

    Args:
        - *profile_path*: The path of a netcdf file with the data needed to
          build the matrix
        - Any other argument accepted by a numpy array (profile_path should be
          the first)
    """
    def __new__(cls, *args, **kwargs):
        profiles_path = args[0]
        args = args[1:]

        with Dataset(profiles_path, 'r') as f:
            tropical = f.variables['tropical'][:]
            midlatitude_summer = f.variables['midlatitude_summer'][:]
            midlatitude_winter = f.variables['midlatitude_winter'][:]
            subartic_summer = f.variables['subartic_summer'][:]
            subartic_winter = f.variables['subartic_winter'][:]
            levels = np.array(f.variables['levels'][:], dtype=np.float32)

        n_levels = levels.size
        assert tropical.size == n_levels
        assert midlatitude_summer.size == n_levels
        assert midlatitude_winter.size == n_levels
        assert subartic_summer.size == n_levels
        assert subartic_winter.size == n_levels

        self = np.ndarray.__new__(cls, (n_levels, 5), *args, **kwargs)

        self[:, 0] = tropical[:]
        self[:, 1] = midlatitude_winter[:]
        self[:, 2] = subartic_winter[:]
        self[:, 3] = midlatitude_summer[:]
        self[:, 4] = subartic_summer[:]

        self.levels = levels

        return self


class Belt(object):
    """
    A Belt is a portion of the Earth, symmetrical respect to the Equator,
    that starts at some latitude and ends at another one.

    The only purpose of this class is that, given a latitude, it is able to
    tell if a point with such a latitude is inside the belt or not.

    Args:
        - *min_lat*: The minimum latitude in the north emisphere where the
          belt starts
        - *max_lat*: The maximum latitude (in the north emisphere) where the
          belt ends
    """
    def __init__(self, min_lat, max_lat):
        if min_lat < 0:
            raise ValueError('min_lat must be positive. Received {}'
                             ''.format(min_lat))
        if max_lat < min_lat:
            raise ValueError('max_lat must be greater than min_lat. Received'
                             'min_lat={} max_lat={}'.format(min_lat, max_lat))
        self.min, self.max = min_lat, max_lat

    def __contains__(self, lat):
        if lat < 0:
            lat *= -1

        if lat >= self.min and lat < self.max:
            return True

        return False


try:
    CLIMAT_MATRIX = OzoneClimatologyMatrix(OZONE_PROFILE_FILES)
except:
    CLIMAT_MATRIX = None

TROPICAL_BELT = Belt(0, 30)
MIDLATITUDE_BELT = Belt(30, 50)
SUBARTIC_BELT = Belt(50, 91) # 91 to be safe if a point has lat 90

SUMMER = (4, 5, 6, 7, 8, 9)



def generate_ozone_profile(profile, levels, month, climat_matrix=CLIMAT_MATRIX,
                           climat_levels=None):
    """
    Generate the ozone component for a profile using only climatological
    considerations.

    Args:
        - *profile*:
        - *levels*: the pressures (in hPa) that define the levels of the
                    ozone profile that will be returned as output
        - *month*: the month when the profile has been chosen
        - *climat_matrix*: an  OzoneClimatologyMatrix with the coefficients
          that define the climatology for the ozone
        - *climat_levels*: if climat_matrix is not an OzoneClimatologyMatrix
          but a simple numpy matrix with the same shape, it is necessary to
          supply an array with the pressure levels the coefficient of the
          numpy matrix refers to.

    Returns:
        A 1D array with the values of the ozone for each level
    """

    if climat_matrix is None:
        raise ValueError('Climatology matrix not specified (and I can not '
                         'find or read {}'.format(OZONE_PROFILE_FILES))

    if climat_levels is None:
        climat_levels = climat_matrix.levels

    min_level = np.min(levels)
    max_level = np.max(levels)
    min_climat_level = np.min(climat_levels)
    max_climat_level = np.max(climat_levels)

    if min_level < min_climat_level:
        raise ValueError('The highest level of the model ({:.2f} hPa) is '
                         'above the heighest level generated by the ozone '
                         'climatology'.format(min_level, min_climat_level))

    if max_level > max_climat_level:
        raise ValueError('The lowest level of the model ({:.2f} hPa) is '
                         'under the lowest level generated by the ozone '
                         'climatology'.format(max_level, max_climat_level))

    lat = profile.lat

    if lat in TROPICAL_BELT:
        log.debug('The point is in the tropical belt')
        if month in SUMMER:
            log.debug('Using the summer profile')
            coeff = np.array([.75, .05, 0, .15, 0.05])
        else:
            log.debug('Using the winter profile')
            coeff = np.array([.75, .15, .05, .05, 0])

    elif lat in MIDLATITUDE_BELT:
        log.debug('The point is in the mid-latitude belt')
        if month in SUMMER:
            log.debug('Using the summer profile')
            coeff = np.array([.05, .05, .15, .60, .15])
        else:
            log.debug('Using the winter profile')
            coeff = np.array([.15, .60, .15, .05, .05])

    elif lat in SUBARTIC_BELT:
        log.debug('The point is in the subartic belt')
        if month in SUMMER:
            log.debug('Using the summer profile')
            coeff = np.array([0, 0, .05, .25, .70])
        else:
            log.debug('Using the winter profile')
            coeff = np.array([0, .25, .70, 0, .05])
    else:
        log.error('No zone suitable for latitude {:.2f}'.format(lat))
        raise ValueError('No zone suitable for latitude {:.2f}'.format(lat))

    # Check that the sum of the coefficients was 1
    assert np.allclose(np.sum(coeff), 1.)

    # Add some random noise to the coefficients (no needed anymore! Disabled!)
    # noise = 1 + np.random.rand(5) * 0.1
    # log.debug('Adding the following noise to the cofficients: '
    #           '({:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f})'.format(*noise))
    # coeff *= noise

    # Renormalize the coefficients
    # coeff /= sum(coeff)

    # log.debug('Using the following weights: '
    #           '({:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f})'.format(*coeff))

    climat_profile = np.dot(climat_matrix, coeff)

    # Create an interpolator in the logarithmic space
    log_ozone_interp = interp1d(
                                np.log(climat_levels),
                                climat_profile,
                                kind='linear',
                                copy=False,
                                bounds_error=True,
                               )

    # Create an interpolator in the usual space (starting from the log one)
    def ozone_interp(x):
        return log_ozone_interp(np.log(x))

    ozone_profile = ozone_interp(np.array(levels))

    # Now we will fix the measure unit multiply by some constants
    mO3 = 48    # molecular mass of o3
    md = 28.966 # molecular mass of air

    return np.log(ozone_profile * mO3 *1e-6 / md)
