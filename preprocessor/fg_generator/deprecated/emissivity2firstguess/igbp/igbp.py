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

from constants import SimpleNamespace

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

log = logging.getLogger(__name__)


def isqrt(n):
    """
    Compute the biggest integer :math:`x` such that :math:`x^2 <= n'.
    
    Args:
        - :math:`n`: An not negative integer number  

    Return:
        The biggest integer :math:`x` such that :math:`x^2 <= n`
    """

    assert n >= 0, "{} must be a positive integer to compute isqrt".format(n)

    if n == 0:
        return 0
    if n == 1:
        return 1
    if n == 2:
        return 1
    
    return int(np.floor(np.sqrt(n)))



class IGBP_CLASSES(SimpleNamespace):
    EVERGREEN_NEEDLE_FORESTS = 1
    EVERGREEN_BROAD_FORESTS = 2
    DECIDUOUS_NEEDLE_FORESTS = 3
    DECIDUOUS_BROAD_FORESTS = 4
    MIXED_FORESTS = 5
    CLOSED_SHRUBSLAND = 6
    OPEN_SHRUBSLAND = 7
    WOODY_SAVANNAS = 8
    SAVANNAS = 9
    GRASSLAND = 10
    WETLAND = 11 
    CROPLAND = 12
    URBAN_AREA = 13
    CROP_MOSAIC = 14
    ANTARTIC_OR_PERMANENT_SNOW = 15
    BARREN_OR_DESERT_LAND = 16
    OCEAN_WATER = 17
    TUNDRA = 18
    FRESH_SNOW = 19
    SEA_ICE = 20


class Igbp(object):
    """
    A wrapper around an IGBP file

    Args:
        - *igbp_file_path*: The path of an IGBP file
    """

    def __init__(self, igbp_file_path):
        raw_data = np.fromfile(igbp_file_path, dtype='>i1')
        l = raw_data.size
        log.debug('Found {} entries in the IGBP file'.format(l))

        # Get the dimension of the IGBP grid
        # The IGBP files should have twice points in latutide than in longitude
        self.y_points = isqrt(l // 2)
        self.x_points = self.y_points * 2

        log.debug('Reshaping IGBP data as a ({}, {}) matrix'
                  ''.format(self.y_points, self.x_points))
        self.data = raw_data.reshape(self.y_points, self.x_points)
        
        # compute the angular size of the IGBP cells
        self.lon_res = 360.0 / self.x_points
        self.lat_res = 180.0 / self.y_points


    def __call__(self, lon, lat):
        """
        Return the value of the IGBP class in a particular point of
        the Earth
        
        Args:
            - *lon*: the longitude of the point
            - *lat*: the latitude of the point
        
        Return:
            The IGBP class of the particular point
        """

        # lon should be normalized betwenn 0 and 360
        lon = lon % 360.
        
        # For the latitude, we will go from 0 (to the north pole)
        # to 180 (south Pole)
        lat = 90. - lat        

        # Get the longitude indices on the self.data matrix
        lon_indx = np.int32(np.floor(lon / self.lon_res))
        
        # The same for the latitude
        lat_indx = np.int32(np.floor(lat / self.lat_res))
        
        # The last problem we have to deal with is the problem of the borders.
        # We divided the earth in cells and each cell can be represented as a
        # square on a plane map. The problem arises if a point is over the
        # border of a cell. The current implementation assign a point to the
        # lower cell if it is on a horizontal border and to the right cell if
        # it is on a vertical border. The problem is that, if we are on the
        # lower border of the lowest cell, we have to assign the point to the
        # lowest cell (that is, the lowest cells contains both the upper border
        # and the lowest border).
        # The following command reduce all the values that are equal to
        # self.num_y_points to self.num_y_points-1. It works both if lat_indx
        # is a scalar or if lat_indx is an array
        np.clip(lat_indx, 0, self.y_points -1)

        return self.data[lat_indx, lon_indx]
