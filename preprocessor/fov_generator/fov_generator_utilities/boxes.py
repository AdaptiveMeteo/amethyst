# This file is part of Iasi2observations.
#
# Iasi2observations is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Iasi2observations is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Iasi2observations. If not, see <http://www.gnu.org/licenses/>.

"""
boxes is the module that contains all the classes related to
select only the FOVs that are in a particular region
"""

import logging

import numpy as np

LOGGER = logging.getLogger(__name__)

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

class Box(object):
    """
    An object that returns True if a point is inside a particular domain (and
    False otherwise). To ask if a point is inside the domain, the "in" operator
    can be used::

        >>> (longitude, latitude) in my_box
        True

    Otherwise, it is possible to use the function contains. The advantage of use
    the latter is that this function is vectorized and, therefore, if latitude
    and longitude are two arrays, the result will be an array
    """
    def __contains__(self, position):
        raise NotImplementedError

    def contains(self, lon, lat):
        """
        Return True if lon and lat are the coordinates of a point inside the
        box
        """
        raise NotImplementedError

class Rectangle(Box):
    """
    A rectangular box

    Args:
        - *lonmin*: the minimum value of the longitude
        - *lonmax*: the maximum value of the longitude
        - *latmin*: the minimum value of the latitude
        - *latmax*: the maximum value of the latitude
    """
    def __init__(self, lonmin, lonmax, latmin, latmax):
        self.lonmin, self.lonmax = lonmin, lonmax
        self.latmin, self.latmax = latmin, latmax

    def __contains__(self, position):
        lon, lat = position
        in_lon = np.logical_and(lon >= self.lonmin, lon <= self.lonmax)
        in_lat = np.logical_and(lat >= self.latmin, lat <= self.latmax)
        return np.logical_and(in_lon, in_lat)

    def contains(self, lon, lat):
        return self.__contains__((lon, lat))
