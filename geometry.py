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

"""
The geometry module contains all the functions that are related with computing
properties of points on a sphere. In particular, if we define the angular
distance among two points on a sphere as the amplitude of the angle that they
define from the center of the sphere, in this module there are several
functions related to compute that quantity or the distance between the points.
"""

import numpy as np

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

def sinsquare_half_angular_distance(lon1, lat1, lon2, lat2):
    """
    Compute the square of the sinus of the half of the angolar distance between
    two (or more, using numpy broadcasting rule) points on a sphere. The main
    advantage compared with the "angular_distance" function is that this
    function should be faster to be computed (because it avoid to call a square
    root function).

    Args:
        - *lon1*: The longitude of the first point
        - *lat1*: The latitude of the first point
        - *lon2*: The longitude of the second point
        - *lon2*: The latitude of the second point

    Returns:
        The square of the sinus of the angular distance. That is a number such
        that the arcsinus of its square root is the angle among the two points.
    """

    # Convert in radiants
    lon1 = lon1 * np.pi / 180.
    lat1 = lat1 * np.pi / 180.
    lon2 = lon2 * np.pi / 180.
    lat2 = lat2 * np.pi / 180.

    lat_dist = np.sin((lat1 - lat2) / 2.)
    lat_dist = np.square(lat_dist)

    lon_dist = np.sin((lon1 - lon2) / 2.)
    lon_dist = np.square(lon_dist)

    return lat_dist + np.cos(lat1) * np.cos(lat2) * lon_dist



def sinus_half_angular_distance(lon1, lat1, lon2, lat2):
    """
    Compute the sinus of the half of the angular distance between two
    (or more, using numpy broadcasting rule) points on a sphere.

    Args:
        - *lon1*: The longitude of the first point
        - *lat1*: The latitude of the first point
        - *lon2*: The longitude of the second point
        - *lon2*: The latitude of the second point

    Returns:
        The sinus of the angular distance. That is a number such that
        its arcsin is the angle among the two points.
    """

    return np.sqrt(sinsquare_half_angular_distance(lon1, lat1, lon2, lat2))



def dist_on_earth(lon1, lat1, lon2, lat2):
    """
    Compute the distance between two (or more, using numpy
    broadcasting rule) points on Earth (approximated as a sphere of radious
    6371 Km).

    Args:
        - *lon1*: The longitude of the first point
        - *lat1*: The latitude of the first point
        - *lon2*: The longitude of the second point
        - *lon2*: The latitude of the second point

    Returns:
        The distance between the two points
    """

    radious = 6371.

    sin_alpha = sinus_half_angular_distance(lon1, lat1, lon2, lat2)

    # Remove numerical errors that could bring sin_alpha over 1
    if isinstance(sin_alpha, np.ndarray):
        sin_alpha[sin_alpha >= 1] = 1 - 1E-9
    elif sin_alpha >= 1:
        sin_alpha = 1 - 1E-9

    return 2 * radious * np.arcsin(sin_alpha)


def min_distance_indx(p_lon, p_lat, lons, lats):
    """
    Given a point P and an array of other points on a sphere, return the index
    of the point which is closest to P. If the other points are in a
    multi-dimensional array, the returned index will be a tuple with a entry
    for each axis.

    Args:
        - *p_lon*: The longitude of the point P
        - *p_lat*: The latitude of the point P
        - *lons*: The longitudes of all the other points
        - *lats*: The latitudes of all the other points

    Returns:
        The index of the closest point
    """

    dists = sinsquare_half_angular_distance(p_lon, p_lat, lons, lats)
    d_min = np.min(dists)

    # min_pos contains one 1D array for each dimension of lons (or lats)
    # The arrays have the same size. For an index i,
    # (min_output[0][i], min_output[1][i], ..., min_output[n][i])
    # is a position of a minimum
    min_pos = np.where(dists == d_min)

    # Get only the first one
    output = [v[0] for v in min_pos]

    return tuple(output)
