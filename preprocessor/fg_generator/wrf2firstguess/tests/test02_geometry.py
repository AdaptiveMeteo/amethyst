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

import numpy as np

from utilities.geometry import sinsquare_half_angular_distance,\
                               sinus_half_angular_distance,\
                               dist_on_earth,\
                               min_distance_indx

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

def sinsquare_angular_distance_test01():
    # Test that among the same points, the angle is 0

    # Generate a random grid of 100 x 100 points
    lons = np.linspace(-180, 180, 100)
    lats = np.linspace(-90, 90, 100)

    for i in range(100):
        for j in range(100):
            sin_alpha_square = sinsquare_half_angular_distance(
                                                               lons[i],
                                                               lats[j],
                                                               lons[i],
                                                               lats[j]
                                                               )
            assert np.abs(sin_alpha_square) < 1e-10


def sinsquare_angular_distance_test02():
    # Test the angle changing only the latitude
    # Generate a grid of 100 x 100 points
    lons = np.linspace(-180, 180, 100)
    lats = np.linspace(-90, 90, 100)

    for lat_change in (1, 5, 10, 30, 60, 90):
        expected_value = np.square(np.sin((lat_change *np.pi / 180.) / 2))
        for i in range(100):
            for j in range(100):
                # Add lat_change to the old latitude
                # The new lat must be between -90 and 90: if we go
                # above the max, subtract instead of add
                new_lat = lats[j] + lat_change

                if new_lat > 90:
                    new_lat -= 2 * lat_change
                sin_alpha_square = sinsquare_half_angular_distance(
                                                                   lons[i],
                                                                   lats[j],
                                                                   lons[i],
                                                                   new_lat
                                                                   )
                assert sin_alpha_square - expected_value < 1e10


def sinsquare_angular_distance_test03():
    # Test the angle changing only the longitude
    # Generate a random grid of 100 x 100 points
    lons = np.linspace(-180, 180, 100)
    lats = np.linspace(-90, 90, 100)

    for lon_change in (1, 5, 10, 30, 100, 150, 200):
        expected_on_equator = np.square(np.sin((lon_change *np.pi / 180.) / 2))
        for i in range(100):
            # Add lat_change to the old latitude
            # The new lon must be between -180 and 180: reduce using the mod
            new_lon = (lons[i] + lon_change + 180.) % 360 - 180.
            for j in range(100):
                scale = np.square(np.cos(lats[j] * np.pi / 180.))
                expected_value = expected_on_equator * scale
                sin_alpha_square = sinsquare_half_angular_distance(
                                                                   lons[i],
                                                                   lats[j],
                                                                   new_lon,
                                                                   lats[j]
                                                                   )
                assert sin_alpha_square - expected_value < 1e10



def sinus_angular_distance_test01():
    # Test that among the same points, the angle is 0

    # Generate a random grid of 100 x 100 points
    lons = np.linspace(-180, 180, 100)
    lats = np.linspace(-90, 90, 100)

    for i in range(100):
        for j in range(100):
            sin_alpha = sinsquare_half_angular_distance(
                                                        lons[i],
                                                        lats[j],
                                                        lons[i],
                                                        lats[j]
                                                        )
            assert np.abs(sin_alpha) < 1e-10


def sinus_angular_distance_test02():
    # Test the angle changing only the latitude
    # Generate a random grid of 100 x 100 points
    lons = np.linspace(-180, 180, 100)
    lats = np.linspace(-90, 90, 100)

    for lat_change in (1, 5, 10, 30, 100):
        expected_value = np.sin((lat_change *np.pi / 180.) / 2)
        for i in range(100):
            for j in range(100):
                # Add lat_change to the old latitude
                # The new lat must be between -90 and 90: if we go
                # subtract instead of add
                new_lat = lats[j] + lat_change
                if new_lat > 90:
                    new_lat -= 2 * lat_change
                sin_alpha = sinus_half_angular_distance(
                                                        lons[i],
                                                        lats[j],
                                                        lons[i],
                                                        new_lat
                                                        )
                assert sin_alpha - expected_value < 1e10


def sinus_angular_distance_test03():
    # Test the angle changing only the longitude
    # Generate a random grid of 100 x 100 points
    lons = np.linspace(-180, 180, 100)
    lats = np.linspace(-90, 90, 100)

    for lon_change in (1, 5, 10, 30, 100, 150, 200):
        expected_on_equator = np.sin((lon_change *np.pi / 180.) / 2)
        for i in range(100):
            # Add lat_change to the old latitude
            # The new lon must be between -180 and 180: reduce using the mod
            new_lon = (lons[i] + lon_change + 180.) % 360 - 180.
            for j in range(100):
                scale = np.cos(lats[j] * np.pi / 180.)
                expected_value = expected_on_equator * scale
                sin_alpha_square = sinsquare_half_angular_distance(
                                                                      lons[i],
                                                                      lats[j],
                                                                      new_lon,
                                                                      lats[j]
                                                                      )
                assert sin_alpha_square - expected_value < 1e10



def dist_on_earth_test01():
    # Test that among the same points, the distance is 0

    # Generate a random grid of 100 x 100 points
    lons = np.linspace(-180, 180, 100)
    lats = np.linspace(-90, 90, 100)

    for i in range(100):
        for j in range(100):
            d = dist_on_earth(
                              lons[i],
                              lats[j],
                              lons[i],
                              lats[j]
                              )
            assert d < 1e-10


def dist_on_earth_test02():
    # Test the angle changing only the latitude
    # Generate a random grid of 100 x 100 points
    lons = np.linspace(-180, 180, 100)
    lats = np.linspace(-90, 90, 100)

    for lat_change in (1, 5, 10, 30, 100):
        expected_value = 6371 * np.pi / 180. * lat_change
        for i in range(100):
            for j in range(100):
                # Add lat_change to the old latitude
                # The new lat must be between -90 and 90: if we go
                # subtract instead of add
                new_lat = lats[j] + lat_change
                if new_lat > 90:
                    new_lat -= 2 * lat_change
                d = dist_on_earth(
                                  lons[i],
                                  lats[j],
                                  lons[i],
                                  new_lat
                                  )
                assert d - expected_value < 1e10


def dist_on_earth_test03():
    # Test the angle changing only the longitude
    # Generate a random grid of 100 x 100 points
    lons = np.linspace(-180, 180, 100)
    lats = np.linspace(-90, 90, 100)

    for lon_change in (1, 5, 10, 30, 100, 150, 200):
        expected_on_equator = 6371 * np.pi / 180. * lon_change
        for i in range(100):
            # Add lat_change to the old latitude
            # The new lon must be between -180 and 180: reduce using the mod
            new_lon = (lons[i] + lon_change + 180.) % 360 - 180.
            for j in range(100):
                scale = np.cos(lats[j] * np.pi / 180.)
                expected_value = expected_on_equator * scale
                d = dist_on_earth(
                                  lons[i],
                                  lats[j],
                                  new_lon,
                                  lats[j]
                                  )
                assert d - expected_value < 1e10


def dist_on_earth_test04():
    # Check if dist_on_earth follows the broadcast rule

    # Generate a random grid of 100 x 100 points
    lons = np.linspace(-180, 180, 100)
    lats = np.linspace(-90, 90, 100)
    d = dist_on_earth(0, 0, lons, lats)

    for i in range(lons.size):
        assert d[i] == dist_on_earth(0, 0, lons[i], lats[i])

    lons_mat, lats_mat = np.meshgrid(lons, lats, indexing='ij')
    
    d_mat = dist_on_earth(0, 0, lons_mat, lats_mat)
    
    for i in range(lons.size):
        for j in range(lats.size):
            assert d_mat[i,j] == dist_on_earth(0, 0, lons[i], lats[j])

    d_mat = dist_on_earth(lons_mat, lats_mat, 0, 90)
    
    for i in range(lons.size):
        for j in range(lats.size):
            assert d_mat[i,j] == dist_on_earth(lons[i], lats[j], 0, 90)

    assert np.allclose(dist_on_earth(lons_mat,lats_mat,lons_mat,lats_mat), 0)


def dist_on_earth_test05():
    # Check dist_on_earth using some real cities
    Rome = (12.50, 41.90)
    Venice = (12.36, 45.44)
    Washington = (-77.02, 38.90)
    Paris = (2.36, 48.86)
    Trieste = (13.8, 45.63)

    def dist_among_cities(c1, c2):
        return dist_on_earth(c1[0], c1[1], c2[0], c2[1])

    RomeVenice = dist_among_cities(Rome, Venice)
    VeniceTrieste = dist_among_cities(Venice, Trieste)
    ParisWashington = dist_among_cities(Paris, Washington)
    ParisVenice = dist_among_cities(Paris, Venice)
    RomeWashington = dist_among_cities(Rome, Washington)

    def rel_error(val, expected):
        return np.abs(val / expected -1)

    assert rel_error(RomeVenice, 393.6) < 1e-02
    assert rel_error(VeniceTrieste, 114.7) < 1e-02
    assert rel_error(ParisWashington, 6165.3) < 1e-02
    assert rel_error(ParisVenice, 844.85) < 1e-02
    assert rel_error(RomeWashington, 7216.5) < 1e-02



def min_distance_indx_test01():
    # Check the minimum distance between some cities
    Rome = (12.50, 41.90)
    Venice = (12.36, 45.44)
    Washington = (-77.02, 38.90)
    Paris = (2.36, 48.86)
    Trieste = (13.8, 45.63)

    cities = [Rome, Venice, Washington, Paris]
    lons = np.array([c[0] for c in cities])
    lats = np.array([c[1] for c in cities])

    indx = min_distance_indx(Trieste[0], Trieste[1], lons, lats)

    # The closest city must be Venice
    assert indx == (1,)


def min_distance_indx_test02():
    # Check thath the function works even with multidimensional arrays
    Rome = (12.50, 41.90)
    Venice = (12.36, 45.44)
    Washington = (-77.02, 38.90)
    Paris = (2.36, 48.86)
    Trieste = (13.8, 45.63)

    cities = [Rome, Venice, Washington, Paris]
    lons = np.array([c[0] for c in cities])
    lats = np.array([c[1] for c in cities])
    
    lons = lons.reshape(2,2)
    lats = lats.reshape(2,2)

    indx = min_distance_indx(Trieste[0], Trieste[1], lons, lats)

    # The closest city must be Venice
    assert indx == (0,1)



def min_distance_indx_test02():
    # Check that there are no problems with more than one minimum
    Rome = (12.50, 41.90)
    Venice = (12.36, 45.44)
    Washington = (-77.02, 38.90)
    Paris = (2.36, 48.86)
    Trieste = (13.8, 45.63)

    cities = [Rome, Venice, Washington, Paris]
    cities *= 2

    lons = np.array([c[0] for c in cities])
    lats = np.array([c[1] for c in cities])

    indx = min_distance_indx(Trieste[0], Trieste[1], lons, lats)

    # The closest city must be Venice (one of the two)
    assert indx == (1,) or indx == (6,)
