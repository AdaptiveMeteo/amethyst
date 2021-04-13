# This file is part of Cris2observations.
#
# Cris2observations is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Cris2observations is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Cris2observation. If not, see <http://www.gnu.org/licenses/>.

"""
scris_wrapper is the module that contains the class SCrISWrapper, which is a
wrapper to handle the SCrIS files containing the data measured by the CrIS
sounder
"""

import h5py
import logging

import numpy as np
from numpy.lib.stride_tricks import as_strided

from utilities.array_reshapers import array_twist

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

LOGGER = logging.getLogger(__name__)

class SCrISWrapper(object):
    """
    A wrapper around the SCrIS files. A SCrISWrapper reads all the information
    inside a specified GCRSO file and loads it into memory. After that, it
    exposes all the data through its interface.

    :ivar str path: The path of the original SCrIS file
    :ivar radiances: The radiances measured by the CrIS instrument
    :ivar wavenumbers: The wavelenght of the channels of the sounder
    """

    def __init__(self, scris_path):
        self.path = scris_path
        with h5py.File(self.path, 'r') as scris_file:
            geo_all = scris_file['All_Data/CrIS-SDR_All']

            LOGGER.debug('Reading ES_RealLW')
            reallw = np.array(geo_all['ES_RealLW'][:], dtype=np.float32)
            LOGGER.debug('Its shape is {}'.format(reallw.shape))
            LOGGER.debug('Twisting ES_RealLW')
            reallw = array_twist(reallw)
            LOGGER.debug('Now its shape is {}'.format(reallw.shape))
            LOGGER.debug('Removing the first and the last two channels')
            reallw = reallw[:, :, 2:-2]

            LOGGER.debug('Reading ES_RealMW')
            realmw = np.array(geo_all['ES_RealMW'][:], dtype=np.float32)
            LOGGER.debug('Its shape is {}'.format(realmw.shape))
            LOGGER.debug('Twisting ES_RealMW')
            realmw = array_twist(realmw)
            LOGGER.debug('Now its shape is {}'.format(realmw.shape))
            LOGGER.debug('Removing the first and the last two channels')
            realmw = realmw[:, :, 2:-2]

            LOGGER.debug('Reading ES_RealSW')
            realsw = np.array(geo_all['ES_RealSW'][:], dtype=np.float32)
            LOGGER.debug('Its shape is {}'.format(realsw.shape))
            LOGGER.debug('Twisting ES_RealSW')
            realsw = array_twist(realsw)
            LOGGER.debug('Now its shape is {}'.format(realsw.shape))
            LOGGER.debug('Removing the first and the last two channels')
            realsw = realsw[:, :, 2:-2]

            LOGGER.debug('Gluing them together to generate the radiances')
            self.radiances = np.concatenate((reallw, realmw, realsw), axis=-1)
            self.n_scanlines = self.radiances.shape[0]
            self.scanline_length = self.radiances.shape[1]

        LOGGER.debug('Generating wavenumbers')
        low_wavenumbers = np.linspace(650 - 0.625*2, 1095 + 0.625*2, 717)
        medium_wavenumbers = np.linspace(1210 - 1.25*2, 1750 + 1.25*2, 437)
        height_wavenumbers = np.linspace(2155 - 2.50*2, 2550 + 2.50*2, 163)
        # Remove the first and the last two channels
        low_wavenumbers = low_wavenumbers[2:-2]
        medium_wavenumbers = medium_wavenumbers[2:-2]
        height_wavenumbers = height_wavenumbers[2:-2]

        self.wavenumbers = np.concatenate((
                                           low_wavenumbers,
                                           medium_wavenumbers,
                                           height_wavenumbers
                                          ))

        LOGGER.debug('Generating FOV angles')
        # These are the angles for the detector 5, the central one. The
        # sensors are sorted as the following:
        #
        #   7  8  9
        #   4  5  6
        #   1  2  3
        #
        # For the detector 2 and 8, the angle is the same of the detector
        # 5. For the 4 and the 6, we have to split by two the angular distance
        # between two position of the detector 5.
        central_fovangles = [48.3409, 45.0068, 41.6727, 38.3386, 35.0045,
                             31.6708, 28.3371, 25.0035, 21.6698, 18.3361,
                             15.0031, 11.6701, 8.3371, 5.0040, 1.6709,
                             -1.6615, -4.9939, -8.3263, -11.6588, -14.9912,
                             -18.3235, -21.6557, -24.9880, -28.3203, -31.6525,
                             -34.9844, -38.3163, -41.6482, -44.9801, -48.3119]
        central_fovangles = np.array(central_fovangles)

        # Now we need to build the angles for all the fov
        fovangles = np.zeros((3*central_fovangles.size), dtype=np.float32)
        # starting from the second element, the right values every three step
        # is the central one
        fovangles[1::3] = central_fovangles
        # Compute the angular distance between two different points of the
        # central detector
        angle_diff = central_fovangles[1:] - central_fovangles[:-1]
        # Compute the angle for the detector 6: we need to add to the
        # position of the detector 5, one third of the angular distance
        # among the current position of the detector 5 and its next one
        # We will not fill the last position of the array because there
        # is no "next position" for the detector 5
        fovangles[2:-1:3] = central_fovangles[:-1] + angle_diff/3
        # Compute the angle for the detector 4: we need to remove to the
        # position of the detector 5, one third of the angular distance
        # among the current position of the detector 5 and its previous
        # one. We will not fill the first position of the array because
        # there is no "previous position" for the detector 5
        fovangles[3:-1:3] = central_fovangles[1:] - angle_diff/3

        # Now we take care of the first and of the last position
        fovangles[0] = 2 * fovangles[1] - fovangles[2]
        fovangles[-1] = 2 * fovangles[-2] - fovangles[-3]

        # The following is a trick to repeat the values of
        # fovangles for each scanline. It is equivalent to
        # np.repeat(fovangles, self.n_scanlines, axes=0)
        self.fovangles = as_strided(
                                    fovangles,
                                    strides=(0,) + fovangles.strides,
                                    shape=(self.n_scanlines,)+fovangles.shape,
                                    )


    def valid_radiances(self):
        """
        Return an array of booleans whose entries are False for the indices
        that represent radiances with negative values
        """
        # Create an array that reports for each entry if it is positive or not
        good_values = self.radiances > 0
        # Multiply the values along the last axis to known if inside the
        # radiance there is at least one negative element
        valid_radiances = np.prod(good_values, axis=-1)
        return np.bool_(valid_radiances)


    def generate_index_map(self):
        """
        An index map is a vector that "counts" the observations. Indeed, in
        the original SCRIS file every variable is saved in a 3d array where
        the first index is the scanline, the second is the position of the
        instrument inside the current scanline and the third is the number
        of the decoder.
        """
        # With "n_of_scanlines" here we refer to the scanlines of the satellite
        # that 3 times less than the vertical lines that we generated before
        n_scanlines = self.n_scanlines // 3
        # The same for the scanline_length, because we have 3 detectors in each
        # line
        scanline_length = self.scanline_length // 3

        numobs = n_scanlines * scanline_length * 9

        # This vector simply counts the observations. Be aware that if you compute
        # i % 9 for each entries of this vector, you have the last index of its
        # position on the array. If you instead do (i // 9) % scanline_length you
        # have the second index and so on
        original_position = np.arange(0, numobs).reshape(n_scanlines,
                                                         scanline_length,
                                                         9
                                                         )

        # let see where the array_twist moves the values
        grid2d_position = array_twist(original_position)

        # Now we prepare an array where we will store the old indices
        index_map = np.empty(grid2d_position.shape[:2] + (3,), dtype=np.int32)

        # And we use the trick explained before to get the old indices
        index_map[:, :, 2] = grid2d_position % 9
        index_map[:, :, 1] = (grid2d_position // 9) % scanline_length
        index_map[:, :, 0] = (grid2d_position // 9) // scanline_length

        return index_map
