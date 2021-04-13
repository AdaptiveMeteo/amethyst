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

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

class UnivariateMultiSpline(object):
    """
    This class extends the UnivariateSpline class. Indeed, it can handle
    interpolation in one dimension for multidimensional array.
    
    Args:
        - *x*: a 1D array for the position
        - *y*: a ND array such that the first dimension is the same than
          the x dimension
    """
    
    def __init__(self, x, y, *args, **kwargs):
        # Create a place where all the UnivariateSpline will be saved
        # (or, to be more precise, their pointers)
        self.interp = np.empty_like(y[0], dtype=np.object)

        # Prepare a iterator for each elements of y (beside the first axis)
        it = np.nditer(y[0], flags=['multi_index'])
        
        # Now run on y, creating a UnivariateSpline for each element
        while not it.finished:
            # Create a slice that take only the first axis of y for a
            # specidic position selected by it.multi_index
            y_indx = (Ellipsis, ) + it.multi_index
            self.interp[it.multi_index] = UnivariateSpline(
                                                           x,
                                                           y[y_indx],
                                                           *args,
                                                           **kwargs
                                                           )
            it.iternext()


    def __call__(self, x):
        # Prepare the space for the output
        output = np.empty_like(self.interp, dtype=np.float32)
        
        # Now iterate on that vector
        it = np.nditer(output, flags=['multi_index'])
        
        # Evaluate the interpolation on x point
        while not it.finished:
            output[it.multi_index] = self.interp[it.multi_index](x)
            it.iternext()
        
        return output



class TempExtrapolator(object):
    """
    A Temperature Extrapolator is an object which aims to extend the range
    of a specific temperature profile over the most external part of the 
    atmosphere.

    To be initialized, such an object requires some coefficients:
    The earth is divided in a few zones based on the latitude and for each
    zone, we need some coefficients for some specific levels of the profile.
    In particular, for each zone and each level of pressure, we need the
    coefficients of the linear approximation of the dipendence of the
    temperature from the temperature measured on lower levels.

    Therefore, the NetCDF file to initialize such an object must contain the
    following variables:

    - *zones*: a 1D variable with the latitudes of the different zones
    - *extrapolated_levels*: a 1D array with the pressure in hPa of the
      levels that must be extrapolated
    - *known_levels*: a 1D array with the pressure in hPa of the levels
      where the temperature should be known to use the extrapolator
    - *coefficients*: A 3D array with the coefficients for the extrapolation

    The extrapolation works in the following way: for a point on the earth
    with latitude :math:`z_i` where :math:`z` is the array of the zones we
    can compute the temperature at the pressure :math:`e_j` (where :math:`e`
    is the array of the extrapolated levels) using the relation

    .. math::

       T(e_j) = c_{ij0} + c_{ij1} \cdot T(k_1) + c_{ij2} \cdot T(k_2) + \ldots 
                + c_{ijn} T(K_n)

    where :math:`(k_1, k_2, \ldots, k_n)` is the array of the levels where the
    temperature is known, :math:`c` is the array of the coefficients and 
    :math:`T` is the function that returns the temperature for a given
    pressure.

    Arg:
        - *coeff_file*: the path of a NetCDF file with the needed coefficients
    """

    def __init__(self, coeff_file):
        with Dataset(coeff_file, 'r') as f:
            # The variable with the latitude of the
            # zone
            self.zones = f.variables['zones'][:]

            # This is the pressure levels that we can compute and
            # for which we have the coefficients
            ext_levels = f.variables['extrapolated_levels'][:]
            self.ext_levels = np.array(ext_levels, dtype=np.float32)

            # Now we save the minimum and the maximum of the extrapolated
            # levels. These values define the range where the temperature
            # profiles will be extrapolated
            self.min_level = np.min(self.ext_levels)
            self.max_level = np.max(self.ext_levels)

            # These are the levels that we need to compute the 
            # extrapolation
            known_levels = f.variables['known_levels'][:]
            self.known_levels = np.array(known_levels, dtype=np.float32)

            # These are the coefficients
            coeff = f.variables['coefficients'][:]
            coeff = np.array(coeff, dtype=np.float32)

            # To get the values for each latitude, we interpolate
            # all the coefficients over the zones
            degree = min(self.zones.size - 1, 3)
            self.interp = UnivariateMultiSpline(self.zones, coeff, k=degree)


    def extrapolate_temperature(self, temp_profile, lat):
        """
        Given a temperature profile of a point and its latitude, extrapolate
        the values of the temperature over the most external part of the 
        atmosphere

        Args:
            - *profile*: a callable object. Calling it with the value of
              the pressure in a specific point (in hPa), it must return
              the value of the temperature in that point
            - *lat*: the latitude of the point

        Return:
            - *extrapolate_profile*: a callable objet. Calling it with the
              value of the pressure in a specific point (in mb), it must return
              the value of the temperature in that point.

              The allowed values are the one in the range defined by
              extrapolate_profile.min and extrapolate_profile.max
        """
        # We approximate the north and the south emispheres as symmetrical
        lat = np.abs(lat)

        # Get the coefficients for the current location
        coeff = self.interp(lat)

        # Now we get the temperature on the known points
        t_on_known_levels = np.empty( coeff.shape[-1])

        # For the first value we put 1 because in that value the coefficient
        # is the mean and must not be changed
        t_on_known_levels[0] = 1

        for i in range(1, t_on_known_levels.size):
            t_on_known_levels[i] = temp_profile(self.known_levels[i-1])

        # To obtain the extrapolated temperature, we can use a matmult
        t_on_extrapolated = np.dot(coeff, t_on_known_levels)

        # Create an interpolator for the temperature which is
        # linear with the logaritm of the pressure.
        log_temp_extrap = interp1d(
                                   np.log(self.ext_levels),
                                   t_on_extrapolated,
                                   kind='linear',
                                   copy=False,
                                   bounds_error=True,
                                   )

        # Finally, we can return the extrapolator as a callable object
        def temp_extrap(x):
            return np.float32(log_temp_extrap(np.log(x)))

        # Save the minimum and the maximum of the values that the
        # extrapolator can handle
        temp_extrap.min_level = self.min_level
        temp_extrap.max_level = self.max_level

        return temp_extrap
