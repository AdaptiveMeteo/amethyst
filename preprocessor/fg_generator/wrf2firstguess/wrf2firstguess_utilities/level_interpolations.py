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

import numpy as np
from scipy.interpolate import interp1d
from supersmoother import SuperSmoother
from amethyst_config import preprocessor_vars

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

# The minimum value of the water vapour reached in the highest levels of the
# atmosphere
MIN_WATER_VAPOUR = preprocessor_vars["first_guess"]["min_water_vapor"]

# The pressure (in hPa) where such value is reached
P_MIN_WATER_VAPOUR = preprocessor_vars["first_guess"]["p_min_water_vapor"]

# The pressure (in hPa) where the profile will start to be smoothed
# If you put here the same value of P_MIN_WATER_VAPOUR, no smoothing
# will be performed
WATER_VAPOUR_SMOOTH_AFTER = preprocessor_vars["first_guess"]["water_vapor_smooth_after"]

log = logging.getLogger(__name__)

def interp_temperature_over_levels(profile, levels, temp_extrap):
    """
    Given a WrfProfile and some levels, return an interpolated profile of the
    temperature over the levels. Use a
    :meth:`utilities.temp_extrapolator.TempExtrapolator` to compute the values
    of the levels that are too hight for the model and try to be smooth in
    the junction.

    For the values that are too low for the model, uses the lowest value.

    To obtain a smooth profile from the extrapolator, SuperSmoother_ is
    applied on the values returned by the temp extrapolator over the levels
    that are in its range.

    .. _SuperSmoother: https://pypi.python.org/pypi/supersmoother

    If some levels are both in the range covered by the extrapolator and in the
    one covered by the model profile, the transition is performed linearly
    respect to the logarithm of the levels.

    Args:
        - *profile*: A WrfProfile
        - *levels*: An array with the values of the pressure in
          hPa ordered in descending order.
        - *temp_extrap*: A TempExtrapolator to generate the values of the
          temperature for the highest levels

    Return:
        A 1D numpy array with the value of the temperature for each level
    """

    # Create an interpolator for the temperature which is
    # linear with the logaritm of the pressure. We need to
    # change the order of the entries so that profile.pressure
    # will be ordered in ascending order
    log_temp_interp = interp1d(
                               np.log(profile.pressure_levels[::-1]),
                               profile.temperature[::-1],
                               kind='linear',
                               copy=False,
                               bounds_error=True,
                               )

    # Create a callable object that, for each value of the pressure in the
    # interval between min(profile.pressure) and max(profile.pressure),
    # returns the temperature
    min_level = np.min(profile.pressure_levels)
    max_level = np.max(profile.pressure_levels)
    def temp_interp(x):
        if x < min_level:
            log.debug('Trying to read the temperature from the model for '
                      'the level {} which is too high; using {} instead'
                      ''.format(x, min_level))
            x = min_level
        if x > max_level:
            log.debug('Trying to read the temperature from the model for '
                      'the level {} which is too low; using {} instead'
                      ''.format(x, max_level))
            x = max_level
        return log_temp_interp(np.log(x))

    # Now we can use the t_extr to create a profile of the
    # temperature for the highest levels
    temp_extrap = temp_extrap.extrapolate_temperature(temp_interp,
                                                      profile.lat,
                                                      )

    # Prepare the temperature profile that will be the output of this function
    temperature = np.empty(levels.shape, dtype=np.float32)

    # The next step is to glue together the values of the extrapolator and
    # of the real profile. For this reason, we divide the levels in different
    # zones
    model_levels = profile.pressure_levels
    start_model = np.max(model_levels)
    end_model = np.min(model_levels)
    start_extrapolator = temp_extrap.max_level
    end_extrapolator = temp_extrap.min_level

    # We need to check that the two zones have no hole in the middle
    in_the_middle = np.where(np.logical_and(levels < end_model,
                                            levels > start_extrapolator))[0]
    if in_the_middle.size != 0:
        log.error('There are some levels that are among the temperature '
                  'profile and the extrapolator. No reasonable temperature '
                  'can be computed for these levels.')
        raise ValueError('Levels among temperature profile and extrapolator '
                         'are not allowed')

    # The zone among the ground and the starting of the model
    under_model = np.where(levels > start_model)[0]
    # which can be filled with the value of the temperature in the
    # first cell of the model
    if under_model.size != 0:
        log.debug('Filling with the first model profile value ({:.2f} degrees)'
                  ' the indices of the levels that are under the first level '
                  'of the model (from {} to {})'.format(profile.temperature[0],
                                                        under_model[0],
                                                        under_model[-1]))
        temperature[under_model] = profile.temperature[0]

    # If the model can cover all the levels also in the upper sky,
    # we can complete the profile using only the model
    if levels[-1] >= end_model:
        log.debug('The model covers all the highest levels. Not using '
                  'the extrapolator')
        # The zone covered by the model
        covered_by_model = np.where(np.logical_and(levels <= start_model,
                                                   levels >= end_model))[0]
        log.debug('Using model for indices from {} to {}'
                  ''.format(covered_by_model[0], covered_by_model[-1]))
        for i in covered_by_model:
            temperature[i] = temp_interp(levels[i])
        return temperature

    # If this is not the case, we need to use the extrapolator!
    # Check which are the levels that the extrapolator can handle
    covered_by_extrap = np.where(np.logical_and(levels <= start_extrapolator,
                                                levels >= end_extrapolator))[0]
    log.debug('Extrapolator covers indices from {} to {}'
              ''.format(covered_by_extrap[0], covered_by_extrap[-1]))

    # Even if there are some levels that are too hight for the extrapolator, it
    # will handle them anyway
    handled_by_extrap = np.where(levels<=start_extrapolator)[0]

    # Allocate the space for the values generated by the extrapolator
    extrap_values = np.empty_like(handled_by_extrap, dtype=np.float32)

    # Fill the values that are completely covered by the extrapolator
    for i, j in enumerate(covered_by_extrap):
        extrap_values[i] = temp_extrap(levels[j])

    # Now the procedure if there are some levels that are too hight even for
    # the extrapolator
    if levels[-1] < end_extrapolator:
        log.warning('The temperature of the top levels can not be correctly '
                    'generated because they are outside the range of the '
                    'extrapolator. They will be filled with the highest '
                    'level of the extrapolator')
        # Fill the last part of extrap_values with the last value of the
        # extrapolator
        extrap_values[len(covered_by_extrap):] = temp_extrap(end_extrapolator)

    # Smooth the values of the extrapolator
    sm = SuperSmoother()
    sm.fit(np.log(levels[handled_by_extrap]), extrap_values)
    extrap_values2 = sm.predict(np.log(levels[handled_by_extrap]))
    max_change = np.max(np.abs(extrap_values2 - extrap_values))
    log.debug('Smoothing the extrapolated function; max change is '
              '{:.2f} degrees'.format(max_change))
    # Comment this line to avoid smoothing
    extrap_values = np.array(extrap_values2, dtype=np.float32)

    # Now we have to deal with the problem of the overlapping between the
    # profile of the model and the extrapolator. In principle, we could use
    # just the values from the profile but, in that case, the function could
    # not be smooth when the profile ends and the extrapolator starts. Instead,
    # we decided to move linearly (in the logarithmic space) from the values
    # of the profile to the values of the extrapolator.
    overlapping = np.where(np.logical_and(levels >= end_model,
                                          levels <= start_extrapolator))[0]
    if overlapping.size > 0:
        overlap_levels = levels[overlapping]
        start_overlapping = overlap_levels[0]
        end_overlapping = overlap_levels[-1]
        log.debug('The model and the extrapolator overlaps between level {} '
                  'and level {} ({:.2f} hPa and {:.2f} hPa)'
                  ''.format(overlapping[0], overlapping[-1],
                            start_overlapping, end_overlapping))

        # The values of the extrapolator that overlaps
        overlap_extrap_values = extrap_values[:len(overlapping)]

        # The values of the profile that overlaps
        overlap_prof_values = np.empty_like(overlap_extrap_values)
        for i, j in enumerate(overlapping):
            overlap_prof_values[i] = temp_interp(levels[j])

        # We now compute the range of the overlapping vector in the log space
        rng = np.log(start_overlapping) - np.log(end_overlapping)

        # Now we compute the coefficients for the transition. This is a vector
        # that is 1 on the level of start_overlapping and is 0 on the level of
        # end_overlapping. Moreover, if x is a level such that log(x) is the 
        # midpoint of log(end_overlapping) and log(start_overlapping), then
        # transition is 0.5 on x
        trans = (np.log(overlap_levels) - np.log(end_overlapping)) / rng

        # Save the temperature for the overlapping levels
        temperature[overlapping] = trans * overlap_prof_values + \
                                   (1-trans) * overlap_extrap_values

    # Now we must save the values of the profile of the model before the
    # range of the extrapolator
    covered_only_by_model = np.where(np.logical_and(levels <= start_model,
                                              levels > start_extrapolator))[0]
    log.debug('Using values of the model between indices {} and {} ({:.2f} '
              'hPa and {:.2f} hPa)'.format(covered_only_by_model[0],
                                            covered_only_by_model[-1],
                                            levels[covered_only_by_model[0]],
                                            levels[covered_only_by_model[-1]]))

    for i in covered_only_by_model:
        temperature[i] = temp_interp(levels[i])

    # Finally, the values that are covered only by the extrapolator, which are
    # the ones just after the overlapping section    
    log.debug('Using values of the extrapolator between indices {} and {} '
              '({:.2f} hPa and {:.2f} hPa)'
              ''.format(overlapping[-1] + 1,
                        len(temperature) - 1,
                        levels[overlapping[-1]+1],
                        levels[len(temperature)-1]))

    temperature[overlapping[-1] + 1:] = extrap_values[len(overlapping) : ]

    return temperature


def interp_water_vapour_over_levels(
                                    profile,
                                    levels,
                                    fixed_vapour_point=MIN_WATER_VAPOUR,
                                    fixed_vapour_reached = P_MIN_WATER_VAPOUR,
                                    smooth_after = WATER_VAPOUR_SMOOTH_AFTER,
                                    min_value = 0.00001
                                    ):
    """
    Given a WrfProfile and some levels, return an interpolated profile of the
    water vapour over the levels. Use a fixed value if there are some levels
    that are too hight for the model and be smooth in the junction.

    We expect that, in the highest layers of the atmosphere, the water vapour
    is fixed at a specific value (usually, 0.003). Therefore, if the model do
    not reach the highest levels that we have as input, we can use this
    information to try to generate good approximations of the real values.

    The values of the levels that are above the pressure where the water vapour
    is fixed will be set to a specific value.

    The values of the levels that are under the lowest level of the model will
    be filled with the lowest value of the model.

    To ensure that the profile is smooth in the top levels (where we put a
    fixed value) we smooth with SuperSmoother_ all the levels above a specified
    pressure value.

    .. _SuperSmoother: https://pypi.python.org/pypi/supersmoother

    Args:
        - *profile*: A WrfProfile
        - *levels*: An array with the values of the pressure in
          hPa ordered in descending order.
        - *fixed_vapour_point* (optional): the fixed value that the returned
          profile will reach in the highest part of the atmosphere
        - *fixed_vapour_reached* (optional): the pressure in (hPa) when the
          fixed value will be reached. Every information in the model after
          that point will be discarded
        - *smooth_after* (optional): the pressure (in hPa) after which the
          profile will be smoothed
        - *min_value* (optiona): if min_value is not None, all the values that
          are smaller than min_value are substituted with min_value. This
          ensures that there are, for example, no negative values related to
          numerical instability

    Return:
        A 1D numpy array with the value of the water vapour for each level
    """
    # Prepare the array where the values will be saved
    water_vapour = np.empty_like(levels)

    # Find the best indices to approximate when the fixed point is reached
    # The best index is the one with the smaller distance in the log space
    dist = np.abs(np.log(levels) - np.log(fixed_vapour_reached) )
    fixed_index = np.argmin(dist)

    log.debug('Setting all the levels of the water vapour above the index '
              '{} ({:.2f} hPa) to {:.3f}'.format(fixed_index,
                                                 levels[fixed_index],
                                                 fixed_vapour_point))
    water_vapour[fixed_index:] = fixed_vapour_point

    end_model = np.min(profile.pressure_levels)
    start_model = np.max(profile.pressure_levels)

    if end_model > levels[fixed_index]:
        log.error('Unable to interpolate the water vapour: the model reaches '
                  '{:.3f} hPa but we need to reach {:.3f} hPa (the fixed '
                  'point).'.format(end_model, levels[fixed_index]))
        raise ValueError('The fixed point is too hight or the model highest '
                         'level is too low')

    # Compute the first index after which we are inside the model
    lowest_model_index = np.where(levels <= start_model)[0][0]

    # Create an interpolator for the water vapour which is
    # linear with the logaritm of the pressure. We need to
    # change the order of the entries so that profile.pressure
    # will be ordered in ascending order
    log_wv_interp = interp1d(
                             np.log(profile.pressure_levels[::-1]),
                             profile.water_vapour[::-1],
                             kind='linear',
                             copy=False,
                             bounds_error=True,
                             )
    # To use the log_wv_interp, we need the value of the levels in log space
    levels_log = np.log(levels[lowest_model_index : fixed_index])

    log.debug('Filling levels from {} to {} with the values of the model '
              'profile interpolated in the log space (from {:.3f} to {:.3f} '
              'hPa)'.format(lowest_model_index, fixed_index-1,
                            levels[lowest_model_index], levels[fixed_index-1]))
    water_vapour[lowest_model_index:fixed_index] = log_wv_interp(levels_log)

    if lowest_model_index != 0:
        log.debug('Filling levels from 0 to {} with the first value of the '
                  'model profile (from {:.3f} to {:.3f} hPa)'
                  ''.format(lowest_model_index -1, levels[0],
                            levels[lowest_model_index -1]))
        water_vapour[:lowest_model_index] = water_vapour[lowest_model_index]

    # Now we ensure that all the values are greater than min_value
    if min_value is not None:
        water_vapour[water_vapour < min_value] = min_value

    # Finally we use the SuperSmoother library to generate a smooth profile
    # (expecially before the fixed point)
    # Find the best indices to approximate when we must start the smoother
    # Again, the best index is the one with the smaller distance in the log
    # space
    dist = np.abs(np.log(levels) - np.log(smooth_after) )
    fixed_index = np.argmin(dist)

    # Check the first index that must be smoothed
    start_smooth_index = np.where(levels <= smooth_after)[0][0]
    log.debug('The smoother will be executed starting from index {}'
              ''.format(start_smooth_index))

    # Smooth the values in the log space
    #print('start_smooth_index: {}'.format(start_smooth_index))
    #print('water_vapour[start_smooth_index:]: {}'.format(water_vapour[start_smooth_index:]))
    sm = SuperSmoother()
    sm.fit(
           np.log(levels[start_smooth_index:]),
           water_vapour[start_smooth_index:]
           )

    smooth_values = sm.predict(np.log(levels[start_smooth_index:]))

    old_values = water_vapour[start_smooth_index:]
    max_change = np.max(np.abs(smooth_values - old_values))
    log.debug('Smoothing the water vapour profile; max change is '
              '{:f}'.format(max_change))

    # Comment this line to avoid smoothing
    water_vapour[start_smooth_index:] = smooth_values

    return water_vapour

