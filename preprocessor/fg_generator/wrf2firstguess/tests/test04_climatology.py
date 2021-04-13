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

from nose.tools import raises
from os.path import isfile
import numpy as np


from utilities.climatology import OZONE_PROFILE_FILES, CLIMAT_MATRIX, Belt,\
                                  generate_ozone_profile

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

def ozone_climatology_file_test01():
    assert isfile(OZONE_PROFILE_FILES), 'Unable to find the file {} with the '\
                                        'coefficients for the ozone '\
                                        'climatology'.format(OZONE_PROFILE_FILES)

def ozone_climatology_file_test02():
    assert CLIMAT_MATRIX is not None, 'Unable to use the coefficients of file'\
                                      ' {} to generate an ozone climatology'\
                                      ' matrix'.format(OZONE_PROFILE_FILES)

def ozone_climatology_file_test03():
    # Check if the matrix is well formatted
    assert len(CLIMAT_MATRIX.shape) == 2
    assert CLIMAT_MATRIX.shape[-1] == 5
    
    assert len(CLIMAT_MATRIX.levels.shape) == 1
    assert CLIMAT_MATRIX.levels.size == CLIMAT_MATRIX.shape[0]


def Belt_test01():
    b = Belt(10, 40)
    assert 10 in b
    assert 40 not in b
    assert 20 in b
    assert 30 in b
    assert 23.123 in b
    assert 0 not in b
    assert -10 in b
    assert -40 not in b
    assert -20 in b
    assert -30 in b
    assert -23.123 in b

@raises(ValueError)
def Belt_test02():
    b = Belt(-10, 10)

@raises(ValueError)
def Belt_test03():
    b = Belt(10, 8)

def generate_ozone_profile_test01():
    # Create an object that I will use as a placeholder
    # for a profile
    class Object(object): pass
    p = Object()

    min_levels = np.min(CLIMAT_MATRIX.levels)
    max_levels = np.max(CLIMAT_MATRIX.levels)
    
    levels = np.linspace(min_levels, max_levels, 100)

    for i in range(-180, 181):
        p.lat = i / 2.
        for month in range(0, 12):
            generate_ozone_profile(p, levels, month)

@raises(ValueError)
def generate_ozone_profile_test02():
    # Test that returns error if the levels are too hight
    # Create an object that I will use as a placeholder
    # for a profile
    class Object(object): pass
    p = Object()
    p.lat = 25

    min_levels = np.min(CLIMAT_MATRIX.levels)
    max_levels = np.max(CLIMAT_MATRIX.levels)
    
    levels = np.linspace(min_levels / 2., max_levels, 100)

    generate_ozone_profile(p, levels, 9)

@raises(ValueError)
def generate_ozone_profile_test03():
    # Test that returns error if the levels are too low
    # Create an object that I will use as a placeholder
    # for a profile
    class Object(object): pass
    p = Object()
    p.lat = 60

    min_levels = np.min(CLIMAT_MATRIX.levels)
    max_levels = np.max(CLIMAT_MATRIX.levels)
    
    levels = np.linspace(min_levels, max_levels+100, 100)

    generate_ozone_profile(p, levels, 6)

@raises(ValueError)
def generate_ozone_profile_test04():
    # Test that returns error if the latitude is meaningless
    # Create an object that I will use as a placeholder
    # for a profile
    class Object(object): pass
    p = Object()
    p.lat = 91

    min_levels = np.min(CLIMAT_MATRIX.levels)
    max_levels = np.max(CLIMAT_MATRIX.levels)
    
    levels = np.linspace(min_levels, max_levels, 100)

    generate_ozone_profile(p, levels, 1)
    

@raises(ValueError)
def generate_ozone_profile_test04():
    # Test that returns error if the latitude is meaningless
    # Create an object that I will use as a placeholder
    # for a profile
    class Object(object): pass
    p = Object()
    p.lat = -91

    min_levels = np.min(CLIMAT_MATRIX.levels)
    max_levels = np.max(CLIMAT_MATRIX.levels)
    
    levels = np.linspace(min_levels, max_levels, 100)

    generate_ozone_profile(p, levels, 1)

@raises(ValueError)
def generate_ozone_profile_test04():
    # Test that returns error if the CLIMAT_MATRIX is None
    class Object(object): pass
    p = Object()
    p.lat = -91

    min_levels = np.min(CLIMAT_MATRIX.levels)
    max_levels = np.max(CLIMAT_MATRIX.levels)
    
    levels = np.linspace(min_levels, max_levels, 100)

    generate_ozone_profile(p, levels, 1, climat_matrix=None)
    
