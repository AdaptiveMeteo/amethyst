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

from utilities.split_by_time import split_by_time

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

def split_by_time_test01():
    # Create an array for all the days of january
    values = np.arange('2016-01', '2016-02', dtype='datetime64[D]')
    # Split the points: regrup the ones that are closest to the first
    # day of January and the ones that are closest to the first day of
    # February. Until 15, the nearest day is 2016-01-01, afterwards it
    # is 2016-02-01
    references = np.array(['2016-01', '2016-02'], dtype='datetime64[D]')
    a = split_by_time(values, references)

    # Check the splitting
    assert len(a) == 2
    assert a[0] == (0,0,16)
    assert a[1] == (1,16,31)

def split_by_time_test02():
    # Create an array for all the days of april
    values = np.arange('2016-04', '2016-05', dtype='datetime64[D]')
    # Split the points: this test looks like the one split_by_time_test01
    # but now the 2016-04-15 it is in the middle
    references = np.array(['2016-04', '2016-05'], dtype='datetime64[D]')
    a = split_by_time(values, references)

    # Check the splitting
    assert len(a) == 2
    assert a[0] == (0,0,16) or a[0] == (0,0,15)
    assert a[1] == (1,16,30) or a[1] == (1,15,30)
    assert a[0][2] == a[1][1]


def split_by_time_test03():
    # Create an array with a lot of different dates
    values1 = np.arange('2016-01', '2016-12', dtype='datetime64[D]')
    values1_s = np.array(values1, dtype='datetime64[s]')
    values2 = np.arange('2016-12-15T00:00:00Z', '2016-12-18T23:59:59Z',
                        dtype='datetime64[s]')
    values = np.concatenate((values1_s, values2))
    
    # Reorder the date in a random way
    np.random.seed(1)
    np.random.shuffle(values)

    references = np.array(['2016-02', '2016-06', '2017-01-10'],
                          dtype='datetime64[D]')

    
    a = split_by_time(values, references)
    
    # Check that all the dates has been counted
    check = 0
    for ref, start, stop in a:
        check += stop - start
    assert check==values.size
    
    # Check that, for all the dates, the ref is the right one
    for ref, start, stop in a:
        ref_val = references[ref]
        for i in range(start, stop):
            val = values[i]
            t_dist = np.abs(ref_val - val)
            for k in references:
                assert t_dist <= np.abs(k-val)
