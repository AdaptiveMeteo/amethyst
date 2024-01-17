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

from __future__ import print_function

from sys import stderr

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"


def test_import_numpy():
    import numpy

def test_import_scipy():
    import numpy

def test_import_argparse():
    import argparse

def test_import_netcdf():
    try:
        import netCDF4
    except:
        print("\nThis package requires netCDF4 to be installed. You can install"
              " it by executing:\n  # pip install netCDF4\non your terminal",
              file=stderr)
        raise

def test_import_supersmoother():
    try:
        import supersmoother
    except:
        print("\nThis package requires supersmoother to be installed. You can install"
              " it by executing:\n  # pip install supersmoother\non your terminal",
              file=stderr)
        raise
