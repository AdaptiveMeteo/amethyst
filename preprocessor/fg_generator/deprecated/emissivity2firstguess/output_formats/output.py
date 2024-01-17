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

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

class Output(object):
    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def save(self, pos, land_or_water, bias, functions, covariance):
        raise NotImplementedError
