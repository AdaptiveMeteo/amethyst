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
import numpy as np

from preprocessor.fg_generator.wrf2firstguess.wrf2firstguess_utilities.first_guess import AtmosphericFirstGuess

__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__copyright__ = "Copyright 2016, eXact-lab and Paolo Antonelli"
__credits__ = ["Stefano Piani", "Paolo Antonelli"]
__license__ = "GPL"
__maintainer__ = "Stefano Piani"
__email__ = "stefano.piani@exact-lab.it"

@raises(NotImplementedError)
def AtmosphericFirstGuess_test01():
    fg = AtmosphericFirstGuess()

def AtmosphericFirstGuess_test02():
    # A shortcut for the name of the class
    FG = AtmosphericFirstGuess

    # Create an AtmosphericFirstGuess object avoiding the init function
    fg = FG.__new__(FG)    
    fg.open()
    fg.close()
    
    with fg:
        pass

@raises(NotImplementedError)
def AtmosphericFirstGuess_test03():
    # A shortcut for the name of the class
    FG = AtmosphericFirstGuess

    # Create an AtmosphericFirstGuess object avoiding the init function
    fg = FG.__new__(FG)    
    
    fg.pressure_levels

@raises(NotImplementedError)
def AtmosphericFirstGuess_test04():
    # A shortcut for the name of the class
    FG = AtmosphericFirstGuess

    # Create an AtmosphericFirstGuess object avoiding the init function
    fg = FG.__new__(FG)    
    
    fg.temperature

@raises(NotImplementedError)
def AtmosphericFirstGuess_test05():
    # A shortcut for the name of the class
    FG = AtmosphericFirstGuess

    # Create an AtmosphericFirstGuess object avoiding the init function
    fg = FG.__new__(FG)    
    
    fg.water_vapour

@raises(NotImplementedError)
def AtmosphericFirstGuess_test06():
    # A shortcut for the name of the class
    FG = AtmosphericFirstGuess

    # Create an AtmosphericFirstGuess object avoiding the init function
    fg = FG.__new__(FG)    
    
    fg.ozone
