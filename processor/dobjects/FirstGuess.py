#!/usr/bin/env python
"""
Class to retrieve data from a First Guess Atmospheric state.
"""

import numpy as np
from netCDF4 import Dataset
import amethyst_config

class FirstGuess(object):
    """
    This class is a wrapper around a netCDF data file
    """
    #PaoloA 31-032021
    EIGENVALUES_LAND=amethyst_config.processor_vars['eigenforland']
    EIGENVALUES_SEA=amethyst_config.processor_vars['eigenforsea']
    CO2 = amethyst_config.processor_vars['constant_co2']['value']
    
    def __init__(self, datafile, co2=CO2, eigen_land=EIGENVALUES_LAND, 
                 eigen_sea=EIGENVALUES_SEA, var_selection = None):
        """
        Initialize the FirstGuess object
        """
        self.df = Dataset(datafile, mode='r')
        self.eigen_land = eigen_land
        self.eigen_sea = eigen_sea

        self.levels = len(self.df.groups['atmospheric_components'].dimensions['number_of_atmospheric_levels'])
        self.on_land = np.array(
            self.df.groups['surface_components'].variables['LandOrWater'][:],
            dtype=np.bool,
        )
        
        # Use all variables if there is no selection
        if var_selection == None:
            var_selection = [ True for x in range(5)]
        
        # Build Level Selection
        self.LEVEL_SELECTION = [ x+1 for x in range(len(var_selection[1:])*self.levels) if var_selection[1:][x//self.levels] ]
        if var_selection[0]:
            self.LEVEL_SELECTION.append( self.LEVEL_SELECTION[-1] + 1)

        self.co2 = co2
        self.eigen_land = eigen_land
        self.eigen_sea = eigen_sea

    def eigenvalues(self, obs):
        if self.on_land[obs]:
            return self.eigen_land
        else:
            return self.eigen_sea

    def xdim(self, obs):
        return np.array([self.levels, self.levels, self.levels, self.levels, 1, self.eigenvalues(obs)],
                        dtype=np.int32)

    def varindx(self, obs):
        #PaoloA 15082017
        levsel = self.LEVEL_SELECTION + list(range(326, 326 + self.eigenvalues(obs)))
        levsel_array = np.array(
            levsel,
            dtype=np.int32
        )
        return levsel_array - 1

    def state_vector(self, obs):
        """
        Get one of the first guesses in the file [0:obsnum]
        """
        size = np.sum(self.xdim(obs))

        p = np.empty(self.levels, dtype=np.float64)
        x0 = np.empty(size, dtype=np.float64)
        xa=x0

        p[:] = self.df.groups['atmospheric_components'].variables['p'][obs,:]

        x0[0: self.levels] = self.df.groups['atmospheric_components'].variables['T'][obs,:]
        x0[self.levels:self.levels*2]   = np.log(self.df.groups['atmospheric_components'].variables['q'][obs,:])
        x0[self.levels*2:self.levels*3] = self.co2
        x0[self.levels*3:self.levels*4] = np.log(self.df.groups['atmospheric_components'].variables['O3'][obs,:])
        x0[self.levels*4] = self.df.groups['atmospheric_components'].variables['skT'][obs]
        x0[self.levels*4 + 1:self.levels*4 + 1 + self.eigenvalues(obs)] = 0
        
        return [p, x0, xa]


class firstguess(object):
    """
    This class is a wrapper around a FirstGuess object.
    It is used to read just an observation and keep it in memory
    until the object is deleted.
    """
    def __init__(self, fguess, obs):
        """
        Initialize the object using a FirstGuess and the obs index
        """
        self.xdim = fguess.xdim(obs)
        self.varindx = fguess.varindx(obs)
        [self.p, self.x0, self.xa] = fguess.state_vector(obs)
        self.obs = obs

    def state_vector(self, obs):
        """
        Return the in-memory object
        """
        if obs != self.obs:
            raise IndexError('Object initilized with obs='+repr(self.obs)+
                             ' but obs='+repr(obs)+' requested !')
        return [self.p, self.x0, self.xa]

