#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 11:02:09 2024

@author: paoloscaccia
"""

import logging

import numpy as np
from scipy.interpolate import interp1d
from supersmoother import SuperSmoother

__author__     = 'Paolo Scaccia <paolo.scaccia@adaptivemeteo.com>'
__copyright__  = "Copyright 2023, Adaptive Meteo S.r.l."
__credits__    = ["Paolo Scaccia", "Paolo Antonelli"]
__license__    = "GPL"
__maintainer__ = "Paolo Scaccia"
__email__      = "paolo.scaccia@adaptivemeteo.com"

log = logging.getLogger(__name__)

def merge_profiles(self, climatology_profile, source_profile):
        """
        This function is a tool for merging an additional source profile (i.e. WRF)
        with the climatology.
        
        Parameters
        ----------
        climatology_profile : class Profile
            Climatology profile as defined in preprocessor.utils.climatology
        source_profile : class SourceProfile (example WrfProfile)
            Additional source profile to be merged with the climatology profile.

        Returns
        -------

        """
        source_levels = source_profile.pressure_levels
        climatology_levels = climatology_profile.pressure
        
        start_model   = np.max(source_levels)
        end_model     = np.min(source_levels)
        start_climatology = np.max(climatology_levels)
        end_climatology   = np.min(climatology_levels)

        return
    
    
        