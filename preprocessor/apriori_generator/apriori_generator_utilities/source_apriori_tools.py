#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__     = 'Paolo Scaccia <paolo.scaccia@adaptivemeteo.com>'
__copyright__  = "Copyright 2023, Adaptive Meteo S.r.l"
__credits__    = ["Paolo Scaccia", "Paolo Antonelli"]
__license__    = "GPL"
__version__    = "1.0"
__maintainer__ = "Paolo Scaccia"
__email__      = "paolo.scaccia@adaptivemeteo.com"

import numpy as np
import logging


if __name__ == '__main__':
    log = logging.getLogger()
else:
    log = logging.getLogger(__name__)

class StaticAprioriGridError(Exception):
    """
    This error is raised when the input pressure grid exceeds the 
    one contained in the static apriori source, used as reference
    for the interpolation.
    """

    pass

def compute_LERP_jacobian(ref_grid, input_grid):
    """
        This function compute the Linear Operator associated
        with the interpolation of the input grid (input_grid)
        above a reference grid (ref_grid).
        
        NOTE: extrapolation is assumed not valid. All values
        outside the data range will be considered interpolated
        with the bound values of the reference grid (as for 
        the scipy function scipy.interpolate.interp1d with 
        the fill_value attribute set to (ref_grid.min(), 
        ref_grid.max()).

    Parameters
    ----------
    ref_grid : np.array
        Reference grid for the interpolation.
        It is assumed sorted in an ascending order.
                                             
    input_grid : np.array
        Reference grid for the interpolation.
        It is assumed sorted in an ascending order.

    Returns
    -------
    jac : np.ndarray
        Jacobian matrix (Linear Operator) associated
        to the given interpolation.

    """
    size_input = input_grid.size
    size_ref   = ref_grid.size
    
    # Raise error if the input grid falls 
    # outside the reference data range
    if ref_grid.min() > input_grid.min() or ref_grid.max() < input_grid.max():
      raise StaticAprioriGridError('Input grid falls outside '
                                   'the reference for the interpolation')

    # Init the Jacobian matrix
    jac = np.zeros((size_input, size_ref))
    
    # Get the index if the reference grid in which 
    # each input value falls.
    bin_indxes = np.digitize(input_grid, ref_grid)

    # Compute dx array (normalization)
    dx         = ref_grid[1:] - ref_grid[:-1]
    
    for i_row, (bin_indx, grid_value) in enumerate(zip(bin_indxes,input_grid)):
        
        edge_indx = bin_indx - 1
        if edge_indx < 0:
            # The input value falls before the reference min
            jac[i_row,0] = 1
        elif edge_indx >= size_ref - 1:
            # The input value falls above the reference max
            jac[i_row,-1] = 1
        else:
            norm = dx[edge_indx]
            # Fill the (i_row)-th Jacobian row
            jac[i_row, edge_indx] = (ref_grid[edge_indx + 1] - grid_value ) / norm
            jac[i_row, edge_indx + 1] = (grid_value - ref_grid[edge_indx]) / norm
            
    return jac

def scale_apriori_covariance(source_apriori, source_grid, fg_grid, J = None):
    
    if J is None:
        # Compute the Matrix associated with the linear interpolation
        J = compute_LERP_jacobian(source_grid, fg_grid)
    # Return the rescaled apriori covariance
    return (J.dot(source_apriori)).dot(J.transpose())
    
    
    
    
