#!/usr/bin/env python

"""
Compute a surface emissivity spectrum from a set of model coefficients and
model functions

Input:
   SEC are the surface emissivity coefficients (vector)
Output:
   SEvalues = computed emissivity values
   SEwn = wavenumber scale of computed emissivity values

Data from File
   SEMF are the surface emissivity model functions
       (matrix with functions in rows)
   SEMFB is the offset (bias) for all the model functions

Note:
 The surface emissivity model functions are defined on their own
 wavenumber scale.
"""
from __future__ import print_function, division

from numpy import array, dot, exp, ones, zeros, interp,seterr
from numpy import bool_ as bool
from netCDF4 import Dataset

seterr(all='warn', over='raise')


class Emissivity(object):
    """
    This class provides a way to compute emissivity values
    """

    def __init__(self, datafile, eigen_land, eigen_sea):
        """Reads Model parameters from file"""
        self.df = Dataset(datafile, mode='r')

        self.obsnum = len(self.df.dimensions['number_of_FOVs'])

        self.SEwn = self.df.groups['surface_components'].variables['ModelWaveNumbers'][:]
        self.on_land = array(
            self.df.groups['surface_components'].variables['LandOrWater'][:],
            dtype=bool,
        )
        self.SEvalues = ones(self.SEwn.size)
        self.SEC = None
        self.eigen_land = eigen_land
        self.eigen_sea = eigen_sea

    def _del__(self):
        self.df.close()

    def eigenvalues(self, obs):
        if self.on_land[obs]:
            return self.eigen_land
        else:
            return self.eigen_sea

    def get(self, obs, SEC):
        """Proxy for getting emissivities"""
        if self.obsnum != 0:
            self.SEMF = self.df.groups['surface_components'].variables['ModelFunctions'][obs,:self.eigenvalues(obs), Ellipsis]
            self.SEMFB = self.df.groups['surface_components'].variables['ModelFunctionsBias'][obs, Ellipsis]
        self.compute_from_model(SEC)
        self.SEC = SEC
        return [self.SEwn, self.SEvalues]

    def get_jacobian(self, wvn, jcb, SEC=None):
        """Proxy for getting Jacobians diagonal elements"""
        # Note: The jacobian we want is the lblrtm surface jacobian
        # times each SEMF interpolated to the
        # calculation scale.
        if SEC is not None:
            self.SEC = SEC
        self.compute_from_model(self.SEC)
        interp_EmissVal = interp(wvn, self.SEwn, self.SEvalues, 0.999, 0.999)
        w = interp_EmissVal * (1.0 - interp_EmissVal)
        jac = zeros((jcb.size, self.SEC.size))
        for i in range(self.SEC.size):
            interp_SEMF = interp(wvn, self.SEwn, self.SEMF[i], 0.0, 0.0)
            #PaoloA 01-04-2021
            jac[:, i] = interp_SEMF * jcb * w
            #jac[:, i] = interp_SEMF.dot(jcb * w)*interp_SEMF
        return jac

    def compute_from_model(self, SEC):
        """
        Compute Emissivity values from emissivity coefficients
        from solution vector
        """
        # load surface emissivity coefficients from solution vector
        self.SEvalues = dot(self.SEMF.T, SEC) + self.SEMFB
        self.SEvalues = exp(self.SEvalues) / (1.0 + exp(self.SEvalues))
        self.SEC = SEC


class emissivity_model(object):
    """
    Reader is in the init. No file access afterwards.
    """
    def __init__(self, emiss, obs):
        self.SEwn = emiss.SEwn
        self.eigenvalues = emiss.eigenvalues(obs)
        if emiss.obsnum != 0:
            self.SEMF = emiss.df.groups['surface_components'].variables['ModelFunctions'][obs,:self.eigenvalues, Ellipsis]
            self.SEMFB = emiss.df.groups['surface_components'].variables['ModelFunctionsBias'][obs, Ellipsis]
        else:
            self.SEMF = emiss.SEMF
            self.SEMFB = emiss.SEMFB
        self.SEC = None
        self.SEvalues = ones(self.SEwn.size)
        self.obs = obs

    def get(self, obs, SEC):
        """Proxy for getting emissivities"""
        if obs != self.obs:
            raise IndexError('Object initilized with obs='+repr(self.obs)+
                             ' but obs='+repr(obs)+' requested !')
        self.SEC = SEC
        self.compute_from_model(self.SEC)
        return [self.SEwn, self.SEvalues]

    def get_jacobian(self, wvn, jcb, SEC=None):
        """Proxy for getting Jacobians diagonal elements"""
        # Note: The jacobian we want is the lblrtm surface jacobian
        # times each SEMF interpolated to the
        # calculation scale.
        if SEC is not None:
            self.SEC = SEC
        self.compute_from_model(self.SEC)
        interp_EmissVal = interp(wvn, self.SEwn, self.SEvalues, 0.999, 0.999)
        w = interp_EmissVal*(1.0-interp_EmissVal)
        jac = zeros((jcb.size, self.SEC.size))
        for i in range(self.SEC.size):
            interp_SEMF = interp(wvn, self.SEwn, self.SEMF[i], 0.0, 0.0)
            #PaoloA 01-04-2021
            jac[:, i] = interp_SEMF * jcb * w
            #jac[:, i] = interp_SEMF.dot(jcb * w)*interp_SEMF
        return jac

    def compute_from_model(self, SEC):
        """
        Compute Emissivity values from emissivity coefficients
        from solution vector
        """
        # load surface emissivity coefficients from solution vector
        self.SEvalues = dot(self.SEMF.T, SEC) + self.SEMFB
        self.SEvalues = exp(self.SEvalues)/(1.0+exp(self.SEvalues))
        self.SEC = SEC

#
# Unit test of the above class
#
if __name__ == '__main__':

    if __package__ is None:
        raise ImportError('The file "Emissivity.py" is embedded into '
                          'the dobjects package\nTo lauch it, use:\n'
                          '"python -m dobjects.Emissivity" '
                          'from the main directory of this project.')

    em = Emissivity('data/emissivitymodel.nc')
    em.test()
