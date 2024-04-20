#!/usr/bin/env python

"""
Compute the radiance vector F using the selected forward model

x is a state vector containing the variables for which the forward
model should be evaluated

Note: The state vector unit for concentration is the natural log or vmr in ppv.
xdim is a vector of dimensions for each of parameters of the solution
(T, wv, ozone, surface)
"""

# pylint: disable=C0103
# pylint: disable=C0325
# pylint: disable=E0611

import numpy as np


class NotConvergentIteration(Exception):
    """
    This error is raised when it is required to compute an observation
    on which the algorithm does not converge.
    """
    pass


class ForwardModel(object):
    """
    Forward Model Class
    """
    def __init__(self, control):
        """
        Initializes the class. Any relevanty parameter is in the
        control object.
        """
        #  fm.wnF = calculation wavenumber vector
        #  fm.F   = calculation radiance vector
        self.cx = control
        self.xdim = control.xdim
        self.outdata = {}
        self.model = control.model
        self.emiss = control.emiss
        self.wnF = self.model.cwvn
        self.F = np.zeros(len(self.model.cwvn))
        self.Md = 28.966  # Molecular mass of dry air
        self.Mw = 18.016  # Molecular Mass of water
        self.Mc = 44.01   # Molecular Mass of CO2
        self.Mo = 48      # Molecular Mass of Ozone
        self.tds = 0
        self.nlev = self.xdim[0]
        self.tde = self.tds+self.xdim[0]
        self.wvs = self.tde
        self.wve = self.wvs+self.xdim[1]
        self.cos = self.wve
        self.coe = self.cos+self.xdim[2]
        self.ozs = self.coe
        self.oze = self.coe+self.xdim[3]
        self.sks = self.oze
        self.ske = self.sks+self.xdim[4]
        self.ems = self.ske
        self.eme = self.ems+self.xdim[5]
        self.wnK = None
        self.K = None
        self.F = None
    
    def compute_forward(self, xhat):
        """
        Define Input to radiance calculator
          Atmospheric profile data:
             pressure (Nx1) level pressure (mbar)
             tdry     (Nx1) level temperature (K)
             w        (Nx1) level water vapor mass mixing ratio (g/kg)
             vmr      (Nx1) level water vapor volume mixing ratio (ppmv)
             alt      (Nx1) level altitudes (km)
             oz       (Nx1) level ozone (ppmv)
        NOTE: when using pressure boundaries the surface elevation is required
        in the lowest altitude level for use in the hypsometric equation.
        All other altitudes are ignored in the profile level input.
        Compute the jacobian K using the selected model
        """
        # convert from log(vmr in ppv) to vmr in ppmv

        try:
            vmr = np.exp(xhat[self.wvs:self.wve])
        except FloatingPointError:
            raise NotConvergentIteration
        co2_ppmv = xhat[self.cos:self.coe]
        oz_vmr = np.exp(xhat[self.ozs:self.oze])

        indata = {}
        xdim = self.xdim
        Jvar = self.cx.Jvar
        ii = self.cx.indx
        indata['sfgrd'] = self.cx.sfgrd
        indata['emrf']  = np.column_stack([self.cx.emrf, 1 - self.cx.emrf])

        # State vector temeprature is in [K]
        indata['tskin']     = xhat[self.sks:self.ske]
        indata['psf']       =  self.cx.surfacePressure_mb
        # State vector Temperature is in [K]
        indata['temp']      = np.flipud(xhat[self.tds:self.tde])
        # State vector Water vapor is in log(q) where q is in [Kg/Kg]
        indata['h2o']       = np.flipud(vmr)
        # State vector CO2 is in ppmv
        indata['co2']       = np.flipud(1.E-6*(self.Mc/self.Md)*co2_ppmv)
        # State vector Ozone is in log(q) where q is in [Kg/Kg]
        indata['o3']        = np.flipud(oz_vmr)
        indata['pressure']  = np.flipud(self.cx.pressure_grid)
        indata['pobs']      = self.cx.observationPressure_mb
        indata['obsang']    = self.cx.FOVangle
        indata['sunang']    = self.cx.SunAngle
        indata['solzenith'] = self.cx.Solar_zenith_angle
        indata['azangle']   = self.cx.Solar_azimuth_angle
        indata['obslevel']  = self.cx.oss_obslevel
        indata['lat']       = self.cx.fov_latitude

        # Debug
        for k,v in indata.items():
              print(k, v)
              print()

        # Call the selected forward model (OSS)
        self.model.compute(indata, self.outdata)
        # Subselect channels which are used in the inversion
        self.F = self.outdata['y'][ii]
        print('F in ForwardModel.py')
        print(self.F)
        SEflag  = np.size(np.where(Jvar == -2)) > 0
        SKTflag = np.size(np.where(Jvar == -1)) > 0
        Tflag   = np.size(np.where(Jvar == 0)) > 0
        WVflag  = np.size(np.where(Jvar == 1)) > 0
        CO2flag = np.size(np.where(Jvar == 2)) > 0
        O3flag  = np.size(np.where(Jvar == 3)) > 0

        # Set number of channels
        krow = len(ii)

        # Inizialize the K matrix with nans
        jac = np.zeros((krow, sum(xdim[0:6])))
        jac.fill(np.NAN)
        nlev = self.nlev

        if Tflag:
            # Temperature Jacobians are in [K]
            jcb = np.transpose(self.outdata['xkt'][0:nlev, ii])
            jac[:, self.tds:self.tde] = np.fliplr(jcb)
        if WVflag:
            # Convert state vector water vapor from log(vmr)
            # [in log(kg/kg)] to vmr in [kg/kg]
            vmr = np.exp(xhat[self.wvs:self.wve])
            # Go from dR/dq in [Kg/Kg] to dR/dlog(q) by multiplying
            # (dR/dq)*q as dlog(q)=dq/q ...
            w_mat = np.tile(vmr, (krow, 1))
            jcb = np.transpose(self.outdata['xkt'][nlev+2:nlev+2+nlev, ii])
            # Jacobians in log(q)
            jac[:, self.wvs:self.wve] = np.fliplr(jcb)*w_mat
        if CO2flag:
            jcb = np.transpose(self.outdata['xkt'][2*nlev+2:2*nlev+2+nlev, ii])
            # Convert Jacobians in [kg/kg] into jacobians in [ppmv]
            jac[:, self.cos:self.coe] = np.fliplr(jcb)*(self.Mc/self.Md)*1.e-6

        if O3flag:
            # Convert state vector ozone from log(vmr)
            # [in log(kg/kg)] to vmr in [kg/kg]
            vmr = np.exp(xhat[self.ozs:self.oze])
            # Go from dR/dq in [Kg/Kg] to dR/dlog(q) by multiplying
            # (dR/dq)*q as dlog(q)=dq/q ...
            w_mat = np.tile(vmr, (krow, 1))
            jcb = np.transpose(self.outdata['xkt'][3*nlev+2:3*nlev+2+nlev, ii])
            # jacobians in log(q)
            jac[:, self.ozs:self.oze] = np.fliplr(jcb)*w_mat

        if SKTflag:
            # Surface temperature Jacobians in [K]
            jcb = np.transpose(self.outdata['xkt'][nlev:nlev+1, ii])
            jac[:, self.sks:self.ske] = jcb
            #print('SKT\n',jac[:, self.sks:self.ske],jac[:, self.sks:self.ske].shape)
        if SEflag:
            # Surface emissivity Jacobians
            jcb = np.transpose(self.outdata['paxkemrf'][0, :])
            jcb = (self.emiss.get_jacobian(self.wnF, jcb))
            jac[:, self.ems:self.eme] = jcb[ii]
            self.cx.Kse = jcb
            #print('SE\n',jac[:, self.ems:self.eme],jac[:, self.ems:self.eme].shape)
        self.K = np.ascontiguousarray(jac)
        self.wnK = self.wnF[ii]

    def compute_residuals(self):
        """
        Compute the obs minus calc difference, yobs_minus_yhat
        Input :
          R     = the observed radiance
          F     = the calculated radiance
        Output:
          out = radiance difference vector yobs_minus_yhat
        """
        yobs_minus_yhat = self.cx.R - self.F
        return yobs_minus_yhat

#
# Unit test of the above
#
if __name__ == '__main__':

    if __package__ is None:
        raise ImportError('The file "ForwardModel.py" is embedded into '
                          'the mirto package\nTo lauch it, use:\n'
                          '"python -m mirto.ForwardModel" '
                          'from the main directory of this project.')

    from os import path
    from amethyst.ParseConfig import ParseConfig
    from dobjects.FirstGuess import FirstGuess
    from dobjects.SounderFOV import SounderFOV
    from dobjects.Solar import Solar
    from dobjects.Hitran import Hitran
    #PaoloA 14082020
    #from dobjects.AprioriCovariance import AprioriCovariance
    from dobjects.AprioriCovariance import NewAprioriCovariance
    from dobjects.ObservationError import ObservationError
    #PaoloA 01-04-2021
    #from dobjects.Emissivity import Emissivity
    from dobjects.Emissivity import NewEmissivity
    from ossFM_files.ossFM import ossFM
    #
    # OSS init input
    #
    class mirto_config(object):
        """ Helper class to mimic configuration in code_main """
        def __init__(self, model, xmlconf, emiss, fg, fov):
            self.model = model
            self.surfacePressure_mb = xmlconf.constantPressure
            self.observationAltitude_km = xmlconf.geometryObservationAltitude
            self.observationPressure_mb = xmlconf.geometryObservationPressure
            self.EstimateK = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                                      dtype=int)
            self.Jvar = np.array(xmlconf.retrievalSelectedStateVectorVariables)
            self.Iteration_limit = sum(self.EstimateK)
            self.xdim = fg.xdim
            self.indx = fov.indx
            self.emiss = emiss
            [self.pressure_grid, self.xhat, self.a] = fg.state_vector(0)
            ems = (self.xdim[0]+self.xdim[1]+self.xdim[2]+
                   self.xdim[3]+self.xdim[4])
            eme = (self.xdim[0]+self.xdim[1]+self.xdim[2]+
                   self.xdim[3]+self.xdim[4]+self.xdim[5])
            [self.sfgrd, self.emrf] = aemiss.get(0, self.xhat[ems:eme])
            [self.fov_latitude, self.fov_longitude, self.fov_time,
             self.FOVangle, self.wnR, self.R] = fov.fov(0)
            self.SunAngle = 90.0
            self.surfacePressure_mb = self.pressure_grid[0]
            self.observationPressure_mb = self.pressure_grid[-1]

    datapath = 'data'
    asolar = Solar(path.join(datapath, 'solar_irradiances.nc'))
    ahitran = Hitran(path.join(datapath, 'leo.iasi.0.05.nc'))
    aemiss = Emissivity(path.join(datapath, 'emissivitymodel.nc'))

    oss = ossFM(asolar, ahitran)
    axmlconf = ParseConfig('configuration.xml')

    afov = SounderFOV(path.join(datapath, 'fov.nc'))
    afg = FirstGuess(path.join(datapath, 'fg.nc'))
    #PaoloA 14082020
    #apriori = AprioriCovariance(path.join(datapath, 'apriori.nc'))
    apriori = NewAprioriCovariance(path.join(datapath, 'apriori.nc'))
    obserr = ObservationError(path.join(datapath, 'obserr.nc'))
    aemiss = Emissivity(path.join(datapath, 'emissivitymodel.nc'))

    cx = mirto_config(oss, axmlconf, aemiss, afg, afov)
    rad = ForwardModel(cx)
    rad.compute_forward(cx.xhat)
    try:
        from netCDF4 import Dataset
        rootgrp = Dataset('radiance_FM.nc', 'w', format='NETCDF4')
    except ImportError:
        from scipy.io import netcdf
        rootgrp = netcdf.netcdf_file('radiance_FM.nc', mode='w')
    nchan = rootgrp.createDimension('nchan', len(rad.wnK))
    ncwn = rootgrp.createVariable('wn', 'f4', ('nchan',))
    ncra = rootgrp.createVariable('radiance', 'f4', ('nchan',))
    ncra.coordinates = 'wn'
    ncwn[:] = rad.wnK
    ncra[:] = rad.F
    rootgrp.close()
