#!/usr/bin/env python
"""
Class to retrieve a priori Covariances.
"""
import numpy as np
from scipy.linalg import inv
from netCDF4 import Dataset
import amethyst_config

class AprioriCovariance(object):
    #PaoloA 29-03-2021 
    #EIGENVALUES_TCOV = 4
    EIGENVALUES_TCOV = amethyst_config.processor_vars['eigenfortcov']
    #PaoloA 18082017
    #CLEAR_OUTSIDE_DIAGONAL=True
    CLEAR_OUTSIDE_DIAGONAL=amethyst_config.processor_vars['clear_outside_diag']
    #PaoloA 31-03-2021
    EIGENVALUES_LAND = amethyst_config.processor_vars['eigenforland']
    EIGENVALUES_SEA  = amethyst_config.processor_vars['eigenforsea']
    CO2_STD          = amethyst_config.processor_vars['constant_co2_std']['value']
    #PaoloA
    #def __init__(self, datafile, aemiss, eigen_land = 2, eigen_sea = 2):
    #PaoloA 29-03-2021
    #def __init__(self, datafile, aemiss, eigen_land = 5, eigen_sea = 5):
    #PaoloA 31-03-2021
    def __init__(self, datafile, aemiss, eigen_land = EIGENVALUES_LAND, eigen_sea = EIGENVALUES_SEA, var_selection = None, co2_std = CO2_STD):
        """
        Initialize the AprioriCovariance object
        """
        self.df = Dataset(datafile, mode='r')
        self.f = self.df.groups['atmospheric_components'].groups['Covariances']

        self.T = self.f.variables['T']
        self.q = self.f.variables['q']
        self.Tq = self.f.variables['T_q']
        self.O3 = self.f.variables['O3']
        self.co2_std  = co2_std
        self.obsnum   = self.T.shape[0]
        self.n_levels = self.T.shape[1]

        # Use all variables if there is no selection
        if var_selection == None:
            var_selection = [ True for x in range(5)]
        
        # Build Level Selection
        self.LEVEL_SELECTION = [ x+1 for x in range(len(var_selection[1:])*self.n_levels) if var_selection[1:][x//self.n_levels] ]
        if var_selection[0]:
            self.LEVEL_SELECTION.append( self.LEVEL_SELECTION[-1] + 1)           

        self.eigen_land = eigen_land
        self.eigen_sea = eigen_sea

        self.emiss_cov = aemiss.df.groups['surface_components'].variables['ModelCovariance']
        self.on_land = np.array(
            aemiss.df.groups['surface_components'].variables['LandOrWater'][:],
            dtype=np.bool_,
        )

        self.single_apriori = False

    def eigenvalues(self, obs):
        if self.on_land[obs]:
            return self.eigen_land
        else:
            return self.eigen_sea

    def varindx(self, obs):
        eigen = self.eigenvalues(obs)
        lvl_selection = np.array(
            ##PaoloA 16082017
            #NewAprioriCovariance.LEVEL_SELECTION + list(range(246, 246 + eigen))
            #NewAprioriCovariance.LEVEL_SELECTION + list(range(366, 366 + eigen))
            #NewAprioriCovariance.LEVEL_SELECTION + list(range(486, 486 + eigen))
            self.LEVEL_SELECTION + list(range(326, 326 + eigen))
        )
        return lvl_selection - 1

    def __del__(self):
        self.df.close()

    def covariance_matrix(self, obs, W=None):
        """Return (Sa, SaInv) for observation obs.

        Parameters
        ----------
        W : ndarray (n_levels, n_levels), optional
            Log-pressure interpolation matrix for regridding atmospheric
            covariance blocks.  When provided, each block B is replaced by
            W @ B @ W.T before assembly.  Pass None (default) for no regrid.
        """
        eigen = self.eigenvalues(obs)
        size = (4 * self.n_levels) + 1 + eigen

        Sa = np.zeros((size , size), dtype=np.float64)
        SaInv = np.zeros((size, size), dtype=np.float64)

        def _rg(block):
            """Apply regridding if W was supplied."""
            b = np.asarray(block, dtype=np.float64)
            return W @ b @ W.T if W is not None else b

        n = self.n_levels
        Sa[n*0: n*1, n*0: n*1] = _rg(self.T[obs,:])
        Sa[n*1: n*2, n*1: n*2] = _rg(self.q[obs,:])
        Sa[n*2: n*3, n*2: n*3] = np.eye(n, dtype=np.float64) * self.co2_std
        # PaoloS 12/04/2024: added perturbation of order 1e-6 to O3 covariance matrix
        Sa[n*3: n*4, n*3: n*4] = _rg(self.O3[obs,:] + np.eye(n)*1e-6)

        Sa[n*0: n*1, n*1: n*2] = _rg(self.Tq[obs,:])
        Sa[n*1: n*2, n*0: n*1] = Sa[n*0: n*1, n*1: n*2].T
        Sa[-eigen - 1, -eigen - 1] = AprioriCovariance.EIGENVALUES_TCOV
        for j in range(1, eigen + 1):
            Sa[-j, -j] = self.emiss_cov[obs, eigen - j]

        if AprioriCovariance.CLEAR_OUTSIDE_DIAGONAL:
            newSa = np.zeros((size, size), dtype=np.float64)
            newSa[np.diag_indices_from(newSa)] = np.diag(Sa)
            Sa = newSa
            SaInv[np.diag_indices_from(SaInv)] = 1 / np.diag(Sa)
        else:
            #SaInv[self.n_levels * 0: self.n_levels * 2, self.n_levels * 0: self.n_levels * 2] = inv(
            #    Sa[self.n_levels * 0: self.n_levels * 2, self.n_levels * 0: self.n_levels * 2])
            #PaoloA 29-03-2021
            #SaInv[self.n_levels * 0: self.n_levels * 1, self.n_levels * 0: self.n_levels * 1] = inv(
            #    Sa[self.n_levels * 0: self.n_levels * 1, self.n_levels * 0: self.n_levels * 1])
            #SaInv[self.n_levels * 1: self.n_levels * 2, self.n_levels * 1: self.n_levels * 2] = inv(
            #    Sa[self.n_levels * 1: self.n_levels * 2, self.n_levels * 1: self.n_levels * 2])
            #SaInv[self.n_levels * 2: self.n_levels * 3, self.n_levels * 2: self.n_levels * 3] = np.eye(self.n_levels,
            #                                               dtype=np.float64) / 16.
            #SaInv[self.n_levels * 3: self.n_levels * 4, self.n_levels * 3: self.n_levels * 4] = inv(
            #    Sa[self.n_levels * 3: self.n_levels * 4, self.n_levels * 3: self.n_levels * 4])
            #SaInv[-eigen - 1, -eigen - 1] = 1. / NewAprioriCovariance.EIGENVALUES_TCOV
            #for j in range(1, eigen + 1):
            #    SaInv[-j, -j] = 1. / self.emiss_cov[obs, eigen-j]
            SaInv[: self.n_levels*2, : self.n_levels*2] = np.linalg.inv(Sa[: self.n_levels*2, : self.n_levels*2])
            #SaInv[: self.n_levels*2, : self.n_levels*2] = np.linalg.pinv(SaInv[: self.n_levels*2, : self.n_levels*2],rcond=1e-15)
            SaInv[self.n_levels*2: self.n_levels*3, self.n_levels*2: self.n_levels*3] = np.eye(self.n_levels, dtype=np.float64) * (1/self.co2_std)
            #SaInv[self.n_levels*3: self.n_levels*4, self.n_levels*3: self.n_levels*4] = np.linalg.inv(self.O3[obs,:])   
            SaInv[self.n_levels*3: self.n_levels*4, self.n_levels*3: self.n_levels*4] = np.linalg.inv( Sa[self.n_levels*3: self.n_levels*4, self.n_levels*3: self.n_levels*4] )
            SaInv[-eigen - 1, -eigen - 1] = 1/AprioriCovariance.EIGENVALUES_TCOV
            for j in range(1, eigen + 1):
                SaInv[-j, -j] = 1/self.emiss_cov[obs, eigen - j]

            # SaInv=np.linalg.pinv(Sa,rcond=1e-15)

        if __debug__:
            #np.savetxt('/home/mirto/amethyst_test_SaInv.txt',SaInv)
            #np.savetxt('/home/mirto/amethyst_test_Sa.txt',Sa)
            test = np.dot(Sa, SaInv)
            test[np.abs(test) < 1e-4] = 0
            # Test PaoloS: 5/7/24
            #assert np.allclose(test, np.eye(size), rtol=1e-04,
            #                  atol=1e-04)

        selSa = Sa[self.varindx(obs), :][:, self.varindx(obs)]
        #PaoloA 29-03-2021
        selSaInv = SaInv[self.varindx(obs), :][:, self.varindx(obs)]
        #selSaInv = np.linalg.pinv(selSa,rcond=1e-15)

        if __debug__:
            test = np.dot(selSa, selSaInv)
            test[np.abs(test) < 1e-4] = 0
            assert np.allclose(
                test,
                np.eye(len(self.LEVEL_SELECTION) + eigen),
                rtol=2e-04,
                atol=2e-04
            )

        return selSa, selSaInv



class covariance_matrix(object):
    """
    This class is a wrapper around an AprioriCovariance object.
    It is used to read just an observation and keep it in memory
    until the object is deleted.
    """
    def __init__(self, apriori, obs):
        """
        Initialize the object using an AprioriCovariance and the obs index
        """
        self.varindx = apriori.varindx(obs)
        [self.Sa, self.SaInv] = apriori.covariance_matrix(obs)
        self.obs = obs

    def covariance_matrix(self, obs):
        """
        Return the in-memory object
        """
        if obs != self.obs:
            raise IndexError('Object initilized with obs='+repr(self.obs)+
                             ' but obs='+repr(obs)+' requested !')
        return [self.Sa, self.SaInv]

