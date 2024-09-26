#!/usr/bin/env python

"""
Transform retrieved profiles in innovations for the data assimilation.
"""

import numpy as np
from scipy.linalg import svd


class DAresult(object):
    """ Empty object to store result """
    def __init__(self):
        self.K = None
        self.Yret_prime = None
        self.Hprime_rad = None
        self.S = None
        self.R = None
        self.Yret = None
        self.Hret_prime = None
        self.costf = None
        self.delta = None
        #PaoloA
        self.Lambda = None
        self.Sa_ret = None
        self.SaInv_ret = None

class DataStructure(object):
    def __init__(self, nlev):
        self.nlev = nlev


class transform(object):
    """
    Implements retirval transformation for DA
    """
    #Paolo Antonelli
    def __init__(self, obserr_object, nlev, conf_vars,approach = 2):
    #def __init__(self, obserr_object, nlev, approach = 1):
        """ Store the observation error and compute transform """
        self.obs_err = obserr_object
        _, self.Sigma, _ = obserr_object.get_svd()
        self.A = np.ascontiguousarray(1.0/np.sqrt(self.Sigma))
        self.A = self.A.reshape((self.Sigma.shape[0],1))
        self.data_structure = DataStructure(nlev)
        self.approach = approach
        
    def transform_retrievals_for_DA(self, Sa, SaInv, x0, xhat, K, yobs_minus_yhat, profile):
        """ Generate square root of inverse of Sa on VET grid """

        if profile.prlinalg is not None:
            profile.prlinalg.enable()

        #VET = DAresult()
        RET = DAresult()
        
        #PaoloA
        RET.Sa_ret=Sa
        RET.SaInv_ret=SaInv
        
        #We will use just the first 2*nlev values of Sa
        #Sa = Sa[0:2*nlev, 0:2*nlev]
        
        #indx = np.loadtxt(self.chan_selection_file).astype(int)
        nlev = self.data_structure.nlev
        nlev_tr = 81
        indx=np.concatenate((np.arange(0, nlev_tr), np.arange(nlev,nlev+nlev_tr)),  axis=0)

        Sa = Sa[:, indx]
        Sa = Sa[indx, :]
        #Sa[0:nlev_sps, 0:nlev_sps] = Sa[0:nlev_sps, 0:nlev_sps]
        #Sa[nlev_sps:2*nlev_sps] = Sa[nlev_sps:nlev_sps+nlev, nlev_sps:nlev_sps+nlev]

        # Generate square root 		and inverse square root of covariance matrix
        # using singular value decomposition to avoid complex numbers
        SaU, SaD, SaVT = svd(Sa)

        sqrtSaD = np.sqrt(SaD).reshape((SaD.shape[0],1))
        sqrtSaDInv = (1.0/sqrtSaD).reshape((SaD.shape[0],1))

        sqSa = np.dot(SaU, sqrtSaD * SaVT)
        sqSaInv = np.dot(SaU, sqrtSaDInv *SaVT)

        # Generate transofremd radiances according to Migliorini 2011 Eq. 3
        # Yrad = residuals + mat_mult(K,xhat)
        # Yprime_rad = mat_mult(self.A,mat_mult(self.LT,Yrad))


        spectral_indx = self.obs_err.oe_sub_indices

        K = (K[spectral_indx.astype('int64'),:])[:,indx]

        RET.Hprime_rad = self.A * self.obs_err.svd.V_dot(K)
        #print("RET.Hprime_rad shape {}".format(RET.Hprime_rad.shape))

        # Generate Signal to Noise on RET grid for the full state vector
        RET.S = np.dot(RET.Hprime_rad, sqSa)

        # Apply SVD to S on RET grid for the full state vector
        # Get only 20 eigenvectors
        RET.R = 20
        U, Lambda, V = svd(RET.S)

        Lambda = Lambda[0:RET.R]
        U = U[:,0:RET.R]
        V = V[0:RET.R, :]

        #PaoloA
        RET.Lambda = Lambda

        if self.approach == 1:
            RET.Hret_prime = np.dot(U.T, RET.Hprime_rad)
            # Here I reduce the dimension of xhat to make it consistent with K
            
            #Yrad_prime = self.A * self.obs_err.svd.V_dot(yobs_minus_yhat(spectral_indx) + np.dot(K,xhat[:2*nlev]))
            #Paolo Antonelli

            Yrad_prime = self.A.reshape(self.A.size,) * self.obs_err.svd.V_dot(yobs_minus_yhat[spectral_indx.astype('int64')] + np.dot(K,xhat[indx]))
            RET.Yret_prime = np.dot(U.T,Yrad_prime)

        elif self.approach == 2:
            temp = np.dot(V, sqSaInv)
            RET.Hret_prime = Lambda.reshape((Lambda.shape[0],1)) * temp
            diag_from_lambda = (Lambda/(1 + Lambda*Lambda)).reshape((Lambda.shape[0],1))
            diag_from_lambda_inv = 1/diag_from_lambda
            W = np.dot(sqSa ,np.dot(V.T , diag_from_lambda * U.T))
            #PaoloA
            #Yret = xhat[:2*nlev] - x0[:2*nlev] + np.dot(W, np.dot(RET.Hprime_rad, x0[:2*nlev]))
            Yret = xhat[indx] - x0[indx] + np.dot(W, np.dot(RET.Hprime_rad, x0[indx]))
            #PaoloA
            RET.Yret_prime = diag_from_lambda_inv.flatten() * np.dot(temp, Yret)
            #RET.Yret_prime = diag_from_lambda_inv * np.dot(temp, Yret)
        else:
            raise ValueError("Invalid approach for the transformer")

        if profile.prlinalg is not None:
            profile.prlinalg.disable()


        return RET
