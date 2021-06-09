from __future__ import division

from numpy import array, allclose, abs, sum, dot
from numpy.linalg import cholesky, LinAlgError

from netCDF4 import Dataset

#from scipy.sparse import csc_matrix
#from scikits.sparse.cholmod import cholesky as sparse_cholesky

class IasiSVD:
    def __init__(self, U, D, V):
        self.U = U
        self.D = D
        self.V = V

    def U_dot(self, M):
        return dot(self.U, M)

    def dot_U(self, M):
        return dot(M, self.U)

    def V_dot(self, M):
        return dot(self.V, M)

    def dot_V(self, M):
        return dot(M, self.V)

    def D_dot(self, M):
        return dot(self.D, M)

    def dot_D(self, M):
        return dot(M, self.D)


class IasiObservationError(object):
    """
    This class is a wrapper around a netCDF data file
    """
    def __init__(self, data_filename):
        """
        Initialize the ObservationError object
        """
        with Dataset(data_filename, mode='r') as df:
            self.__obs_err = array(df.variables['obs_err'][:])
            self.__inv_obs_err = array(df.variables['inv_obs_err'][:])

            U = array(df.variables['obs_err_U'][:])
            V = U.T
            D = array(df.variables['obs_err_D'][:])
            self.__svd = IasiSVD(U,D,V)
        
        # Check for some kind of properties of the obs_err matrix
        self.symmetric = False
        if allclose(self.__obs_err.T, self.__obs_err):
            self.symmetric = True

        self.sparse = False
        not_zeros = sum(abs(self.__obs_err)>1e-15)
        if not_zeros < self.__obs_err.shape[0] * self.__obs_err.shape[1]/100:
            self.sparse = True
        
        self.positive_definite = True
        try:
            cholesky(self.__obs_err)
        except LinAlgError:
            self.positive_definite = False
        
        # Now prepare the decomposition of the matrix accordlingly to the properties
        # if self.positive_definite and self.sparse and self.symmetric:
        #     self.__sparse_cholesky = sparse_cholesky(csc_matrix(self.__obs_err))
        

    @property
    def obs_err(self):
        """
        Get one of the observation error matrix
        """
        return self.__obs_err

    @property
    def inv_obs_err(self):
        return self.__inv_obs_err
    
    def inv_obs_err_dot(self, M):
        # if self.positive_definite and self.sparse and self.symmetric:
        if False:
            return self.__sparse_cholesky.solve_A(M)
        else:
            return dot(self.inv_obs_err, M)
    
    def dot_inv_obs_err(self, M):
        # if self.positive_definite and self.sparse and self.symmetric:
        if False:
            return self.__sparse_cholesky.solve_A(M.T).T
        else:
            return dot(M, self.inv_obs_err)

    @property
    def svd(self):
        return self.__svd
    
    def get_svd(self):
        return (self.__svd.U, self.__svd.D, self.__svd.V)

class IasiObservationErrorSps(object):
    """
    This class is a wrapper around a netCDF data file
    """
    def __init__(self, data_filename):
        """
        Initialize the ObservationError object
        """
        with Dataset(data_filename, mode='r') as df:
            self.__obs_err = array(df.variables['obs_err'][:])
            self.__inv_obs_err = array(df.variables['inv_obs_err'][:])
            self.__oe_sub_indices = array(df.variables['oe_sub_indices'][:])

            U = array(df.variables['obs_err_U'][:])
            V = U.T
            D = array(df.variables['obs_err_D'][:])
            self.__svd = IasiSVD(U,D,V)

        # Check for some kind of properties of the obs_err matrix
        self.symmetric = False
        if allclose(self.__obs_err.T, self.__obs_err):
            self.symmetric = True

        self.sparse = False
        not_zeros = sum(abs(self.__obs_err)>1e-15)
        if not_zeros < self.__obs_err.shape[0] * self.__obs_err.shape[1]/100:
            self.sparse = True

        self.positive_definite = True
        try:
            cholesky(self.__obs_err)
        except LinAlgError:
            self.positive_definite = False

        # Now prepare the decomposition of the matrix accordlingly to the properties
        # if self.positive_definite and self.sparse and self.symmetric:
        #     self.__sparse_cholesky = sparse_cholesky(csc_matrix(self.__obs_err))



    @ property
    def obs_err(self):
        """
        Get one of the observation error matrix
        """
        return self.__obs_err

    @property
    def inv_obs_err(self):
        return self.__inv_obs_err

    @property
    def oe_sub_indices(self):
        """
        Get the subselction indices relative to the obs_err channels 
        """
        return self.__oe_sub_indices-1

    def inv_obs_err_dot(self, M):
        # if self.positive_definite and self.sparse and self.symmetric:
        if False:
            return self.__sparse_cholesky.solve_A(M)
        else:
            return dot(self.inv_obs_err, M)

    def dot_inv_obs_err(self, M):
        # if self.positive_definite and self.sparse and self.symmetric:
        if False:
            return self.__sparse_cholesky.solve_A(M.T).T
        else:
            return dot(M, self.inv_obs_err)

    @property
    def svd(self):
        return self.__svd

    def get_svd(self):
        return (self.__svd.U, self.__svd.D, self.__svd.V)


