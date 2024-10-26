from numpy           import diag, identity, array, allclose, dot, abs,sum, in1d
from netCDF4         import Dataset
from numpy.linalg    import cholesky, LinAlgError
from amethyst_config import processor_vars


def create_obs_err(data_filename, indx_file = None):
    if indx_file is None:
        return ObservationError(data_filename)
    else:
        from numpy import loadtxt
        return ObservationError(data_filename, indx = loadtxt(indx_file,dtype=int))

def create_obs_err_tr(data_filename, indx_file, indx_tr_file, obs_err_type='cris'):
    from numpy import loadtxt
    return ObservationErrorTR(data_filename, loadtxt(indx_file,dtype=int), loadtxt(indx_tr_file,dtype=int))

class SVD:
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

class ObservationError(object):
    """
    This class is a wrapper around a netCDF data file
    """
    def __init__(self, data_filename, indx = None):
        """
        Initialize the ObservationError object
        """
        with Dataset(data_filename, mode='r') as df:
            if indx is None:
                 self.__obs_err = array(df.variables['obs_err'][:])
                 self.__inv_obs_err = array(df.variables['inv_obs_err'][:])
                 U = array(df.variables['obs_err_U'][:])
                 D = array(df.variables['obs_err_D'][:])
            else:
                 self.__obs_err = array(df.variables['obs_err'][:])[indx[:,None],indx]
                 self.__inv_obs_err = array(df.variables['inv_obs_err'][:])[indx[:,None],indx]
                 U = array(df.variables['obs_err_U'][:])[indx[:,None],indx]
                 D = array(df.variables['obs_err_D'][indx])
            V = U.T
            self.__svd = SVD(U,D,V)
        
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

class ObservationErrorTR(object):
    """
    This class is a wrapper around a netCDF data file
    """
    def __init__(self, data_filename, obserr_indx, obserr_indx_tr):
        """
        Initialize the ObservationError object
        """
        if in1d(obserr_indx_tr, obserr_indx).size != obserr_indx_tr.size:
            raise IndexError('TR indices are not a sub selection of the instrument channel list !')

        
        with Dataset(data_filename, mode='r') as df:
            self.__obs_err = array(df.variables['obs_err'][obserr_indx_tr,obserr_indx_tr])
            self.__inv_obs_err = array(df.variables['inv_obs_err'][obserr_indx_tr,obserr_indx_tr])
            self.__oe_sub_indices = array([  i for i,x in enumerate(obserr_indx)  if x in obserr_indx_tr ])
            U = array(df.variables['obs_err_U'][obserr_indx_tr,obserr_indx_tr])
            D = array(df.variables['obs_err_D'][obserr_indx_tr])

            V = U.T
            self.__svd = SVD(U,D,V)

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
