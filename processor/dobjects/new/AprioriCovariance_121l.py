#!/usr/bin/env python
"""
Class to retrieve a priori Covariances.
"""
import numpy as np
from scipy.linalg import inv
from netCDF4 import Dataset

class AprioriCovariance(object):
    """
    This class is a wrapper around a netCDF data file
    """
    def __init__(self, datafile):
        """
        Initialize the AprioriCovariance object
        """
        self.df = Dataset(datafile, mode='r')

        try:
            self.obsnum = len(self.df.dimensions['obsnum'])
            self.single_apriori = False
        except KeyError:
            self.obsnum = 1
            self.single_apriori = True

        self.variables = self.df.variables
        self.varindx = (self.variables['varindx'][:]).astype('int')-1
        
    def __del__(self):
        self.df.close()

    def covariance_matrix(self, obs):
        """
        Get one of the covariance matrices and its inverse from file
        """
        if self.single_apriori:
            try:
                Sa = self.variables['selSa'][Ellipsis]
                SaInv = self.variables['selSaInv'][Ellipsis]
            except (KeyError, TypeError):
                Sa = self.variables['Sa'][Ellipsis]
                Sa = Sa[self.varindx, :]
                Sa = Sa[:, self.varindx]
                SaInv = self.variables['SaInv'][Ellipsis]
                SaInv = SaInv[self.varindx, :]
                SaInv = SaInv[:, self.varindx]
        else:
            try:
                Sa = self.variables['selSa'][obs, Ellipsis]
                SaInv = self.variables['selSaInv'][obs, Ellipsis]
            except (KeyError, TypeError):
                Sa = self.variables['Sa'][obs, Ellipsis]
                Sa = Sa[self.varindx, :]
                Sa = Sa[:, self.varindx]
                SaInv = self.variables['SaInv'][obs, Ellipsis]
                SaInv = SaInv[self.varindx, :]
                SaInv = SaInv[:, self.varindx]
        return [np.copy(Sa), np.copy(SaInv)]



class NewAprioriCovariance(object):
    EIGENVALUES_TCOV = 9
    #PaoloA 18082017
    #CLEAR_OUTSIDE_DIAGONAL=True
    CLEAR_OUTSIDE_DIAGONAL=False

    #PaoloA 16082017
    #LEVEL_SELECTION = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                       #17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                       #31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
                       #45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
                       #59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72,
                       #73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86,
                       #87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100,
                       #101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
                       #112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122,
                       #184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,
                       #195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205,
                       #206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216,
                       #217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227,
                       #228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238,
                       #239, 240, 241, 242, 243, 244, 245]
#    LEVEL_SELECTION = [1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16, 
#                       17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32, 
#                       33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48, 
#                       49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64, 
#                       65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80, 
#                       81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96, 
#                       97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 
#                       113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 
#                       129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 
#                       145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 
#                       161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 
#                       177, 178, 179, 180, 181, 182, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 
#                       284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 
#                       300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 
#                       316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 
#                       332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 
#                       348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363,
#                       364, 365] 
    LEVEL_SELECTION = [1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16, 
                        17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32, 
                        33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48, 
                        49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64, 
                        65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80, 
                        81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96, 
                        97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 
                       113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 
                       129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 
                       145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 
                       161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 
                       177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 
                       193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 
                       209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 
                       225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 
                       241, 242, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 
                       378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 
                       394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 
                       410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 
                       426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 
                       442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 
                       458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 
                       474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485]

    #PaoloA
    #def __init__(self, datafile, aemiss, eigen_land = 2, eigen_sea = 2):
    def __init__(self, datafile, aemiss, eigen_land = 5, eigen_sea = 2):
        """
        Initialize the AprioriCovariance object
        """
        self.df = Dataset(datafile, mode='r')
        self.f = self.df.groups['atmospheric_components'].groups['Covariances']

        self.T = self.f.variables['T']
        self.q = self.f.variables['q']
        self.Tq = self.f.variables['T_q']
        self.O3 = self.f.variables['O3']
        self.obsnum = self.T.shape[0]
        self.n_levels = self.T.shape[1]

        self.eigen_land = eigen_land
        self.eigen_sea = eigen_sea

        self.emiss_cov = aemiss.df.groups['surface_components'].variables['ModelCovariance']
        self.on_land = np.array(
            aemiss.df.groups['surface_components'].variables['LandOrWater'][:],
            dtype=np.bool,
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
            NewAprioriCovariance.LEVEL_SELECTION + list(range(486, 486 + eigen))
        )
        return lvl_selection - 1

    def __del__(self):
        self.df.close()

    def covariance_matrix(self, obs):
        eigen = self.eigenvalues(obs)
        size = (4 * self.n_levels) + 1 + eigen
        Sa = np.zeros((size , size), dtype=np.float64)
        SaInv = np.zeros((size, size), dtype=np.float64)

        Sa[self.n_levels*0: self.n_levels*1, self.n_levels*0: self.n_levels*1] = self.T[obs,:]
        Sa[self.n_levels*1: self.n_levels*2, self.n_levels*1: self.n_levels*2] = self.q[obs,:]
        Sa[self.n_levels*2: self.n_levels*3, self.n_levels*2: self.n_levels*3] = np.eye(self.n_levels, dtype=np.float64) * 16
        Sa[self.n_levels*3: self.n_levels*4, self.n_levels*3: self.n_levels*4] = self.O3[obs,:]
        Sa[self.n_levels*0: self.n_levels*1, self.n_levels*1: self.n_levels*2] = self.Tq[obs,:]
        Sa[self.n_levels*1: self.n_levels*2, self.n_levels*0: self.n_levels*1] = self.Tq[obs,:].T
        Sa[-eigen - 1, -eigen - 1] = NewAprioriCovariance.EIGENVALUES_TCOV
        for j in range(1, eigen + 1):
            Sa[-j, -j] = self.emiss_cov[obs, eigen - j]

        if NewAprioriCovariance.CLEAR_OUTSIDE_DIAGONAL:
            newSa = np.zeros((size, size), dtype=np.float64)
            newSa[np.diag_indices_from(newSa)] = np.diag(Sa)
            Sa = newSa
            SaInv[np.diag_indices_from(SaInv)] = 1 / np.diag(Sa)
        else:
            #SaInv[self.n_levels * 0: self.n_levels * 2, self.n_levels * 0: self.n_levels * 2] = inv(
            #    Sa[self.n_levels * 0: self.n_levels * 2, self.n_levels * 0: self.n_levels * 2])
            SaInv[self.n_levels * 0: self.n_levels * 1, self.n_levels * 0: self.n_levels * 1] = inv(
                Sa[self.n_levels * 0: self.n_levels * 1, self.n_levels * 0: self.n_levels * 1])
            SaInv[self.n_levels * 1: self.n_levels * 2, self.n_levels * 1: self.n_levels * 2] = inv(
                Sa[self.n_levels * 1: self.n_levels * 2, self.n_levels * 1: self.n_levels * 2])
            SaInv[self.n_levels * 2: self.n_levels * 3, self.n_levels * 2: self.n_levels * 3] = np.eye(self.n_levels,
                                                           dtype=np.float64) / 16.
            SaInv[self.n_levels * 3: self.n_levels * 4, self.n_levels * 3: self.n_levels * 4] = inv(
                Sa[self.n_levels * 3: self.n_levels * 4, self.n_levels * 3: self.n_levels * 4])
            SaInv[-eigen - 1, -eigen - 1] = 1. / NewAprioriCovariance.EIGENVALUES_TCOV
            for j in range(1, eigen + 1):
                SaInv[-j, -j] = 1. / self.emiss_cov[obs, eigen-j]

        if __debug__:
            test = np.dot(Sa, SaInv)
            test[np.abs(test) < 1e-4] = 0
            assert np.allclose(test, np.eye(size), rtol=1e-04,
                               atol=1e-04)

        selSa = Sa[self.varindx(obs), :][:, self.varindx(obs)]
        selSaInv = SaInv[self.varindx(obs), :][:, self.varindx(obs)]

        if __debug__:
            test = np.dot(selSa, selSaInv)
            test[np.abs(test) < 1e-4] = 0
            assert np.allclose(
                test,
                np.eye(len(NewAprioriCovariance.LEVEL_SELECTION) + eigen),
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

