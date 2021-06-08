#!/usr/bin/env python
from __future__ import print_function

"""
This class export an object which contains data from a netCDF file.
File expected content is:

      float cWvn(nchan) ;
      short molidfix(nmolfix) ;
      short molid(nmol) ;
      float pref(nlev) ;
      float tmptab(nlayod, ntmpod) ;
      float wvptab(nlayod, nmoltab) ;
      float coef(nf, nchmax) ;
        coef:_FillValue = -9999.f ;
      int ichmap(nf, nchmax) ;
        ichmap:_FillValue = -1 ;
      short nch(nf) ;
      int isels(nf) ;
      double vwvn(nf) ;
      short imols(nf, mxmols) ;
        imols:_FillValue = -1s ;
      float kfix(nf, ntmpod, nlayod) ;
        kfix:_FillValue = -9999.f ;
      float dkh2o(nf, ntmpod, nlayod) ;
        dkh2o:_FillValue = -9999.f ;
      float kh2o(nf, ntmpod, nlayod) ;
        kh2o:_FillValue = -9999.f ;
      float kvar(nf, ntmpod, nlayod, mxmols) ;
        kvar:_FillValue = -9999.f ;

"""


from numpy import array, shape, empty
from netCDF4 import Dataset

class Hitran(object):
    """Loads hitran pre-computed tables in memory"""

    def __init__(self, datafile):
        """Opens the data file and get basic header informations"""
        data_file = Dataset(datafile, mode='r')

        self.__nchan  = len(data_file.dimensions['nchan'])
        self.__nlayod = len(data_file.dimensions['nlayod'])
        self.__ntmpod = len(data_file.dimensions['ntmpod'])
        self.__nf     = len(data_file.dimensions['nf'])
        self.__nmol   = len(data_file.dimensions['nmol'])
        self.__mxmols = len(data_file.dimensions['nmol'])
        self.__nmoltab= len(data_file.dimensions['nmoltab'])

        self.__data = dict()
        for name in data_file.variables:
            self.__data[name] = array(data_file.variables[name][:])

        data_file.close()


    def get(self, name, part=None):
        """Reads one of the file variable or a subset of it"""
        if part is None:
            return self.__data[name]
        else:
            output = empty([len(part[0]), len(part[1])], dtype=int)
            for i, ind_i in enumerate(part[0]):
                for j, ind_j in enumerate(part[1]):
                    output[i, j] = self.__data[name][ind_i, ind_j]
            return output

    def nchan(self):
        """Reads nchan dimension"""
        return self.__nchan

    def nlayod(self):
        """Reads nlayod dimension"""
        return self.__nlayod

    def ntmpod(self):
        """Reads ntmpod dimension"""
        return self.__ntmpod

    def nf(self):
        """Reads nf dimension"""
        return self.__nf

    def nmol(self):
        """Reads nmol dimension"""
        return self.__nmol

    def mxmols(self):
        """Reads mxmols dimension"""
        return self.__mxmols

    def nmoltab(self):
        """Reads nmoltab dimension"""
        return self.__nmoltab

    def test(self):
        """Test of the class"""
        inhmap = self.get('ichmap', [range(0, 20), range(0, 10)])
        print(shape(inhmap))
        inhmap = self.get('ichmap')
        print(shape(inhmap))
        print(self.nchan())

#
# Unit test of the above class
#
if __name__ == '__main__':
    import sys
    import time

    if __package__ is None:
        raise ImportError('The file "Hitran.py" is embedded into '
                          'the dobjects package\nTo lauch it, use:\n'
                          '"python -m dobjects.Hitran" '
                          'from the main directory of this project.')

    if sys.argv != 2:
       file_to_open = 'data/leo.iasi.0.05.nc'
    else:
       file_to_open = sys.argv[1]
    pre_read = time.time()
    hit = Hitran(file_to_open)
    print('File opened in {:.3f} seconds'.format(time.time()-pre_read))
    hit.test()
