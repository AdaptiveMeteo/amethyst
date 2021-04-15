#!/usr/bin/env python
#
# This script launches the software for the Amethyst Preprocessor.
#

"""
AMETHYST Preprocessor 

This script is an open source software written to launch the Amethyst preprocessor
taking  a wrf model as a first guess for Mirto, an open source software for
elaborating meteorological interferometer data

"""

import logging
from argparse         import ArgumentParser
from sys              import exit as sysexit
from traceback        import format_exc
from amethyst_config  import preprocessor_vars, common_vars

__author__     = 'Paolo Scaccia <paolo.scaccia@adaptivemeteo.com>'
__copyright__  = "Copyright 2021, Adaptive Meteo S.r.l"
__credits__    = ["Paolo Scaccia", "Paolo Antonelli"]
__license__    = "GPL"
__version__    = "1.0"
__maintainer__ = "Paolo Scaccia"
__email__      = "paolo.scaccia@adaptivemeteo.com"

if __name__ == '__main__':
    log = logging.getLogger()
else:
    log = logging.getLogger(__name__)

def launch_preprocessing(argv):
    """
        Preprocessor Launcher written to be called from shell and from other python scripts.        
    ----------
    argv : class ArgumentParser
        Class containing .

    Returns
    -------
    None.

    """
    
    if argv.instrument == 'cris':
        
        
        
    else:
        
    
    
    
    return

def main():
    
    v_levels = ['debug', 'info', 'warning']
    
    parser = ArgumentParser()
    parser.add_argument('--output','-o',type=str,required=True,
                        help='Output directory containing all NETCDF preprocessed data')
    parser.add_argument('--overpass',default = None,type=str,
                        help='Overpass date and time YYYYMMDD_HHMMSS')
    parser.add_argument('--wrf', '-w', type=str, default=common_vars['wrfdir'],
                        help='Where are the Wrf model data.')
    parser.add_argument('--gcrso', type=str, default=preprocessor_vars['gcrso'],
                        help='')
    parser.add_argument('--scris', type=str, default=preprocessor_vars['scris'],
                        help='')
    parser.add_argument('--iasi_native', type=str, default=preprocessor_vars['iasinat'],
                        help='')
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='The level of verbosity of the software')
    parser.add_argument('--levels', '-l', type=int, default=common_vars['levels'],
                        help='The total number of levels')
    parser.add_argument('--instrument', '-m', type=str, required=True,choices=['cris','info'],
                        help='Instrument type')
    parser.add_argument('--latmin', type=str, default = preprocessor_vars['geobox']['latmin'],
                        help='Min Latitude')
    parser.add_argument('--latmax', type=str, default = preprocessor_vars['geobox']['latmax'],
                        help='Max Latitude')
    parser.add_argument('--lonmin', type=str, default = preprocessor_vars['geobox']['lonmin'],
                        help='Min Longitude')
    parser.add_argument('--lonmax', type=str, default = preprocessor_vars['geobox']['lonmax'],
                        help='Max Longitude')

    argv = parser.parse_args()
    
    # Call launcher
    launch_preprocessing(argv)

    return
if __name__ == '__main__':
    sysexit(main())