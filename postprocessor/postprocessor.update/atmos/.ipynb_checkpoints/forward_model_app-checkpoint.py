#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 14:00:35 2021

@author: paolo
"""
import numpy as np
import json
from os.path import isfile
from os import system, path
import sys

def read_JSON_scene(file, fm_version):
    
    indoss = json.load(open(file, 'r'))
    for key in indoss.keys():
        if type(indoss[key]) == list:
            indoss[key] = np.array(indoss[key])
    if fm_version == 2:     
        indoss['emrf'] = indoss['emrf'][0,:]
    return [ indoss ]

def read_fg(fg, fm_version, emissivity_mode = 'default'):
    from common.data_reader.mirto_wrapper import MIRTO_fov_Wrapper, MIRTO_fg_Wrapper
    
    indoss = {}
    fg_data = MIRTO_fg_Wrapper(fg)
    fg_data.to_standard()
    if not isfile(fg.replace('fg.nc','fov.nc')):
        sys.exit("Missing fov.nc file in {}".format(fg.replace('fg.nc','')))
    fov_data = MIRTO_fov_Wrapper(fg.replace('fg.nc','fov.nc'))
    fov_data.to_standard()
    
    nfovs = fg_data.number_of_FOVs
    indata = []
    
    for ifov in range(nfovs):
        indoss = {}
        indoss['sfgrd'] = fg_data.surface_components_ModelWaveNumbers
    
        # Compute Surface Emissivity Coefficient
        if emissivity_mode == 'default':
            coef = np.zeros( fg_data.surface_components_ModelFunctions.shape[1]  )
        elif emissivity_mode == 'randomize':
            coef = np.array([ np.random.normal(0,np.sqrt(eig)) 
                             for eig in fg_data.surface_components_ModelCovariance[ifov] ] )
    
        tmp  = fg_data.surface_components_ModelFunctions[ifov].T.dot(coef) + \
                          fg_data.surface_components_ModelFunctionsBias[ifov]
        emrf = np.exp(tmp)/(1.0+np.exp(tmp))
        indoss['emrf']              = np.column_stack([ emrf, 1 - emrf])
        indoss['tskin']             = fg_data.atmospheric_components_skT[ifov]
        indoss['temp']              = np.flipud(fg_data.air_temperature[ifov]) 
        indoss['h2o']               = 1e-3*(np.flipud( fg_data.air_water_vapor_mr[ifov] ) ) # kg/kg
        indoss['o3']                = np.exp(np.flipud( fg_data.atmospheric_components_O3[ifov]) ) # kg/kg
        indoss['pressure']          = np.flipud(fg_data.air_pressure[ifov])
        indoss['constant_co2']      = 1.E-6*(44.01 /28.966)*405  # Costant 405 ppmv CO2 
        indoss['obslevel']          = -1
        indoss['obsang']            = fov_data.fov_angle[ifov]
        indoss['solzenith']         = fov_data.solar_zenith_angle[ifov]
        indoss['azangle']           = fov_data.solar_azimuth_angle[ifov]
        indoss['lat']               = fg_data.latitude[ifov]
        indoss['psf']               = fg_data.atmospheric_components_sp[ifov]
        indoss['pobs']              = indoss['pressure'][0]
        
        if fm_version == 2:
            indoss['sunang']          = 90.0
            indoss['emrf'] = emrf
        
        """
        if fm_version == 3:
             
        elif fm_version == 2:
             sys.exit("Not implemented yet!")
        else: 
             sys.exit("Wrong FM version!")
        """   
        indata.append(indoss)
    return indata
   
    
def save_output(fm, outdata, outfile, nl, nspe, version):
        
    # Delete old .nc output
    if isfile(outfile):
        system('rm {}'.format(outfile))
    try:
        from netCDF4 import Dataset
        rootgrp = Dataset(outfile, 'w', format='NETCDF4')
    except ImportError:
        from scipy.io import netcdf
        rootgrp = netcdf.netcdf_file(outfile, mode='w')
    
    nchand = rootgrp.createDimension('channels', fm.nchan)
    nsped = rootgrp.createDimension('nspe', fm.nfs)
    npargd = rootgrp.createDimension('levels', nl)
    nfovs  = rootgrp.createDimension('number_of_FOVs', len(outdata))
    jj = rootgrp.createDimension('jj', 2)
    jj = rootgrp.createDimension('sfc', 1)
    ncx      = rootgrp.createVariable('x', 'f4', ('number_of_FOVs','channels'))
    ncy      = rootgrp.createVariable('y', 'f4', ('number_of_FOVs','channels'))
    
    """
    nckt     = rootgrp.createVariable('Kt', 'f4', ('number_of_FOVs','channels', 'levels'))
    nckwv    = rootgrp.createVariable('Kwv', 'f4', ('number_of_FOVs','channels', 'levels'))
    nckco2   = rootgrp.createVariable('Kco2', 'f4', ('number_of_FOVs','channels', 'levels'))
    ncko3    = rootgrp.createVariable('Ko3', 'f4', ('number_of_FOVs','channels', 'levels'))
    nckskt   = rootgrp.createVariable('Kskt', 'f4', ('number_of_FOVs','channels', 'sfc'))
    ncksp    = rootgrp.createVariable('Ksp', 'f4', ('number_of_FOVs','channels', 'sfc'))
    ncxkemrf = rootgrp.createVariable('xkemrf', 'f4',
                                    ('number_of_FOVs','nspe', 'channels', 'jj'))
    """
    for ifov,fov_output in enumerate(outdata):
        ncx[ifov,:] = fm.cwvn
        ncy[ifov,:] = fov_output['y']
        """
        nckt[ifov,:]         = np.transpose(fov_output['xkt'][0:nl, :])
        nckskt[ifov,:,:]     = np.transpose(fov_output['xkt'][nl, :])
        ncksp[ifov,:,:]      = np.transpose(fov_output['xkt'][nl+1, :])
        nckwv[ifov,:,:]      = np.transpose(fov_output['xkt'][nl+2:nl+2+nl, :])
        nckco2[ifov,:,:]     = np.transpose(fov_output['xkt'][2*nl+2:2*nl+2+nl, :])
        ncko3[ifov,:,:]      = np.transpose(fov_output['xkt'][3*nl+2:3*nl+2+nl, :])
        
        if version == 2:
            ncxkemrf[ifov,:,:,:] = fov_output['xkemrf'].transpose()
        else:
            ncxkemrf[ifov,:,:,:] = fov_output['xkemrf']
        """

    rootgrp.close()
    print("Saved output test in {}".format(outfile))

    return
    
def generate_spectra(file, hitran, solar, fm, outfile, fm_version = 3):
        
    
    # Output data (empty first call)
    
    # Read Test Case
    if '.json' in path.basename(file):
        indoss = read_JSON_scene(file, fm_version)
    elif 'fg.nc' == path.basename(file):
        indoss = read_fg(file,fm_version, emissivity_mode = 'randomize')
    else:
        sys.exit("Wrong input file format!")
    
    outdata = []
    for ossinput in indoss: # Iterate over number of FOVs
        outoss = {}
        fm.compute(ossinput, outoss)
        outdata.append(outoss)


    nspe = indoss[0]['sfgrd'].size
    nlev = len(indoss[0]['temp'])

    # Write output netCDF4
    save_output(fm, outdata, outfile, nlev, nspe, fm_version)
        
    return
    
 
def parser():
    import argparse


    parser=argparse.ArgumentParser()
    parser.add_argument('-fmv', 
                        '--fm_version', 
                        type    = int, 
                        choices = [2,3],
                        default = 3,
                        help    = 'OSS version')
    parser.add_argument('-i', 
                        '--input_file', 
                        type    =str, 
                        required = True,
                        help    = 'Input File')
    parser.add_argument('-o',
                        '--output_file', 
                        type=str,
                        required = True,
                        help='Output file')
    parser.add_argument('--hitran', 
                        type=str,
                        required = False,
                        default = 'data/leo.iasi.0.05.nc',
                        help='Coefficient file')
    parser.add_argument('-s',
                        '--solar_irradiance', 
                        type=str,
                        required = False,
                        default = 'data/solar_irradiances.nc',
                        help='Solar Irradiance file')
    parser.add_argument('-fm',
                        '--forward_model_lib', 
                        type=str,
                        required = False,
                        default = 'mirto/iasi/mirto_oss_v3/ossFM_files/ossFM',
                        help='OSSFM Library Path')


    return parser.parse_args()

############################################################################## MAIN
if __name__ == '__main__':
    from dobjects.Hitran import Hitran   # Check what Hitran is imported
    from dobjects.Solar import Solar     # Check what Solar  is imported

    args = parser()
    
    # Import FM Libraries
    FM_IMPORT_CMD = "from {} import ossFM".format(
        args.forward_model_lib.replace('/','.').replace('.py','').replace('..','')
        )
    exec(FM_IMPORT_CMD)
    
    solar = Solar(args.solar_irradiance)
    hitran = Hitran(args.hitran)
    fm     = ossFM(solar,hitran)
    generate_spectra(args.input_file, hitran, solar, fm, args.output_file, fm_version = args.fm_version)
    
