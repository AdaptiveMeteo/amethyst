import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import argparse
import amethyst_config

'''
Generate the TR Obesrvation Error for non diagonal matrix.

- obserr_channelsel_paolo_cnes.nc: is the CNES non diagonal matrix for IASI ObsErr
- obserr_diag_res_new_1_3_sps.nc: is the diagonal observation error for IASI already reduced in number of 
  channels for TR (sps) calculation. This file contains the second subselection indices: oe_sub_indices
- Output: obserr_channelsel_paolo_cnes_sps.nc. This is the non diagonal file for IASI noise to be used in the
  TR (SPS) calculation.

'''

def convert_obs_err_matrix(nc_fid,nc_fid_tr,outfile):
    
    # Ncdump     
    #nc_attrs, nc_dims, nc_vars         = read_nc(nc_fid)
    #nc_tr_attr, nc_tr_dims, nc_tr_vars = read_nc(nc_fid_tr)
    
    
    oe_sub_indices = np.loadtxt(amethyst_config.processor_vars["tr_chan_list"])
    
    #Read Observation Error matrix used in the inversion
    try:
        in_obs_err =  nc_fid.variables['obs_err'][:]
    except:
        in_obs_err =  nc_fid.variables['obserr'][:]
        
    
    #Remove channels which are picking high in the atmosphere
    in_oe_sub_indices = nc_fid_tr.variables['oe_sub_indices'][:]
    sub_obs_err = in_obs_err[np.ix_(oe_sub_indices.astype(int)-1,oe_sub_indices.astype(int)-1)]
    
    #Calculate Inverse
    inv_obs_err_tr = np.linalg.inv(sub_obs_err)
    
    #Get Observation Error SVD
    sub_obs_err_U, sub_obs_err_D, sub_obs_err_V = np.linalg.svd(sub_obs_err, full_matrices=True)
    
    #Check if decomposition is ok
    np.allclose(sub_obs_err, np.dot(sub_obs_err_U[:, :len(sub_obs_err_D)] * sub_obs_err_D, sub_obs_err_V))
    
    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    w_nc_fid = Dataset(outfile, 'w', format='NETCDF4')
    w_nc_fid.description = "The TR subselection of the Inversion Observation Error Covariance "
    
    selchannels = w_nc_fid.createDimension('selchannels',len(in_oe_sub_indices))
    
    obs_err        = w_nc_fid.createVariable('obs_err','f4',('selchannels','selchannels'))
    inv_obs_err    = w_nc_fid.createVariable('inv_obs_err','f4',('selchannels','selchannels'))
    obs_err_U      = w_nc_fid.createVariable('obs_err_U','f4',('selchannels','selchannels'))
    obs_err_D      = w_nc_fid.createVariable('obs_err_D','f4',('selchannels'))
    oe_sub_indices = w_nc_fid.createVariable('oe_sub_indices','f4',('selchannels'))

    obs_err[:,:] = sub_obs_err
    inv_obs_err[:,:] = inv_obs_err_tr
    obs_err_U[:,:] = sub_obs_err_U
    obs_err_D[:] = sub_obs_err_D
    oe_sub_indices[:] = in_oe_sub_indices
    
    w_nc_fid.close()


def read_nc(ncfid):
    # NetCDF global attributes
    nc_attrs = ncfid.ncattrs()
    for nc_attr in nc_attrs:
        print("Attributes {}".format(nc_attr))
    nc_dims = [dim for dim in ncfid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    for dim in nc_dims:
        print("Name: {}".format(dim))
        print("size: {}".format(len(ncfid.dimensions[dim])))
    # Variable information.
    nc_vars = [var for var in ncfid.variables]  # list of nc variables
    for var in nc_vars:
        if var not in nc_dims:
            print("Name: {}".format(var))
            print("dimensions: {}".format(ncfid.variables[var].dimensions))
    return nc_attrs, nc_dims, nc_vars


if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
                        help="Input file path")
    parser.add_argument('-i_sps', '--input_tr', required=True,
                        help="Input TR file path")   
    parser.add_argument('-o', '--output', required=True,
                        help="Output file path")   
    argv =parser.parse_args()
    nc_f = './obserr_channelsel_paolo_cnes.nc'  # Your filename
    nc_fid = Dataset(argv.input, 'r')  # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class    
    
    nc_fid_tr = Dataset(argv.input_tr, 'r')  # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class
    
    convert_obs_err_matrix(nc_fid,nc_fid_tr,argv.output)
    
    