import logging
import numpy as np
from netCDF4 import Dataset, MFDataset
import sys, os
from datetime import datetime, timedelta
from atmos.mirto_atmos_tools import mr2rh


__author__    = "Paolo Antonelli and Paolo Scaccia"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__   = ["Paolo Antonelli","Paolo Scaccia"]
__license__   = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

if __name__ == '__main__':
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)

   
class dimension(object):
    def __init__(self,size):
        self.__dict__ = {'size':size}
        
class netcdf_structure(object):
    def __init__(self,dataset):
       self.variables  = {}
       self.groups     = {}
       self.dimensions = { dim : dimension(dataset.dimensions[dim].size) for dim in dataset.dimensions   } 

       for variable in list(dataset.variables):
            self.variables[variable] = np.asarray(dataset[variable])
            
       if list(dataset.groups) != []:
           for g in list(dataset.groups):
               self.groups[g] = {}                
               for variable in list(dataset[g].variables):
                   self.groups[g][variable] = np.asarray(dataset[g][variable])
       return

    def combine_dataset(self,dataset):
        for variable in list(dataset.variables):
            self.variables[variable] = np.concatenate([self.variables[variable],np.asarray(dataset[variable])])
        if list(dataset.groups) != []:
            for g in list(dataset.groups):
                for variable in list(dataset[g].variables):
                       self.groups[g][variable] = np.concatenate([self.groups[g][variable], np.asarray(dataset[g][variable])])
        return

class da_structure(object):
    def __init__(self,n_obs,nlev,filetype):
        
        self.year   = np.zeros(n_obs,dtype=np.int)
        self.month  = np.zeros_like(self.year)
        self.day    = np.zeros_like(self.year)
        self.hour   = np.zeros_like(self.year)        
        self.minute = np.zeros_like(self.year)
        self.second = np.zeros_like(self.year)
        self.lat    = np.zeros(n_obs,dtype=np.float64)
        self.lon    = np.zeros_like(self.lat)
        self.field_of_view                      = np.zeros_like(self.lat)
        self.posterior_surface_skin_temperature = np.zeros_like(self.lat)
        self.prior_surface_skin_temperature     = np.zeros_like(self.lat)
        self.sigma_level = np.arange(nlev)
        self.d2     = np.zeros_like(self.lat)
        self.air_pressure = np.zeros((n_obs,nlev))
        
        if filetype == 'TR':
            self.scaled_eigenvalues = None
            self.scaled_observation = None
            self.scaled_observation_operator = None
        elif filetype == 'BASIC':
            self.prior_se_coef_logit     = None
            self.posterior_se_coef_logit = None            
            self.prior_specific_humidity = np.zeros((n_obs,nlev))
            self.posterior_specific_humidity = np.zeros_like(self.prior_specific_humidity)
            self.prior_relative_humidity     = np.zeros_like(self.prior_specific_humidity)
            self.posterior_relative_humidity = np.zeros_like(self.prior_specific_humidity)
            self.prior_ozone     = np.zeros_like(self.prior_specific_humidity)
            self.posterior_ozone = np.zeros_like(self.prior_specific_humidity)
            self.priori_air_temperature = np.zeros_like(self.prior_specific_humidity)
            self.posterior_air_temperature = np.zeros_like(self.prior_specific_humidity)

        return
        
        
        
def read_nc(filepaths):
    """
    %Read in simple netcdf files (with no groups)
    %
    %Paolo Antonelli
    %Copyright  AdaptiveMeteo
    %Mon Dec 12 12:00:04 UTC 2016
    %Edited by Paolo Scaccia Mon Aug 3 15:00:00 UTC 2020
    """
    if type(filepaths) is str:
        if os.path.isfile(filepaths):
            LOGGER.debug("{} found!".format(filepaths))
        else:
            LOGGER.info("MIRTO file {} not found!".format(filepaths))
            LOGGER.error("MIRTO file {} not found!".format(filepaths))
    
            return None
        try:
            data = Dataset(filepaths,"r",dtype=np.float64) # read as a double
            return data
        
        except:
            LOGGER.info("Reading error with file {}".format(filepaths))
            LOGGER.error("Reading error with file {}".format(filepaths))
            return None
    elif type(filepaths) is list:
            # Aggregate netcdf files
            if len(filepaths) != 0:
                for i,filepath in enumerate(filepaths):
                    
                       data = Dataset(filepath,"r") # read as a double
                       if i == 0:
                           aggregated_data = netcdf_structure(data)
                       else:
                           aggregated_data.combine_dataset(data)
            
                return aggregated_data
            else:
                LOGGER.info("No overpassees found!")
    else:
        raise IOError        
        
def clean_dir_list(path,dir_list):
    """
        Remove from the list directories with missing files
    """
    outlist = []

    for dir_name in dir_list:        
        results_file = "/".join([path,dir_name,"/mirto/results.nc"])
        fov_file     = results_file.replace("results","fov")
        fg_file      = results_file.replace("results","fg")

        if os.path.isfile(results_file):
            if os.path.isfile(fov_file):
                if os.path.isfile(fg_file):
                    outlist.append(dir_name)
                else:
                    LOGGER.info("MIRTO file {} not found!".format(fg_file))
            else:
                LOGGER.info("MIRTO file {} not found!".format(fov_file))
        else:
            LOGGER.info("MIRTO file {} not found!".format(results_file))
            
    return outlist

def qc_innovation(good_retr,up2level,fg,results,path,N):
    """
    Returns the indices of the retrievals which depart (at any level) from the first guess
    for less than N*sigma
    """
    
    out_indx = []
    nlev = fg.groups['atmospheric_components'].dimensions['number_of_atmospheric_levels'].size
    sigma = Dataset(apriori_path,"r")
    
    for indx in good_retr:
        sigma_T = sigma['atmospheric_components']['Covariances']['T'][indx,:nlev,:nlev]
        diag_sigma_T = np.diag(sigma_T)
        sigma_logq   = sigma['atmospheric_components']['Covariances']['q'][indx,:nlev,:nlev]
        diag_sigma_logq = np.diag(sigma_logq)
        
        fg_T = fg.grous['atmospheric_components']['T'][indx,:,:nlev]
        fg_q = fg.groups['atmospheric_components']['q'][indx,:nlev,:nlev]
        results_T = results.variables['t'][indx,:nlev,:nlev]
        results_q = results.variables['q'][indx,:nlev,:nlev]
        delta_T = fg_T - results_T
        
        # q in fg.nc e results.nc is in g/kg; q in mirto is log(kg/kg)
        delta_logq = np.log(fg_q/1000) - np.log(results_q/1000)
        # to be completed
        
        if ~np.any( np.abs(delta_T[:up2level]) > 3*np.sqrt(diag_sigma_T[:up2level]) ) and \
           ~np.any( delta_logq[:up2level] > N*np.sqrt(diag_sigma_logq[:up2level])  ):
               out_indx.append(indx)
             
    return np.array(out_indx)

def thin_tr(fov,thin_dist,good_indx):    
    from geometry.earth_geometry import haversine
    
    nrec = good_indx.size        
    lat = fov.variables['Latitude'][good_indx]
    lon = fov.variables['Longitude'][good_indx]
    
    keep = []
    indx = np.arange(nrec)
    rec_indx = np.copy(good_indx)
    
    while indx.size != 0:
        rec_indx = rec_indx[indx]
        tmp_lat = lat[indx]
        tmp_lon = lon[indx]
        keep.append(rec_indx[0])
        d = haversine(lon[0],lat[0],tmp_lon,tmp_lat)
        indx = np.where(d > thin_dist)[0]
    
    LOGGER.info("THINNED RETRIEVALS: {:d} ({:.1f}% of the total FOVs)".format(len(keep),len(keep)*100/fov.dimensions['number_of_FOVs'].size))
    return np.array(keep)

def fill_tr_basic_structure(fov,results,fg,rh,good_indx):
    """
    % function sps_conc, basic_conc]=mirto_fill_sps_basic_structures(fov,results,fg)
    %
    % Initizizlie the sps_conc and basic_conc structures with empty fileds and
    % with the proper date and time.
    %
    % Paolo Antonelli
    % Wed Aug  2 04:00:23 CEST 2017
    """
    nlev = results.dimensions['levels'].size

    tr_conc = da_structure(good_indx.size,nlev,'TR')
    basic_conc = da_structure(good_indx.size,nlev,'BASIC')

    nrec = results.variables['t'][:].shape[1]
    nlev0 = results.variables['p'][:].shape[0]
    
    for i,indx in enumerate(good_indx):
        
        time_unix      = np.double(fov.variables['Time'][indx])
        time_obs       = datetime(1970,1,1) + timedelta(milliseconds=time_unix)
        
        tr_conc.year[i]  = time_obs.year
        tr_conc.month[i] = time_obs.month
        tr_conc.day[i]   = time_obs.day
        tr_conc.minute[i] = time_obs.minute
        tr_conc.second[i] = time_obs.second
        
        basic_conc.year[i]  = time_obs.year
        basic_conc.month[i] = time_obs.month
        basic_conc.day[i]   = time_obs.day
        basic_conc.minute[i] = time_obs.minute
        basic_conc.second[i] = time_obs.second
        
        
    # Concatenate variables for TR (Transformed Retrievals) files
    tr_conc.lat = fov.variables['Latitude'][good_indx]
    tr_conc.lon = fov.variables['Longitude'][good_indx]
    tr_conc.air_pressure  = results.variables['p'][good_indx,:]
    tr_conc.field_of_view = fov.variables['FOV_angle'][good_indx]
    tr_conc.d2            = results.variables['d2'][good_indx]
    tr_conc.posterior_surface_skin_temperature = results.variables['skt'][good_indx]
    tr_conc.prior_surface_skin_temperature     = fg.groups['atmospheric_components']['skT'][good_indx]
    tr_conc.scaled_eigenvalues = results.variables['DA_Lambda'][good_indx,:]
    tr_conc.scaled_observation = results.variables['DA_Yret'][good_indx,:]
    mnel = results.dimensions['mnel'].size
    tr_conc.scaled_observation_operator = np.zeros((good_indx.size,mnel,2*nlev) )
    tr_conc.scaled_observation_operator[:,:,:nlev] = results.variables['DA_Hret'][good_indx,:,:nlev]
    tr_conc.scaled_observation_operator[:,:,nlev:] = results.variables['DA_Hret'][good_indx,:,nlev:]
    
    # Concatenate variables for BASIC files
    
    # Retrievals (FOV) Elements
    basic_conc.lat = fov.variables['Latitude'][good_indx]
    basic_conc.lon = fov.variables['Longitude'][good_indx]
    basic_conc.field_of_view = fov.variables['FOV_angle'][good_indx]
    basic_conc.d2            = results.variables['d2'][good_indx]
    
    # Surface Elements
    basic_conc.posterior_surface_skin_temperature = results.variables['skt'][good_indx]
    basic_conc.prior_surface_skin_temperature     = fg.groups['atmospheric_components']['skT'][good_indx]
    basic_conc.prior_se_coef_logit = np.zeros_like(results.variables['ems_coeff'][good_indx,:])
    basic_conc.posterior_se_coef_logit = results.variables['ems_coeff'][good_indx,:]
    
    # Atmospheric Elements
    basic_conc.air_pressure  = results.variables['p'][good_indx,:]
    basic_conc.prior_air_temperature     = fg.groups['atmospheric_components']['T'][good_indx,:nlev]
    basic_conc.posterior_air_temperature = results.variables['t'][good_indx,:nlev]
    basic_conc.prior_specific_humidity   = fg.groups['atmospheric_components']['q'][good_indx,:nlev]/1000
    basic_conc.posterior_specific_humidity   = results.variables['q'][good_indx,:nlev]/1000
    basic_conc.prior_relative_humidity       = mr2rh(fg.groups['atmospheric_components']['p'][good_indx,:nlev],
                                                     fg.groups['atmospheric_components']['T'][good_indx,:nlev],
                                                     fg.groups['atmospheric_components']['q'][good_indx,:nlev],
                                                     273.15)[0]
    basic_conc.posterior_relative_humidity   = mr2rh(results.variables['p'][good_indx,:nlev],
                                                     results.variables['t'][good_indx,:nlev],
                                                     results.variables['q'][good_indx,:nlev],
                                                     273.15)[0]
    basic_conc.priori_ozone = fg.groups['atmospheric_components']['O3'][good_indx,:nlev]
    basic_conc.posterior_ozone = results.variables['o3'][good_indx,:nlev]
    
    
    return [tr_conc , basic_conc]



def da_write_tr_ncfile(tr_conc,astime_str,outpath):
    
    """    
    %function da_write_sps_ncfile(basic_conc,astime_str)
    %
    %Writes tr data into WRFDA input format (netcdf file)
    %
    % Paolo Antonelli
    % Mon Aug  8 01:54:11 CEST 2016
    %
    % -- Python version -- 
    % Paolo Scaccia
    % Wed Aug  5 16:00:10 CEST 2020
    """
    
    n_obs = tr_conc.lat.size
    n_tr_eig = tr_conc.scaled_observation.shape[1]
    n_obs_statev_dim = tr_conc.scaled_observation_operator.shape[-1]
    n_levels = tr_conc.sigma_level.size
    
        
    """
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Write Input File in netcdf format
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    """
    
    tr_file = "{}_TR_FULL_XX.nc".format(astime_str)

    with Dataset(outpath +'/' + tr_file,"w",format="NETCDF4") as ncfile:

        ncfile.createDimension('number_of_valid_retrievals',n_obs)
        ncfile.createDimension('number_of_sigma_levels',n_levels)
        ncfile.createDimension('number_of_eigenvalues',n_tr_eig)
        ncfile.createDimension('size_of_reduced_state_vector',n_obs_statev_dim)

        record_number = ncfile.createVariable('record_number',np.int,('number_of_valid_retrievals',))
        record_number[:] = np.arange(n_obs)  # starting point 0 or 1 ???
        record_number.long_name     = 'Record number of the original orbit'
        record_number.standard_name = 'rid'
        record_number.valid_min     = 0
        record_number.valid_max     = 10000

        year = ncfile.createVariable('year',np.int,('number_of_valid_retrievals',))
        year[:] = tr_conc.year[:]
        month = ncfile.createVariable('month',np.int,('number_of_valid_retrievals',))
        month[:] = tr_conc.month[:]
        day = ncfile.createVariable('day',np.int,('number_of_valid_retrievals',))
        day[:] = tr_conc.day[:]
        hour = ncfile.createVariable('hour',np.int,('number_of_valid_retrievals',))
        hour[:] = tr_conc.hour[:]
        minute = ncfile.createVariable('minute',np.int,('number_of_valid_retrievals',))
        minute[:] = tr_conc.minute[:]
        second = ncfile.createVariable('second',np.int,('number_of_valid_retrievals',))
        second[:] = tr_conc.second[:]
        
        lat = ncfile.createVariable('latitude',np.float64,('number_of_valid_retrievals',))
        lat[:] = tr_conc.lat[:]
        lat.long_name     = 'Latitude'
        lat.standard_name = 'Latitude'
        lat.units         = 'degree'
        lat.valid_min     = -90
        lat.valid_max     = 90
        
        lon = ncfile.createVariable('longitude',np.float64,('number_of_valid_retrievals',))
        lon[:] = tr_conc.lon[:]
        lon.long_name     = 'Longitude'
        lon.standard_name = 'Longitude'
        lon.units         = 'degree'
        lon.valid_min     = 1
        lon.valid_max     = 360

        fovs = ncfile.createVariable('field_of_view',np.float64,('number_of_valid_retrievals',))
        fovs[:] = tr_conc.field_of_view[:]
        
        d2 = ncfile.createVariable('d2',np.float64,('number_of_valid_retrievals',))
        d2[:] = tr_conc.d2[:]
        
        sigma = ncfile.createVariable('sigma_level',np.int,('number_of_sigma_levels',))
        sigma[:] = tr_conc.sigma_level[:]
        
        pres = ncfile.createVariable('air_pressure',np.float64,('number_of_sigma_levels','number_of_valid_retrievals'))
        pres.units = 'Pa'
        pres[:,:] = tr_conc.air_pressure.transpose()

        pres.long_name = 'atmospheric pressure at sigma levels'
        pres.standard_name = 'pressure'
        pres.valid_min     = 0
        pres.valid_max     = 150000
        
        post_skt = ncfile.createVariable('posterior_surface_skin_temperature',np.float64,('number_of_valid_retrievals',))
        post_skt[:] = tr_conc.posterior_surface_skin_temperature[:]
        post_skt.long_name     = 'Posterior Surface Temperature'
        post_skt.standard_name = 'Surface Temperature'
        post_skt.units         = 'K'
        post_skt.valid_min     = 0
        post_skt.valid_max     = 360

        prior_skt = ncfile.createVariable('prior_surface_skin_temperature',np.float64,('number_of_valid_retrievals',))
        prior_skt[:] = tr_conc.prior_surface_skin_temperature[:]
        prior_skt.long_name     = 'Prior Surface Temperature'
        prior_skt.standard_name = 'Surface Temperature'
        prior_skt.units         = 'K'
        prior_skt.valid_min     = 0
        prior_skt.valid_max     = 360
        
        if type(tr_conc.scaled_eigenvalues) != type(None):
            scaled_eig = ncfile.createVariable('scaled_eigenvalues',np.float64,('number_of_eigenvalues','number_of_valid_retrievals'))
            scaled_eig[:,:] = tr_conc.scaled_eigenvalues.transpose()
            scaled_eig.long_name     = 'Scaled tq eigenvalues for all eigenvectors'
            scaled_eig.standard_name = 'Scaled tqproduct'
            scaled_eig.units         = '-'
        
        scaled_obs = ncfile.createVariable('scaled_observation',np.float64,('number_of_eigenvalues','number_of_valid_retrievals'))
        
        scaled_obs[:,:] = tr_conc.scaled_observation[:,:].transpose()
        scaled_obs.long_name     = 'Scaled tq products for all eigenvectors'
        scaled_obs.standard_name = 'Scaled tqproduct'
        scaled_obs.units         = '-'
        
        scaled_obs_operator = ncfile.createVariable('scaled_observation_operator',np.float64,('size_of_reduced_state_vector','number_of_eigenvalues','number_of_valid_retrievals'))
        scaled_obs_operator[:,:] = tr_conc.scaled_observation_operator.transpose()
        scaled_obs_operator.long_name = 'Observation operator for tq products for all eigenvectors'
        scaled_obs_operator.standard_name = 'Observation operator for tq products'
        scaled_obs_operator.units         = '-'
    
        ncfile.creation_date = datetime.utcnow().isoformat()[:-7]
        ncfile.created = "Paolo Antonelli"
        ncfile.company = "AdaptiveMeteo"
        
        LOGGER.info("{} FOVs written in {}".format(n_obs,outpath + '/' + tr_file))
    return


def da_write_basic_ncfile(basic_conc,astime_str,outpath):
    """    
    %function da_write_sps_ncfile(basic_conc,astime_str)
    %
    % Writes  data into WRFDA input format (netcdf file)
    %
    % Paolo Antonelli
    % Mon Aug  8 01:54:11 CEST 2016
    %
    % -- Python version -- 
    % Paolo Scaccia
    % Wed Aug  5 16:00:10 CEST 2020
    """
    
    n_obs      = basic_conc.lat.size
    n_ems_coef = basic_conc.prior_se_coef_logit.shape[1]
    n_levels   = basic_conc.sigma_level.size
    
    """
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Write Input File in netcdf format
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    """

    basic_file = "{}_BASIC_FULL_XX.nc".format(astime_str)
    
    with Dataset(outpath +'/' + basic_file,"w",format="NETCDF4") as ncfile:

        ncfile.createDimension('number_of_valid_retrievals',n_obs)
        ncfile.createDimension('number_of_sigma_levels',n_levels)
        ncfile.createDimension('number_of_surface_emissivity_elements',n_ems_coef)

        record_number = ncfile.createVariable('record_number',np.int,('number_of_valid_retrievals',))
        record_number[:] = np.arange(n_obs) + 1
        record_number.long_name     = 'Record number of the original orbit'
        record_number.standard_name = 'rid'
        record_number.valid_min     = 0
        record_number.valid_max     = 10000

        year = ncfile.createVariable('year',np.int,('number_of_valid_retrievals',))
        year[:] = basic_conc.year[:]
        month = ncfile.createVariable('month',np.int,('number_of_valid_retrievals',))
        month[:] = basic_conc.month[:]
        day = ncfile.createVariable('day',np.int,('number_of_valid_retrievals',))
        day[:] = basic_conc.day[:]
        hour = ncfile.createVariable('hour',np.int,('number_of_valid_retrievals',))
        hour[:] = basic_conc.hour[:]
        minute = ncfile.createVariable('minute',np.int,('number_of_valid_retrievals',))
        minute[:] = basic_conc.minute[:]
        second = ncfile.createVariable('second',np.int,('number_of_valid_retrievals',))
        second[:] = basic_conc.second[:]
        
        lat = ncfile.createVariable('latitude',np.float64,('number_of_valid_retrievals',))
        lat[:] = basic_conc.lat[:]
        lat.long_name     = 'Latitude'
        lat.standard_name = 'Latitude'
        lat.units         = 'degree'
        lat.valid_min     = -90
        lat.valid_max     = 90
        
        lon = ncfile.createVariable('longitude',np.float64,('number_of_valid_retrievals',))
        lon[:] = basic_conc.lon[:]
        lon.long_name     = 'Longitude'
        lon.standard_name = 'Longitude'
        lon.units         = 'degree'
        lon.valid_min     = 1
        lon.valid_max     = 360

        fovs = ncfile.createVariable('field_of_view',np.float64,('number_of_valid_retrievals',))
        fovs[:] = basic_conc.field_of_view[:]
        
        d2 = ncfile.createVariable('d2',np.float64,('number_of_valid_retrievals',))
        d2[:] = basic_conc.d2[:]
        
        sigma = ncfile.createVariable('sigma_level',np.int,('number_of_sigma_levels',))
        sigma[:] = basic_conc.sigma_level[:]
        
        pres = ncfile.createVariable('air_pressure',np.float64,('number_of_sigma_levels','number_of_valid_retrievals'))
        pres.units = 'hPa'
        pres[:,:] = basic_conc.air_pressure.transpose()
        pres.long_name = 'atmospheric pressure at sigma levels'
        pres.standard_name = 'pressure'
        pres.valid_min     = 0
        pres.valid_max     = 1500
        
        post_skt = ncfile.createVariable('posterior_surface_skin_temperature',np.float64,('number_of_valid_retrievals',))
        post_skt[:] = basic_conc.posterior_surface_skin_temperature[:]
        post_skt.long_name     = 'Posterior Surface Temperature'
        post_skt.standard_name = 'Surface Temperature'
        post_skt.units         = 'K'
        post_skt.valid_min     = 0
        post_skt.valid_max     = 360

        prior_skt = ncfile.createVariable('prior_surface_skin_temperature',np.float64,('number_of_valid_retrievals',))
        prior_skt[:] = basic_conc.prior_surface_skin_temperature[:]
        prior_skt.long_name     = 'Prior Surface Temperature'
        prior_skt.standard_name = 'Surface Temperature'
        prior_skt.units         = 'K'
        post_skt.valid_min     = 0
        post_skt.valid_max     = 360

        prior_air_t = ncfile.createVariable('prior_air_temperature',np.float64,('number_of_sigma_levels','number_of_valid_retrievals'))
        prior_air_t[:] = basic_conc.prior_air_temperature.transpose()
        prior_air_t.long_name     = 'Prior Atmospheric Temperature at sigma level'
        prior_air_t.standard_name = 'prior temperature'
        prior_air_t.units         = 'K'
        prior_air_t.valid_min     = 0
        prior_air_t.valid_max     = 350
        
        post_air_t = ncfile.createVariable('posterior_air_temperature',np.float64,('number_of_sigma_levels','number_of_valid_retrievals'))
        post_air_t[:] = basic_conc.posterior_air_temperature.transpose()
        post_skt.long_name     = 'Posterior Atmospheric Temperature at sigma level'
        post_skt.standard_name = 'posterior temperature'
        post_skt.units         = 'K'
        post_skt.valid_min     = 0
        post_skt.valid_max     = 350

        prior_spec_hum = ncfile.createVariable('prior_specific_humidity',np.float64,('number_of_sigma_levels','number_of_valid_retrievals'))
        prior_spec_hum[:] = basic_conc.prior_specific_humidity.transpose()
        prior_spec_hum.long_name     = 'Prior atmospheric specific humidity at sigma level'
        prior_spec_hum.standard_name = 'prior specific humidity'
        prior_spec_hum.units         = 'kg/kg'
        prior_spec_hum.valid_min     = 0
        prior_spec_hum.valid_max     = 0.06

        post_spec_hum = ncfile.createVariable('posterior_specific_humidity',np.float64,('number_of_sigma_levels','number_of_valid_retrievals'))
        post_spec_hum[:] = basic_conc.posterior_specific_humidity.transpose()
        post_spec_hum.long_name     = 'Posterior atmospheric specific humidity at sigma level'
        post_spec_hum.standard_name = 'posterior specific humidity'
        post_spec_hum.units         = 'kg/kg'
        post_spec_hum.valid_min     = 0
        post_spec_hum.valid_max     = 0.06

        prior_rel_hum = ncfile.createVariable('prior_relative_humidity',np.float64,('number_of_sigma_levels','number_of_valid_retrievals'))
        prior_rel_hum[:] = basic_conc.prior_relative_humidity.transpose()
        prior_rel_hum.long_name     = 'Prior atmospheric relative humidity at sigma level'
        prior_rel_hum.standard_name = 'prior relative humidity'
        prior_rel_hum.units         = 'percent'
        prior_rel_hum.valid_min     = 0
        prior_rel_hum.valid_max     = 100

        post_rel_hum = ncfile.createVariable('posterior_relative_humidity',np.float64,('number_of_sigma_levels','number_of_valid_retrievals'))
        post_rel_hum[:] = basic_conc.posterior_relative_humidity.transpose()
        post_rel_hum.long_name     = 'Posterior atmospheric relative humidity at sigma level'
        post_rel_hum.standard_name = 'posterior relative humidity'
        post_rel_hum.units         = 'percent'
        post_rel_hum.valid_min     = 0
        post_rel_hum.valid_max     = 100

        prior_ozone = ncfile.createVariable('prior_ozone',np.float64,('number_of_sigma_levels','number_of_valid_retrievals'))
        prior_ozone[:] = basic_conc.prior_ozone.transpose()
        prior_ozone.long_name     = 'Prior atmospheric ozone at sigma level'
        prior_ozone.standard_name = 'prior ozone'
        prior_ozone.units         = 'kg/kg'
        prior_ozone.valid_min     = 0
        prior_ozone.valid_max     = 10
    
        posterior_ozone = ncfile.createVariable('posterior_ozone',np.float64,('number_of_sigma_levels','number_of_valid_retrievals'))
        posterior_ozone[:] = basic_conc.posterior_ozone.transpose()
        posterior_ozone.long_name     = 'Posterior atmospheric ozone at sigma level'
        posterior_ozone.standard_name = 'posterior ozone'
        posterior_ozone.units         = 'kg/kg'
        posterior_ozone.valid_min     = 0
        posterior_ozone.valid_max     = 10
    
        prior_se_coef_logit = ncfile.createVariable('prior_se_coef_logit',np.float64,('number_of_surface_emissivity_elements','number_of_valid_retrievals'))
        prior_se_coef_logit[:] = basic_conc.prior_se_coef_logit.transpose()
        prior_se_coef_logit.long_name  = 'Prior pc coefficient of se in logit'
        prior_ozone.standard_name = 'prior se_ceof_logit'
        prior_ozone.units         = '-'
        prior_ozone.valid_min     = np.nan
        prior_ozone.valid_max     = np.nan
        
        posterior_se_coef_logit = ncfile.createVariable('posterior_se_coef_logit',np.float64,('number_of_surface_emissivity_elements','number_of_valid_retrievals'))
        posterior_se_coef_logit[:] = basic_conc.posterior_se_coef_logit.transpose()
        posterior_se_coef_logit.long_name  = 'Posterior pc coefficient of se in logit'
        posterior_ozone.standard_name = 'posterior se_ceof_logit'
        posterior_ozone.units         = '-'
        posterior_ozone.valid_min     = np.nan
        posterior_ozone.valid_max     = np.nan
                
        ncfile.creation_date = datetime.utcnow().isoformat()[:-7]
        ncfile.created = "Paolo Antonelli"
        ncfile.company = "AdaptiveMeteo"
        
        LOGGER.info("{} FOVs written in {}".format(n_obs,outpath + '/' + basic_file))

    return

def generate_da_input_files(path,apriori_path,outdir,astime,delta):
    """
    generate_da_input_files read MIRTO files and write results in DA format
    Paolo Antonelli
    Sat Sep  2 03:25:22 UTC 2017
    
    Copyright AtaptiveMeteo S.r.l.
    
    INPUT
        path - Path to MIRTO files
        astime - Assimilation time (datetime object)
        delta  - Time interval for the assimilation (float with minutes)
    """
    
    if path[-1] == "/":
        path = path[:-1]
        
    # Define time interval centered around astime
    delta_min = timedelta(minutes=delta)
    time_max = astime + delta_min
    time_min = astime - delta_min
    
    # Select overpasses with time in the defined interval  
    dir_list = [ dir_name for dir_name in os.listdir(path) if '.' not in dir_name]
    dir_list = [ dir_name for dir_name in dir_list if datetime.strptime(dir_name,"%Y%m%d_%H%M%S") <= time_max
                                                        and datetime.strptime(dir_name,"%Y%m%d_%H%M%S") >= time_min ]
    n_overpasses = len(dir_list)    
    if n_overpasses != 0:
        LOGGER.info("{} overpasses found".format(n_overpasses))
    else:
        LOGGER.info("No overpasses found")
        return
    
    if n_overpasses == 1: 
        # Case with only one overpass
        dir_name = dir_list[0]
        results_file = "/".join([path,dir_name,"mirto/results.nc"])
        fov_file = "/".join([path,dir_name,"mirto/fov.nc"])
        fg_file = "/".join([path,dir_name,"mirto/fg.nc"])
        
        LOGGER.info("Overpass {} found in time window".format(dir_name))

        results = read_nc(results_file)
        fov     = read_nc(fov_file)
        fg      = read_nc(fg_file)
        
    else:
        # Case with multi overpasses

        dir_list = clean_dir_list(path,dir_list)
        
        results_files = [ "/".join([path,name,"mirto/results.nc"])for name in dir_list]
        fov_files     = [ "/".join([path,name,"mirto/fov.nc"])for name in dir_list]
        fg_files = [ "/".join([path,name,"mirto/fg.nc"])for name in dir_list]

        results = read_nc(results_files)
        fov = read_nc(fov_files)
        fg  = read_nc(fg_files)
        

    if results != None and fov != None and fg != None:
        
        # Compute RH and QC        
        rh = mr2rh(results.variables['p'][:],results.variables['t'][:],results.variables['q'][:],273.15)[0]
        qc = np.zeros_like(results.variables['d2'][:],dtype=bool)
        good_indx = np.where( (results.variables['d2'][:] < 5) & (np.amax(rh,axis=1) <= 100))[0]
        LOGGER.info("SUCCESSFUL RETRIEVALS: {:d} ({:.1f}%)".format(good_indx.size,good_indx.size*100/qc.size))
        
        #good_indx = qc_innovation(good_indx,50,fg,results,apriori_path,3)
        #LOGGER.info("OPTIMAL RETRIEVALS: {:d} ({:.1f}%)".format(good_indx.size,good_indx.size*100/qc.size))
        qc[good_indx] = True
        
        # Thin TR
        thin_val = 80
        LOGGER.info("Thinning at {:d} km".format(thin_val))
        good_indx = thin_tr(fov,thin_val,good_indx)
        
        # Prepare out
        if good_indx.size!=0:
    
            [tr_conc, basic_conc] = fill_tr_basic_structure(fov,results,fg,rh,good_indx)
            
            LOGGER.info("Writing TR file...")
            da_write_tr_ncfile(tr_conc,astime.strftime("%Y%m%d_%H%M%S"),outdir)
            LOGGER.info("Writing BASIC file...")
            da_write_basic_ncfile(basic_conc,astime.strftime("%Y%m%d_%H%M%S"),outdir)
        else:
            LOGGER.info("No retrievals after the thinning")
    return 



###### MAIN
    

def config_logger(verbosity):
    LOGGER.setLevel(verbosity)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - '
                                  '%(funcName)s: %(message)s',
                                  datefmt='%m/%d/%Y %H:%M:%S')

    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(verbosity)
    streamhandler.setFormatter(formatter)
    LOGGER.addHandler(streamhandler)
    
def test():
    from datetime import datetime
    
    config_logger(getattr(logging, "INFO"))
    outdir = "."
    # Test single file
    date = datetime(2020,8,5,12)
    path = os.getcwd()
    apriori_path = './apriori.nc'
    generate_da_input_files(path,apriori_path,outdir,date,100)
    
    return


def prepare_parser():
    from argparse import ArgumentParser

    v_levels = ['debug', 'info', 'warning']
    dyns     = ['syn','asyn','map_syn','map_syn_langevin',
                'int_syn','map_asyn','int_syn_field',
                'map_syn_langevin_RK','int_syn_milstein']
    
    parser = ArgumentParser()
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='The level of verbosity of the software')
    parser.add_argument('--minutes','-m', default=75, type=int,
                        help='Time window (in minutes) for the aggregation')
    parser.add_argument('--outdir','-o',default='.',type=str,
                        help='Output directory for TR and BASIC file')
    parser.add_argument('--apriori_path','-a',default='.',type=str,
                        help='Directory for the apriori nc file')
    parser.add_argument('--workdir','-w',default='.',type=str,
                        help='Working directory with mirto data')
    parser.add_argument('--source_iasi','-si',default='.',type=str,
                        help='Directory with iasi data')
    parser.add_argument('--source_cris','-sc',default='.',type=str,
                        help='Directory with cris data')
    parser.add_argument('--time','-t',default=datetime.utcnow().strftime("%Y%m%d_%H%M%S"),type=str,
                        help='Directory with iasi data')
   
    return parser.parse_args()

def main():
    
    # Config logger
    config_logger(getattr(logging, "INFO"))
    
    # Parse Arguments
    argv = prepare_parser()    

    # Assign arguments
    outdir = argv.outdir
    path   = argv.workdir
    source_iasi = argv.source_iasi
    source_cris = argv.source_cris
    date = datetime.strptime(argv.time,"%Y%m%d_%H%M%S")
    apriori_path = argv.apriori_path
    
    # Link source to workdir
    os.system("ln -sf {}/2* {}".format(source_iasi,path))
    os.system("ln -sf {}/2* {}".format(source_cris,path))
    
    # Generate DA files    
    generate_da_input_files(path,apriori_path,outdir,date,argv.minutes)
    
if __name__ == '__main__':
    main()
    
