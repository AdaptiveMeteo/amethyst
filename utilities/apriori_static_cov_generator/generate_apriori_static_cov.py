from netCDF4 import Dataset
import numpy as np
#from numpy.random import Generator, PCG64
import numpy.matlib
import glob, os
import logging
import xarray as xr

__author__    = "Paolo Antonelli and Paolo Scaccia"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__   = ["Paolo Antonelli","Paolo Scaccia"]
__license__   = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

def moving_avg(x, n):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[n:] - cumsum[:-n]) / float(n)


if __name__ == '__main__':
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)

LONGITUDE = 'Longitude'
LATITUDE = 'Latitude'
NATMLEVELS = 'number_of_atmospheric_levels'
ATMGROUP = 'atmospheric_components'
COVGROUP = 'Covariances'
NFOVS = 'number_of_FOVs'

np.set_printoptions(threshold=np.inf)

# Paolos 15/01/24
area      = 'arctic'
data_path ="/galaxy/data/Arctic_22/validation/test0001"  
cut_pressure_lev = True
n_lev_cut        = 70

#
os.chdir(data_path)
allfiles = glob.glob("*.nc")
nfiles = len(allfiles)
fovs = np.zeros([nfiles])



for index, file in enumerate(allfiles):

  #print(file)

  # Using readline()
  nc_f = file

  nc_fid = Dataset(nc_f, 'r')

  #Read fovs_levels in hPa
  fovs[index] = nc_fid.dimensions['n_fovs'].size       
  c = nc_fid.dimensions['n_lev']

nfovs = np.int(np.sum(fovs))
nlev = c.size

print('NFOVS: {}'.format(nfovs)) 
print('NLEV: {}'.format(nlev)) 

fovs_levels = np.matrix(np.zeros([nlev,nfovs]))

pre_dt = np.ma.zeros([nlev,nfovs])
post_dt = np.ma.zeros([nlev,nfovs])

pre_drh = np.ma.zeros([nlev,nfovs])
post_drh = np.ma.zeros([nlev,nfovs])

pre_dwv = np.ma.zeros([nlev,nfovs])
post_dwv = np.ma.zeros([nlev,nfovs])

pre_dlogwv = np.ma.zeros([nlev,nfovs])
post_dlogwv = np.ma.zeros([nlev,nfovs])

pre_dtdlogwv = np.ma.zeros([2*nlev,nfovs])
post_dtdlogwv = np.ma.zeros([2*nlev,nfovs])

cum_sum = 0

# Paolos 15/01/24
all_pressures = None

for index, file in enumerate(allfiles):

  print(file)

  # Using readline()
  nc_f = file

  nc_fid = Dataset(nc_f, 'r')

  #Read fovs_levels in hPa
  p = nc_fid.variables['fovs_levels'][:]
  #print('Pressure grid: {}'.format(p)) 
  print('PMAX index ',index,np.ma.max(p))
  # Paolos 15/01/24
  # Store pressure grid
  if all_pressures is None:
     all_pressures = np.copy(p)
  else:
     all_pressures = np.vstack((all_pressures,p))
  
  first_ele = cum_sum
  last_ele =  cum_sum + np.int(fovs[index])
  n_ele = np.int(fovs[index])
  fovs_levels[:,first_ele:last_ele] = p.reshape(nlev,n_ele)

  #
  # TEMPERATURE
  #
  #Read Temperature for Prior Posterior and Sonde
  prior_temp = nc_fid.variables['prior_temp'][:]
  post_temp = nc_fid.variables['post_temp'][:]
  sonde_temp = nc_fid.variables['sonde_temp'][:]

  #Get Temperature differences
  prior_delta_temp = prior_temp - sonde_temp
  post_delta_temp = post_temp - sonde_temp

  prior_delta_temp.mask = np.copy(sonde_temp.mask)
  post_delta_temp.mask = np.copy(sonde_temp.mask)

  pre_dt[:,first_ele:last_ele] = np.transpose(prior_delta_temp)
  post_dt[:,first_ele:last_ele] = np.transpose(post_delta_temp)

  #
  # RELATIVE HUMIDITY
  #
  #Read RH for Prior Posterior and Sonde
  prior_rh = nc_fid.variables['prior_rh'][:]
  post_rh = nc_fid.variables['post_rh'][:]
  sonde_rh = nc_fid.variables['sonde_rh'][:]

  #Get RH  differences
  prior_delta_rh = prior_rh - sonde_rh
  post_delta_rh = post_rh - sonde_rh

  prior_delta_rh.mask = np.copy(sonde_rh.mask)
  post_delta_rh.mask = np.copy(sonde_rh.mask)

  pre_drh[:,first_ele:last_ele] = np.transpose(prior_delta_rh)
  post_drh[:,first_ele:last_ele] = np.transpose(post_delta_rh)


  #
  # WATER VAPOUR 
  #
  #Read WV for Prior Posterior and Sonde
  prior_wv = nc_fid.variables['prior_water_vapour'][:]
  post_wv = nc_fid.variables['post_water_vapour'][:]
  sonde_wv = nc_fid.variables['sonde_water_vapour'][:]

  #Get WV  differences
  prior_delta_wv = prior_wv - sonde_wv
  post_delta_wv = post_wv - sonde_wv

  prior_delta_log_wv = np.log(prior_wv/1000) - np.log(sonde_wv/1000)
  post_delta_log_wv = np.log(post_wv/1000) - np.log(sonde_wv/1000)

  prior_delta_wv.mask = np.copy(sonde_wv.mask)
  post_delta_wv.mask = np.copy(sonde_wv.mask)

  prior_delta_log_wv.mask = np.copy(sonde_wv.mask)
  post_delta_log_wv.mask = np.copy(sonde_wv.mask)

  pre_dwv[:,first_ele:last_ele] = np.transpose(prior_delta_wv)
  post_dwv[:,first_ele:last_ele] = np.transpose(post_delta_wv)

  pre_dlogwv[:,first_ele:last_ele] = np.transpose(prior_delta_log_wv)
  post_dlogwv[:,first_ele:last_ele] = np.transpose(post_delta_log_wv)

  pre_dtdlogw = np.ma.append(pre_dt,pre_dlogwv,axis=0);
  post_dtdlogw = np.ma.append(post_dt,post_dlogwv,axis=0);

  cum_sum = cum_sum + np.int(fovs[index])

    

#Calculate mean and std of dt and dlogq over the levels  

mean_dt = np.ma.mean(pre_dt,axis=1)
std_dt = np.ma.std(pre_dt,axis=1)

max_m_dt = np.ma.max(mean_dt)
min_m_dt = np.ma.min(mean_dt)

#if abs(max_m_dt) > abs(min_m_dt):
#    ext_m_dt = max_m_dt
#else:
#    ext_m_dt = min_m_dt

#ext_s_dt = np.ma.max(std_dt)

midx = np.ma.notmasked_edges(mean_dt)
sidx = np.ma.notmasked_edges(std_dt)

#Fill mean and std missing values above last level of rawinsonde for T
ext_m_dt = mean_dt[midx[1]]
ext_s_dt = std_dt[sidx[1]]

print('mean_dt: {}'.format(mean_dt))
print('ext_m_dt: {}'.format(ext_m_dt))
print('ext_s_dt: {}'.format(ext_s_dt))

mean_dt[midx[1]]=ext_m_dt*2.0
std_dt[sidx[1]]=ext_s_dt*2.0
for i in np.arange(len(mean_dt)-midx[1]-1):
    mean_dt[midx[1]+i+1] = np.ma.filled(mean_dt[midx[1]+i+1].astype(float), mean_dt[midx[1]+i]*1.05)
for i in np.arange(len(std_dt)-sidx[1]-1):
    std_dt[sidx[1]+i+1] = np.ma.filled(std_dt[sidx[1]+i+1].astype(float), std_dt[sidx[1]+i]*1.05)

for i in np.arange(len(mean_dt)-1):
    if abs(mean_dt[i+1])>abs(mean_dt[i]*1.3) and abs(mean_dt[i+1])>abs(mean_dt[i+2]*1.3):
        mean_dt[i+1] = (mean_dt[i]+mean_dt[i+2])/2.0
        print('FLAG MEAN {}'.format(mean_dt[i+1]))

for i in np.arange(len(std_dt)-1):
    if std_dt[i+1]>std_dt[i]*1.3 and std_dt[i+1]>std_dt[i+2]*1.3:
        std_dt[i+1] = (std_dt[i]+std_dt[i+2])/2.0
        print('FLAG STD {}'.format(std_dt[i+1]))

#print(mean_dt)
#print(std_dt)
print('mean_dt: {}'.format(mean_dt))
print('std_dt: {}'.format(std_dt))

mean_dlogwv = np.ma.mean(pre_dlogwv,axis=1)
std_dlogwv = np.ma.std(pre_dlogwv,axis=1)

max_m_dlogwv = np.ma.max(mean_dlogwv)
min_m_dlogwv = np.ma.min(mean_dlogwv)

print('mean_dlogwv: {}'.format(mean_dlogwv))
print('std_dlogwv: {}'.format(std_dlogwv))

#if abs(max_m_dlogwv) > abs(min_m_dlogwv):
#    ext_m_dlogwv = max_m_dlogwv
#else:
#    ext_m_dlogwv = min_m_dlogwv
#
#
#ext_s_dlogwv = np.ma.max(std_dlogwv)


midx = np.ma.notmasked_edges(mean_dlogwv)
sidx = np.ma.notmasked_edges(std_dlogwv)

#Fill mean and std missing values above last level of rawinsonde for LOGQ
ext_m_dlogwv = mean_dlogwv[midx[-1]]
ext_s_dlogwv = std_dlogwv[midx[-1]]


#Fill mean and std missing values above last level of rawinsonde for LOGQ
mean_dlogwv[midx[1]]=ext_m_dlogwv*2.0
std_dlogwv[sidx[1]]=ext_s_dlogwv*2.0
for i in np.arange(len(mean_dlogwv)-midx[1]-1):
    mean_dlogwv[midx[1]+i+1] = np.ma.filled(mean_dlogwv[midx[1]+i+1].astype(float), mean_dlogwv[midx[1]+i]*1.05)
for i in np.arange(len(std_dlogwv)-sidx[1]-1):
    std_dlogwv[sidx[1]+i+1] = np.ma.filled(std_dlogwv[sidx[1]+i+1].astype(float), std_dlogwv[sidx[1]+i]*1.05)

for i in np.arange(len(mean_dlogwv)-1):
    if abs(mean_dlogwv[i+1])>abs(mean_dlogwv[i]*1.3) and abs(mean_dlogwv[i+1])>abs(mean_dlogwv[i+2]*1.3):
        mean_dlogwv[i+1] = (mean_dlogwv[i]+mean_dlogwv[i+2])/2.0
        print('FLAG MEAN {}'.format(mean_dlogwv[i+1]))

for i in np.arange(len(std_dlogwv)-1):
    if std_dlogwv[i+1]>std_dlogwv[i]*1.3 and std_dlogwv[i+1]>std_dlogwv[i+2]*1.3:
        std_dlogwv[i+1] = (std_dlogwv[i]+std_dlogwv[i+2])/2.0
        print('FLAG STD {}'.format(std_dlogwv[i+1]))


#Now use extrapolated mean and std to fill the difference above last available level of rawinsonde
#Values are filled with random gaussian valus with given mean and std
print('Sahpe of pre_dt: {}'.format(np.shape(pre_dt)))
print('Sahpe of mean_dt: {}'.format(np.shape(mean_dt)))


for i in np.arange(len(mean_dt)):
    idx = np.where(pre_dt[i,:].mask==True)
    for j in idx[0]:
        #print('j: {}'.format(j))
        pre_dt[i,j] =  np.random.normal(mean_dt[i],std_dt[i],1) 
        if pre_dt[i,j] > 10.0:
            pre_dt[i,j] = pre_dt[i,j]/2.0
            #print('FLAG pre_dt[{}{}]: {}'.format(i,j,pre_dt[i,j]))
for i in np.arange(len(mean_dlogwv)):
    idx = np.where(pre_dlogwv[i,:].mask==True)
    #print('idx: {}'.format(idx))
    for j in idx[0]:
        pre_dlogwv[i,j] = np.random.normal(mean_dlogwv[i],std_dlogwv[i],1)

fac=2.0
for i in np.arange(nlev-2): 
    for j in np.arange(nfovs):
        if abs(pre_dt[i+1,j])>abs(pre_dt[i,j]*fac) and abs(pre_dt[i+1,j])>abs(pre_dt[i+2,j]*fac):
            print('PRE T: {}'.format(pre_dt[i+1,j]))
            pre_dt[i+1,j]=(pre_dt[i,j]+pre_dt[i+2,j])/2
            print('FIXING T: {}'.format(pre_dt[i+1,j]))
        if abs(pre_dlogwv[i+1,j])>abs(pre_dlogwv[i,j]*fac) and abs(pre_dlogwv[i+1,j])>abs(pre_dlogwv[i+2,j]*fac):
            print('PRE LOGWV: {}'.format( pre_dlogwv[i+1,j]))
            pre_dlogwv[i+1,j]=(pre_dlogwv[i,j]+pre_dlogwv[i+2,j])/2
            print('FIXING LOGWV: {}'.format( pre_dlogwv[i+1,j]))


#Append DLOGWV matrix to DT matrix so that the covariance can include the cross terms  
#the new matrix to be used to generate the covaraince is pre_dtdlogwv
pre_dtdlogw = np.ma.append(pre_dt,pre_dlogwv,axis=0); 

#Get the cross terms T LOGWV
C_pre = np.ma.cov(pre_dtdlogw)

C_T = C_pre[:nlev,:nlev]
print("C_T: {}".format(np.diagonal(C_T)))
C_LOGQ = C_pre[nlev:,nlev:]
print("C_LOGQ: {}".format(np.diagonal(C_LOGQ)))

C_T_q = C_pre[:nlev,nlev:]

test = numpy.allclose(C_T, C_T.T, rtol=1e-05, atol=1e-08)
print('Simmetric test: {}'.format(test))
det = np.linalg.det(C_T) 
print('Determinant C(T): {}'.format(det))
test = numpy.allclose(C_LOGQ, C_LOGQ.T, rtol=1e-05, atol=1e-08)
print('Simmetric test: {}'.format(test))
det = np.linalg.det(C_LOGQ) 
print('Determinant C(C_LOGQ): {}'.format(det))
test = numpy.allclose(C_T_q, C_T_q.T, rtol=1e-05, atol=1e-08)
print('Simmetric test: {}'.format(test))
det = np.linalg.det(C_T_q) 
print('Determinant C(C_T_q): {}'.format(det))




#print(np.diag(C_T))
#print(np.diag(C_LOGQ))
#print(np.diag(C_T_q))

print('Reading Ozone Covariance')
ozone_static_covariance = '/home/mirto/amethyst/ancillary/atmosphere/static_ozone.nc'
#ozone_static_covariance = '/work/cris/mirto_lsr/fixed/ozone.nc'
with Dataset(ozone_static_covariance,'r') as ozone_f:
    ozone=ozone_f.variables['ozone'][:]
    #ozone_f.close()

#Write output file
mode = 'w'
enable_cmp = False
cmp_level = 4


os.chdir("/home/adaptive")
with Dataset(area+'_apriori.nc', mode) as output_f:

    output_f.set_auto_mask(False)

    numobs = 1 
    output_f.createDimension(NFOVS, numobs)

    lats_var = output_f.createVariable(LATITUDE,
                                      'f4',
                                      (NFOVS,),
                                      zlib=True,
                                      fill_value=1e9
                                      )
    lats_var[:] = 25.0 


    lons_var = output_f.createVariable(LONGITUDE,
                                      'f4',
                                      (NFOVS,),
                                      zlib=True,
                                      fill_value=1e9
                                      )
    lons_var[:] = -160.0 

    #create the atmospheric group
    atm = output_f.createGroup(ATMGROUP)
    if cut_pressure_lev:
        output_f.createDimension(NATMLEVELS, n_lev_cut)
    else:
        output_f.createDimension(NATMLEVELS, nlev)

    # PaoloS 15/01/24
    #save average pressure grid
    first_dim_chunk_size = min(50, numobs)
    nc_pressure = atm.createVariable( 'p', 
                                      'f4',
                                       (NATMLEVELS,),
                                      fill_value = 0.,
                                      zlib=enable_cmp,
                                      complevel=cmp_level)
    all_pressures = np.ma.masked_invalid(all_pressures)
    all_pressures = np.ma.masked_where(all_pressures > 1e4,all_pressures)
    if cut_pressure_lev:
        nc_pressure[:] = np.ma.mean(all_pressures,axis=0)[:n_lev_cut]
    else:
        nc_pressure[:] = np.ma.mean(all_pressures,axis=0)
    nc_pressure[0]    = 1100 # Fixed first level
    nc_pressure.units = 'hPa'

    #create a COVGROUP inside the atmospheric one
    cov_group = atm.createGroup(COVGROUP)


    mols = ['T', 'q', 'O3', 'T_q']

    output_tables = {}
    for mol in mols:
       print('Saving table {} in the output file'.format(mol))
       table_name = mol
       table_type = 'f4'  # floating point, 32 bytes
       table_dims = (NFOVS, NATMLEVELS, NATMLEVELS)
       first_dim_chunk_size = min(50, numobs)
       if cut_pressure_lev:
         chunks = (first_dim_chunk_size, n_lev_cut, n_lev_cut)
       else:
         chunks = (first_dim_chunk_size, nlev, nlev)
       new_table = cov_group.createVariable(
                                            table_name,
                                            table_type,
                                            table_dims,
                                            fill_value=1e9,
                                            zlib=enable_cmp,
                                            chunksizes=chunks,
                                            complevel=cmp_level
                                            )
       output_tables[mol] = new_table

    if cut_pressure_lev:
        output_tables['T'][0,:,:] = C_T[:n_lev_cut,:n_lev_cut]
        output_tables['q'][0,:,:] = C_LOGQ[:n_lev_cut,:n_lev_cut]
        output_tables['O3'][0,:,:] =  ozone[:n_lev_cut,:n_lev_cut]
        output_tables['T_q'][0,:,:] =  C_T_q[:n_lev_cut,:n_lev_cut]
    else:
        output_tables['T'][0,:,:] = C_T
        output_tables['q'][0,:,:] = C_LOGQ
        output_tables['O3'][0,:,:] =  ozone
        output_tables['T_q'][0,:,:] =  C_T_q

    #output_f.close()

# %%
