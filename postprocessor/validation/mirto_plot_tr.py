import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from argparse import ArgumentParser
import sys,os
import matplotlib.patches as mpatches

DELTA_TIME  = 75
TIF_N_HOURS = 4

def extract_date(string):
    try:
        return datetime.strptime(string[:15],"%Y%m%d_%H%M%S")
    except:
        datestring = "_".join(string.split('_')[4:6])
        return datetime.strptime(datestring,"%Y%m%d_%H%M%S")
        

def temporal_check(filelist,filedate):
    
    keep = []
    for file_index,file in enumerate(filelist):
        if abs((extract_date(file) - filedate).total_seconds())//60 < DELTA_TIME:
            keep.append(file_index)
    keep = np.array(keep)
    if keep.size != 0:
        return filelist[keep]
    else: 
        return []

def plot_da_dehaan_test(basdset,trdset,argv,filetype):

    size_of_state_vector = trdset.dimensions['size_of_reduced_state_vector'].size
    n_coef   = size_of_state_vector//2
    eigen = trdset['scaled_eigenvalues'][:]
    n_eigen  = np.where( eigen <= 1 )
    m_eigen = eigen.mean(axis=0)
    
    #Read Observation operator from TR file 
    DA_Hret = trdset['scaled_observation_operator'][:,:,:]
    #print('DA_Hret T: {}'.format(DA_Hret[0,0,0:81]))
    #print('DA_Hret WV: {}'.format(DA_Hret[0,0,81:]))
    #print('DA_Hret shape: {}'.format(DA_Hret.shape))
    
    #Read and concatenate  T and WV (log(q)) first guesses from BASIC file into x0 and xhat
    x0_T=basdset['prior_air_temperature'][:,:]
    x0_WV=np.log(basdset['prior_specific_humidity'][:,:])
    #print('x0_T shape: {}'.format(x0_T.shape))
    #print('x0_WV shape: {}'.format(x0_WV.shape))
    
    x0 = np.concatenate((x0_T,x0_WV), axis=1)
    #print('x0 shape: {}'.format(x0.shape))
    
    xhat_T=basdset['posterior_air_temperature'][:,:]
    xhat_WV=np.log(basdset['posterior_specific_humidity'][:,:])
    xhat = np.concatenate((xhat_T,xhat_WV), axis=1)
    #print('xhat: {}'.format(xhat))
    #print('xhat shape: {}'.format(xhat.shape))
    
    #Read Transformed Retrievals from TR file
    DA_Yret = trdset['scaled_observation'][:,:]
    #print('DA_Yret: {}'.format(DA_Yret))
    #print('DA_Yret shape: {}'.format(DA_Yret.shape))
    
    nfovs = DA_Yret.shape[0]
    #print('nfovs shape: {}'.format(nfovs))
    
    #Calculate De Haan quantities
    xhatdelta_list = []
    x0delta_list = []
    
    a = DA_Yret[0,:]-np.matmul(np.squeeze(DA_Hret[0,:,:]),xhat[0,:])
    #print('a: {}'.format(a))
     
    for i in range(nfovs):
        xhatdelta_list.append(DA_Yret[i,:]-np.matmul(np.squeeze(DA_Hret[i,:,:]),xhat[i,:]))
        x0delta_list.append(DA_Yret[i,:]-np.matmul(np.squeeze(DA_Hret[i,:,:]),x0[i,:]))
    
    xhatdelta = np.array(xhatdelta_list)
    x0delta = np.array(x0delta_list)
    #print('xhatdelta: {}'.format(xhatdelta))
    #print('x0delta: {}'.format(x0delta))
    #print('xhatdelta shape: {}'.format(xhatdelta.shape))
    #print('x0delta shape: {}'.format(x0delta.shape))
    
    xhat_res=np.sqrt(np.diag(np.cov(xhatdelta.T)))
    x0_res=np.sqrt(np.diag(np.cov(x0delta.T)))
    #print('xhat_res shape: {}'.format(xhat_res.shape))
    #print('x0_res shape: {}'.format(x0_res.shape))
    
    DA_Lambda = trdset['scaled_eigenvalues'][:,:]
    #print('DA_Lambda shape: {}'.format(DA_Lambda.shape))
    
    tmp = np.square(DA_Lambda)
    #print('tmp shape: {}'.format(tmp.shape))
    
    xhatlambda=1./np.sqrt(1+np.square(DA_Lambda))
    x0lambda=np.sqrt(1+np.square(DA_Lambda))
    #print('xhatlambda shape: {}'.format(xhatlambda.shape))
    #print('x0lambda shape: {}'.format(x0lambda.shape))
    
    xhat_mean_lambda=1./np.sqrt(1+np.square(DA_Lambda.mean(axis=0)))
    x0_mean_lambda=np.sqrt(1+np.square(DA_Lambda.mean(axis=0)))
    #print('xhat_mean_lambda shape: {}'.format(xhat_mean_lambda.shape))
    #print('x0_mean_lambda shape: {}'.format(x0_mean_lambda.shape))
    
    #PLOT XHAT
    
    plt.figure(figsize=(10,8))
    plt.clf()
    
    plt.title('{} De Haan Test Xhat'.format(filetype))
    n_eigen  = DA_Lambda.shape[1]
    #print('n_eigen: {}'.format(n_eigen))
    
    x = range(1,n_eigen+1)
    
    y = xhatlambda
    #print('y: {}'.format(y))
    #print('y shape: {}'.format(y.shape))
    
    y_m = abs(y).mean(axis=0)
    #print('y_m: {}'.format(y_m))
    #print('y_m shape: {}'.format(y_m.shape))
    
    y_std = abs(y).std(axis=0)
    #print('y_std: {}'.format(y_std))
    #print('y_std shape: {}'.format(y_std.shape))
    
    plt.plot(x,y_m,label=r'$\mu \left( {(1+\lambda^2)}^{- \frac{1}{2}} \right)$',marker='o',ms=3)
    plt.fill_between(x, y_m-y_std, y_m+y_std, color='grey', alpha=0.5,label=r'$\sigma \left( {(1+\lambda^2)}^{- \frac{1}{2}} \right)$')
    y = xhat_res
    #print('y: {}'.format(y))
    #print('y shape: {}'.format(y.shape))
    
    plt.plot(x,y,marker='o',ms=3,label=r"$\sigma \left( {y'}_{ret} - {H'}_{ret}\hat{x} \right)$")
    plt.yscale('log')
    plt.xticks(np.arange(1,n_eigen+1))
    
    plt.xlabel('Eigenvector Number',size=15)
    plt.grid()
    plt.grid(which='minor',ls='-',alpha=0.2)
    plt.legend(prop={'size':16})
    overpass = "".join( os.path.basename(argv.input).split('_')[:2])
    outfile = "/".join(  [ argv.outdir , 'TR_' + argv.world_area + '_' + overpass + '_de_haan_xhat.png'])

    if argv.eps:
        outfile = outfile.replace('.png','.eps')

    try:
        plt.savefig(outfile,dpi=300)
        print("Saved plot {}".format(outfile))
    except:
        print("Error in saving plot {}".format(outfile))

    plt.clf()

    #PLOT X0
    
    plt.figure(figsize=(10,8))
    plt.clf()
    
    plt.title('{} De Haan Test X0'.format(filetype))
    n_eigen  = DA_Lambda.shape[1]
    #print('n_eigen: {}'.format(n_eigen))
    
    x = range(1,n_eigen+1)
    
    y = x0lambda
    #print('y: {}'.format(y))
    #print('y shape: {}'.format(y.shape))
    
    y_m = abs(y).mean(axis=0)
    #print('y_m: {}'.format(y_m))
    #print('y_m shape: {}'.format(y_m.shape))
    
    y_std = abs(y).std(axis=0)
    #print('y_std: {}'.format(y_std))
    #print('y_std shape: {}'.format(y_std.shape))
    
    plt.plot(x,y_m,marker='o',ms=3,label=r'$\mu \left( {(1+\lambda^2)}^{- \frac{1}{2}} \right)$')
    plt.fill_between(x, y_m-y_std, y_m+y_std, color='grey', alpha=0.5,label=r'$\sigma \left( {(1+\lambda^2)}^{- \frac{1}{2}} \right)$')

    y = x0_res
    #print('y: {}'.format(y))
    #print('y shape: {}'.format(y.shape))
    
    plt.plot(x,y,marker='o',ms=3,label=r"$\sigma \left( {y'}_{ret} - {H'}_{ret} x_{0} \right)$")
    plt.yscale('log')
    plt.xticks(np.arange(1,n_eigen+1))
    
    plt.xlabel('Eigenvector Number',size=15)
    plt.grid()
    plt.grid(which='minor',ls='-',alpha=0.2)
    plt.legend(prop={'size':16})

    overpass = "".join( os.path.basename(argv.input).split('_')[:2])
    outfile = "/".join(  [ argv.outdir , 'TR_' + argv.world_area + '_' + overpass + '_de_haan_x0.png'])

    if argv.eps:
        outfile = outfile.replace('.png','.eps')

    try:
        plt.savefig(outfile,dpi=300)
        print("Saved plot {}".format(outfile))
    except:
        print("Error in saving plot {}".format(outfile))

    plt.clf()

    return

def plot_tr_eigenvalues(dataset,argv,filetype):

    eigen = dataset['scaled_eigenvalues'][:,:]
    #print('EIGENVALUES: {}'.format(eigen))
    m_eigen = eigen.mean(axis=0)
    n_eigen = eigen.shape[1]

    #Plot Eigenvalues
    print('n_eigen: {}'.format(n_eigen))
    plt.figure(figsize=(10,6))

    x = range(1,n_eigen+1)
    y = eigen

    plt.clf()
    plt.title('S2N eigenvalues')
    leg = []
    y_m = abs(y).mean(axis=0)
    y_std = abs(y).std(axis=0)
    plt.plot(x,y_m,marker='o',ms=3,label=r'$ \mu $')

    dif = y_m - y_std
    sum = y_m + y_std
    plt.fill_between(x,dif,sum, color='grey', alpha=0.5,label=r'$ \mu \pm \sigma$')

    #plt.plot(x,dif,c='red',ls='-')
    #plt.plot(x,sum,c='red',ls='-')
    #plt.yscale('log')
    plt.xticks(np.arange(1,n_eigen+1))

    #plt.ylim(10,1015)
    plt.ylabel('Eigenvalues',size=16)
    plt.xlabel('Eigenvector Number',size=15)
    plt.yscale('log')
    plt.grid()
    plt.legend(prop={'size':16})

    overpass = "".join( os.path.basename(argv.input).split('_')[:2])
    outfile = "/".join(  [ argv.outdir , 'EIGEN_' + argv.world_area + '_' + overpass + '_values.png'])

    if argv.eps:
        outfile = outfile.replace('.png','.eps')

    try:
        plt.savefig(outfile,dpi=300)
        print("Saved plot {}".format(outfile))
    except:
        print("Error in saving plot {}".format(outfile))

    plt.clf()

    return

    
        
def plot_transformed_retrievals(dataset,argv,filetype):

    eigen = dataset['scaled_eigenvalues'][:]
    m_eigen = eigen.mean(axis=0)


    # Plot temperature profile
    y = dataset['scaled_observation'][:,:]
    n_eigen  = y.shape[1] 
    print('n_eigen: {}'.format(n_eigen))
    plt.figure(figsize=(10,6))

    x = range(1,n_eigen+1)

    plt.clf()
    plt.title('Transformed Retrievals')
    leg = []
    y_m = abs(y).mean(axis=0)
    y_std = abs(y).std(axis=0)
    plt.plot(x,y_m,marker='o',ms=3,label=r'$ \mu $')
    
    # Debug plot
    #for i,ll in enumerate(y_m):
    #      print(i,x[i],ll,y_std[i],(y_m-y_std)[i],(y_m+y_std)[i])
    dif = y_m - y_std
    sum = y_m + y_std
    plt.fill_between(x,dif,sum, color='grey', alpha=0.5,label=r'$ \mu \pm \sigma$')

    #plt.plot(x,dif,c='red',ls='-')
    #plt.plot(x,sum,c='red',ls='-')
    #plt.yscale('log')
    plt.xticks(np.arange(1,n_eigen+1))

    #plt.ylim(10,1015)
    plt.ylabel('TR',size=16)
    plt.xlabel('Eigenvector Number',size=15)
    plt.grid()
    plt.legend(prop={'size':16})

    overpass = "".join( os.path.basename(argv.input).split('_')[:2])
    outfile = "/".join(  [ argv.outdir , 'TR_' + argv.world_area + '_' + overpass + '_values.png'])

    if argv.eps:
        outfile = outfile.replace('.png','.eps')

    try:
        plt.savefig(outfile,dpi=300)
        print("Saved plot {}".format(outfile))
    except:
        print("Error in saving plot {}".format(outfile))

    plt.clf()

    return


def plot_scaled_operator(dataset,argv,filetype):

    size_of_state_vector = dataset.dimensions['size_of_reduced_state_vector'].size
    n_coef   = size_of_state_vector//2
    eigen = dataset['scaled_eigenvalues'][:]
    n_eigen  = np.where( eigen <= 1 )
    m_eigen = eigen.mean(axis=0)
    
    # Plot temperature profile
    x = np.ma.masked_invalid(dataset['scaled_observation_operator'][:,:,:n_coef])
    for ind in zip(n_eigen[0],n_eigen[1]):
        x[ind,:].mask = [ True for l in range(n_coef)]
    plt.figure(figsize=(10,8))

    y = dataset['air_pressure'][:,:]
    plt.clf()
    plt.title('{} Observation Operator for Temperature'.format(filetype))
    leg = []
    for i in range(x.shape[1]):
        sleg = 'EIG({}): {:.1f}'.format(i+1,m_eigen[i])
        plt.plot(x[:,i,:].mean(axis=0),y.mean(axis=0),label=i+1)
        if m_eigen[i] >= 1: 
          leg.append(sleg)
    plt.yscale('log')

    plt.legend(leg)

    #plt.ylim(10,1015)
    plt.ylabel('Pressure (hPa)',size=15)
    plt.xlabel('Average Observation Operator for Water Vapor',size=15)
    plt.grid()
    plt.gca().invert_yaxis()

    overpass = "".join( os.path.basename(argv.input).split('_')[:2])
    outfile = "/".join(  [ argv.outdir , 'TR_' + argv.world_area + '_' + overpass + '_obs_operator_temp.png'])

    if argv.eps:
        outfile = outfile.replace('.png','.eps')

    try:
        plt.savefig(outfile,dpi=300)
        print("Saved plot {}".format(outfile))
    except:
        print("Error in saving plot {}".format(outfile))

    plt.clf()

    plt.figure(figsize=(10,8))
    plt.title('{} Observation Operator for Water Vapour'.format(filetype))
    x = np.ma.masked_invalid(dataset['scaled_observation_operator'][:,:,n_coef:])
    for ind in zip(n_eigen[0],n_eigen[1]):
        x[ind,:].mask = [ True for l in range(n_coef)]
    for i in range(x.shape[1]):
        plt.plot(x[:,i,:].mean(axis=0),y.mean(axis=0),label=i+1)
    plt.yscale('log')
    plt.yscale('log')
    plt.ylabel('Pressure (hPa)',size=15)
    plt.xlabel('Average Observation Operator',size=15)
    plt.grid()
    plt.gca().invert_yaxis()
    plt.legend(leg)
    overpass = "".join( os.path.basename(argv.input).split('_')[:2])
    outfile = "/".join(  [ argv.outdir , 'TR_' + argv.world_area + '_' + overpass + '_obs_operator_wv.png'])

    if argv.eps:
        outfile = outfile.replace('.png','.eps')

    try:
        plt.savefig(outfile,dpi=300)
        print("Saved plot {}".format(outfile))
    except:
        print("Error in saving plot {}".format(outfile))
    
    return

def plot_coordinates(lat,lon,colour,title):
    from mpl_toolkits.basemap import Basemap

    lat_max = lat.max()
    lon_max = lon.max()
    lon_min = lon.min()
    lat_min = lat.min()
    
    x,y = np.meshgrid(lon,lat)  
    m = Basemap(projection='merc',lat_ts=20,
                llcrnrlat=lat_min,urcrnrlat=lat_max,\
                llcrnrlon=lon_min,urcrnrlon=lon_max)
    
    m.drawcoastlines()
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='aqua')
    x,y = m(x,y)
    plt.scatter(x, y, 1,c=colour,marker='+')    
    
    plt.title(title)
    
def sounder_mapper_on_imager_background(tiffile,tr_file,basic_file,argv):
    
    from osgeo import gdal
    import cartopy.crs as ccrs
    try:
        ds=gdal.Open(tiffile) # Open GeoTiff File
        if ds == None:
            raise IOError
        tif_background = True
    except Exception as e:
        log.info("Cannot open geoTiff file: {}".format(e))                 
        tif_background = False

    if tif_background:
        print("Geo Description:")
        print(ds.GetDescription())
    
        data = ds.ReadAsArray()
        gt   = ds.GetGeoTransform()
        proj = ds.GetProjection()
        print('Geotif projection:')
        print('{}'.format(proj))
        print("")

    
    #Proj4 string used to generate VIIRS truecolor image
    #proj4str='+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'   
 
    #Cartopy CRS definition based on proj4 parameters used to generate VIIRS truecolor image
    #defined in: /home/oper/cris/cspp/polar2grid_v_2_3/lib/python3.7/site-packages/polar2grid/grids/grids.conf
    
    #For the Arctic:
    #polar_alaska, proj4, +proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lat_ts=60.0 +lon_0=-150 +units=m, 400,-400
    #crs = ccrs.Stereographic(central_latitude=90.0, central_longitude=-150.0,
    #                 false_easting=0.0, false_northing=0.0,
    #                 true_scale_latitude=60, globe=None)
    
        
    # Define Projection
    if argv.world_area == 'arctic':
        print('AREA: ARCTIC')
        myccrs=  ccrs.Stereographic(central_latitude=90.0, central_longitude=-150.0,
                                    false_easting=0.0,     false_northing=0.0,
                                    true_scale_latitude=60)
    elif argv.world_area == 'pacific':
        print('AREA: PACIFIC')
        globe = ccrs.Globe(semimajor_axis=6378137, flattening=1/298.257223563)
        myccrs=ccrs.PlateCarree(central_longitude=-160.0, globe=globe)           

    # Draw Map wit coast lines
    plt.figure(figsize=(15, 15))
    ax = plt.axes(projection=myccrs)
    gc=ax.coastlines(resolution='10m',linewidth=.75,color='orange')
    #gl=ax.gridlines(linewidth=.5,color='black')
    gl=ax.gridlines(draw_labels=True, linewidth=.75, color='black')

    gl.bottom_labels = False
    gl.left_labels = False

    tr_lon = tr_file['longitude'][:] 
    tr_lat = tr_file['latitude'][:]
      
    basic_lon = basic_file['longitude']
    basic_lat = basic_file['latitude']

    if tif_background:

        # Define Map extention
        extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
                  gt[3] + ds.RasterYSize * gt[5], gt[3])
        print('Geografic extent: {}'.format(extent))
    
        if argv.world_area == 'pacific':
            a , b = myccrs.transform_point(extent[0], extent[2], ccrs.Geodetic())
            c , d = myccrs.transform_point(extent[1], extent[3], ccrs.Geodetic())
            extent = (a, c, b, d)
            print('Projected extent: {}'.format(extent))
    
    
        print('Data size in pixels: {}'.format(data.shape))
        if data.ndim == 3:
            img = ax.imshow(data[:3, :, :].transpose((1, 2, 0)),extent=extent,origin='upper')
        else:
            img = ax.imshow(data,extent=extent,cmap='gray',origin='upper')
        
    else:
        
        xmin, xmax = basic_lon.min(), basic_lon.max()
        ymin, ymax = basic_lat.min(), basic_lat.max()

        ax.set_extent([ xmin + 12.5, xmax-12.5, ymin+12.5,ymax-12.5],crs=ccrs.PlateCarree())
        ax.background_img()

    
    #Get FOV parameters using satellite altitude and instrument apreture
    sat_ape=[0.963, 0.897] # degrees
    sat_alt = 824; # km
    earth_radius = 6371 # km
    tand = lambda x: np.tan(np.radians(x))

    cmap = None
    normalize = None
    fov_instrument_color = { 'iasi' : 'tab:blue',
                             'cris' : 'salmon'}
    # Plot FOVs ellipses
    for fov_index, xlon in enumerate(tr_lon):
       
        # Compute axis and angle
        # Double check if axis are semi-axis or full-axis
        L1=( tand(sat_ape[0]+tr_file['field_of_view'][fov_index]) - tand(tr_file['field_of_view'][fov_index]))*sat_alt
        L2= tand(sat_ape[1])*sat_alt
                     
        if not tr_file['azimuth_angle'][fov_index].mask:
           angle = -90 - tr_file['azimuth_angle'][fov_index]
        else:
           angle = -90 
        
        #Transform lat and lon in cartesian coordinates
        projx1, projy1 = myccrs.transform_point(xlon, tr_lat[fov_index], ccrs.Geodetic())
        
        # Define Ellipse Face and Edge Colour
        color = None
        edgecolor = fov_instrument_color[ tr_file['instrument'][fov_index]  ]
        facecolor = edgecolor
        zorder    = 1 
        # Draw the patch
        ax.add_patch(      mpatches.Ellipse(   xy=[projx1, projy1], 
                                               width=L1*1000, 
                                               height=L2*1000, 
                                               angle=angle, 
                                               edgecolor= edgecolor,
                                               alpha=1,
                                               facecolor= facecolor, 
                                               transform=myccrs,
                                               linewidth=0.4,
                                               zorder=zorder)
                     )  
        
    # Plot FOVs ellipses
    for fov_index, xlon in enumerate(basic_lon):
       
        # Compute axis and angle
        # Double check if axis are semi-axis or full-axis
        L1=( tand(sat_ape[0]+basic_file['field_of_view'][fov_index]) - tand(basic_file['field_of_view'][fov_index]))*sat_alt
        L2= tand(sat_ape[1])*sat_alt
                     
        if not basic_file['azimuth_angle'][fov_index].mask:
           angle = -90 - basic_file['azimuth_angle'][fov_index]
        else:
           angle = -90 

        #Transform lat and lon in cartesian coordinates
        projx1, projy1 = myccrs.transform_point(xlon, basic_lat[fov_index], ccrs.Geodetic())
        
        # Define Ellipse Face and Edge Colour
        color = None
        edgecolor = fov_instrument_color[ basic_file['instrument'][fov_index]  ]
        facecolor = edgecolor
        zorder    = 1 
        # Draw the patch
        ax.add_patch(      mpatches.Ellipse(   xy=[projx1, projy1], 
                                               width=L1*1000, 
                                               height=L2*1000, 
                                               angle=angle, 
                                               edgecolor= edgecolor,
                                               alpha=1,
                                               facecolor= facecolor, 
                                               transform=myccrs,
                                               linewidth=0.4,
                                               zorder=zorder)
                     ) 
        """
        if fov_index not in good_indx and plot_var != None:
            ax.scatter(projx1, projy1,marker='X',c='red',s=2,zorder=11)
        """
        
    #<----- cycle end
    c = [ mpatches.Circle((0.5, 0.5), radius = 0.25, facecolor='tab:blue', edgecolor="none" ) ,
          mpatches.Circle((0.5, 0.5), radius = 0.25, facecolor='salmon',   edgecolor="none" )
         ]
    plt.legend(c,['IASI','CrIS'],prop={'size':16})
    overpass = "".join( os.path.basename(argv.input).split('_')[:2])
    outfile = "/".join(  [ argv.outdir , 'TR_' + argv.world_area + '_' + overpass + '_mapplot.png'])
    
    if argv.eps:
        outfile = outfile.replace('.png','.eps')
    
    #Note that plt.show opens a new figure. If image has to be saved, savefig command has to be issued before plt.show.
    #plt.legend(prop={'size':17})
    try:
        plt.savefig(outfile,dpi=300)
        print("Map saved: {}".format(outfile))
    except Exception as e:
        print("Error in saving map plot: {}".format(e))

    return

if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', '-i', default='./',
                        help='Input TR or BASIC file')    
    parser.add_argument('--tiff', '-t', default=None, required = False,
                        help='Input TIFF file')
    parser.add_argument('-wa','--world_area',type=str, required=True,choices=['arctic','pacific'],
                        help='World Area')
    parser.add_argument('-o','--outdir',type=str, required=False,default='./',
                        help='Output dir')
    parser.add_argument('--eps',type=bool, required=False,default=False,
                        help='EPS Mode')
    
    argv = parser.parse_args()
    
    if 'TR' in argv.input:
        tr_file = Dataset(argv.input,'r')
        basic_file = Dataset(argv.input.replace('TR','BASIC'),'r')
    else:
        basic_file = Dataset(argv.input,'r')
        tr_file    = Dataset(argv.input.replace('BASIC','TR'), 'r' )
    
    sounder_mapper_on_imager_background(argv.tiff, tr_file, basic_file, argv)
    plot_scaled_operator(tr_file,argv, 'TR')
    plot_transformed_retrievals(tr_file,argv, 'TR')
    plot_da_dehaan_test(basic_file,tr_file,argv,'TR')
    plot_tr_eigenvalues(tr_file,argv,'TR')
    #plot_scaled_operator(basic_file,argv.outdir,'BASIC')

    
    # Plot with blue markers TR or BASIC FOVs
    #plot_coordinates(data['latitude'][:],data['longitude'][:],'blue','DATA')
    
    
    
