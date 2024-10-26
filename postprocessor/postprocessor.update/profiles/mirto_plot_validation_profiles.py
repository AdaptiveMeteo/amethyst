import matplotlib.pyplot as plt
from matplotlib import cm
from data_reader.read_netcdf import netCDFReader
import glob
import logging
import os , sys
from matplotlib import ticker
from matplotlib.ticker import FormatStrFormatter
import numpy as np

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


def mirto_plot_validation_profiles(x1,x2,x3,y,var,sup_title_str=None):

   fig, ax = plt.subplots(1,2)
   fig.set_size_inches(16.5, 10.5)

   if var in 'T':
      xlabel_str_0 = 'Air Temperature [K]'
      xlabel_str_1 = 'Delta Air Temperature [K]'
      ylabel_str = 'p [hPa]'
      title_str =  'Temperature'
      ylim = [1014, 0.1]
   elif var in 'RH':
      xlabel_str_0 = 'Air Relative Humidity [%]'
      xlabel_str_1 = 'Delta Air Relative Humidity [%]'
      ylabel_str = 'p [hPa]'
      title_str =  'Relative Humidity'
      ylim = [1014, 200]
   elif var in 'WV':
      xlabel_str_0 = 'Air Water Vapor MR [g/kg]'
      xlabel_str_1 = 'Delta Air Water Vapor MR [g/kg]'
      ylabel_str = 'p [hPa]'
      title_str =  'Water Vapor MR'
      ylim = [1014, 200]
   
   ax[0].plot(x1, y, 'b+', x2, y, 'ro', x3, y,'gx', linestyle='solid')
   ax[0].set(xlabel = xlabel_str_0, ylabel = ylabel_str)
   ax[0].set_title(title_str, pad = 12, fontsize = 14)
   ax[0].set_ylim(ylim)
   
   lgnd = ax[0].legend(['wrf fg', 'mirto ret', 'sonde'], loc='upper right', shadow=True)
   ax[0].grid(True)
   ax[0].set(yscale = 'log')

   d1 = x1 - x3 
   d2 = x2 - x3
   
   ax[1].plot(d1, y, 'b+', d2, y, 'ro', linestyle='solid')
   ax[1].set(xlabel = xlabel_str_1, ylabel = ylabel_str)
   ax[1].set_title(title_str, pad = 12, fontsize = 14)
   lgnd = ax[0].legend(['wrf fg', 'mirto ret', 'sonde'], loc='upper right', shadow=True)
   ax[1].set_ylim(ylim)
   lgnd = ax[1].legend(['wrf fg - sonde', 'mirto ret - sonde'], loc='upper right', shadow=True)
   ax[1].grid(True)
   ax[1].set(yscale = 'log')
   fig.suptitle(sup_title_str, fontsize = 14)

   #plt.show()

   return fig

def plot_wrf_validation_profiles(x1,x2,y,var,sup_title_str=None,rms=None,plot_mode='aggregated'):

   fig, ax = plt.subplots(1,2)
   fig.set_size_inches(16.5, 10.5)
   if var in 'T':
      xlabel_str_0 = 'Air Temperature [K]'
      xlabel_str_1 = 'Delta Air Temperature [K]'
      ylabel_str = 'p [hPa]'
      title_str =  'Temperature'
      ylim = [1014, 40]
      xlim = [-3,3]
   elif var in 'RH':
      xlabel_str_0 = 'Air Relative Humidity [%]'
      xlabel_str_1 = 'Delta Air Relative Humidity [%]'
      ylabel_str = 'p [hPa]'
      title_str =  'Relative Humidity'
      ylim = [1014, 200]
      xlim = [-30,30]

   elif var in 'WV':
      xlabel_str_0 = 'Air Water Vapor MR [g/kg]'
      xlabel_str_1 = 'Delta Air Water Vapor MR [g/kg]'
      ylabel_str = 'p [hPa]'
      title_str =  'Water Vapor MR'
      ylim = [1014, 200]
      xlim = [-3,3]
      
   ygrid = np.arange(100,1100,100)      

   if plot_mode == 'aggregated':   
       x1_mean = x1.mean(axis=0)
       x2_mean = x2.mean(axis=0)
       x1_std = x1.std(axis=0)
       x2_std = x2.std(axis=0)
       second_color = 'black'

   else:
       x1_mean = x1.reshape(x1.size)
       x2_mean = x2.reshape(x2.size)
       second_color = 'black'
   ax[0].plot(x1_mean, y,'b+',label = 'WRF ' + var,  
              linestyle='solid', marker='+')

   ax[0].plot(x2_mean, y,second_color,
              label = 'Sonde ' + var, linestyle = 'solid',marker='+')

   if plot_mode == 'aggregated':          
       ax[0].fill_betweenx(y, x1_mean - x1_std, x1_mean +x1_std,color='tab:blue',alpha=0.3,label=r'$ \mu \pm \sigma $' + '  (WRF)') 
       ax[0].fill_betweenx(y, x2_mean - x2_std, x2_mean +x1_std,color='red',alpha=0.3,label=r'$ \mu \pm \sigma $' + '  (Sonde)') 
   #ax[0].xaxis.ticklabel_format(axis='x',style='sci')   

   ax[0].set(yscale = 'log')   
   ax[0].set(xlabel = xlabel_str_0, ylabel = ylabel_str)
   ax[0].set_title(title_str, pad = 12, fontsize = 14)
   ax[0].set_ylim(ylim)
   lgnd = ax[0].legend(loc='upper right', shadow=True,prop={'size':20})
   ax[0].grid(True,which='both')
   ax[0].grid(True,which='both',c='grey',ls='--',alpha=0.5)
   #ax[0].set_yticklabels(ygrid[::-1].astype(str))
   ax[0].set_yticks(ygrid)
   ax[0].set_yticklabels(ygrid.astype(str))


   d = x1 - x2 
   if plot_mode == 'aggregated':   
       mean_d = d.mean(axis=0)
       std_d  = d.std(axis=0)
       rcolor = 'black'

   else:
       rcolor = 'red'
       mean_d = d.reshape(x1.size)
   ax[1].plot(mean_d, y, rcolor, linestyle='solid',marker='+',label='WRF - Sonde')

   if plot_mode == 'aggregated':   
       ax[1].fill_betweenx(y, mean_d - std_d , mean_d +std_d,color='grey',alpha=0.5,label=r'$ \mu \pm \sigma $') 
       
   ax[1].set(xlabel = xlabel_str_1, ylabel = ylabel_str)
   ax[1].set_title(title_str, pad = 12, fontsize = 14)
   lgnd = ax[0].legend(loc='upper right', shadow=True)
   ax[1].set_ylim(ylim)
   ax[1].grid(True,which='both')
   ax[1].grid(True,which='both',c='grey',ls='--',alpha=0.5)
   ax[1].set(yscale = 'log')
   
   if rms is not None:
       if plot_mode == 'aggregated':   
           ax[1].plot(rms.mean(axis=0),y,marker='o',label='RMS',c='red')
   fig.suptitle(sup_title_str, fontsize = 14)
   ax[1].set_yticks(ygrid)
   ax[1].set_yticklabels(ygrid.astype(str))
   ax[1].set_xlim(xlim)

   lgnd = ax[1].legend(loc='upper right', shadow=True,prop={'size':14})
       
   return fig

#PaoloA 31Aug2022
#def plot_wrf_comparison_profile_diff(ctrl_dif,oper_dif,foper_dif,y,var,sup_title_str=None):
def plot_wrf_comparison_profile_diff(oper_dif,foper_dif,y,var,sup_title_str=None):

   fig, ax = plt.subplots(1,1)
   fig.set_size_inches(16.5, 10.5)

   if var in 'T':
      xlabel_str_0 = 'Air Temperature [K]'
      xlabel_str_1 = 'Delta Air Temperature [K]'
      ylabel_str = 'p [hPa]'
      title_str =  'Temperature'
      ylim = [1014, 40]
      xlim = [-3,3]
      xgrid = np.arange(-3,4,1)

   elif var in 'RH':
      xlabel_str_0 = 'Air Relative Humidity [%]'
      xlabel_str_1 = 'Delta Air Relative Humidity [%]'
      ylabel_str = 'p [hPa]'
      title_str =  'Relative Humidity'
      ylim = [1014, 200]
      xlim = [-30,30]
      xgrid = np.arange(-30,40,10)

   elif var in 'WV':
      xlabel_str_0 = 'Air Water Vapor MR [g/kg]'
      xlabel_str_1 = 'Delta Air Water Vapor MR [g/kg]'
      ylabel_str = 'p [hPa]'
      title_str =  'Water Vapor MR'
      ylim = [1014, 200]
      xlim = [-3,3]
      xgrid = np.arange(-3,4,1)
   ygrid = np.arange(100,1100,100)      

   #PaoloA 31Aug2022
   # Plot CNTRL Diff
   #ctrl_std  = ctrl_dif.std(axis=0)
   #ctrl_mean = ctrl_dif.mean(axis=0)
   #ax.plot(ctrl_mean, y,'black',label = 'CNTRL - MW',  
   #                    linestyle='solid')

   #ax.fill_betweenx(y, ctrl_mean-ctrl_std , ctrl_mean+ctrl_std,
   #                    color='grey',alpha=0.5,label=r'$ \mu \pm \sigma $')
   
   # Plot OPER Diff
   oper_std  = oper_dif.std(axis=0)
   oper_mean = oper_dif.mean(axis=0)
   #PaoloA 31Aug2022
   #ax.plot(oper_mean, y,'blue',label = 'CNTRL - TR',  
   ax.plot(oper_mean, y,'blue',label = 'CONV - TR',  
                       linestyle='solid')
   #ax.fill_betweenx(y, oper_mean-oper_std , oper_mean+oper_std,
   #                   color='tab:blue',alpha=0.5,label=r'$ \mu \pm \sigma $')

   # Plot F-OPER Diff
   foper_std  = foper_dif.std(axis=0)
   foper_mean = foper_dif.mean(axis=0)
   #PaoloA 31Aug2022
   #ax.plot(foper_mean, y,'red',label = 'CNTRL - (TR+MW)',  
   ax.plot(foper_mean, y,'red',label = 'CONV - (TR+MW)',  
                       linestyle='solid')
   #ax.fill_betweenx(y, foper_mean-foper_std , foper_mean+foper_std,
   #                    color='salmon',alpha=0.5,label=r'$ \mu \pm \sigma $')

   ax.set(xlabel = xlabel_str_1, ylabel = ylabel_str)
   
   ax.set_xlabel(xlabel_str_1,size=15)
   ax.set_ylabel(ylabel_str,size=15)
   
   ax.set_title(title_str, pad = 12, fontsize = 14)
   ax.set_ylim(ylim)
   ax.grid(True,which='both')
   ax.grid(True,which='both',c='grey',ls='--',alpha=0.5)
   ax.set(yscale = 'log')
   

   fig.suptitle(sup_title_str, fontsize = 14)
   ax.set_yticks(ygrid)
   ax.set_yticklabels(ygrid.astype(str))
   ax.set_xticks(xgrid)
   ax.set_xticklabels(xgrid.astype(str))
   ax.set_xlim(xlim)

   plt.legend(loc='upper right', shadow=True,prop={'size':14})
       
   return fig

#PaoloA 31Aug2022
#def plot_wrf_comparison_profiles(conv,mw,oper,foper,y,var,sup_title_str=None):
def plot_wrf_comparison_profiles(conv,oper,foper,y,var,sup_title_str=None):

   fig, ax = plt.subplots(1,1)
   fig.set_size_inches(16.5, 10.5)

   if var in 'T':
      xlabel_str_0 = 'Air Temperature [K]'
      xlabel_str_1 = 'Delta Air Temperature [K]'
      ylabel_str = 'p [hPa]'
      title_str =  'Temperature'
      ylim = [1014, 40]
      xlim = [0,3]
      xgrid = np.arange(-3,4,1)

   elif var in 'RH':
      xlabel_str_0 = 'Air Relative Humidity [%]'
      xlabel_str_1 = 'Delta Air Relative Humidity [%]'
      ylabel_str = 'p [hPa]'
      title_str =  'Relative Humidity'
      ylim = [1014, 200]
      xlim = [0,30]
      xgrid = np.arange(-30,40,10)

   elif var in 'WV':
      xlabel_str_0 = 'Air Water Vapor MR [g/kg]'
      xlabel_str_1 = 'Delta Air Water Vapor MR [g/kg]'
      ylabel_str = 'p [hPa]'
      title_str =  'Water Vapor MR'
      ylim = [1014, 200]
      xlim = [0,3]
      xgrid = np.arange(-3,4,1)
   ygrid = np.arange(100,1100,100)

   # Plot CNTRL 
   conv_std  = conv.std(axis=0)
   conv_mean = conv.mean(axis=0)
   ax.plot(conv_mean, y,'red',label = 'CONV',
                       linestyle='solid')
   #ax.fill_betweenx(y, ctrl_mean-ctrl_std , ctrl_mean+ctrl_std,
   #                    color='grey',alpha=0.5,label=r'$ \mu \pm \sigma $')

   # Plot MW , 
   #PaoloA 31Aug2022
   #mw_std  = mw.std(axis=0)
   #mw_mean = mw.mean(axis=0)
   #ax.plot(mw_mean, y,'cyan',label = 'CONV+MW',
   #                    linestyle='solid')


   #ax.fill_betweenx(y, ctrl_mean-ctrl_std , ctrl_mean+ctrl_std,
   #                    color='grey',alpha=0.5,label=r'$ \mu \pm \sigma $')

   # Plot OPER 
   oper_std  = oper.std(axis=0)
   oper_mean = oper.mean(axis=0)
   ax.plot(oper_mean, y,'blue',label = 'CONV+TR',
                       linestyle='solid')
   #ax.fill_betweenx(y, oper_mean-oper_std , oper_mean+oper_std,
   #                   color='tab:blue',alpha=0.5,label=r'$ \mu \pm \sigma $')

   # Plot F-OPER 
   foper_std  = foper.std(axis=0)
   foper_mean = foper.mean(axis=0)
   ax.plot(foper_mean, y,'green',label = 'CONV+TR+MW',
                       linestyle='solid')
   #ax.fill_betweenx(y, foper_mean-foper_std , foper_mean+foper_std,
   #                    color='salmon',alpha=0.5,label=r'$ \mu \pm \sigma $')

   ax.set(xlabel = xlabel_str_1, ylabel = ylabel_str)

   ax.set_xlabel(xlabel_str_1,size=15)
   ax.set_ylabel(ylabel_str,size=15)

   ax.set_title(title_str, pad = 12, fontsize = 14)
   ax.set_ylim(ylim)
   ax.grid(True,which='both')
   ax.grid(True,which='both',c='grey',ls='--',alpha=0.5)
   ax.set(yscale = 'log')


   fig.suptitle(sup_title_str, fontsize = 14)
   ax.set_yticks(ygrid)
   ax.set_yticklabels(ygrid.astype(str))
   ax.set_xticks(xgrid)
   ax.set_xticklabels(xgrid.astype(str))
   ax.set_xlim(xlim)

   plt.legend(loc='upper right', shadow=True,prop={'size':14})

   return fig
    
def prepare_parser():
    from argparse import ArgumentParser

    
    parser = ArgumentParser()
    parser.add_argument('--outdir','-o',default='/data/mirto/db/aggregated',type=str,
                        required=True, help='Output directory for TR and BASIC file')
    parser.add_argument('--input','-i',default='plot_overpass',type=str,required=True,
                        help='Plot mode')
    parser.add_argument('--filetype','-ft',
                        type=str,default='aggregated',choices=['aggregated','single'],
                        help='Validation file type')
    parser.add_argument('--mode','-m',
                        type=str,default='wrf',choices=['wrf','retrieval','wrf_comparison'],
                        help='System type')
    parser.add_argument('--forecast_hour',
                    type=int,default=-1, help='WRF Forecast hour')
    parser.add_argument('--wrf_control_dir',required=False,default=None, type=str,
                        help='Dir with Control files for wrf_comparison mode'  )
    parser.add_argument('--wrf_oper_dir',required=False,default=None, type=str,
                        help='Dir with Oper files for wrf_comparison mode'  )
    parser.add_argument('--wrf_foper_dir',required=False,default=None, type=str,
                        help='Dir with F-oper files for wrf_comparison mode'  )
    parser.add_argument('--wrf_conv_dir',required=False,default=None, type=str,
                        help='Dir with Conventional files for wrf_comparison mode'  )

    return parser.parse_args()


def main():
   import os
    
   argv = prepare_parser()

   if argv.mode == 'retrieval':

       ncfiles = glob.glob(argv.input + '/*.nc')
    
       for file in ncfiles:
    
          a = netCDFReader(file)
    
          title_str = os.path.basename(file).replace('_MIRTO_validation.nc','').replace('_',' ')
    
          y = a.pressure
          x1 = a.prior_temp
          x2 = a.post_temp
          x3 = a.sonde_temp
          mirto_plot_validation_profiles(x1,x2,x3,y,'T',title_str) 
          
          ofile = file.replace(os.path.basename(file),'T_'+os.path.basename(file).replace('nc','png')) 
          plt.savefig(ofile, transparent = True)
    
          plt.close()
    
          x1 = a.prior_rh
          x2 = a.post_rh
          x3 = a.sonde_rh
          mirto_plot_validation_profiles(x1,x2,x3,y,'RH',title_str,rms=rms_rh)
    
          ofile = file.replace(os.path.basename(file),'RH_'+os.path.basename(file).replace('nc','png'))
          plt.savefig(ofile, transparent = True)
    
          plt.close()
          
          x1 = a.prior_water_vapour
          x2 = a.post_water_vapour
          x3 = a.sonde_water_vapour
    
          mirto_plot_validation_profiles(x1,x2,x3,y,'WV',title_str)
    
          ofile = file.replace(os.path.basename(file),'WV_'+os.path.basename(file).replace('nc','png'))
          plt.savefig(ofile, transparent = True)    
          plt.close()
          
   elif argv.mode == 'wrf':
       
       if argv.filetype == 'aggregated':
           filetype = 'aggregated' 
       elif argv.filetype == 'single' and argv.forecast_hour == 0:
           filetype = 'analysis' 
       elif argv.filetype == 'single' and argv.forecast_hour != 0:
           filetype = 'forecast' 

       if argv.forecast_hour == 0:
           forecast_hour = '_analysis_'  
       elif argv.forecast_hour > 0:
           forecast_hour = '+{}H'.format(argv.forecast_hour)
       else:
           forecast_hour = '_all_forecasts_'  
       
       print('Path: {}'.format(argv.input))
       ncfiles = [ argv.input+'/'+file for file in os.listdir(argv.input) if '.nc' in file \
                                                                       and filetype in file \
                                                                       and forecast_hour in file]
       print('Ncfiles: {}'.format(ncfiles))
           
       if len(ncfiles) == 0:
           print("No {} file found in {}".format(argv.filetype,argv.input))
           
       for file in ncfiles:
           
           print("Plotting validation file {}".format(os.path.basename(file)))
           data = netCDFReader(file)

           title_str = os.path.basename(file).replace('_',' ').replace('.nc','')
           if argv.filetype == 'aggregated':
               title_str = title_str + ' (statistics over {} cases)'.format(data.n_cicles)
               y = data.sonde_levels.mean(axis=0)
               plot_mode = 'aggregated'
           else:
               y =  data.sonde_levels
               plot_mode = 'single'
               
           x1 = data.temp
           x2 = data.sonde_temp
           rms_temp = data.rms_profile_temp

           plot_wrf_validation_profiles(x1,x2,y,'T',title_str, rms=rms_temp,plot_mode=plot_mode) 
            
           ofile = argv.outdir +'/'+ os.path.basename(file).replace('.nc','.png').replace('Validation','Validation_T')
           plt.savefig(ofile, transparent = True)
           print("Saved file {}".format(ofile))
           plt.close()
           """
           from atmos.mirto_atmos_tools import mr2rh
           x1 = mr2rh(data.wrf_levels*1e-2,
                      data.temp,
                      data.water_vapour)[0]
           """
           x1 = data.rh
           x2 = data.sonde_rh
           rms_rh = data.rms_profile_rh

           ofile = ofile.replace('T','RH')
           plot_wrf_validation_profiles(x1,x2,y,'RH',title_str,rms=rms_rh,plot_mode=plot_mode)     
           plt.savefig(ofile, transparent = True)   
           print("Saved file {}".format(ofile))
           plt.close()
          
           x1 = data.water_vapour
           x2 = data.sonde_water_vapour
           rms_wv = data.rms_profile_water_vapour
           ofile = ofile.replace('RH','WV')

           plot_wrf_validation_profiles(x1,x2,y,'WV',title_str,rms=rms_wv,plot_mode=plot_mode)     
           plt.savefig(ofile, transparent = True)   
           print("Saved file {}".format(ofile))
           plt.close()
    
   elif argv.mode == 'wrf_comparison':
           
           filetype = 'aggregated' 
        
           if argv.forecast_hour == 0:
                forecast_hour = '_analysis_'  
           elif argv.forecast_hour > 0:
                forecast_hour = '+{}H'.format(argv.forecast_hour)
           else:
                forecast_hour = '_all_forecasts_'  
           
           # Read all aggregated files
           nc_files = []
           #PaoloA 31Aug2022
           #for path in [argv.wrf_control_dir, argv.wrf_oper_dir, argv.wrf_foper_dir,argv.wrf_conv_dir]:
           for path in [argv.wrf_oper_dir, argv.wrf_foper_dir,argv.wrf_conv_dir]:
                        file = [ path+'/'+file for file in os.listdir(path) \
                                                         if '.nc' in file \
                                                         and filetype in file \
                                                         and forecast_hour in file]
                        if len(file) == 0:
                            sys.exit("Aggregated file for {} forecast hours not found "
                                     "in {}".format(argv.forecast_hour,path))
                        else:
                            print("Reading ",file[0],'...')
                            nc_files.append(netCDFReader(file[0]))
                       
           #PaoloA 31Aug2022
           #control, oper, foper, conv = nc_files
           oper, foper, conv = nc_files
           print('NCfiles: {}'.format(nc_files))
           
           #PaoloA 31Aug2022
           #common_forecasts = np.intersect1d(control.ID, 
           #                                  np.intersect1d(oper.ID,
           #                                                 np.intersect1d(foper.ID,conv.ID)
           common_forecasts = np.intersect1d(oper.ID, 
                                             np.intersect1d(foper.ID,conv.ID)
                                             )
           if len(common_forecasts) == 0:
               sys.exit('No common forecast!')
           else:
               print("!!!  Found {} common forecasts   !!!".format(len(common_forecasts)))
           
           #PaoloA 31Aug2022
           #control_indx = []
           oper_indx = []
           foper_indx = []
           conv_indx = []
           
           for common_id in common_forecasts:
               #PaoloA 31Aug2022
               #control_indx.append(np.where(control.ID == common_id)[0][0])
               oper_indx.append(np.where(oper.ID == common_id)[0][0])
               foper_indx.append(np.where(foper.ID == common_id)[0][0])
               conv_indx.append(np.where(conv.ID == common_id)[0][0])
           #PaoloA 31Aug2022
           #control_indx=np.array(control_indx)
           oper_indx=np.array(oper_indx)
           foper_indx=np.array(foper_indx)
           conv_indx=np.array(conv_indx)

           y = conv.sonde_levels.mean(axis=0)
           title = os.path.basename(file[0]).replace('Validation',
                                                     'Validation_MAE_Comparison_').replace('aggregated.nc',
                                                                                           '').replace('_',
                                                                                                       ' ')
           title += ' ( {} cases) '.format(conv.n_cicles)
           
           for var in ['temp','rh','water_vapour']:
               rms_string = 'rms_profile_'+var
               #PaoloA 31Aug2022
               #ctrl_rms_dif  = conv.__dict__[rms_string][conv_indx,:] - control.__dict__[rms_string][control_indx,:]
               oper_rms_dif  = conv.__dict__[rms_string][conv_indx,:] - oper.__dict__[rms_string][oper_indx,:]
               foper_rms_dif = conv.__dict__[rms_string][conv_indx,:]  - foper.__dict__[rms_string][foper_indx,:]

               conv_rms  = conv.__dict__[rms_string][conv_indx,:]
               #PaoloA 31Aug2022
               #mw_rms  =  control.__dict__[rms_string][control_indx,:]
               oper_rms  = oper.__dict__[rms_string][oper_indx,:]
               foper_rms = foper.__dict__[rms_string][foper_indx,:]

               var_name = {'temp':'T','rh':'RH','water_vapour':'WV'}.get(var)
               
               #PaoloA 31Aug2022
               #plot_wrf_comparison_profile_diff(ctrl_rms_dif,
               plot_wrf_comparison_profile_diff(oper_rms_dif,
                                            foper_rms_dif,
                                            y,
                                            var_name,
                                            sup_title_str = title)
               
               ofile = argv.outdir +'/' + os.path.basename(file[0]).replace('Validation','Validation_MAE_Comparison_'+var_name).replace('_aggregated.nc','.png')
               if '+-1' in ofile:
                   ofile = ofile.replace('forecast+-1H','all_forecasts')
                   
               plt.savefig(ofile, transparent = True)   
               print("Saved file ",ofile)


               #PaoloA 31Aug2022
               plot_wrf_comparison_profiles(conv_rms,
               #                             mw_rms,
                                            oper_rms,
                                            foper_rms,
                                            y,
                                            var_name,
                                            sup_title_str = title)

               ofile = argv.outdir +'/' + os.path.basename(file[0]).replace('Validation','Validation_Single_MAE_Comparison_'+var_name).replace('_aggregated.nc','.png')
               if '+-1' in ofile:
                   ofile = ofile.replace('forecast+-1H','all_forecasts')

               plt.savefig(ofile, transparent = True)
               print("Saved file ",ofile)
           

if __name__ == "__main__":
    main()

