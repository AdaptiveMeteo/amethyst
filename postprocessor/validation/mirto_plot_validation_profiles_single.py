import matplotlib.pyplot as plt
from matplotlib import cm
from utilities.data_reader.read_netcdf import netCDFReader
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


def plot_validation_profiles(x1,x2,x3,y,var,sup_title_str=None):

   fig, ax = plt.subplots(1,2)
   fig.set_size_inches(16.5, 10.5)

   if var in 'T':
      xlabel_str_0 = 'Air Temperature [K]'
      xlabel_str_1 = 'Delta Air Temperature [K]'
      ylabel_str = 'p [hPa]'
      title_str =  'Temperature'
      ylim = [1014, 200]
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
   
   lgnd = ax[0].legend(['wrf fg', 'amethyst ret', 'sonde'], loc='upper right', shadow=True)
   ax[0].grid(True)
   ax[0].set(yscale = 'log')

   d1 = x1 - x3 
   d2 = x2 - x3
   
   ax[1].plot(d1, y, 'b+', d2, y, 'ro', linestyle='solid')
   ax[1].set(xlabel = xlabel_str_1, ylabel = ylabel_str)
   ax[1].set_title(title_str, pad = 12, fontsize = 14)
   lgnd = ax[0].legend(['wrf fg', 'amethyst ret', 'sonde'], loc='upper right', shadow=True)
   ax[1].set_ylim(ylim)
   lgnd = ax[1].legend(['wrf fg - sonde', 'amethyst ret - sonde'], loc='upper right', shadow=True)
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
      ylim = [1014, 200]
      xlim = [-5.0,5.0]
   elif var in 'RH':
      xlabel_str_0 = 'Air Relative Humidity [%]'
      xlabel_str_1 = 'Delta Air Relative Humidity [%]'
      ylabel_str = 'p [hPa]'
      title_str =  'Relative Humidity'
      ylim = [1014, 200]
      xlim = [-35,35]

   elif var in 'WV':
      xlabel_str_0 = 'Air Water Vapor MR [g/kg]'
      xlabel_str_1 = 'Delta Air Water Vapor MR [g/kg]'
      ylabel_str = 'p [hPa]'
      title_str =  'Water Vapor MR'
      ylim = [1014, 200]
      xlim = [-3.5,3.5]
      
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
           ax[1].plot(rms,y,marker='o',label='RMS',c='red')
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
      ylim = [1014, 200]
      xlim = [-5,5]
      xgrid = np.arange(-5,6,1)

   elif var in 'RH':
      xlabel_str_0 = 'Air Relative Humidity [%]'
      xlabel_str_1 = 'Delta Air Relative Humidity [%]'
      ylabel_str = 'p [hPa]'
      title_str =  'Relative Humidity'
      ylim = [1014, 200]
      xlim = [-35,35]
      xgrid = np.arange(-30,40,10)

   elif var in 'WV':
      xlabel_str_0 = 'Air Water Vapor MR [g/kg]'
      xlabel_str_1 = 'Delta Air Water Vapor MR [g/kg]'
      ylabel_str = 'p [hPa]'
      title_str =  'Water Vapor MR'
      ylim = [1014, 200]
      xlim = [-3.5,3.5]
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
   ax.plot(oper_dif, y,'blue',label = 'CONV - TR',  
                       linestyle='solid')
   #ax.fill_betweenx(y, oper_mean-oper_std , oper_mean+oper_std,
   #                   color='tab:blue',alpha=0.5,label=r'$ \mu \pm \sigma $')

   # Plot F-OPER Diff
   foper_std  = foper_dif.std(axis=0)
   foper_mean = foper_dif.mean(axis=0)
   #PaoloA 31Aug2022
   #ax.plot(foper_mean, y,'red',label = 'CNTRL - (TR+MW)',  
   ax.plot(foper_dif, y,'red',label = 'CONV - (TR+MW)',  
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
def plot_wrf_comparison_profiles(conv,oper,foper,y,var,sup_title_str=None):

   fig, ax = plt.subplots(1,1)
   fig.set_size_inches(10.5, 16.5)

   if var in 'T':
      xlabel_str_0 = 'Air Temperature [K]'
      xlabel_str_1 = 'Delta Air Temperature [K]'
      ylabel_str = 'p [hPa]'
      title_str =  'Temperature'
      ylim = [1014, 200]
      xlim = [-5-5-5,5.0]
      xgrid = np.arange(-5,6,1)

   elif var in 'RH':
      xlabel_str_0 = 'Air Relative Humidity [%]'
      xlabel_str_1 = 'Delta Air Relative Humidity [%]'
      ylabel_str = 'p [hPa]'
      title_str =  'Relative Humidity'
      ylim = [1014, 200]
      #Paolo
      xlim = [5,35]
      xgrid = np.arange(-30,40,10)

   elif var in 'WV':
      xlabel_str_0 = 'Air Water Vapor MR [g/kg]'
      xlabel_str_1 = 'Delta Air Water Vapor MR [g/kg]'
      ylabel_str = 'p [hPa]'
      title_str =  'Water Vapor MR'
      ylim = [1014, 200]
      xlim = [0,3.5]
      xgrid = np.arange(-3,4,1)
   ygrid = np.arange(100,1100,100)

   # Plot CNTRL 
   conv_std  = conv.std(axis=0)
   conv_mean = conv.mean(axis=0)
   ax.plot(conv, y,'red',label = 'CONV',
                       linestyle='solid')

   # Plot OPER 
   oper_std  = oper.std(axis=0)
   oper_mean = oper.mean(axis=0)
   ax.plot(oper, y,'blue',label = 'CONV+TR',
                       linestyle='solid')

   # Plot F-OPER 
   foper_std  = foper.std(axis=0)
   foper_mean = foper.mean(axis=0)
   ax.plot(foper, y,'green',label = 'CONV+TR+MW',
                       linestyle='solid')

   ax.set(xlabel = xlabel_str_1, ylabel = ylabel_str)

   ax.set_xlabel(xlabel_str_1,size=15)
   ax.set_ylabel(ylabel_str,size=15)

   ax.set_title(title_str, pad = 12, fontsize = 18)
   ax.set_ylim(ylim)
   ax.grid(True,which='both')
   ax.grid(True,which='both',c='grey',ls='--',alpha=0.5)
   ax.set(yscale = 'log')


   fig.suptitle(sup_title_str, fontsize = 18)
   ax.set_yticks(ygrid)
   ax.set_yticklabels(ygrid.astype(str), fontsize = 16)
   ax.set_xticks(xgrid)
   ax.set_xticklabels(xgrid.astype(str), fontsize = 16)
   ax.set_xlim(xlim)

   plt.legend(loc='upper right', shadow=True,prop={'size':16})

   return fig
    
#PaoloA 15Apr2022
def plot_wrf_comparison_profiles_per(conv,oper,foper,y,var,sup_title_str=None):

   fig, ax = plt.subplots(1,1)
   fig.set_size_inches(5.5, 16.5)

   if var in 'T':
      xlabel_str_0 = 'Air Temperature [K]'
      xlabel_str_1 = 'Delta Air Temperature [K]'
      ylabel_str = 'p [hPa]'
      title_str =  'Temperature'
      ylim = [1014, 200]
      xlim = [-30.0,30.0]
      xgrid = np.arange(-30,30,10)

   elif var in 'RH':
      xlabel_str_0 = 'Air Relative Humidity [%]'
      xlabel_str_1 = 'Delta Air Relative Humidity [%]'
      ylabel_str = 'p [hPa]'
      title_str =  'Relative Humidity'
      ylim = [1014, 200]
      xlim = [-30.0,30.0]
      xgrid = np.arange(-30,30,10)

   elif var in 'WV':
      xlabel_str_0 = 'Air Water Vapor MR [g/kg]'
      xlabel_str_1 = 'Delta Air Water Vapor MR [g/kg]'
      ylabel_str = 'p [hPa]'
      title_str =  'Water Vapor MR'
      ylim = [1014, 200]
      xlim = [-30.0,30.0]
      xgrid = np.arange(-30,30,10)
   ygrid = np.arange(100,1100,100)

   # Plot OPER 
   ax.plot(oper, y,'blue',label = '(TR-CONV)/CONV',
                       linestyle='solid')

   # Plot F-OPER 
   ax.plot(foper, y,'green',label = '(TR+MW-CONV)/CONV',
                       linestyle='solid')

   ax.set(xlabel = xlabel_str_1, ylabel = ylabel_str)

   ax.set_xlabel(xlabel_str_1,size=15)
   ax.set_ylabel(ylabel_str,size=15)

   ax.set_title(title_str, pad = 12, fontsize = 18)
   ax.set_ylim(ylim)
   ax.grid(True,which='both')
   ax.grid(True,which='both',c='grey',ls='--',alpha=0.5)
   ax.set(yscale = 'log')


   fig.suptitle(sup_title_str, fontsize = 18)
   ax.set_yticks(ygrid)
   ax.set_yticklabels(ygrid.astype(str), fontsize = 16)
   ax.set_xticks(xgrid)
   ax.set_xticklabels(xgrid.astype(str), fontsize = 16)
   ax.set_xlim(xlim)

   plt.legend(loc='upper right', shadow=True,prop={'size':16})

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
    parser.add_argument('--format','-f',default='png', type=str,
                        help='Output file format'  )
    parser.add_argument('--errortype','-et',default='rmse', type=str,choices=['rmse','rmsebc','rmseper'],
                        help='Output file format'  )

    return parser.parse_args()


def main():
   import os
    
   argv = prepare_parser()

   if argv.mode == 'retrieval':

       ncfiles = glob.glob(argv.input + '/*.nc')
    
       for file in ncfiles:
    
          a = netCDFReader(file)
    
          title_str = os.path.basename(file).replace('_amethyst_validation.nc','').replace('_',' ')
    
          y = a.pressure
          x1 = a.prior_temp
          x2 = a.post_temp
          x3 = a.sonde_temp
          plot_validation_profiles(x1,x2,x3,y,'T',title_str) 
          
          ofile = file.replace(os.path.basename(file),'T_'+os.path.basename(file).replace('nc',argv.format)) 
          plt.savefig(ofile, transparent = True)
    
          plt.close()
    
          x1 = a.prior_rh
          x2 = a.post_rh
          x3 = a.sonde_rh
          plot_validation_profiles(x1,x2,x3,y,'RH',title_str,rms=rms_rh)
    
          ofile = file.replace(os.path.basename(file),'RH_'+os.path.basename(file).replace('nc',argv.format))
          plt.savefig(ofile, transparent = True)
    
          plt.close()
          
          x1 = a.prior_water_vapour
          x2 = a.post_water_vapour
          x3 = a.sonde_water_vapour
    
          plot_validation_profiles(x1,x2,x3,y,'WV',title_str)
    
          ofile = file.replace(os.path.basename(file),'WV_'+os.path.basename(file).replace('nc',argv.format))
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
               
           #x1 = np.reshape(data.temp, (-1, len(y)))
           #x2 = np.reshape(data.sonde_temp, (-1, len(y)))
           x1 = data.temp
           x2 = data.sonde_temp
           #Paolo
           if argv.errortype == 'rmse':
               rms_temp = np.sqrt(np.mean(data.delta_profile_temp**2,axis=0)) 
           else:
               rms_temp = np.sqrt(np.mean(data.delta_profile_temp**2,axis=0)-(np.mean(data.delta_profile_temp,axis=0))**2) 

           print('T {}'.format(data.delta_profile_temp,axis=0))
           print('Size of rms T: {}'.format(rms_temp.size))
           print('mean(T) {}'.format(np.mean(data.delta_profile_temp,axis=0)))
           print('std(T) {}'.format(np.std(data.delta_profile_temp,axis=0)))
           print('rms(T) {}'.format(rms_temp))

           plot_wrf_validation_profiles(x1,x2,y,'T',title_str, rms=rms_temp,plot_mode=plot_mode) 
            
           ofile = argv.outdir +'/'+ os.path.basename(file).replace('.nc','.'+argv.format).replace('Validation','Validation_T')
           plt.savefig(ofile, transparent = True)
           print("Saved file {}".format(ofile))
           plt.close()

           x1 = data.rh
           x2 = data.sonde_rh
           #Paolo
           if argv.errortype == 'rmse':
               rms_rh = np.sqrt(np.mean(data.delta_profile_rh**2,axis=0))
           else:
               rms_rh = np.sqrt(np.mean(data.delta_profile_rh**2,axis=0)-(np.mean(data.delta_profile_rh,axis=0))**2)

           ofile = ofile.replace('T','RH')
           plot_wrf_validation_profiles(x1,x2,y,'RH',title_str,rms=rms_rh,plot_mode=plot_mode)     
           plt.savefig(ofile, transparent = True)   
           print("Saved file {}".format(ofile))
           plt.close()
          
           x1 = data.water_vapour
           x2 = data.sonde_water_vapour
           #Paolo
           if argv.errortype == 'rmse':
               rms_wv = np.sqrt(np.mean(data.delta_profile_water_vapour**2,axis=0)) 
           else:
               rms_wv = np.sqrt(np.mean(data.delta_profile_water_vapour**2,axis=0)-(np.mean(data.delta_profile_water_vapour,axis=0))**2) 


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
                                                     'Validation_RMS_Comparison_').replace('aggregated.nc',
                                                                                           '').replace('_',
                                                                                                       ' ')
           title += ' ( {} cases) '.format(conv.n_cicles)
           
           for var in ['temp','rh','water_vapour']:
               rms_string = 'delta_profile_'+var
               #PaoloA 31Aug2022
               #np.sqrt(np.mean(data.delta_profile_water_vapour**2,axis=0))
               #ctrl_rms_dif  = conv.__dict__[rms_string][conv_indx,:] - control.__dict__[rms_string][control_indx,:]
               #Paolo
               if argv.errortype == 'rmse':
                   oper_rms_dif  = np.sqrt(np.mean(conv.__dict__[rms_string][:,:]**2,axis=0))-np.sqrt(np.mean(oper.__dict__[rms_string][:,:]**2,axis=0))
                   foper_rms_dif = np.sqrt(np.mean(conv.__dict__[rms_string][:,:]**2,axis=0))-np.sqrt(np.mean(foper.__dict__[rms_string][:,:]**2,axis=0)) 
               else: 
                   oper_rms_dif  = np.sqrt(np.mean(conv.__dict__[rms_string][:,:]**2,axis=0)-(np.mean(conv.__dict__[rms_string][:,:],axis=0))**2)-np.sqrt(np.mean(oper.__dict__[rms_string][:,:]**2,axis=0)-(np.mean(oper.__dict__[rms_string][:,:],axis=0))**2)
                   foper_rms_dif = np.sqrt(np.mean(conv.__dict__[rms_string][:,:]**2,axis=0)-(np.mean(conv.__dict__[rms_string][:,:],axis=0))**2)-np.sqrt(np.mean(foper.__dict__[rms_string][:,:]**2,axis=0)-(np.mean(foper.__dict__[rms_string][:,:],axis=0))**2)

               #Paolo
               if argv.errortype == 'rmse':
                   conv_rms  = np.sqrt(np.mean(conv.__dict__[rms_string][:,:]**2,axis=0))
               else:
                   conv_rms  = np.sqrt(np.mean(conv.__dict__[rms_string][:,:]**2,axis=0)-(np.mean(conv.__dict__[rms_string][:,:],axis=0))**2)
               #PaoloA 31Aug2022
               #mw_rms  =  control.__dict__[rms_string][control_indx,:]
               #Paolo
               if argv.errortype == 'rmse':
                   oper_rms  = np.sqrt(np.mean(oper.__dict__[rms_string][:,:]**2,axis=0))
                   foper_rms = np.sqrt(np.mean(foper.__dict__[rms_string][:,:]**2,axis=0)) 
               elif argv.errortype == 'rmseper':
                   conv_rms  = np.sqrt(np.mean(conv.__dict__[rms_string][:,:]**2,axis=0)-(np.mean(conv.__dict__[rms_string][:,:],axis=0))**2)
                   oper_rms  = np.sqrt(np.mean(oper.__dict__[rms_string][:,:]**2,axis=0)-(np.mean(oper.__dict__[rms_string][:,:],axis=0))**2)
                   foper_rms = np.sqrt(np.mean(foper.__dict__[rms_string][:,:]**2,axis=0)-(np.mean(foper.__dict__[rms_string][:,:],axis=0))**2) 
                   oper_rms_p  = 100.0*(oper_rms-conv_rms)/conv_rms
                   foper_rms_p = 100.0*(foper_rms-conv_rms)/conv_rms
               else:
                   oper_rms  = np.sqrt(np.mean(oper.__dict__[rms_string][:,:]**2,axis=0)-(np.mean(oper.__dict__[rms_string][:,:],axis=0))**2)
                   foper_rms = np.sqrt(np.mean(foper.__dict__[rms_string][:,:]**2,axis=0)-(np.mean(foper.__dict__[rms_string][:,:],axis=0))**2) 


               var_name = {'temp':'T','rh':'RH','water_vapour':'WV'}.get(var)

               print('Shape of conv_rms: {}'.format(conv_rms.shape))
               print('Shape of oper_rms: {}'.format(oper_rms.shape))
               print('Shape of foper_rms: {}'.format(foper_rms.shape))
               print('Shape of y: {}'.format(y.shape))
               
               #PaoloA 31Aug2022
               #plot_wrf_comparison_profile_diff(ctrl_rms_dif,
               plot_wrf_comparison_profile_diff(oper_rms_dif,
                                            foper_rms_dif,
                                            y,
                                            var_name,
                                            sup_title_str = title)
               
               ofile = argv.outdir +'/' + os.path.basename(file[0]).replace('Validation','Validation_' + argv.errortype + '_Comparison_'+var_name).replace('_aggregated.nc','.'+argv.format)
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

               ofile = argv.outdir +'/' + os.path.basename(file[0]).replace('Validation','Validation_Single_' + argv.errortype + '_Comparison_'+var_name).replace('_aggregated.nc','.'+argv.format)
               if '+-1' in ofile:
                   ofile = ofile.replace('forecast+-1H','all_forecasts')

               plt.savefig(ofile, transparent = True)
               print("Saved file ",ofile)

               #PaoloA 15Apr2022
               if argv.errortype == 'rmseper':
                   plot_wrf_comparison_profiles_per(conv_rms,
                                            oper_rms_p,
                                            foper_rms_p,
                                            y,
                                            var_name,
                                            sup_title_str = title)

                   ofile = argv.outdir +'/' + os.path.basename(file[0]).replace('Validation','Validation_Single_' + argv.errortype + '_Comparison_'+var_name).replace('_aggregated.nc','.'+argv.format)
                   if '+-1' in ofile:
                       ofile = ofile.replace('forecast+-1H','all_forecasts')

                   print("OFILE: {}".format(ofile))

                   plt.savefig(ofile, transparent = True)
                   print("Saved file ",ofile)

           

if __name__ == "__main__":
    main()


