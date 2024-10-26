#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 15:32:02 2021

@author: paolo
"""

__author__    = "Paolo Antonelli and Paolo Scaccia"
__copyright__ = "Copyright 2021, AdaptiveMeteo S.r.l."
__credits__   = ["Paolo Antonelli","Paolo Scaccia"]
__license__   = "--"
__version__ = "0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

from data_reader.mirto_wrapper import MIRTO_joined_Wrapper
import sys
from argparse import ArgumentParser
import numpy as np

######  CLASS FOR EMISSIVITY  #################################################

class Emissivity(object):   
    """
    This class provides a way to compute emissivity values
    """

    def __init__(self, fg):
        """Reads Model parameters from file"""
        from netCDF4 import Dataset

        self.df = Dataset(fg, mode='r')
        self.obsnum = len(self.df.dimensions['number_of_FOVs'])
        self.SEC = None
        self.SEwn = self.df['surface_components']['ModelWaveNumbers'][:]
        
    def __del__(self):
        self.df.close()

    def get(self, SEC):
        """Proxy for getting emissivities"""
        
        self.n_coef = SEC.shape[1]
        self.SEMF   = self.df['surface_components']['ModelFunctions'][:,:self.n_coef,Ellipsis]
        self.SEMFB  = self.df['surface_components']['ModelFunctionsBias'][:,Ellipsis]
        self.compute_from_model(SEC)
        self.SEC = SEC

        return [self.SEwn, self.SEvalues]

    def compute_from_model(self, SEC):
        """
        Compute Emissivity values from emissivity coefficients
        from solution vector
        """        
        SEC[SEC.mask] = 0
        # load surface emissivity coefficients from solution vector
        self.SEvalues = np.zeros_like(self.SEMFB)
        for i in range(self.obsnum):            
            self.SEvalues[i,:] = np.dot(self.SEMF[i,:,:].T, SEC[i,:]) + self.SEMFB[i,:]

        self.SEvalues = np.exp(self.SEvalues) / (1.0 + np.exp(self.SEvalues))
        # self.SEvalues = 1.0 / (1.0 + np.exp(-self.SEvalues))
        
        self.SEC = SEC
        
        return

########################################################################################  END  ### 

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]

units = {  'air_temperature'       : '° C',
           'air_relative_humidity' : '%',
           'air_water_vapor_mr'    : 'g/Kg'
           }

def compute_diff(ref_data,input_data,outdir, ozone = False, emiss_file = None, mode = 'dif',overpass=None):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines
    
    plt_title = "Average Difference" if mode in ["dif","fg_comparison",'fg_vs_wrf'] else "Single Retrieval"
    
    # <------------  Insert FOV CONTROL
 
    # Prepare y axes
    y = ref_data.air_pressure.mean(axis=0)
    y_surf = 1015.0
    y_top = 10
    ygrid = np.arange(1000.0, 100.0, -100.0) 
    ygrid = np.append(ygrid,[100.0, 10.0, 1.0, 0.1])
    yrange = [y_surf, y_top]
    i_surf, val_surf = find_nearest( ygrid, y_surf )
    i_top, val_top   = find_nearest( ygrid, y_top )
    
    """
    i2_top, val2_top   = find_nearest(y,y_top)
    i2_surf, val2_surf = find_nearest(y,y_surf)

    def grid2lev(pressure_grid):
            out = np.zeros_like(pressure_grid)
            if len(pressure_grid.shape) > 1:
                for i,praw in enumerate( pressure_grid):
                    for j,pvalue in enumerate(praw):
                        out[i,j] =  find_nearest(y,pvalue)[0]
            else:
                for i,pvalue in enumerate(pressure_grid):
                    out[i] =  find_nearest(y,pvalue)[0]
            
            return out
    def lev2grid(levels):

            out = np.zeros_like(levels,dtype=int)
            if len(levels.shape) > 1:    
                for i,nraw in enumerate( levels):
                    for j,n in enumerate(nraw):
                        if int(n) > y.size:
                            out[i,j] = y.size
                        else:
                            out[i,j] =  find_nearest(ygrid,y[int(n)])[1]
            else:
                for i,n in enumerate(levels):
                    if int(n) > y.size:
                            out[i] = y.size
                    else:
                        out[i] =  find_nearest(ygrid,y[int(n)])[1]
 
            return out
                
    """  
    for var in ['air_water_vapor_mr','air_temperature','air_relative_humidity']:
        
        delta = ref_data.__dict__[var] - input_data.__dict__[var]
        delta = delta.mean(axis=0)
        std   = delta.std(axis=0)
                
        fig = plt.figure()
        fig.set_size_inches(16.5, 10.5)

        plt.title("{} - {}".format(plt_title,var.replace('_',' ')))
        plt.plot(delta,y,label='Mean',marker='s',ms=2)
        
        legend_patch = mpatches.Patch(color='#539ecd', label=r'$ 1\sigma $  Interval',alpha = 0.3)
        legend_line = mlines.Line2D([], [], color='tab:blue', marker='s',
                          markersize=3, label='Mean')
        
        plt.fill_betweenx(y,delta - std, delta + std, color='#539ecd',alpha=0.3)
        ax = plt.gca()
        
        ax.set_ylim(yrange)
        plt.ylabel('Pressure [hPa]',size=13)        
        ax.set_yticks(ygrid[np.arange(i_surf,i_top+1)])
        ax.set_yticklabels(ygrid[np.arange(i_surf,i_top+1)].astype(str))
        
        plt.yscale('log')   
        for i in np.arange(0,y.size,5):
            value = y[i]
            if value < val_surf and value > val_top:
                plt.axhline(y=value,c='grey',linewidth=0.5)
                plt.text(delta.max() + std.max(), value,c='grey',s='lev {}'.format(i),verticalalignment='top')
            
        """
        ax2 = ax.secondary_yaxis('right',functions=(grid2lev,lev2grid))

        #ax2.set_yticks(ygrid[np.arange(i_surf,i_top+1)])
        ax2.set_yticklabels(np.arange(i2_surf, i2_top).astype(str) )

        ax2.set_ylabel('Level')
        ax2.grid(True)                       
        ax2.set_ylim([val2_surf,val2_top])        
        ax2.set_yticks(y[i2_surf:i2_top])
        
        ax2.grid(True)
        """
        #ax.grid(True,ls='dashed',c='black',linewidth=1)
        plt.gca().legend(handles=[legend_line,legend_patch])

        plt.xlabel(r'$ \Delta $ ( {} )'.format(units[var]),size=13)
        
        filename = '{}_{}.png'.format(mode,var)
        if overpass != None:
            filename = filename.replace('.png','_{}.png'.format(overpass))
        plt.ylabel('Pressure (hPa)',size=13)
        plt.savefig(outdir + "/" + filename,dpi=300)
        plt.clf()
        
        print("Saved plot {}".format(outdir + '/' + filename))
        
    if emiss_file is not None:
        
        fig = plt.figure()
        fig.set_size_inches(16.5, 10.5)
        ax = plt.gca()
        ax.grid(True)
        emiss_ref   = Emissivity(emiss_file[0])
        ref_sfgrd, ref_emrf = emiss_ref.get(ref_data.ems_coeff) if "ems_coeff" in ref_data.__dict__.keys() else emiss_ref.get(np.zeros_like(input_data.ems_coeff))

        if mode != 'single':
            emiss_input = Emissivity(emiss_file[1])
            in_sfgrd, in_emrf   = emiss_input.get(input_data.ems_coeff)

        
        plt.title("{} - Emissivity".format(plt_title))

            
        dif = (ref_emrf - in_emrf) if mode != 'single' else ref_emrf

        std = dif.std(axis=0)
        dif = dif.mean(axis=0)
        
        plt.plot(ref_sfgrd,dif,label='Mean',marker='s',ms=2,c='red')
        plt.xlabel(r'Wave Number ($cm^{-1}$)',size=13)
        plt.ylabel('Emissivity',size=13)
        
        plt.fill_between(ref_sfgrd,dif + std, dif - std, color='salmon',alpha=0.3)
        
        legend_patch = mpatches.Patch(color='salmon', label=r'$ 1\sigma $  Interval',alpha = 0.3)
        legend_line = mlines.Line2D([], [], color='red', marker='s',
                          markersize=3, label='Mean')

        #plt.yscale('log')
        #plt.gca().invert_yaxis()

        #plt.xlabel(r'$ \Delta $ ( {} )'.format(units[var]),size=13)
        #plt.ylabel('Pressure (hPa)',size=13)
        plt.gca().legend(handles=[legend_line,legend_patch])
        filename = '{}_emissivity.png'.format(mode)
        if overpass != None:
            filename = filename.replace('.png','_{}.png'.format(overpass))

        plt.savefig(outdir + '/' + filename,dpi=300)
        print("Saved plot {}".format(outdir + '/' + filename))

        
    return   
    
def residual_analisis(ref_residuals, input_residuals):
    
    
    
    
    
    return

    
def compute_diff_with_wrf(ref_data,wrf_data,outdir, ozone = False, emiss_file = None, mode = 'dif',overpass=None):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines
    
    plt_title = "Average Difference" if mode in ["dif","fg_comparison",'fg_vs_wrf'] else "Single Retrieval"
    
 
    # Prepare y axes
    y = ref_data.air_pressure.mean(axis=0)
    y_surf = 1015.0
    y_top = 10
    ygrid = np.arange(1000.0, 100.0, -100.0) 
    ygrid = np.append(ygrid,[100.0, 10.0, 1.0, 0.1])
    yrange = [y_surf, y_top]
    i_surf, val_surf = find_nearest( ygrid, y_surf )
    i_top, val_top   = find_nearest( ygrid, y_top )
    
    WRF_VAR = {'air_water_vapor_mr'    : 'water_vapour',
               'air_temperature'       : 'temperature',
               'air_relative_humidity' : 'rh'}
    for var in ['air_water_vapor_mr','air_temperature','air_relative_humidity']:
        
        
        delta = ref_data.__dict__[var] - wrf_data.profiles[0].__dict__[WRF_VAR[var]]
        delta = delta.mean(axis=0)
        std   = delta.std(axis=0)
                
        fig = plt.figure()
        fig.set_size_inches(16.5, 10.5)

        plt.title("{} - {}".format(plt_title,var.replace('_',' ')))
        plt.plot(delta,y,label='Mean',marker='s',ms=2)
        
        legend_patch = mpatches.Patch(color='#539ecd', label=r'$ 1\sigma $  Interval',alpha = 0.3)
        legend_line = mlines.Line2D([], [], color='tab:blue', marker='s',
                          markersize=3, label='Mean')
        
        plt.fill_betweenx(y,delta - std, delta + std, color='#539ecd',alpha=0.3)
        ax = plt.gca()
        
        ax.set_ylim(yrange)
        plt.ylabel('Pressure [hPa]',size=13)        
        ax.set_yticks(ygrid[np.arange(i_surf,i_top+1)])
        ax.set_yticklabels(ygrid[np.arange(i_surf,i_top+1)].astype(str))
        
        plt.yscale('log')   
        for i in np.arange(0,y.size,5):
            value = y[i]
            if value < val_surf and value > val_top:
                plt.axhline(y=value,c='grey',linewidth=0.5)
                plt.text(delta.max() + std.max(), value,c='grey',s='lev {}'.format(i),verticalalignment='top')
            
        plt.gca().legend(handles=[legend_line,legend_patch])

        plt.xlabel(r'$ \Delta $ ( {} )'.format(units[var]),size=13)
        
        filename = 'WRF_fg_{}_{}.png'.format(mode,var)
        if overpass != None:
            filename = filename.replace('.png','_{}.png'.format(overpass))
        plt.ylabel('Pressure (hPa)',size=13)
        plt.savefig(outdir + "/" + filename,dpi=300)
        plt.clf()
        
        print("Saved plot {}".format(outdir + '/' + filename))
        
    if emiss_file is not None:
        sys.exit("Not yet implemented!")
        fig = plt.figure()
        fig.set_size_inches(16.5, 10.5)
        ax = plt.gca()
        ax.grid(True)
        emiss_ref   = Emissivity(emiss_file[0])
        ref_sfgrd, ref_emrf = emiss_ref.get(ref_data.ems_coeff) if "ems_coeff" in ref_data.__dict__.keys() else emiss_ref.get(np.zeros_like(input_data.ems_coeff))

        if mode != 'single':
            emiss_input = Emissivity(emiss_file[1])
            in_sfgrd, in_emrf   = emiss_input.get(input_data.ems_coeff)

        
        plt.title("{} - Emissivity".format(plt_title))

            
        dif = (ref_emrf - in_emrf) if mode != 'single' else ref_emrf

        std = dif.std(axis=0)
        dif = dif.mean(axis=0)
        
        plt.plot(ref_sfgrd,dif,label='Mean',marker='s',ms=2,c='red')
        plt.xlabel(r'Wave Number ($cm^{-1}$)',size=13)
        plt.ylabel('Emissivity',size=13)
        
        plt.fill_between(ref_sfgrd,dif + std, dif - std, color='salmon',alpha=0.3)
        
        legend_patch = mpatches.Patch(color='salmon', label=r'$ 1\sigma $  Interval',alpha = 0.3)
        legend_line = mlines.Line2D([], [], color='red', marker='s',
                          markersize=3, label='Mean')

        #plt.yscale('log')
        #plt.gca().invert_yaxis()

        #plt.xlabel(r'$ \Delta $ ( {} )'.format(units[var]),size=13)
        #plt.ylabel('Pressure (hPa)',size=13)
        plt.gca().legend(handles=[legend_line,legend_patch])
        filename = 'WRF_FG_{}_emissivity.png'.format(mode)
        if overpass != None:
            filename = filename.replace('.png','_{}.png'.format(overpass))

        plt.savefig(outdir + '/' + filename,dpi=300)
        print("Saved plot {}".format(outdir + '/' + filename))

        
    return   
def main():
        
    parser = ArgumentParser()
    parser.add_argument('--mode', '-m', required=True,choices=['dif','single','fg_comparison','fg_vs_wrf'],
                        help='Software Mode')
    parser.add_argument('--reference','-r', required=False,
                        help='Directory with reference results')
    parser.add_argument('--input','-i', required=True,
                        help='Directory with input results')
    parser.add_argument('--outdir','-o', default='./',
                        help='Output directory')
    parser.add_argument('--emiss','-e',default=False,type=bool,
                        help='Emissivity Model')
    parser.add_argument('-ov','--overpass', required=False,type=str,
                        help="Overpass name")
    parser.add_argument('--wrf', required=False,type=str,
                        help="WRF output")

    
    argv = parser.parse_args()
        
    if argv.mode == 'dif':
        print("##########################")
        print("      Difference Mode   ")
        print("##########################")
        
        if argv.reference is None:
            sys.exit("I need a path for the reference file!")
            
        ref   = MIRTO_joined_Wrapper(argv.reference+'/results.nc')
        input = MIRTO_joined_Wrapper(argv.input+'/results.nc')
        
        if not argv.emiss:        
            compute_diff(ref,input,argv.outdir,overpass=argv.overpass)
        else:
            compute_diff(ref,input,argv.outdir, emiss_file=[argv.reference +'/fg.nc',
                                                            argv.input + '/fg.nc'],overpass=argv.overpass)
    elif argv.mode == 'single':
        print("##########################")
        print("     Single Result Mode   ")
        print("##########################")

        input = MIRTO_joined_Wrapper(argv.input+'/results.nc')
        ref   = MIRTO_joined_Wrapper(argv.input+'/results.nc')
        # Set to zero  all reference fields
        for k in ref.__dict__.keys():
            ref.__dict__[k] = np.zeros_like(input.__dict__[k])
        if not argv.emiss:        
            compute_diff(input,ref,argv.outdir,mode=argv.mode,overpass=argv.overpass)
        else:
            compute_diff(input,ref,argv.outdir, emiss_file=[argv.input + '/fg.nc',
                                                            argv.input + '/fg.nc'],mode = argv.mode,overpass=argv.overpass)
    elif argv.mode == 'fg_comparison':
        from data_reader.mirto_wrapper import MIRTO_fg_Wrapper

        print("#####################################")
        print("    Comparison with First Guess  ")
        print("#####################################")
        ref_file = argv.reference if argv.reference != None else argv.input + "/fg.nc"
        input = MIRTO_joined_Wrapper(argv.input+'/results.nc')
        ref   = MIRTO_fg_Wrapper(ref_file)
        ref.compute_rh()
        ref.to_standard()
        if not argv.emiss:
            # The input order is right: do not change!
            compute_diff(ref,input,argv.outdir,mode=argv.mode,overpass=argv.overpass)
        else:
            compute_diff(ref,input,argv.outdir, emiss_file=[argv.input + '/fg.nc', ref_file],mode = argv.mode,overpass=argv.overpass)
    elif argv.mode == 'fg_vs_wrf':
        from data_reader.mirto_wrapper import MIRTO_fg_Wrapper
        from data_reader.input_reader import InputData
        
        
        print("############################################")
        print("    Comparison First Guess vs WRF Output   #")
        print("############################################")
        
        ref_file = argv.input if 'fg.nc' in argv.input \
                                  else argv.input + "/fg.nc"
        wrf_data = InputData(argv.wrf)
        fg_ref   = MIRTO_fg_Wrapper(ref_file)
        fg_ref.compute_rh()
        fg_ref.to_standard()
        
        compute_diff_with_wrf(fg_ref,wrf_data,argv.outdir,mode=argv.mode,overpass=argv.overpass)

        
    else:
        sys.exit("Not yet implemented!")
        
    return
    
    
    
    
    
    
if __name__ == '__main__':
    sys.exit(main())
    
    
