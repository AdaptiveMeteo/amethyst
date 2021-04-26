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
from sys              import exit as sysexit
import sys
from datetime         import datetime, timedelta
from amethyst_config  import preprocessor_vars, common_vars
import os
import glob
import numpy as np
import subprocess

__author__     = 'Paolo Scaccia <paolo.scaccia@adaptivemeteo.com>'
__copyright__  = "Copyright 2021, Adaptive Meteo S.r.l"
__credits__    = ["Paolo Scaccia", "Paolo Antonelli"]
__license__    = "GPL"
__version__    = "1.0"
__maintainer__ = "Paolo Scaccia"
__email__      = "paolo.scaccia@adaptivemeteo.com"

if __name__ == '__main__':
    log = logging.getLogger()
    AMETHYST_PATH = [ x for x in sys.path if os.path.basename(x) == 'amethyst' ][0]
else:
    log = logging.getLogger(__name__)


def get_wrf_dir_time(list_date,base_date):

    #http://stackoverflow.com/a/17249529/846892
    sat_fmt = '%Y%m%d_%H%M%S'
    wrf_fmt = '%Y%m%d%H'

    d = datetime.strptime(base_date, sat_fmt)
    wrf_d = [datetime.strptime(x, wrf_fmt) for x in list_date]
    wrf_dt = [ (x - d).total_seconds() for x in wrf_d ]
    wrf_dt_neg =  [x for x in wrf_dt if x  <= 0 ]
    min_dt = np.sort( wrf_dt_neg )[-1]
    min_index = np.where(wrf_dt == min_dt)[0][0]
    
    return list_date[min_index]

def get_wrf_file_time(list_file_date,base_date):
    #http://stackoverflow.com/a/17249529/846892
    sat_fmt = '%Y%m%d_%H%M%S'
    wrf_f_fmt = '%Y-%m-%d_%H:%M:%S'
    d = datetime.strptime(base_date, sat_fmt)
    wrf_d = [datetime.strptime(x[11:], wrf_f_fmt) for x in list_file_date]
    #print('wrf_d: {}'.format(wrf_d))
    wrf_dt = [ (x - d).total_seconds() for x in wrf_d ]
    abs_wrf_dt = [ abs(x) for x in wrf_dt ]
    #print('abs_wrf_dt: {}'.format(abs_wrf_dt))
    #wrf_dt_neg =  [x for x in wrf_dt if x  <= 0 ]
    #min_dt = np.sort( abs_wrf_dt )[-1]
    #min_dt = np.sort( abs_wrf_dt )
    #print('min_dt: {}'.format(min_dt))
    min_index = np.argmin(abs_wrf_dt)
    #print('min_index: {}'.format(min_index))

    return list_file_date[min_index]

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
    
    # Launcher Support function
    def if_exists_runs(name,function,mode):
        flag = os.path.isdir(name) if mode=='dir' else os.path.isfile(name)
        if flag:
            return function(name)    
        else:
            sysexit('{} {} does not exist. Exiting'.format("Directory" if mode == "dir" else "File",name))
    get_wrf_times = lambda x : [ glob.glob(x + '/2*'), [ os.path.split(d)[-1] for d in glob.glob(x + '/2*') ] ]

    sat_ov_times                = if_exists_runs(argv.l1dir, os.listdir,'dir')
    wrf_fct_dirs, wrf_fct_times = if_exists_runs(argv.wrfdir,get_wrf_times,'dir')  
    sat_ov_times.sort(reverse=True)
    
    if argv.overpass == None:
        sat_pass_date = sat_ov_times[0]
        log.info('Last DB dir with L1 files is {}'.format(sat_pass_date))
    else:
        sat_pass_date = argv.overpass
        log.info('Selected DB dir with L1 files is {}'.format(sat_pass_date))    
    
        # Check if all the L1 data have been downloaded or are still downloading
    ready = "/".join([ argv.l1dir, sat_pass_date,'ready4processing'])
    if not os.path.isfile(ready):
            sysexit('{}/{} exists but not all needed data have been downloaded yet!\n'
                    '... retrying in 3 min\n'   
                    '... exiting for now'.format(argv.l1dir,sat_pass_date))
        
    # Check if this passage has already been processed, if so exit;
    ret_status = [sat_pass_date + '.retrieval.' + s for s in ['success', 'failed', 'processing']]
    if any([os.path.isfile(argv.logdir + '/' + f) for f in ret_status]):
       sysexit('Execution of Main stops here as overpass has already been sumbmitted to processing.\nExiting')
    else:
       open(argv.logdir + '/' + ret_status[0] , 'a').close()
 
    # sat_ov_time = datetime.strptime(sat_pass_date, "%Y%m%d_%H%M%S")
       
    #Find colosest wrf output (previous to overpass) to be used as retrieval first guess 
    sel_wrf_run=get_wrf_dir_time(wrf_fct_times,sat_pass_date)
    log.info('Selected WRF run: {}'.format(sel_wrf_run))

    #Get available wrf forecast  times
    get_wrf_file_times = lambda x: [ glob.glob(x + '/wrf/wrfout_d01*') , [os.path.split(f)[-1] for f in glob.glob(x + '/wrf/wrfout_d01*') ]]
    wrf_fct_files, wrf_fct_file_times = if_exists_runs("/".join([argv.wrfdir,sel_wrf_run]), get_wrf_file_times, 'dir')                       
    sel_wrf_run_file = get_wrf_file_time(wrf_fct_file_times, sat_pass_date)
    log.info('Selected WRF Run File: {}'.format(sel_wrf_run_file))

    # Build Selected WRF Filename
    wrf_fmt   = '%Y%m%d%H'
    wrf_f_fmt = '%Y-%m-%d_%H:%M:%S'
    wrf_o_fmt = 'wrfout_d01_%Y-%m-%d_%H:%M:%S'
    wrffilename = datetime.strptime(sel_wrf_run_file[11:], wrf_f_fmt).strftime(wrf_o_fmt)
    wrffile = "/".join([argv.wrfdir,sel_wrf_run,'wrf', wrffilename])
            
    if not os.path.isfile(   wrffile.replace(os.path.basename(wrffile),'DONE') ):
        prev_wrf_dir = datetime.strptime(sel_wrf_run, wrf_fmt)-timedelta(hours=6)
        sel_wrf_run = prev_wrf_dir.strftime(wrf_fmt)
        wrffile = "/".join([argv.wrfdir,sel_wrf_run,'wrf',wrffilename])

    if not os.path.isfile(wrffile):
         sysexit('WRF file to be used as retrieval FG: {} NOT FOUND\n'
                 '... exiting for now\n'.format(wrffilename))

    log.info('WRF file to be used as retrieval FG: {}'.format(wrffilename))   

    if argv.cluster_mode:
        #######################
        #                     # 
        #   CLUSTER MODE ON   #
        #                     #
        #######################
        
        #Write and submit shell script that runs mirto 
        sysexit("To be implemented...") # Still implementinig submit function...
        """
        cmd =  + L1DIR + sat_pass_date + ' ' + wrffile + ' ' + SHMDIR + sat_pass_date + ' v3 > ' + LOGDIR + instr + '_logfile-`date +\%Y\%m\%d\%H\%M\%S`.log'

        submit(cmd, runtime=240, cores=16, threads=1, 
               queue_name='ops',    
               directory=RUNDIR,
               shm_directory=SHMDIR,
               jobscript='run_mirto_' + instr + '_' + sat_pass_date + '.sh', 
               output=LOGDIR + '$PBS_JOBNAME.OUT', 
               error=LOGDIR + '$PBS_JOBNAME.ERR')
        """
    else:
        ########################
        #                      # 
        #   CLUSTER MODE OFF   #
        #                      #
        ########################


        
        #Build command cascade for qsub
        if argv.instrument == 'cris':
            
            
            # Calls to CrIS preprocessor scripts
            cmd_cascade = [  "python {}/preprocessor/fov_generator/cris/cloudmask/viirscris2cm.py {} {} "
                                                           "--outfile {}/cloudmask.nc -v info "
                                                           "--lonmin {} --lonmax {} --latmin {} "
                                                           "--latmax {}".format(AMETHYST_PATH,argv.gcrso, argv.scris,argv.output,
                                                                                argv.lonmin,argv.lonmax,argv.latmin,argv.latmax),
                                                           
                             "python {}/preprocessor/fov_generator/cris/cris2observations.py {} {} {}/fov.nc "
                                                           "-cmf {}/cloudmask.nc -cmt {} -v info -m {}/geo_indices.nc "
                                                           "--lonmin {} --lonmax {} --latmin {}  "
                                                           "--latmax {}".format(AMETHYST_PATH,argv.gcrso,argv.scris,argv.output,
                                                                                argv.output,argv.cmt,argv.output,
                                                                                argv.lonmin,argv.lonmax,argv.latmin,argv.latmax),
                                                    
                             "python {}/preprocessor/fg_generator/wrf2firstguess/wrf2firstguess.py --input {} "
                                                           " {}/fov.nc {}/fg.nc -v info --levels 81".format(AMETHYST_PATH,wrffile,argv.output,argv.output,
                                                                                argv.lonmin,argv.lonmax,argv.latmin,argv.latmax),

    
                             "python {}/preprocessor/fg_generator/emissivity2firstguess/emiss2firstguess.py {}/fov.nc"
                                                           " {}/fg.nc -v info".format(AMETHYST_PATH,argv.output,argv.output),

                                                           
                             "python {}/preprocessor/apriori_generator/covtable2firstguesscov.py {}/fov.nc "
                                                           " {}/apriori.nc -v info --compression 9".format(AMETHYST_PATH,argv.output,argv.output)
                                                           
                             ]

            # Logger Printouts
            logger_cascade = ["Generating CloudMask for CrIS...",
                              "Generating CrIS observations...",
                              "Generating atmospheric first guess...",
                              "Generating surface first guess...",
                              "Generating atmospheric first guess covariance..."]

            # Run cascade
            for cmd, printout in zip(cmd_cascade,logger_cascade):

                log.debug(cmd)  # Debug printout
                    
                log.info("______________________________________________")
                log.info(printout)
                subprocess.run(cmd.replace('  ',' ').split())
                
        else:
            # Iasi
            sysexit("To be added...")
            
    return

    
def preproccessor_parser():
    """
        This function parses all the arguments necessary for the preprocessor.
        Default variables are taken from amethyst/amethyst_config.py

    Returns
    -------
        argv : ArgumentParser

    """
    from argparse         import ArgumentParser
    
    
    v_levels = ['debug', 'info', 'warning']
    
    parser = ArgumentParser()
    parser.add_argument('--output','-o',type=str,required=True,
                        help='Output directory containing all NETCDF preprocessed data')
    parser.add_argument('--overpass',default = None,type=str,
                        help='Overpass date and time YYYYMMDD_HHMMSS')
    parser.add_argument('--wrfdir', '-w', type=str, default=common_vars['wrfdir'],
                        help='Where are the Wrf model data.')
    parser.add_argument('--l1dir', type=str, default=preprocessor_vars['l1dir'],
                        help='Where the L1 products are.')
    parser.add_argument('--gcrso', type=str, default=preprocessor_vars['gcrso'],
                        help='')
    parser.add_argument('--scris', type=str, default=preprocessor_vars['scris'],
                        help='')
    parser.add_argument('--iasidir', type=str, default=preprocessor_vars['iasidir'],
                        help='')
    parser.add_argument('--logdir', type=str, default=preprocessor_vars['logdir'],
                        help='Directory with log files')
    parser.add_argument('--verbose', '-v', choices=v_levels, default='info',
                        help='The level of verbosity of the software')
    parser.add_argument('--levels', '-l', type=int, default=common_vars['levels'],
                        help='The total number of levels')
    parser.add_argument('--instrument', type=str, required=True,choices=['cris','info'],
                        help='Instrument type')
    parser.add_argument('--latmin', type=str, default = common_vars['geobox']['latmin'],
                        help='Min Latitude')
    parser.add_argument('--latmax', type=str, default = common_vars['geobox']['latmax'],
                        help='Max Latitude')
    parser.add_argument('--lonmin', type=str, default = common_vars['geobox']['lonmin'],
                        help='Min Longitude')
    parser.add_argument('--lonmax', type=str, default = common_vars['geobox']['lonmax'],
                        help='Max Longitude')
    parser.add_argument('--cluster_mode', type=bool, default = False,
                        help='Cluster Mode Flag')
    parser.add_argument('--cmt', type=str, default = None,
                    help='Max Longitude')
    argv = parser.parse_args()
    
    # If missing set Cloudmask threshold 
    if argv.cmt == None:
        argv.cmt = preprocessor_vars['{}_cmt'.format(argv.instrument)]   
    
    return argv


def main():

    # Read inline arguments if the script is called from terminal    
    argv = preproccessor_parser()
         
    # Prepare the log class
    verbosity = getattr(logging, argv.verbose.upper())
    log.setLevel(verbosity)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - '
                                  '%(funcName)s: %(message)s',
                                  datefmt='%m/%d/%Y %H:%M:%S')

    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(verbosity)
    streamhandler.setFormatter(formatter)
    log.addHandler(streamhandler)

    # Call launcher
    launch_preprocessing(argv)

    return

if __name__ == '__main__':
    sysexit(main())