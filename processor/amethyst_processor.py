#!/usr/bin/env python

from __future__ import print_function, division

import argparse
from os import path
import sys
import time

# Get python major version
py_version = int(sys.version_info[0])

# Parallel libraries
from multiprocessing import Process, Queue # @UnresolvedImport
if py_version > 2:
    from queue import Empty    # @UnresolvedImport @UnusedImport
else:
    from Queue import Empty    # @UnresolvedImport @Reimport

import traceback
from ossfm.ossFM import ossFM
from dobjects.Solar import Solar
from dobjects.Hitran import Hitran
from dobjects.SounderFOV import FOVCount
from dobjects.ObservationError import create_obs_err
from dobjects.ObservationError import create_obs_err_sps
from amethyst.ForwardModel import NotConvergentIteration
from main.amethyst_code_main import core
from main.reader import reader_f
from main.scriba import scriba_f, ProgressBar

import cProfile, pstats

import amethyst_config

class logger():
    def __init__(self, verbose_level, file_log=sys.stdout):
        self.__verbose_level = int(verbose_level)
        self.output_str = ""
        self.__file_log = file_log

    def log(self, txt_val, verbose_val, print_now=True, end='\n'):
        if verbose_val <= self.__verbose_level:
            self.output_str += txt_val
            self.output_str += end
        if print_now:
            print(self.output_str, end="", file=self.__file_log)
            if self.__file_log == sys.stdout:
                sys.stdout.flush()
            self.output_str = ""

    def close(self):
        if self.__file_log != sys.stdout and self.__file_log != sys.stderr:
            self.__file_log.close()

def read_log_file(arg_passed):
    if arg_passed == 'stderr':
        return sys.stderr
    elif arg_passed == 'stdout':
        return sys.stdout
    else:
        try:
            return open(arg_passed, 'w')
        except:
            raise ValueError('Invalid file name')

def invert_process(inverter, nlev, proc_num, to_compute, to_write,
                   reader_is_alive, L, transformer,
                   output_vars, log_data, start_process ):

    proc_str = '{0:03d}'.format(proc_num+1)
    obs_computed_by_me = list()
    Myprofile = type('Myprofile', (object,), {})
    profile = Myprofile()
    profile.prall = None
    profile.prlinalg = None
    if VERBOSE == 5:
        profile.prall = cProfile.Profile()
    if VERBOSE == 6:
        profile.prlinalg = cProfile.Profile()

    first_time = True
    while reader_is_alive.empty() or first_time:
        first_time = False
        while not to_compute.empty():
            try:
                obs, data = to_compute.get(block=True, timeout=3)
            except:
                break
            L.log('%%%%%%%%%%%%%%%%%% PROCESSING OBS' + str(obs) +
                  ' %%%%%%%%%%%%%%%%%%', 4, print_now = False)
            obs_computed_by_me.append(obs)
            start = time.time()
            fov, fg, apriori, aemiss = data
            output = dict()
            if profile.prall is not None:
                profile.prall.enable( )
            try:
                time_invert=time.time()
                solution = inverter.invert(fov, fg, apriori, aemiss,
                                           obs, log=L, profile=profile)
                L.log('Elapsed Time in the inverter for OBS '+str(obs)+' : '+
                      repr(time.time()-time_invert)+' s', 3, print_now=False)
                time_output=time.time()
                if output_vars["pressure"]:
                    output['pressure'] = inverter.cx.pressure_grid
                if output_vars["temperature"]:
                    output['temperature'] = solution.xhat[0:nlev]
                if output_vars["water_vapor"]:
                    output['water_vapor'] = solution.xhat[nlev:2*nlev]
                if output_vars["ozone"]:
                    output['ozone'] = solution.xhat[2*nlev:3*nlev]
                if output_vars["surface_temperature"]:
                    output['surface_temperature'] = solution.xhat[3*nlev]
                if 'surface_emissivity_coefficients' in output_vars:
                    output['surface_emissivity_coefficients'] = solution.xhat[3*nlev+1:]
                if 'd2' in output_vars:
                    output['d2'] = solution.d2
                if 'jacobian' in output_vars:
                    output['jacobian'] = solution.jacobian
                if 'residuals' in output_vars:                   
                    output['residuals'] = solution.residuals
                #PaoloA 12112018
                if 'fgresiduals' in output_vars:                   
                    output['fgresiduals'] = solution.fgresiduals
                #PaoloA
                if 'Sa'  in output_vars:
                    output['Sa_ret'] = solution.Sa_ret
                    output['SaInv_ret'] = solution.SaInv_ret

                L.log('Elapsed Time in the output for OBS '+str(obs)+' : '+
                      repr(time.time()-time_output)+' s', 3, print_now=False)
                # Apply retrieval transfromation
                # residulas are passed as argument to allow for
                # a second option in calculating tranformed retrievals
                time_transform=time.time()
                DA = transformer.transform_retrievals_for_DA(
                              solution.Sa_ret,
                              solution.SaInv_ret,
                              solution.xa,
                              solution.xhat,
                              solution.fm.K,
                              inverter.yobs_minus_yhat,
                              profile)
                L.log('Elapsed Time in the transform for OBS '+str(obs)+' : '+
                      repr(time.time()-time_transform)+' s', 3, print_now=False)

                #These values should be added to the output netcdf file
                if 'da_r' in output_vars:
                    output['da_r'] = DA.R #optimal number of elements to be assimilated
                if 'transformed_retrievals' in output_vars:
                    output['transformed_retrievals'] = DA.Yret_prime
                if 'sn_eigenvalues' in output_vars:
                    output['sn_eigenvalues'] = DA.Lambda
                if 'assimilation_operator' in output_vars:
                    output['assimilation_operator'] = DA.Hret_prime
                    

                L.log('Elapsed Time computing (inv + trans) OBS ' + str(obs) + ' : '+
                      repr(time.time()-start)+' s', 3, print_now=False)
            except NotConvergentIteration:
                # traceback.print_exc(file=sys.stderr)
                print('%%%%%%%%%%%%%%%%%% Inversion failed for FOV in OBS'+
                      str(obs) + ' %%%%%%%%%%%%%%%%%%', file=sys.stderr)
                L.log('Wasted Time computing OBS ' + str(obs) + ': '+
                      repr(time.time()-start)+' s', 3, print_now=False)
                output = None
            if profile.prall is not None:
                profile.prall.disable( )
            to_write.put([obs,output])

    # This will alert if a process for some strange reason exits from
    # the main cycle when there are others jobs to be computed.
    # This could happen, for example, if to_compute.get reach a timeout
    if not to_compute.empty():
        print('Process ' + proc_str + ' anomaly ended after ' +
              str(time.time()-start_process) + ' seconds', file=sys.stderr)
        return

    if profile.prall is not None:
        try:
            f1 = open('./profile_process'+ proc_str +'.dat', 'w+')
            ps = pstats.Stats(profile.prall, stream=f1)
            ps.strip_dirs().sort_stats('cumulative').print_stats()
            ps.strip_dirs().sort_stats('time').print_stats()
            f1.close( )
        except:
            pass
    if profile.prlinalg is not None:
        try:
            f2 = open('./linalg_profile_process' + proc_str + '.dat', 'w+')
            ps = pstats.Stats(profile.prlinalg, stream=f2)
            ps.strip_dirs().sort_stats('cumulative').print_stats()
            ps.strip_dirs().sort_stats('time').print_stats()
            f2.close()
        except:
            pass

    if len(obs_computed_by_me) < 10:
        L.log('Process ' + proc_str + ' computed the following observations:',
               3, print_now = False)
        L.log(str(obs_computed_by_me), 3, print_now=False)
    else:
        L.log('Process ' + proc_str + ' computed ' +
               str(len(obs_computed_by_me)) +
               ' observations', 3, print_now = False)
    log_data.put(L.output_str) 



def processor(log_file, numobs, process_number, startobs, verbose,
              output_file = amethyst_config.processor_vars["outfile"]["filename"]):

 
    start_process = time.time()

    L = logger(VERBOSE, LOG_FILE)
    
    #
    # OSS init input
    #
    L.log('Reading OSS init input... ', 1, end='')

    workingDir  = amethyst_config.common_vars['wrkdir']
    co2_std     = amethyst_config.processor_vars['constant_co2_std']
    eigen_land  = amethyst_config.processor_vars['eigenforland']
    eigen_sea   = amethyst_config.processor_vars['eigenforsea']

    asolar  = Solar(amethyst_config.processor_vars["constant_solar_irradiance_file"])
    ahitran = Hitran(amethyst_config.processor_vars["od_file"])
    obs_err = create_obs_err(amethyst_config.processors_vars["noise_file"])
    obs_err_sps = create_obs_err_sps(amethyst_config.processors_vars["noise_file_sps"])

    oss_time=time.time()

    L.log('Done in ' +str(oss_time-start_process)+' seconds', 1)

    # This part of code must be repeated for each input profile in data
    # directory. Must find a way to have names here. Probably the errors
    # also can be preloaded.
    L.log('Creating an inverter... ', 1, end='')
    oss = ossFM(asolar, ahitran)
    inverter = core( oss, obs_err )
    inverter_time=time.time()
    L.log('Done in ' +str(inverter_time-oss_time)+' seconds', 1)
    # <-----
    # Check which observations should be computed
    allobs = FOVCount(amethyst_config.processor_vars["fg_file"])
    obsnum = list(range(STARTOBS, allobs))
    if NUMOBS > 0:
        obsnum = list(range(STARTOBS, min(STARTOBS+NUMOBS, allobs)))
    L.log('Processing Range '+str(obsnum[0])+":"+str(obsnum[-1]), 1)

    #create a reader
    L.log('Creating a data reader... ', 1, end='')

    to_compute = Queue(20)
    reader_is_alive = Queue()
    
    reader = Process(target=reader_f,
                     args=[to_compute, 
                           obs_err_sps, 
                           obsnum, 
                           L, 
                           eigen_land, 
                           eigen_sea, 
                           co2_std, 
                           amethyst_config.processor_vars,
                           inverter.cx.variable_selection]
                     )
    reader.start()
    reader_time=time.time()
    L.log('Done in ' +str(reader_time-inverter_time)+' seconds', 1)

    # Wait for the data
    L.log('Loading data objects... ', 1, end='')
    while to_compute.empty():
        time.sleep(0.05)
    nlev, nselstate, selchannels, transformer = to_compute.get()
    
    L.log('Done in ' +str(time.time()-reader_time)+' seconds', 1)
  
    # Create a queue where the output should be stored
    log_data = Queue()
    # Create a queue for the progress bar
    progress_bar_queue = Queue()
    # Create a progress bar to show the progress
    pb = ProgressBar(len(obsnum), 0, progress_bar_queue)
    
    #create a scriba
    L.log('Creating the scriba process... ', 1, end='')

    to_write    = Queue()
    stop_scriba = Queue()

    eigen_max = max(eigen_sea, eigen_land)

    outfile = workingDir + '/' + OUTPUT_FILE

    scriba = Process(target = scriba_f,
                     args=[to_write, stop_scriba, outfile, obsnum, nlev,
                           eigen_max, nselstate + eigen_max, selchannels,
                           STARTOBS, amethyst_config.processor_vars["output"], progress_bar_queue])
    scriba.start()
    L.log('Done!', 1)

    # Create the processes needed for the computation
    process_list = []

    L.log('Starting computation... ',1, end='')
    for i in range(min(PROCESSES_NUMBER,len(obsnum))):
        p = Process(target=invert_process,
                    args=[inverter,nlev, i, to_compute, to_write, reader_is_alive,
                          logger(VERBOSE, LOG_FILE), transformer,
                          amethyst_config.processor_vars["output"] , log_data, start_process])
        process_list.append((i, p))
        p.start()
    L.log('Running!',1)
   
    start_computing_processes=time.time()
    while True:
        time.sleep(1)
        pb.update()
        pb.draw()

        if any([p.is_alive() for (i, p) in process_list]):
            dead_processes = [p for (i, p) in process_list if not p.is_alive()]
            n_dead_processes = len(dead_processes)
            if len(dead_processes) > 0  and to_compute.empty():
                print(str(n_dead_processes) +
                      " processes have crashed", file=sys.stderr)
                process_list = [(i, p) for (i, p) in process_list if p.is_alive()]
                print("Now there are " + str(len(process_list)) +
                      " processes running", file=sys.stderr)
        else:
            print("All the compute processes are dead!!!", file=sys.stderr)
            break
        if not reader.is_alive():
            reader_is_alive.put(False)
            break


    while any([p.is_alive() for (i, p) in process_list]):
        if not log_data.empty():
            process_log = log_data.get()
            if process_log != "":
                L.log(log_data.get(), 2)

    for i,p in process_list:
        p.join()

    strange_error = False
    if not to_compute.empty():
        strange_error = True

    stop_scriba.put(True)
    scriba.join()

    pb.update()
    pb.draw()

    end_process = time.time()
    total_computing_processes=  end_process -start_computing_processes
    L.log('\nTotal running time of processes : ' + str(total_computing_processes), 1)
    total_time = end_process - start_process
    L.log('\nTotal running time: ' + str(total_time), 1)


    L.close()

    err_comp = list()

    if strange_error:
        while not to_compute.empty():
            try:
                not_comp = to_compute.get(block=False)
                err_comp.append(not_comp)
            except Empty:
                pass
        if len(err_comp) == 0:
            print('A strange behaviour has been detected! Some '
                  'observations could have been skipped!', file=sys.stderr)
        elif len(err_comp) < 20:
            print('The following observations were NOT executed'+
                  ' due to some strange errors:', file=sys.stderr)
            print(str(err_comp), file=sys.stderr)
        else:
            print(str(len(err_comp)) + ' observations were NOT executed'
                  ' due to some strange errors', file=sys.stderr)
        sys.exit(1)

    return


if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', default='None',
                        help="Define which debug variables must be saved")
    parser.add_argument('-l', '--log_file', default='stdout',
                        help="Define where print the log file")
    parser.add_argument('-n', '--numobs', type=int, default=-1,
                        help="Define the max number of observations to be computed")
    parser.add_argument('-p', '--processes', type=int, default=1,
                        help="Define the number of COMPUTING processes")
    parser.add_argument('-s', '--startobs', type=int, default=0,
                        help="Define the number of the first observation to be computed")
    parser.add_argument('-o', '--output', default=amethyst_config.processor_vars["outfile"]["filename"],
                        help="Define the output filename")
    parser.add_argument('-v', '--verbose', type=int, default=2,
                        help="Define the level of verbosity")

    # INPUT 4 FUNCTION
    LOG_FILE         = read_log_file(parser.parse_args().log_file)
    OUTPUT_FILE      = parser.parse_args().output
    NUMOBS           = parser.parse_args().numobs
    PROCESSES_NUMBER = parser.parse_args().processes
    STARTOBS         = parser.parse_args().startobs
    VERBOSE          = parser.parse_args().verbose

    # Launch processor
    processor(LOG_FILE, NUMOBS, PROCESSES_NUMBER, STARTOBS, VERBOSE,
              output_file = OUTPUT_FILE)
