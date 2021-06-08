from __future__ import print_function, division
from netCDF4 import Dataset
from numpy import exp
import traceback
import sys

# Get python major version
py_version = int(sys.version_info[0])
if py_version > 2:
    from queue import Empty    # @UnresolvedImport @UnusedImport
else:
    from Queue import Empty    # @UnresolvedImport @Reimport


class ProgressBar(object):
    def __init__(self, tot_values, start_value, update_queue):
        self.__tot_val = tot_values
        self.__val = start_value
        self.__update_queue = update_queue

    def update(self):
        while not self.__update_queue.empty():
            self.__val = self.__update_queue.get()
    
    def draw(self):
        percentage = self.__val / self.__tot_val *100
        perc_int = (self.__val*100) // self.__tot_val
        bar_position = perc_int // 2
        mod = perc_int % 2
        if mod == 0 and perc_int != 100:
            bar_symbol = '_'
        elif mod == 1:
            if py_version == 3:
                bar_symbol = '\u258C'
            else:
                bar_symbol = '+'
        else:
            bar_symbol = ''
        if py_version == 3:
            bar_string =  '\u2588' * bar_position
        else:
            bar_string =  '*' * bar_position
        bar_string += bar_symbol
        bar_string += '_' * (50-bar_position-1)

        print('\r', end='')
        if py_version == 3:
            sys.stdout.buffer.write(bar_string.encode('UTF-8'))
        else:
            print(bar_string, end='')
        print('  {:.2f}%     '.format(percentage), end = '')
        sys.stdout.flush()
    


def scriba_f(to_write, stop_now, root_file, obsnum, nlev, nem, nselstate,
             selchannels, STARTOBS, output_vars, progress_bar_queue):

    with Dataset(root_file, 'w', format='NETCDF4') as rootgrp:
        """ Create the structure of the file test.nc """
        rootgrp.createDimension('obsnum', len(obsnum))
        rootgrp.createDimension('levels', nlev) 
        #PaoloA
        #rootgrp.createDimension('twice_levels', 2*nlev)
        #rootgrp.createDimension('twice_levels', 242)
        rootgrp.createDimension('twice_levels', 162)
        rootgrp.createDimension('coefficients', nem)
        rootgrp.createDimension('mnel', 20)
        rootgrp.createDimension('selchannels', selchannels)
        rootgrp.createDimension('numselstatevar', nselstate)
        v_tp = 'f4'
        v_tp2 = 'f8'
        dims = ('obsnum','levels')
        f_v = 1.0E20
        nco  =  rootgrp.createVariable('obs', 'i4', dims[0], -1)
        if 'pressure' in output_vars:
            ncp  =  rootgrp.createVariable('p',   v_tp, dims,    f_v)
        if 'temperature' in output_vars:
            nct  =  rootgrp.createVariable('t',   v_tp2, dims,    f_v)
        if 'water_vapor' in output_vars:
            ncq  =  rootgrp.createVariable('q',   v_tp2, dims,    f_v)
        if 'ozone' in output_vars:
            nco3  = rootgrp.createVariable('o3',  v_tp, dims,    f_v)
        if 'surface_temperature' in output_vars:
            ncskt = rootgrp.createVariable('skt', v_tp, dims[0], f_v)
        if 'd2' in output_vars:
            ncd2  = rootgrp.createVariable('d2',  v_tp, dims[0], f_v)
        if 'surface_emissivity_coefficients' in output_vars:
            ncems = rootgrp.createVariable('ems_coeff', v_tp,
                                           ('obsnum','coefficients'), f_v)
        if 'jacobian' in output_vars:
            jacobian = rootgrp.createVariable('jacobian', v_tp,
                            ('obsnum', 'selchannels', 'numselstatevar'), f_v)
        if 'residuals' in output_vars:
            residuals = rootgrp.createVariable('residuals', v_tp, 
                                               ('obsnum', 'selchannels'), f_v)
        #PaoloA 12112018
        if 'fgresiduals' in output_vars:
            fgresiduals = rootgrp.createVariable('fgresiduals', v_tp, 
                                               ('obsnum', 'selchannels'), f_v)
        if 'da_r' in output_vars:
            DA_R  = rootgrp.createVariable('DA_R', 'i2', dims[0], 0)
        if 'transformed_retrievals' in output_vars:
            DA_Yret=rootgrp.createVariable('DA_Yret', v_tp2, ('obsnum', 'mnel'), f_v)
        if 'assimilation_operator' in output_vars:
            DA_Hret=rootgrp.createVariable('DA_Hret', v_tp2, 
                                           ('obsnum', 'mnel', 'twice_levels'), f_v)
        #PaoloA
        if 'sn_eigenvalues' in output_vars:
            DA_Lambda  = rootgrp.createVariable('DA_Lambda', v_tp2, ('obsnum', 'mnel'), f_v)
        if 'Sa' in output_vars:
            Sa_ret = rootgrp.createVariable('Sa_ret',  v_tp2, ('obsnum', 'numselstatevar', 'numselstatevar'), f_v)
            SaInv_ret = rootgrp.createVariable('SaInv_ret',  v_tp2, ('obsnum', 'numselstatevar', 'numselstatevar'), f_v)

        rootgrp.sync()

        printed = 0
        while stop_now.empty():
            while not to_write.empty():
                get_output = False
                try:
                    output = to_write.get(block=True, timeout=1)
                    get_output = True
                except Empty:
                    pass
                if get_output:
                    obs = output[0]
                    nco[obs-STARTOBS] = obs
                    sol = output[1]
                    
                    # If the observation is not convergent, do not
                    # write anything
                    if sol is None:
                        continue
                    try:
                        if 'pressure' in output_vars:
                            ncp[obs-STARTOBS, :] = sol['pressure']
                        if 'temperature' in output_vars:
                            nct[obs-STARTOBS, :] = sol['temperature']
                        if 'water_vapor' in output_vars:
                            ncq[obs-STARTOBS, :] = exp(sol['water_vapor'])*1e3
                        if 'ozone' in output_vars:
                            nco3[obs-STARTOBS, :] = exp(sol['ozone'])*1e3
                        if 'surface_temperature' in output_vars:
                            ncskt[obs-STARTOBS] = sol['surface_temperature']
                        if 'surface_emissivity_coefficients' in output_vars:
                            len_sc = sol['surface_emissivity_coefficients'].size
                            ncems[obs-STARTOBS, :len_sc] = sol['surface_emissivity_coefficients']
                        if 'd2' in output_vars:
                            ncd2[obs-STARTOBS] = sol['d2']
                        if 'da_r' in output_vars:
                            DA_R[obs-STARTOBS] = sol['da_r']
                        if 'transformed_retrievals' in output_vars:
                            DA_Yret[obs-STARTOBS, :] = sol['transformed_retrievals']
                        #PaoloA
                        if 'sn_eigenvalues' in output_vars:
                            DA_Lambda[obs-STARTOBS, :] = sol['sn_eigenvalues']
                        if 'assimilation_operator' in output_vars:
                            DA_Hret[obs-STARTOBS, :, :] = sol['assimilation_operator']
                        if 'jacobian' in output_vars:
                            #PaoloA 14112018
                            #print(" shape {}".format(jacobian[obs-STARTOBS, :, :].shape))
                            jacobian[obs-STARTOBS, :, :] = sol['jacobian']
                        if 'residuals' in output_vars:
                            residuals[obs-STARTOBS, :] = sol['residuals']
                        #PaoloA 12112018
                        if 'fgresiduals' in output_vars:
                            fgresiduals[obs-STARTOBS, :] = sol['fgresiduals']
                        if 'Sa' in output_vars:
                            Sa_ret[obs-STARTOBS, :, :] = sol['Sa_ret']
                            SaInv_ret[obs-STARTOBS, :, :] = sol['SaInv_ret']
                    except:
                        traceback.print_exc(file=sys.stdout)
                        print('Error writing observation ' + str(obs) + ' into output file!') 
                        raise RuntimeError('Exiting...')
                    rootgrp.sync()
                    printed += 1
                    progress_bar_queue.put(printed)
