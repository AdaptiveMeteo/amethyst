#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:19:59 2020

@author: paolos
"""
from datetime import datetime, timedelta
import os
from argparse import ArgumentParser

def parse_gmtco_date(filename,time):
    if time == 'start':
        return datetime.strptime('_'.join([ filename.split('_')[2].replace('d',''),filename.split('_')[3][:-1].replace('t','') ]), '%Y%m%d_%H%M%S')
    elif time == 'end':
        return datetime.strptime('_'.join([ filename.split('_')[2].replace('d',''),filename.split('_')[4][:-1].replace('e','') ]), '%Y%m%d_%H%M%S')
    else:
        return None

def parse_jrr_date(filename,time):
    if time == 'start':
        return datetime.strptime( filename.split('_')[3].replace('s','')[:-1],  '%Y%m%d%H%M%S')    
    elif time == 'end':
        return datetime.strptime( filename.split('_')[4].replace('e','')[:-1], '%Y%m%d%H%M%S')
    else:
        return None
    
def select_viirs_file(cris_start_date, cris_end_date, viirs_path, outdir):
    for file in os.listdir(viirs_path):
        if 'GMTCO' in file:
            filename = os.path.basename(file)
            start_time = parse_gmtco_date(filename,'start') 
            end_time   = parse_gmtco_date(filename, 'end')
            
            if end_time < start_time:
                end_time += timedelta( days = 1)

            if not (start_time < cris_start_date and end_time < cris_start_date ) and \
               not (start_time > cris_end_date and end_time > cris_end_date ):
                   
                try:
                    os.system("cp {} {}/".format(viirs_path +'/'+filename,outdir))
                    print("VIIRS file {} copied in {}".format(filename,outdir))
                except:
                    print("Cannot copy VIIRS file {}".format(filename))
                    pass
                viirs_cm_file = [ vfile for vfile in os.listdir(viirs_path) if 'JRR' in vfile ]
                for cm_file in viirs_cm_file:
                    cm_filename = os.path.basename(cm_file)
                    
                    # Parse starting and ending date of the cloudmask file
                    start_cm_time =  parse_jrr_date(cm_filename, 'start')
                    end_cm_time   =  parse_jrr_date(cm_filename, 'end')
                    if end_cm_time < start_cm_time:
                        end_cm_time += timedelta( days = 1)

                    # If CM file fall inside the GEO VIIRS temporal window add it to the overpass dir
                    if not (start_cm_time < start_time and end_cm_time < end_time) and \
                       not (start_cm_time > end_time and end_cm_time > end_time):
                        os.system("cp {} {}/".format( viirs_path +'/'+cm_file, outdir ))
                        print("Copied Cloud Mask file {} in {}".format( cm_filename, outdir ))
                                
    return

def check_jrr_continuity(start_dates, end_dates):
    out_flag = False
    for i, end_date in enumerate(end_dates[:-1]):
        if abs((end_date - start_dates[i+1]).total_seconds()//60) > 1:
            out_flag = True
    return out_flag
        

def main():
    import numpy as np
    parser = ArgumentParser()
    parser.add_argument('--viirs_dir', '-vd',type=str,
                        help='Path to viirs files')
    parser.add_argument('--cris_dir','-cd',type=str,
                        help='Path to cris files')
    parser.add_argument('--workdir','-w',type=str,
                        help='Workdir')
    argv = parser.parse_args()
    
    # Untar all viirs file
    print('Unpacking VIIRS Cloud mask files...')
    tar_files = [ file for file in os.listdir(argv.viirs_dir) if '.tar' in file]
    for tar_file in tar_files:
        os.system("tar xfv {} -C {}".format(argv.viirs_dir+'/'+tar_file,argv.viirs_dir))
        os.system("rm {}".format(argv.viirs_dir+'/'+tar_file))
    
    # Prepare overpass directories
    cris_native_files = [ file for file in os.listdir(argv.cris_dir) if 'GCRSO' in file]
    for native in cris_native_files:
        overpass_name = '_'.join(os.path.basename(native).split('_')[2:4]).replace('d','').replace('t','')[:-1]
        start = datetime.strptime(overpass_name,'%Y%m%d_%H%M%S')
        end   = datetime.strptime( '_'.join([os.path.basename(native).split('_')[2] ,os.path.basename(native).split('_')[4] ]).replace('d','').replace('e','')[:-1],'%Y%m%d_%H%M%S')
        if end < start:
            end += timedelta( days = 1)
        
        try:
            os.mkdir(argv.workdir +'/' + overpass_name)
        except FileExistsError:
            pass
        os.system("cp {} {}".format(argv.cris_dir +'/'+ native,argv.workdir + '/' + overpass_name))
        os.system("cp {}*.h5 {}".format(argv.cris_dir + '/'+native.replace('GCRSO','SCRIF')[:-35],argv.workdir + '/' + overpass_name))

        print('Created dir {}'.format(overpass_name))   
        select_viirs_file(start, end, argv.viirs_dir, argv.workdir + '/' + overpass_name)
    
    print()    
    print("########################################")
    print("Validating results...")
    print()
    # Check 4 good overpass dirs
    
    good_ov = 0
    missing_fov_ov = 0
    total_ov = len(os.listdir(argv.workdir))
    
    for overpass_dir in os.listdir(argv.workdir):

         start_overpass = datetime.strptime(overpass_dir,'%Y%m%d_%H%M%S')
         gcrso_filename = [name for name in os.listdir(argv.workdir+'/'+overpass_dir) if 'GCRSO' in name][0]
         end_overpass   = datetime.strptime( '_'.join([os.path.basename(gcrso_filename).split('_')[2] ,os.path.basename(gcrso_filename).split('_')[4] ]).replace('d','').replace('e','')[:-1],'%Y%m%d_%H%M%S')
         
         gmtco_files = [ gmtco_file for gmtco_file in os.listdir(argv.workdir+'/'+overpass_dir) if 'GMTCO' in gmtco_file ]
         gmtco_files.sort()

         jrr_list = [ name for name in os.listdir('/'.join([argv.workdir,overpass_dir])) if 'JRR' in name ]
         jrr_list.sort()    
         cm_start_dates = np.array( [ parse_jrr_date(os.path.basename(name),'start') for name in jrr_list  ])
         cm_end_dates   = np.array( [ parse_jrr_date(os.path.basename(name),'end') for name in jrr_list ])
         
         skip_bool = False
         
         if parse_gmtco_date(gmtco_files[0],'start') > start_overpass or parse_gmtco_date(gmtco_files[-1],'end') < end_overpass:

             # Check if gmtco cover the overpass temporal window
             print('Wrong GMTCO covarage in overpass {}'.format(overpass_dir))
             missing_fov_ov += 1
             # Uncomment to block overpass with missing FOVs for the CM
             skip_bool = True
         else:

            for geo_viirs_file in gmtco_files:
                 
                 filename = os.path.basename(geo_viirs_file)
                 start_time = parse_gmtco_date(filename,'start') 
                 end_time   = parse_gmtco_date(filename,'end')
                 if end_time < start_time:
                     end_time += timedelta( days = 1)

                 # Check if jrr files cover the gmtco temporal window                
                 min_start_dif = 999
                 min_end_dif = 999
                 for i,jrr_starting_date in enumerate(cm_start_dates):
                        dt_start = abs((jrr_starting_date - start_time).total_seconds()//60) 
                        dt_end   = abs((cm_end_dates[i] - end_time).total_seconds()//60) 
                        if  dt_start < min_start_dif:
                            min_start_dif = dt_start
                        if  dt_end < min_end_dif:
                            min_end_dif = dt_end                            
                 if min_start_dif > 5 or min_end_dif > 5:
                     print('Wrong JRR covarage in overpass {}'.format(overpass_dir))
                     skip_bool = True
                 else:
                     indx = np.where( (cm_start_dates >= start_time) & (cm_end_dates <= end_time) )[0]
                     if len(indx) != 4 :
                         # Check if the number of JRR files for gmtco is bigger than 2
                         print('Wrong number of JRR in overpass {} (missing FOVs for the Cloud Mask)'.format(overpass_dir))
                         skip_bool = True
                         
         
         if skip_bool:
             print('Overpass {} not ready for CM processing'.format(overpass_dir))
         else:
             # Check JRR files continuity
             skip_bool = check_jrr_continuity(cm_start_dates, cm_end_dates)
             if skip_bool:
                 print('JRR files in overpass {} are not contiguous'.format(overpass_dir))
                 print('Overpass {} not ready for CM processing'.format(overpass_dir))
             else:                
                 print('Overpass {} READY for CM processing'.format(overpass_dir))

                 os.system("touch {}/ready4cmprocessing".format(argv.workdir + '/' + os.path.basename(overpass_dir)))
                 good_ov +=1    
                         
    print()   

    print('Overpasses with missing FOVs: {} out of  {}'.format(missing_fov_ov,total_ov))
    print('Total number of ready overpasses: {} out of  {}'.format(good_ov,total_ov))

if __name__ == '__main__':
    main()  