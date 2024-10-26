import sys, getopt
import os
import glob
from datetime import datetime, timedelta
#import datetime
import subprocess


def main():


   # --- Dir definitions
   BASEDIR = '/data/mirto/oper/'
   PYDIR = '/data/mirto/adaptivemeteo_local/validation/'
   LOGDIR = BASEDIR + 'log/validation/'
   INDIR = '/data/mirto/db/combined/'
   OUTDIR = '/data/mirto/db/validation_results/'

   ctime = datetime.utcnow()

   if ctime.hour > 15:
      timestampStr = ctime.strftime("%Y-%m-%dT12:00")
   else:
      timestampStr = ctime.strftime("%Y-%m-%dT00:00")

   #go to input directory
   os.chdir("/data/mirto/db/combined/")

   #clean old inactive links in input directory
   print("Cleaning inactive links in  /data/mirto/db/combined/")
   os.system("find -L /data/mirto/db/combined/ -maxdepth 1 -type l -delete")

   #update links in input directory
   print("Updating links in  /data/mirto/db/combined/")
   os.system("ln -sf /data/mirto/db/cris/2* /data/mirto/db/combined/.")
   os.system("ln -sf /data/mirto/db/iasi/2* /data/mirto/db/combined/.")

   #go to input directory
   os.chdir("/data/mirto/oper/wrkdir/validation/")

   #Build command to be run from qsub without using CM
   cmd = 'python3 ' + PYDIR + 'validate_retrievals.py -r ' + INDIR + ' -s PHTO -d ' + timestampStr + ' -w ' + OUTDIR + ' > ' + LOGDIR + 'mirto_validation_logfile-`date +\%Y\%m\%d\%H\%M\%S`.log'

   print(cmd)

   #Write and submit shell script that runs mirto 
   submit(cmd, runtime=5, cores=1, threads=1, queue_name='ops', jobscript='validate_retrievals_to.sh', output= LOGDIR + '$PBS_JOBNAME.OUT', error= LOGDIR + '$PBS_JOBNAME.ERR')

   #Build command to be run from qsub without using CM
   cmd = 'python3 ' + PYDIR + 'validate_retrievals.py -r ' + INDIR + ' -s PHLI -d ' + timestampStr + ' -w ' + OUTDIR + ' > ' + LOGDIR + 'mirto_validation_logfile-`date +\%Y\%m\%d\%H\%M\%S`.log'

   print(cmd)

   #Write and submit shell script that runs mirto 
   submit(cmd, runtime=5, cores=1, threads=1, queue_name='ops', jobscript='validate_retrievals_li.sh', output= LOGDIR + '$PBS_JOBNAME.OUT', error= LOGDIR + '$PBS_JOBNAME.ERR')



def submit(command, runtime, cores, threads, queue_name, directory='', modules='',jobscript='jobscript',
    output='/dev/null', error='/dev/null'):
    """
    Function to submit a job to the Queueing System - with jobscript file
    Parameters are:
    command:   The command/program you want executed together with any parameters.
               Must use full path unless the directory is given and program is there. 
    directory: Working directory - where should your program run, place of your data.
               If not specified, uses current directory.
    modules:   String of space separated modules needed for the run.
    runtime:   Time in minutes set aside for execution of the job.
    cores:     How many cores are used for the job.
    ram:       How much memory in GB is used for the job.
    group:     Accounting - which group pays for the compute.
    jobscript: Standard name for the jobscript that needs to be made.
               You should number your jobscripts if you submit more than one.
    output:    Output file of your job.
    error:     Error file of your job.
    """
    runtime = int(runtime)
    cores = int(cores)
    #ram = int(ram)
    if cores > 32:
        print("Can't use more than 32 cores on a node")
        sys.exit(1)
    #if ram > 64:
        #print("Can't use more than 64 GB on a node")
        #sys.exit(1)
    if runtime < 1:
        print("Must allocate at least 1 minute runtime")
        sys.exit(1)
    minutes = runtime % 60
    hours = int(runtime/60)
    walltime = "{:d}:{:02d}:00".format(hours, minutes)
    if directory == '':
        directory = os.getcwd()

    # Making a jobscript
    script = '#!/bin/sh\n'
    script += '#PBS -e ' + error + ' -o ' + output + '\n'
    script += '#PBS -d ' + directory + '\n'
    script += '#PBS -l nodes=1:ppn=' + str(cores) + '\n'
    script += '#PBS -l walltime=' + walltime + '\n'
    script += 'export OMP_NUM_THREADS='  + str(threads) + '\n'
    script += 'echo "script running on host `hostname`"' + '\n'

    script += 'source activate gdal' + '\n'

    script += 'export PYTHONPATH=/data/mirto/adaptivemeteo_local/adaptive_tools' + '\n'


    if modules != '':
        script += 'module load ' + modules + '\n'
    script += command + '\n'
    if not jobscript.startswith('/'):
        jobscript = directory + '/' + jobscript
    with open(jobscript, 'wt') as jobfile:
        jobfile.write(script)
    print(script)
    print("job = subprocess.run([qsub -q ops, jobscript],stdout=subprocess.PIPE, universal_newlines=True)")

    # The submit
    job = subprocess.run(['qsub', '-q', queue_name, jobscript],stdout=subprocess.PIPE, universal_newlines=True)
    jobid = job.stdout.split('.')[0]
    return jobid

if __name__ == '__main__':
    main()

