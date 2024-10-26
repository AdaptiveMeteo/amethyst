import datetime
import os
import numpy as np

2022082000
2022082012

ini_date = datetime.datetime(2022, 8, 8, 0, 0, 0)
fin_date = datetime.datetime(2022, 8, 8, 12, 0, 0)

#blackdates = [datetime.datetime(2020,11,20,0,0,0),datetime.datetime(2020,11,23,3,0,0),datetime.datetime(2020,11,26,6,0,0),datetime.datetime(2020,11,29,9,0,0),datetime.datetime(2020,12,2,12,0,0),datetime.datetime(2020,12,5,15,0,0),datetime.datetime(2020,12,8,18,0,0),datetime.datetime(2020,12,11,21,0,0),datetime.datetime(2020,11,15,0,0,0),datetime.datetime(2020,12,5,18,0,0),datetime.datetime(2020,12,5,21,0,0),datetime.datetime(2020,12,6,6,0,0),datetime.datetime(2020,12,6,9,0,0),datetime.datetime(2020,12,6,18,0,0),datetime.datetime(2020,12,6,21,0,0),datetime.datetime(2020,12,7,6,0,0),datetime.datetime(2020,12,7,9,0,0),datetime.datetime(2020,12,7,18,0,0),datetime.datetime(2020,12,7,21,0,0),datetime.datetime(2020,12,8,6,0,0),datetime.datetime(2020,12,8,9,0,0)] 

WRFDIR='/mnt/miwa/amethyst/amethyst_test_data/input_data/wrf/XXXXXXXXXX/'
OUTDIR='/mnt/miwa/amethyst/amethyst_test_data/input_data/wrf/'

PYDIR='/home/amethyst/amethyst/postprocessor/validation/'
PYPLOTDIR='/home/amethyst/amethyst/postprocessor/profiles/'

date = ini_date
while date  < fin_date :

   #if date in blackdates:
   #     date += datetime.timedelta(hours=3)
   #     continue
   print(date.strftime('%Y%m%d_%H%M%S'))

   WRFTIME = date.strftime('%Y%m%d%H')
   WRFASSDIR =  WRFDIR.replace('XXXXXXXXXX',WRFTIME)

   cmd = 'python ' + PYDIR + 'validate_wrf_single.py --wrf_dir ' + WRFASSDIR + ' -o ' + OUTDIR
   print(cmd)

   os.system(cmd)

   date += datetime.timedelta(hours=3)


dt = 0
while dt  < 7 :

   cmd = 'python ' + PYDIR + 'validate_wrf_single.py --wrf_dir ' + OUTDIR + ' -o ' + OUTDIR + ' --forecast_h ' + str(dt) + ' --only_agg True'
   print(cmd)

   os.system(cmd)

   cmd = 'python ' + PYPLOTDIR + 'mirto_plot_validation_profiles_single.py --i ' + OUTDIR + ' --o ' + OUTDIR + ' --m wrf --forecast_hour ' + str(dt)
   print(cmd)

   os.system(cmd)

   dt += 1
