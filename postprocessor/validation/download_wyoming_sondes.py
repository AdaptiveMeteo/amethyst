import numpy as np
from os import system
import sys
import logging

__author__ = "Paolo Scaccia, Paolo Antonelli and Stefano Piani"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__ = ["Paolo Scaccia", "Paolo Antonelli"]
__license__ = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

LOGGER = logging.getLogger(__name__)

class raob_data(object):
    
    def __init__(self,MATRIX,DATE,LAT,LON,ELEV):

        CONV_TEMP             = 273.15   # Convert temperature from Celsius to Kelvin
        self.levels           = np.array(MATRIX[:,0])
        self.h                = np.array(MATRIX[:,1])
        self.temp             = np.array(MATRIX[:,2]) + CONV_TEMP
        self.dwp_temp         = np.array(MATRIX[:,3]) + CONV_TEMP
        self.relh             = np.array(MATRIX[:,4])
        self.mixr             = np.array(MATRIX[:,5])
        self.wind_dir         = np.array(MATRIX[:,6])
        self.wind_speed       = np.array(MATRIX[:,7])
        self.pot_temp         = np.array(MATRIX[:,8])
        self.eq_pot_temp      = np.array(MATRIX[:,9])
        self.virtual_pot_temp = np.array(MATRIX[:,10])
        self.lat              = np.float(LAT)
        self.lon              = np.float(LON)
        self.elev          = np.float(ELEV)
        self.date          = '20'+DATE[:2]+'-'+DATE[2:4]+'-'+DATE[4:6]+'T'+DATE[7:9]+':'+DATE[9:]

def parse_raob_data(string):
        """
        INPUT 
            string: url with university of wyoming website
        OUTPUT
            sondes: a raob_data class with rowinsondes
        """
        from urllib.request import urlopen
        import re


        # parse txt from website and devide sections
        try:
            sections = repr(urlopen(string.replace("\'","")).read()).split('<PRE>')
        except Exception as error:
            LOGGER.info("Rowinsonde not found for station: "+string.split('STNM=')[-1].split('&')[0])
            LOGGER.debug('Error: {}'.format(error))
            return None
        
        # clean rows and define a matrix from the second section
        try:
            data_rows = sections[1].split('\\n')[5:-1]
            data_rows = np.array([ [x for x in row.replace('       ',' 0').split(' ') if x != ''] for row in data_rows],dtype=float)
        except:
            LOGGER.info("Rowinsonde not found for station: "+string.split('STNM=')[-1].split('&')[0])
            return None
        # extract station info from the last section
        info_rows = sections[2].split('\\n')
        
        
        try:
            return raob_data(data_rows,
                      re.findall("\d+/\d+",info_rows[3])[0],      # DATE
                      re.findall("[+-]?\d+.\d+",info_rows[4])[0], # LAT
                     re.findall("[+-]?\d+.\d+",info_rows[5])[0], # LON
                     re.findall("\d+.\d+",info_rows[6])[0])      # ELEV
        except:
            return raob_data(data_rows,
                      re.findall("\d+/\d+",info_rows[2])[0],      # DATE
                      re.findall("[+-]?\d+.\d+",info_rows[3])[0], # LAT
                     re.findall("[+-]?\d+.\d+",info_rows[4])[0], # LON
                     re.findall("\d+.\d+",info_rows[5])[0])      # ELEV

def download_sondes(station,date,*args):
    """
    INPUT:
        - station: class 'station_class' defined in read_stations_list.py
        - date:    date in datetime format
    OUTPUT:
        - raob_ICAO_YYYYTHH00_UTC.txt : raob file inside outdir (current directory if not passed)

    """
    # assign outdir if passed otherwise defined as './'

    outdir = '' if len(args) == 0 else args[0]+"/"

    year  = date.strftime('%Y')
    month = date.strftime('%m')
    day   = date.strftime('%d')

    # select closest hour
    allowed_hours = np.array([0,12,24])
    hour = date.hour if (date.hour in allowed_hours) else allowed_hours[int(np.around(date.hour/12))]
    if hour == 24:
        hour = 0
    hour = '%02d' % hour

    # select region
    region = station.region
    if region == None:
        region = ''
    string = "http://weather.uwyo.edu/cgi-bin/sounding?region="+region+"f&TYPE=TEXT%3ALIST&YEAR="+year+"&MONTH="+month+"&FROM="+day+hour+"&TO="+day+hour+"&STNM="+station.synop+"&REPLOT=1"

    if region == '':
        string = string.replace('region=','')

    # if specified save raob file
    if outdir != '':
        # clean outdir string
        outdir = outdir.replace('//','/').replace('././','./')

        # compose string for output filename
        date_string = date.date().isoformat()+'T'+hour+'00'
        output_file = outdir+"raob_"+station.icao+"_"+date_string+"_UTC.txt"

        # download raob file
        cmd = "wget -O "+output_file+" "+string
        system(cmd)

    return parse_raob_data(string)


def main(argv):
    """
    Main function: used if the script is called by shell and not imported.
    It saves raob_file inside outdir if '-o' or '--outdir' is passed,
    otherwise the file is saved in the current directory.

    INPUT:   - (optional) outdir 
         - date :   date with format "YYYY-MM-DDTHH:MM" (T is only a letter)
         - station: station code (ICAO , IATA or SYNOP)

    """
    from read_stations_list import stations
    from datetime import datetime
    import getopt


    outdir  = ''
    date    = ''
    station = ''

    opts, args = getopt.getopt(argv,"o:d:s:",["outdir=","date=","station="])

    for opt,arg in opts:
        if opt in ('-o','--outdir'):
            outdir = arg
        elif opt in ('-d','--date'):
            date = arg
        elif opt in ('-s','--station'):
            station = arg

    if date != '':
        date = datetime.strptime(date,'%Y-%m-%dT%H:%M')
    else:
        sys.exit("I need a date: (-d or --date=) YYYY-MM-DDThh:mm" )

    if station != '':
        station = stations.get(station)
    else:
        sys.exit("I need a station: (-s or --station=) ICAO, IATA or SYNOP" )

    raob = download_sondes(station,date,outdir) if outdir != '' else download_sondes(station,date)
    print()
    print("PRES HGHT TEMP DWPT RELH")
    for i in range(raob.temp.size):
        print("{} {} {} {} {}".format(raob.levels[i],raob.h[i],raob.temp[i],raob.dwp_temp[i],raob.relh[i]))

if __name__=="__main__":
    main(sys.argv[1:])

