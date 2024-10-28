import logging
import sys
import numpy as np
from os import system
from os.path import exists
sys.path.append('/home/amethyst/amethyst/postprocessor')
from validation.regions_dictionary import regions

__author__ = "Paolo Scaccia, Paolo Antonelli and Stefano Piani"
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__ = ["Paolo Scaccia", "Paolo Antonelli"]
__license__ = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Scaccia"
__email__ = "paolo.scaccia@adaptivemeteo.com"
__status__ = "Development"

LOGGER = logging.getLogger(__name__)

def download_stations_list():
    """
    Downloads stations.txt from url and cleans the file (removes 
    commments and blank lines).

    TO DO: add workdir option

    """
    url = 'http://weather.rap.ucar.edu/surface/stations.txt'
    system('wget '+url)
    system("awk '$1 !~ /!/ {print}' stations.txt > stations_tmp.txt")
    system("awk '{if (length==83) {print} }' stations_tmp.txt > stations.txt")
    system("rm -f stations_tmp.txt stations.txt.*")
    
    # Correct station list (PaoloS 27-07-2021)
    system("sed -i '/PHTO/ s/43N/72N/' stations.txt")   # Hilo Lat
    system("sed -i '/PHTO/ s/03W/05W/' stations.txt")   # Hilo Lon
    system("sed -i '/PHLI/ s/58N/99N/' stations.txt")   # Lihue Lat
    system("sed -i '/PHLI/ s/20W/34W/' stations.txt")   # Lihue Lon
    

class station_class(object):
    """
    Defines a class for a single station containing every relevant field
    """
    def __init__(self,CD,STATION,ICAO,IATA,SYNOP,LAT,LON,ELEV,COUNTRY):
        self.cd        = CD
        self.name      = STATION
        self.icao      = ICAO
        self.iata      = IATA
        self.synop     = SYNOP
        self.lat       = LAT
        self.lon       = LON
        self.elev      = ELEV
        self.country   = COUNTRY
        self.region    = regions.get(SYNOP)
        self.lat_unit  = 'degrees-north'
        self.lon_unit  = 'degrees-east'
        self.region    = regions.get(SYNOP)

class stations_list(object):
    """
    It creates a class containg the list of all station. It also allows 
    the search for a particular station using ICAO, IATA or SYNOP.
    """

    def __init__(self):
        
        stations_file = open(__file__.replace('read_stations_list.py','stations.txt'),'r')
        rows = [x for x in stations_file.read().split('\n')[:-1] if x[71]=='X']
        stations_file.close()
        
        self.CD      =  []
        self.STATION =  []
        self.ICAO    =  []
        self.IATA    =  []
        self.SYNOP   =  []
        self.LAT     =  []
        self.LON     =  []
        self.ELEV    =  []
        self.COUNTRY =  []
        self.missing_stations = False
        
        # Fill inner data structure with rows read from stations.txt
        for row in rows:
            self.CD.append(row[:2])
            self.STATION.append(row[3:19])
            self.ICAO.append(row[20:24])
            self.IATA.append(row[26:29])
            self.SYNOP.append(row[32:37])
            if row[44] == 'N':
                self.LAT.append(float(".".join(row[39:44].split(' '))))
            else:
                self.LAT.append(-1*float(".".join(row[39:43].split(' '))))
            if row[53] == 'E':
                self.LON.append(float(".".join(row[47:52].split(' '))))
            else:
                self.LON.append(-1*float(".".join(row[47:53].split(' '))))
            self.ELEV.append(int(row[57:61]))
            self.COUNTRY.append(row[-2:])
            
        # Read missing stations from missing_stations.txt
        self.read_missing_stations() 
        
    def get(self,string):

        # check length and assign relative field for the search
        if len(string) == 3:
            search_list = self.IATA
        elif len(string) == 4:
            search_list = self.ICAO
        elif len(string) == 5:
            search_list = self.SYNOP
        else:
            raise ValueError("Can't find a station with "+string+": Unknown field")

        # IMPROVE THIS PART!
        index = np.where(np.array(search_list) == string)[0]

        # check if it's in the list, if True return the result as a 
        # station class
        if not index.size and not self.missing_stations:
                         """
                            Create missing stations dictionary 
                         """
                         LOGGER.info('Station {} is not in stations.txt!'.format(string))

                         station = self.read_missing_stations(string)
                         
                         if station == None:
                             sys.exit('Station ' + string + ' is not in any stations list!\n'
                                                  'Check if the station has rowinsondes: '
                                                 'http://weather.rap.ucar.edu/surface/stations.txt')
                         else:
                             LOGGER.info("Station {} found in missing_stations.txt".format(string))
                             return station                           

        elif index.size:
            index = index[0]

            station = station_class(self.CD[index],
                        self.STATION[index],
                        self.ICAO[index],
                        self.IATA[index],
                        self.SYNOP[index],
                        self.LAT[index],
                        self.LON[index],
                        self.ELEV[index],
                        self.COUNTRY[index])

            return station

        else:
            sys.exit('Station ' + string + ' is not in any stations list!\n'
                     'Check if the station has rowinsondes: '
                    'http://weather.rap.ucar.edu/surface/stations.txt')
    
    def read_missing_stations(self,*string):
            LOGGER.info("Creating dictionary for missing stations...")
   
            missing_stations_file = open(__file__.replace('read_stations_list.py','missing_stations.txt'),'r')   
            rows = [x for x in missing_stations_file.read().split('\n')[:-1] ]
            missing_stations_file.close()

            station = None
            for row in rows:
              
                # If given SYNOP return the station class
                if len(string) > 0:
                    if row[33:38] == string[0]:
                           station = station_class(row[:2],
                                                   row[3:19],
                                                   row[20:24],
                                                   row[26:29],
                                                   row[33:38],
                                                   float(".".join(row[39:46].split(' '))),
                                                   float(".".join(row[48:55].split(' '))),
                                                   row[57:61],
                                                   row[-2:] )
                           if row[46] == 'S':
                               station.lat *= -1
                           if row[55] == 'W':
                               station.lon *= -1      
                           station.region = regions.get(station.synop)
   
            
                self.CD.append(row[:2])
                self.STATION.append(row[3:19])
                self.ICAO.append(row[20:24])
                self.IATA.append(row[26:29])
                self.SYNOP.append(row[33:38])
                if row[46] == 'N':
                    self.LAT.append(float(".".join(row[39:46].split(' '))))
                else:
                    self.LAT.append(-1*float(".".join(row[39:46].split(' '))))
                if row[55] == 'E':
                    self.LON.append(float(".".join(row[48:55].split(' '))))
                else:
                    self.LON.append(-1*float(".".join(row[48:55].split(' '))))
                try:
                    self.ELEV.append(int(row[57:61]))
                except:
                    self.ELEV.append('--')                                 
                self.COUNTRY.append(row[-2:])
                
            self.missing_stations = True

            return station
        
    def get_nearest_station(self,lat,lon,R = 75):
        """
            Return nearest station to given point

        Parameters
        ----------
            lat : float
            lon : float
            R:    float for the maximum  radius for the search area

        Returns
            station: station_class 
            dist   : float, distance betwen FOV and station
            
        """
        from geometry.earth_geometry import haversine
        
        # Read missing stations if not already done
        if not self.missing_stations:
            self.read_missing_stations()
        
        dists = np.array([ haversine(stat_lon, self.LAT[i], lon, lat) for i,stat_lon in enumerate(self.LON)  ])
        dists = np.ma.masked_where( dists > R , dists)

        if dists.mask.all():
            LOGGER.debug("No station found in an area within {:.2f} km radius".format(R))           
            return None, np.nan
            
        index = np.argmin(dists)
        # print('Minimum distance found is {:.3f} km'.format(dists[index]))
        
        return station_class(self.CD[index],
                             self.STATION[index],
                             self.ICAO[index],
                             self.IATA[index],
                             self.SYNOP[index],
                             self.LAT[index],
                             self.LON[index],
                             self.ELEV[index],
                             self.COUNTRY[index]), dists[index]


def read_missing_stations_coord():
    """
        NEVER CALL THIS FUNCTION IF NOT NECESSARY
    """
    
    from datetime import datetime
    from download_wyoming_sondes import download_sondes

    ref_date = datetime.strptime('2020-11-02T12:00','%Y-%m-%dT%H:%M')
    
    ref_file = open('missing_stations.txt','r')
    ref_rows = [x for x in ref_file.read().split('\n')[:-1] ]

    station_file = open('stations.txt','r')
    rows = [x for x in station_file.read().split('\n')[:-1] if x[71]=='X']
    stations = stations_list(rows)

    station_file.close()
    ref_file.close()
    
    out_rows = []
    out_file = open('missing_stations_2.txt','w')

    for i,row in enumerate(ref_rows):
            try:

                raob = download_sondes( stations.get(row[32:37]) , ref_date )
                lat = raob.lat
                lon = raob.lon
                
                lat_str = repr(abs(lat)) + 'N'
                lon_str = repr(abs(lon)) + 'E'

                if lon < 0:
                    lon_str[-1] = 'W'
                if lat < 0:
                    lat_str[-1] = 'S'
                    
                lon_str = lon_str.replace('.',' ')
                lat_str = lat_str.replace('.',' ')

                tmp_str =  row[:39] + lat_str + '  ' + lon_str + row[54:]
                if len(tmp_str) < 82:
                    tmp_str =  row[:9] + ' ' + row[9:39] + lat_str + '  ' + lon_str + row[54:]
                               
                out_rows.append( tmp_str )
                if len(ref_rows[i]) < 82:
                    out_rows[i] = out_rows[-1][:9] + ' ' + out_rows[-1][9:]

                out_file.write((out_rows[-1]))
                out_file.write('\n')

            except:
                    LOGGER.info("Skipped station {}!".format(row[32:37]))
                    out_file.write('# {}'.format(row))
                    out_file.write('\n')

    out_file.close()                
    return

# Main function ############################################à
def main():
    # check if stations.txt exists and in case download it
    if exists('stations.txt') is False:
        download_stations_list()
    
    # define class with a list of all stations
    stations = stations_list()
    return

if __name__ == '__main__':
    main()

