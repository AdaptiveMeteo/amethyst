from datetime import datetime
from download_wyoming_sondes import download_sondes
from read_stations_list import stations
import sys

arg = "16245" if len(sys.argv) == 1 else sys.argv[1]
date = datetime(2018,2,11,10)
#date = datetime.now()
station = stations.get(arg)
sondes = download_sondes(station,date)

sonde_coor = (sondes.lon,sondes.lat)

from tools.earth_geometry import point_in_polygon

domain = [(-52.345181, 57.992209),
	  (-52.345181, 74.70592),
	  (0.185736, 74.70592),
	  (0.185736, 57.992209),
	  (-52.345181, 57.992209)]

if point_in_polygon(sonde_coor,domain) :
	print('The station '+station.name+' is inside the domain:\n')
else:
	print('The station '+station.name+' is not inside the domain:\n')
for x in domain:
	print(x,sep=',')
print('\nRowinsonde coordinates (lon,lat): ',sonde_coor)
