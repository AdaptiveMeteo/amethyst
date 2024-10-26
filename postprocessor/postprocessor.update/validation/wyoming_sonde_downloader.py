# Project: Wyoming Sonde Downloader
# Copyright AdaptiveMeteo
# Author: Paolo Antonelli
# Date: October 15, 2024
# License: GPL-3.0

import requests
import os
import logging
from datetime import datetime
import re
import numpy as np
from geometry.earth_geometry import haversine

LOG = logging.getLogger(__name__)

class StationClass:
    def __init__(self, latitude, longitude, SYNOP, STATION):
        self.latitude = latitude
        self.longitude = longitude
        self.SYNOP = SYNOP
        self.STATION = STATION


class WyomingSondeDownloader:
    def __init__(self, output_dir, missing_stations_file=None):
        """
        Initializes the WyomingSondeDownloader class with configurable variables.

        Parameters:
        output_dir (str): Directory where downloaded data will be saved.
        missing_stations_file (str): Path to the file containing missing stations information.
        """
        self.base_url = "http://weather.uwyo.edu/cgi-bin/sounding"
        self.output_dir = output_dir
        self.missing_stations_file = missing_stations_file
        self.missing_stations = []
        self.CD = []
        self.STATION = []
        self.ICAO = []
        self.IATA = []
        self.SYNOP = []
        self.LAT = []
        self.LON = []
        self.ELEV = []
        self.COUNTRY = []
        self.stations = []

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        LOG.info("Initialized WyomingSondeDownloader with base_url: %s and output_dir: %s", self.base_url, output_dir)

        # Read missing stations from missing_stations.txt if provided
        if self.missing_stations_file:
            self.read_missing_stations()

    def read_missing_stations(self, *string):
        """
        Reads the list of missing stations from a file and populates class variables.
        """
        LOG.info("Creating dictionary for missing stations...")

        try:
            missing_stations_file = open(self.missing_stations_file, 'r')
            rows = [x for x in missing_stations_file.read().split('\n') if x]
            missing_stations_file.close()
        except Exception as e:
            LOG.error("An error occurred while reading the missing stations list from %s: %s", self.missing_stations_file, str(e))
            return None

        station = None
        for row in rows:
            if len(string) > 0 and row[33:38] == string[0]:
                station = self.create_station_class(row)
                if row[46] == 'S':
                    station.lat *= -1
                if row[55] == 'W':
                    station.lon *= -1

            self.CD.append(row[:2])
            self.STATION.append(row[3:19])
            self.ICAO.append(row[20:24])
            self.IATA.append(row[26:29])
            self.SYNOP.append(row[33:38])
            if row[46] == 'N':
                self.LAT.append(float(".".join(row[39:46].split(' '))))
            else:
                self.LAT.append(-1 * float(".".join(row[39:46].split(' '))))
            if row[55] == 'E':
                self.LON.append(float(".".join(row[48:55].split(' '))))
            else:
                self.LON.append(-1 * float(".".join(row[48:55].split(' '))))
            try:
                self.ELEV.append(int(row[57:61]))
            except ValueError:
                self.ELEV.append('--')
            self.COUNTRY.append(row[-2:])

        self.missing_stations = []
        return station

    def create_station_class(self, row):
        """
        Creates a station class object from a row of data.
        """
        station = Station(
            row[:2],
            row[3:19],
            row[20:24],
            row[26:29],
            row[33:38],
            float(".".join(row[39:46].split(' '))),
            float(".".join(row[48:55].split(' '))),
            row[57:61],
            row[-2:]
        )
        return station

    def download_sonde(self, station, date):
        """
        Downloads a rawinsonde file for a given station and date.

        Parameters:
        station (str): Station identifier.
        date (datetime): Date of the rawinsonde observation.

        Returns:
        str: Path to the downloaded file or None if download failed.
        """
        if station in self.missing_stations or station in [s.SYNOP for s in self.missing_stations]:
            LOG.warning("Station %s is marked as missing. Skipping download.", station)
            return None

        date_str = date.strftime("%Y%m%d")
        hour_str = date.strftime("%H")
        url = (f"{self.base_url}?region=naconf&TYPE=TEXT%3ALIST&YEAR={date.year}&MONTH={date.month:02d}&"
               f"FROM={date.day:02d}{hour_str}&TO={date.day:02d}{hour_str}&STNM={station}")
        output_file = os.path.join(self.output_dir, f"{date_str}_{hour_str}_{station}.html")
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            with open(output_file, 'w') as file:
                file.write(response.text)
            LOG.info("Successfully downloaded file: %s", output_file)
            return output_file
        except requests.RequestException as e:
            LOG.error("Failed to download file from %s: %s", url, str(e))
        except Exception as e:
            LOG.error("An error occurred while saving the file: %s", str(e))
        
        return None

    def parse_sonde_file(self, file_path):
        """
        Parses a downloaded rawinsonde file to extract relevant data.

        Parameters:
        file_path (str): Path to the downloaded sonde file.

        Returns:
        dict: A dictionary containing parsed data or None if parsing failed.
        """
        try:
            with open(file_path, 'r') as file:
                content = file.read()
                data = {
                    'latitude': None,
                    'longitude': None,
                    'pressure': [],
                    'altitude': [],
                    'temperature': [],
                    'dewpoint': [],
                    'relative_humidity': [],
                    'wind_direction': [],
                    'wind_speed': []
                }

                # Extract latitude and longitude using regex
                lat_match = re.search(r'Station latitude:\s*([\d.]+)', content)
                lon_match = re.search(r'Station longitude:\s*(-?[\d.]+)', content)
                if lat_match:
                    data['latitude'] = float(lat_match.group(1))
                if lon_match:
                    data['longitude'] = float(lon_match.group(1))

                # Find the data section using regex
                data_section_match = re.search(r'PRES.*?\n-+\n(.*?)\n\s*</PRE>', content, re.DOTALL)
                if not data_section_match:
                    LOG.error("Could not find the data section in the file.")
                    return None

                data_section = data_section_match.group(1)
                lines = data_section.splitlines()

                # Extract atmospheric data from each line
                for line in lines:
                    parts = line.split()
                    if len(parts) >= 11:
                        data['pressure'].append(float(parts[0]))
                        data['altitude'].append(float(parts[1]))
                        data['temperature'].append(float(parts[2]))
                        data['dewpoint'].append(float(parts[3]))
                        data['relative_humidity'].append(float(parts[4]))
                        data['wind_direction'].append(float(parts[6]))
                        data['wind_speed'].append(float(parts[7]))

                if data['latitude'] is not None and data['longitude'] is not None:
                    LOG.info("Parsed data - Latitude: %.2f, Longitude: %.2f", data['latitude'], data['longitude'])
                if data['pressure']:
                    LOG.info("Parsed data - Number of pressure levels: %d", len(data['pressure']))
                    LOG.info("Parsed data - Pressure range: %.2f hPa to %.2f hPa", min(data['pressure']), max(data['pressure']))
                    LOG.info("Parsed data - Altitude range: %.2f m to %.2f m", min(data['altitude']), max(data['altitude']))
                    LOG.info("Parsed data - Temperature range: %.2f C to %.2f C", min(data['temperature']), max(data['temperature']))
                    LOG.info("Parsed data - Dewpoint range: %.2f C to %.2f C", min(data['dewpoint']), max(data['dewpoint']))
                    LOG.info("Parsed data - Wind speed range: %.2f knot to %.2f knot", min(data['wind_speed']), max(data['wind_speed']))
                else:
                    LOG.warning("No pressure data found in the file.")

                return data
        except Exception as e:
            LOG.error("An error occurred while parsing the file: %s", str(e))
        
        return None

    def read_stations_list(self):

    #Reads the list of stations from the stations.txt file and stores them in self.stations.

        try:
            with open('stations.txt', 'r') as station_file:
                rows = station_file.readlines()
                for row in rows:
                    if len(row) >= 61:  # Ensure the row contains enough information
                        # Extract latitude and longitude strings
                        lat_str = row[39:46].strip()
                        lon_str = row[48:55].strip()
    
                        # Remove the cardinal directions and convert to float
                        lat = float(lat_str[:-1].replace(' ', '.'))
                        lon = float(lon_str[:-1].replace(' ', '.'))
    
                        # Apply negative sign if necessary
                        if row[46] == 'S':
                            lat *= -1
                        if row[55] == 'W':
                            lon *= -1
    
                        # Create station object
                        station = StationClass(
                            latitude=lat,
                            longitude=lon,
                            SYNOP=row[33:38].strip(),
                            STATION=row[3:19].strip()
                        )
                        self.stations.append(station)
            LOG.info("Successfully loaded %d stations from stations.txt", len(self.stations))
        except FileNotFoundError:
            LOG.error("Stations file 'stations.txt' not found. Please provide the correct file.")
        except ValueError as e:
            LOG.error("Error parsing stations.txt: %s", e)

    # Get list of stations within the bounding box
    def get_stations_within_area(self, min_lat, max_lat, min_lon, max_lon):
        """
        Returns a list of stations within the specified bounding box.

        Parameters:
        min_lat (float): Minimum latitude of the area.
        max_lat (float): Maximum latitude of the area.
        min_lon (float): Minimum longitude of the area.
        max_lon (float): Maximum longitude of the area.

        Returns:
        list: List of station objects that are within the bounding box.
        """
        stations_in_area = [
            station for station in self.stations
            if min_lat <= station.latitude <= max_lat and min_lon <= station.longitude <= max_lon
        ]
        return stations_in_area

    def get_nearest_station(self, lat, lon, R=75):
        """
        Return nearest station to given point

        Parameters
        ----------
            lat : float
            lon : float
            R:    float for the maximum radius for the search area

        Returns
            station: station_class
            dist   : float, distance between FOV and station
        """
        # Read missing stations if not already done
        if not self.missing_stations:
            self.read_missing_stations()

        dists = np.array([haversine(stat_lon, self.LAT[i], lon, lat) for i, stat_lon in enumerate(self.LON)])
        dists = np.ma.masked_where(dists > R, dists)

        if dists.mask.all():
            LOG.debug("No station found in an area within {:.2f} km radius".format(R))
            return None, np.nan

        index = np.argmin(dists)
        return self.create_station_class_from_index(index), dists[index]

    def create_station_class_from_index(self, index):
        """
        Creates a station class object from an index of the class variables.
        """
        station = Station(
            self.CD[index],
            self.STATION[index],
            self.ICAO[index],
            self.IATA[index],
            self.SYNOP[index],
            self.LAT[index],
            self.LON[index],
            self.ELEV[index],
            self.COUNTRY[index]
        )
        return station

    def read_missing_stations_coord(self):
        """
        NEVER CALL THIS FUNCTION IF NOT NECESSARY
        Reads missing stations coordinates and updates the missing stations file.
        """
        from download_wyoming_sondes import download_sondes

        ref_date = datetime.strptime('2020-11-02T12:00', '%Y-%m-%dT%H:%M')

        with open('missing_stations.txt', 'r') as ref_file:
            ref_rows = [x for x in ref_file.read().split('\n') if x]

        with open('stations.txt', 'r') as station_file:
            rows = [x for x in station_file.read().split('\n') if len(x) > 71 and x[71] == 'X']
            stations = self.read_stations_list(rows)

        out_rows = []
        with open('missing_stations_2.txt', 'w') as out_file:
            for i, row in enumerate(ref_rows):
                try:
                    raob = download_sondes(stations.get(row[32:37]), ref_date)
                    lat = raob.lat
                    lon = raob.lon

                    lat_str = repr(abs(lat)) + 'N'
                    lon_str = repr(abs(lon)) + 'E'

                    if lon < 0:
                        lon_str = lon_str[:-1] + 'W'
                    if lat < 0:
                        lat_str = lat_str[:-1] + 'S'

                    lon_str = lon_str.replace('.', ' ')
                    lat_str = lat_str.replace('.', ' ')

                    tmp_str = row[:39] + lat_str + '  ' + lon_str + row[54:]
                    if len(tmp_str) < 82:
                        tmp_str = row[:9] + ' ' + row[9:39] + lat_str + '  ' + lon_str + row[54:]

                    out_rows.append(tmp_str)
                    if len(ref_rows[i]) < 82:
                        out_rows[i] = out_rows[-1][:9] + ' ' + out_rows[-1][9:]

                    out_file.write(out_rows[-1])
                    out_file.write('\n')

                except Exception as e:
                    LOG.info("Skipped station %s!", row[32:37])
                    out_file.write('# {}'.format(row))
                    out_file.write('\n')

# Example usage for testing
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    output_dir = "./sondes"
    missing_stations_file = "/home/amethyst/amethyst/postprocessor/validation/missing_stations.txt"
    downloader = WyomingSondeDownloader(output_dir, missing_stations_file)
    
    # Example station and date for testing
    station = "91285"
    date = datetime(2024, 10, 14, 0)  # October 14, 2024, 00:00 UTC
    downloaded_file = downloader.download_sonde(station, date)
    
    # Parse the file if it was successfully downloaded
    if downloaded_file:
        parsed_data = downloader.parse_sonde_file(downloaded_file)
        
        # Print parsed data in more detail
        if parsed_data:
            LOG.info("Detailed Parsed Data:")
            LOG.info("Latitude: %s", parsed_data['latitude'])
            LOG.info("Longitude: %s", parsed_data['longitude'])
            LOG.info("Pressure Levels: %s", parsed_data['pressure'])
            LOG.info("Altitude Levels: %s", parsed_data['altitude'])
            LOG.info("Temperature Levels: %s", parsed_data['temperature'])
            LOG.info("Dewpoint Levels: %s", parsed_data['dewpoint'])
            LOG.info("Relative Humidity Levels: %s", parsed_data['relative_humidity'])
            LOG.info("Wind Direction: %s", parsed_data['wind_direction'])
            LOG.info("Wind Speed: %s", parsed_data['wind_speed'])

    # Example usage of reading stations list
    stations_list_file = "/home/amethyst/amethyst/postprocessor/validation/stations.txt"
    stations = downloader.read_stations_list(stations_list_file)
    LOG.info("Stations: %s", stations)
