# Project: Validate Amethyst Retrievals
# Copyright AdaptiveMeteo
# Author: Paolo Antonelli
# Date: October 15, 2024
# License: GPL-3.0

import logging
from datetime import datetime, timedelta
from wyoming_sonde_downloader import WyomingSondeDownloader
from amethyst_read_retrievals import AmethystRetrievals
import numpy as np
from scipy.interpolate import interp1d
import os

LOG = logging.getLogger(__name__)

class ValidateAmethystRetrievals:
    def __init__(self, output_dir, retrievals_file, missing_stations_file=None):
        """
        Initializes the ValidateAmethystRetrievals class.

        Parameters:
        output_dir (str): Directory where downloaded data will be saved.
        retrievals_file (str): Path to the retrievals NetCDF file.
        missing_stations_file (str): Path to the file containing missing stations information.
        """
        self.output_dir = output_dir
        self.retrievals_file = retrievals_file
        self.missing_stations_file = missing_stations_file
        
        # Initialize the WyomingSondeDownloader and AmethystRetrievals classes
        self.downloader = WyomingSondeDownloader(output_dir, missing_stations_file)
        self.downloader.read_stations_list()

        # Extract the path of the retrievals file
        retrievals_path = os.path.dirname(retrievals_file)

        self.retrievals = AmethystRetrievals(retrievals_path, files={
            'amethyst_output': 'amethyst_output.nc',
            'fov': 'fov.nc',
            'fg': 'fg.nc',
            'apriori': 'apriori.nc'
        })

        LOG.info("Initialized ValidateAmethystRetrievals with output_dir: %s and retrievals_file: %s", output_dir, retrievals_file)

    def download_and_validate(self):
        """
        Downloads the corresponding rawinsonde data for the area covered by the retrievals
        and validates the retrievals against the sonde measurements.
        """
        # Define the bounding box for the retrievals
        fov_data = self.retrievals.read_fov_file()
        latitudes = fov_data['Latitude']
        longitudes = fov_data['Longitude']
        min_lat, max_lat = min(latitudes), max(latitudes)
        min_lon, max_lon = min(longitudes), max(longitudes)

        LOG.info("Defined bounding box - Min Lat: %.2f, Max Lat: %.2f, Min Lon: %.2f, Max Lon: %.2f", min_lat, max_lat, min_lon, max_lon)

        # Get list of stations within the bounding box
        stations = self.downloader.get_stations_within_area(min_lat, max_lat, min_lon, max_lon)
        if not stations:
            LOG.warning("No stations found within the defined area.")
            return

        LOG.info("Found %d stations within the area.", len(stations))

        # Download sonde data for each station
        downloaded_sonde_data = {}
        for station in stations:
            LOG.info("Downloading sonde data for station: %s", station.STATION)
            sonde_file = self.downloader.download_sonde(station.SYNOP, datetime.utcfromtimestamp(fov_data['Time'][0] / 1000).strftime('%Y%m%d%H'))
            if sonde_file is None:
                LOG.warning("Failed to download sonde data for station: %s", station.SYNOP)
                continue

            if not station.SYNOP:
                LOG.error("Station SYNOP ID is missing, unable to form the correct URL.")
                continue
                continue

            sonde_data = self.downloader.parse_sonde_file(sonde_file)
            if sonde_data is None:
                LOG.warning("Failed to parse sonde data for station: %s", station.SYNOP)
                continue

            downloaded_sonde_data[station.SYNOP] = sonde_data

        # Loop through each FOV and validate with the nearest available station
        for i in range(len(latitudes)):
            fov = {
                'latitude': latitudes[i],
                'longitude': longitudes[i],
                'pressure': fov_data['Pressure'][i, :],
                'temperature': fov_data['Temperature'][i, :],
                'relative_humidity': fov_data['Relative_Humidity'][i, :],
                'water_vapour_mixing_ratio': fov_data['Water_Vapour_Mixing_Ratio'][i, :]
            }
            lat, lon = fov['latitude'], fov['longitude']
            LOG.info("Processing FOV at latitude: %.2f, longitude: %.2f", lat, lon)

            # Find the nearest station with downloaded sonde data
            station, dist = self.downloader.get_nearest_station(lat, lon)
            if station is None or station.SYNOP not in downloaded_sonde_data:
                LOG.warning("No nearby station with sonde data found for FOV at latitude: %.2f, longitude: %.2f", lat, lon)
                continue

            LOG.info("Nearest station: %s, Distance: %.2f km", station.STATION, dist)

            # Perform validation
            sonde_data = downloaded_sonde_data[station.SYNOP]
            self.validate_retrieval_with_sonde(fov, sonde_data)

    def validate_retrieval_with_sonde(self, fov, sonde_data):
        """
        Validates the retrieval data against sonde data.

        Parameters:
        fov (dict): Field of view data from the retrievals file.
        sonde_data (dict): Parsed sonde data.
        """
        LOG.info("Validating retrieval for FOV at latitude: %.2f, longitude: %.2f", fov['latitude'], fov['longitude'])

        retrieval_pressure = fov['pressure']
        # Validate Temperature
        retrieval_temp = fov['temperature']
        sonde_pressure = sonde_data['pressure']
        sonde_temp = sonde_data['temperature']

        # Interpolate sonde temperature to retrieval pressure levels
        interpolation_function = interp1d(sonde_pressure, sonde_temp, bounds_error=False, fill_value='extrapolate')
        interpolated_sonde_temp = interpolation_function(retrieval_pressure)

        # Calculate temperature differences
        temp_diff = retrieval_temp - interpolated_sonde_temp
        mean_temp_diff = np.mean(temp_diff)
        LOG.info("Mean temperature difference: %.2f K", mean_temp_diff)

        # Validate Relative Humidity
        retrieval_rh = fov['relative_humidity']
        sonde_rh = sonde_data['relative_humidity']

        # Interpolate sonde relative humidity to retrieval pressure levels
        interpolation_function_rh = interp1d(sonde_pressure, sonde_rh, bounds_error=False, fill_value='extrapolate')
        interpolated_sonde_rh = interpolation_function_rh(retrieval_pressure)

        # Calculate relative humidity differences
        rh_diff = retrieval_rh - interpolated_sonde_rh
        mean_rh_diff = np.mean(rh_diff)
        LOG.info("Mean relative humidity difference: %.2f %%", mean_rh_diff)

        # Validate Water Vapour Mixing Ratio
        retrieval_wvmr = fov['water_vapour_mixing_ratio']
        sonde_wvmr = sonde_data['water_vapour_mixing_ratio']

        # Interpolate sonde water vapour mixing ratio to retrieval pressure levels
        interpolation_function_wvmr = interp1d(sonde_pressure, sonde_wvmr, bounds_error=False, fill_value='extrapolate')
        interpolated_sonde_wvmr = interpolation_function_wvmr(retrieval_pressure)

        # Calculate water vapour mixing ratio differences
        wvmr_diff = retrieval_wvmr - interpolated_sonde_wvmr
        mean_wvmr_diff = np.mean(wvmr_diff)
        LOG.info("Mean water vapour mixing ratio difference: %.2f g/kg", mean_wvmr_diff)

# Example usage for testing
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    output_dir = "/home/amethyst/amethyst/postprocessor/validation/validation_sondes"
    retrievals_file = "/mnt/miwa/amethyst/amethyst_test_data/input_data/cris/20220806_233400/amethyst/amethyst_output.nc"
    missing_stations_file = "/home/amethyst/amethyst/postprocessor/validation/missing_stations.txt"
    validator = ValidateAmethystRetrievals(output_dir, retrievals_file, missing_stations_file)
    
    validator.download_and_validate()


