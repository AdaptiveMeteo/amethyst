# Project: Simplified Amethyst Validation
# Copyright AdaptiveMeteo
# Author: Paolo Antonelli, with aid of ChatGPT
# Date: October 15, 2024
# License: GPL-3.0

# This project is dedicated to the validation of retrievals using rawinsonde data.
# It aims to provide a clean and consistent framework for comparing retrieval data
# with observational data from rawinsondes, focusing on simplifying the process
# and making it more robust.

import netCDF4 as nc
import numpy as np
import logging
from datetime import datetime

LOG = logging.getLogger(__name__)

class AmethystRetrievals:
    def __init__(self, path, files):
        """
        Initializes the AmethystRetrievals class with configurable variables.

        Parameters:
        path (str): Path to the directory containing the files.
        files (dict): Dictionary containing the file names.
        """
        self.path = path
        self.files = files
        self.file_paths = {name: f"{self.path}/{filename}" for name, filename in files.items()}
        LOG.info("Initialized AmethystRetrievals with path: %s and files: %s", path, files)

    def read_amethyst_output(self):
        """
        Reads in the NetCDF file 'amethyst_output.nc' and extracts its variables.

        Returns:
        dict: A dictionary containing the extracted variables.
        """
        file_path = self.file_paths.get('amethyst_output')
        if not file_path:
            LOG.error("File 'amethyst_output.nc' not provided.")
            return None
        
        try:
            with nc.Dataset(file_path, 'r') as dataset:
                data = {
                    'obs': dataset.variables['obs'][:],
                    'p': dataset.variables['p'][:],
                    't': dataset.variables['t'][:],
                    'q': dataset.variables['q'][:],
                    'o3': dataset.variables['o3'][:],
                    'skt': dataset.variables['skt'][:],
                    'd2': dataset.variables['d2'][:],
                    'ems_coeff': dataset.variables['ems_coeff'][:],
                    'DA_R': dataset.variables['DA_R'][:],
                    'DA_Yret': dataset.variables['DA_Yret'][:],
                    'DA_Hret': dataset.variables['DA_Hret'][:],
                    'DA_Lambda': dataset.variables['DA_Lambda'][:]
                }
                LOG.info("Successfully read the NetCDF file: %s", file_path)
                return data
        except FileNotFoundError:
            LOG.error("File not found: %s", file_path)
        except Exception as e:
            LOG.error("An error occurred while reading the file: %s", str(e))

        return None

    def read_fov_file(self):
        """
        Reads in the NetCDF file 'fov.nc' and extracts its variables.

        Returns:
        dict: A dictionary containing the extracted variables.
        """
        file_path = self.file_paths.get('fov')
        if not file_path:
            LOG.error("File 'fov.nc' not provided.")
            return None
        
        try:
            with nc.Dataset(file_path, 'r') as dataset:
                data = {
                    'Latitude': dataset.variables['Latitude'][:],
                    'Longitude': dataset.variables['Longitude'][:],
                    'Time': dataset.variables['Time'][:],
                    'Radiance': dataset.variables['Radiance'][:],
                    'Wavenumbers': dataset.variables['Wavenumbers'][:],
                    'FOV_angle': dataset.variables['FOV_angle'][:],
                    'Satellite_zenith_angle': dataset.variables['Satellite_zenith_angle'][:],
                    'Solar_azimuth_angle': dataset.variables['Solar_azimuth_angle'][:],
                    'Solar_zenith_angle': dataset.variables['Solar_zenith_angle'][:],
                    'Satellite_azimuth_angle': dataset.variables['Satellite_azimuth_angle'][:]
                }
                LOG.info("Successfully read the NetCDF file: %s", file_path)
                return data
        except FileNotFoundError:
            LOG.error("File not found: %s", file_path)
        except Exception as e:
            LOG.error("An error occurred while reading the file: %s", str(e))

        return None

    def read_fg_file(self):
        """
        Reads in the NetCDF file 'fg.nc' (First Guess file) and extracts its variables.

        Returns:
        dict: A dictionary containing the extracted variables.
        """
        file_path = self.file_paths.get('fg')
        if not file_path:
            LOG.error("File 'fg.nc' not provided.")
            return None
        
        try:
            with nc.Dataset(file_path, 'r') as dataset:
                data = {
                    'Latitude': dataset.variables['Latitude'][:],
                    'Longitude': dataset.variables['Longitude'][:],
                    'Time': dataset.variables['Time'][:],
                    'atmospheric_components': {
                        'p': dataset['atmospheric_components'].variables['p'][:],
                        'T': dataset['atmospheric_components'].variables['T'][:],
                        'q': dataset['atmospheric_components'].variables['q'][:],
                        'O3': dataset['atmospheric_components'].variables['O3'][:],
                        'skT': dataset['atmospheric_components'].variables['skT'][:],
                        'sp': dataset['atmospheric_components'].variables['sp'][:]
                    },
                    'surface_components': {
                        'ModelWaveNumbers': dataset['surface_components'].variables['ModelWaveNumbers'][:],
                        'ModelFunctions': dataset['surface_components'].variables['ModelFunctions'][:],
                        'ModelFunctionsBias': dataset['surface_components'].variables['ModelFunctionsBias'][:],
                        'ModelCovariance': dataset['surface_components'].variables['ModelCovariance'][:],
                        'LandOrWater': dataset['surface_components'].variables['LandOrWater'][:]
                    }
                }
                LOG.info("Successfully read the NetCDF file: %s", file_path)
                return data
        except FileNotFoundError:
            LOG.error("File not found: %s", file_path)
        except Exception as e:
            LOG.error("An error occurred while reading the file: %s", str(e))

        return None

    def read_apriori_file(self):
        """
        Reads in the NetCDF file 'apriori.nc' and extracts its variables.

        Returns:
        dict: A dictionary containing the extracted variables.
        """
        file_path = self.file_paths.get('apriori')
        if not file_path:
            LOG.error("File 'apriori.nc' not provided.")
            return None
        
        try:
            with nc.Dataset(file_path, 'r') as dataset:
                data = {
                    'Latitude': dataset.variables['Latitude'][:],
                    'Longitude': dataset.variables['Longitude'][:],
                    'atmospheric_components': {
                        'Covariances': {
                            'T': dataset['atmospheric_components/Covariances'].variables['T'][:],
                            'q': dataset['atmospheric_components/Covariances'].variables['q'][:],
                            'O3': dataset['atmospheric_components/Covariances'].variables['O3'][:],
                            'T_q': dataset['atmospheric_components/Covariances'].variables['T_q'][:]
                        }
                    }
                }
                LOG.info("Successfully read the NetCDF file: %s", file_path)
                return data
        except FileNotFoundError:
            LOG.error("File not found: %s", file_path)
        except Exception as e:
            LOG.error("An error occurred while reading the file: %s", str(e))

        return None

# Example usage for testing
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    path = "/mnt/miwa/amethyst/amethyst_test_data/input_data/cris/20220806_233400/amethyst"
    files = {
        'amethyst_output': "amethyst_output.nc",
        'fov': "fov.nc",
        'fg': "fg.nc",
        'apriori': "apriori.nc"
    }
    amethyst = AmethystRetrievals(path, files)
    
    amethyst_output_data = amethyst.read_amethyst_output()
    if amethyst_output_data:
        LOG.info("amethyst_output: Number of observations: %d", len(amethyst_output_data['obs']))
        LOG.info("amethyst_output: Mean temperature: %.2f K", np.mean(amethyst_output_data['t']))
        num_d2_less_than_1 = np.sum(amethyst_output_data['d2'] < 1)
        LOG.info("amethyst_output: Number of d2 values less than 1: %d", num_d2_less_than_1)

    fov_data = amethyst.read_fov_file()
    if fov_data:
        LOG.info("fov: Number of FOVs: %d", len(fov_data['Latitude']))
        LOG.info("fov: Mean Radiance value: %.2f", np.mean(fov_data['Radiance']))

    fg_data = amethyst.read_fg_file()
    if fg_data:
        LOG.info("fg: Number of FOVs: %d", len(fg_data['Latitude']))
        LOG.info("fg: Mean Surface Pressure: %.2f hPa", np.mean(fg_data['atmospheric_components']['sp']))

    apriori_data = amethyst.read_apriori_file()
    if apriori_data:
        LOG.info("apriori: Number of FOVs: %d", len(apriori_data['Latitude']))
        LOG.info("apriori: Mean Temperature Covariance (first FOV): %.2f", np.mean(apriori_data['atmospheric_components']['Covariances']['T'][0]))
