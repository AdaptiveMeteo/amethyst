"""
This modules contains all the classes that the scene_analysis script has to
access to the input files.
"""

import re
from netCDF4 import Dataset
from os.path import basename, dirname, join, isfile
import numpy as np

# The following is the regular expression that a file name must match to
# be elaborated as a SPS file
SPS_MASK = r'^L2VDP-(?P<type>AHP|ATP|META|REST)-\d{14}Z' \
           r'-(IASI-M01|IASI-M02|IASI-METOPA|CRIS-NPP|IASI-METOPB)-\d{14}Z'

# The following is the regular expression that a file name must match to
# be elaborated as a WRF output file
WRF_MASK = r'^wrf(out|var_output)_d\d\d_\d\d\d\d-\d\d-\d\d_\d\d:\d\d:\d\d'

# The following is the regular expression that a file name must match to
# be elaborated as a MIRTO output file
MIRTO_MASK = r'(?P<type>results|fov|fg).nc'

# Ignore warning when numpy operates with masked nan or Inf
np.seterr(divide='ignore', invalid='ignore', over='ignore')

class Profile(object):
    def __init__(self, name, pressure, temperature, water_vapour,rh):
        self.name = name
        self.pressure = pressure
        self.temperature = temperature
        self.water_vapour = water_vapour
        self.rh           = rh
        
    def filter(self, filter_array):
        self.pressure = self.pressure[filter_array]
        self.temperature = self.temperature[filter_array]
        self.water_vapour = self.water_vapour[filter_array]
        self.rh           = self.rh[filter_array]


class InputData(object):
    """
    An abstract class that wraps one or more file used as input by the scene
    analysis script

    The constructor of this class returns an object of a suitable subclass
    described in this module
    """

    def __new__(cls, input_file_path):
        input_file_name = basename(input_file_path)

        if re.match(SPS_MASK, input_file_name):
            return super(InputData, cls).__new__(SPSInputData)

        elif re.match(WRF_MASK, input_file_name):
            return super(InputData, cls).__new__(WRFInputData)

        elif re.match(MIRTO_MASK, input_file_name):
            return super(InputData, cls).__new__(MIRTOInputData)

        else:
            raise IOError('The input file type has not been recognized because '
                          'its name ({}) does not honour any known schema'

                          .format(input_file_name))


class SPSInputData(InputData):
    def __init__(self, input_file_path):
        input_file_name = basename(input_file_path)
        input_file_dir = dirname(input_file_path)

        file_match = re.match(SPS_MASK, input_file_name)
        file_type = file_match.group('type')

        meta_file_name = input_file_name.replace(file_type, 'META')
        meta_file_path = join(input_file_dir, meta_file_name)

        ahp_file_name = input_file_name.replace(file_type, 'AHP')
        ahp_file_path = join(input_file_dir, ahp_file_name)

        atp_file_name = input_file_name.replace(file_type, 'ATP')
        atp_file_path = join(input_file_dir, atp_file_name)

        if not isfile(meta_file_path):
            raise IOError('Meta file {} not found!'.format(meta_file_path))

        if not isfile(ahp_file_path):
            raise IOError('AHP file {} not found!'.format(ahp_file_path))

        if not isfile(atp_file_path):
            raise IOError('ATP file {} not found!'.format(atp_file_path))


        with Dataset(meta_file_path, 'r') as meta_file:
                self.longitude = np.array(meta_file.variables['retrieval_geographical_position'][:,0])
                self.latitude  = np.array(meta_file.variables['retrieval_geographical_position'][:,1])
                self.quality   = np.array(meta_file.variables['quality_and_completeness'][:])
                pressure = np.array(meta_file.variables['pressure'])
                try:
                        pa_units = meta_file.variables['pressure'].units
                except:
                        raise ValueError('Pressure units not specified!')
                if pa_units == 'Pa':
                        pressure = pressure*100.  
                elif pa_units != 'hPa':
                        raise ValueError('{} is not a valid unit for the pressure'
                        	         .format(pa_units))

        with Dataset(ahp_file_path, 'r') as ahp_file:
            ahp_vars = ahp_file.variables

            try:
                apriori_units = ahp_vars['apriori_specific_humidity'].units
            except:
                raise ValueError('Apriori specific humidity units not '
                                 'specified!')

            try:
                post_units = ahp_vars['posterior_specific_humidity'].units
            except:
                raise ValueError('Posterior specific humidity units not '
                                 'specified!')

            apriori_q = np.array(ahp_vars['apriori_specific_humidity'])
            if apriori_units == 'kg/kg':
                apriori_q /= 1000.
            elif apriori_units != 'g/kg':
                raise ValueError('Bad units for apriori specific humidity')

            post_q = np.array(ahp_vars['posterior_specific_humidity'])
            if post_units == 'kg/kg':
                post_q /= 1000.
            elif post_units != 'g/kg':
                raise ValueError('Bad units for posterior specific humidity')

        with Dataset(atp_file_path, 'r') as atp_file:
            atp_vars = atp_file.variables

            try:
                apriori_units = atp_vars['apriori_air_temperature'].units
            except:
                raise ValueError('Apriori air temperature units not '
                                 'specified!')
            try:
                post_units = atp_vars['posterior_air_temperature'].units
            except:
                raise ValueError('Posterior air temperature units not '
                                 'specified!')

            if apriori_units != 'K':
                raise ValueError('Apriori air temperature in the ATP file is '
                                 'defined using unknown units: {}'
                                 .format(apriori_units))

            if post_units != 'K':
                raise ValueError('Posterior air temperature in the ATP file is '
                                 'defined using unknown units: {}'
                                 .format(apriori_units))

            prior_temp = np.array(atp_vars['apriori_air_temperature'])
            post_temp = np.array(atp_vars['posterior_air_temperature'])

            profile1 = Profile(
                'apriori',
                pressure,
                prior_temp,
                apriori_q,
            )
            profile2 = Profile(
                'posterior',
                pressure,
                post_temp,
                post_q,
            )

            # Build a string that identifies this input
            self.name = ahp_file_name.replace('L2VDP-AHP-', '')

            self.profiles = [profile1, profile2]


class WRFInputData(InputData):
    """
    A wrapper around a file genertated by a WRF model.

    This object is able to return different parameters of the WRF model as
    numpy arrays.

    Args:
        - *file_name*: The file with the data of the WRF model

    :ivar times: An array of datetime64 objects that are the time step of the
                 model
    :ivar lons: A 2D array with the longitude of the points of the model grid
    :ivar lats: A 2D array with the latitude of the points of the model grid
    """

    def __init__(self, input_file_path):
        from atmos.mirto_atmos_tools import mr2rh

        with Dataset(input_file_path, 'r') as f:
            # Read times
            raw_time = f.variables['Times']
            string_times = []
            for t in raw_time:
                t_str = ''.join(k.decode('ASCII') for k in t)
                string_times.append(t_str)
            # Replace the _ with a T to divide date from time
            string_times = [t.replace('_', 'T') + u'Z' for t in
                            string_times]
            self.times = np.array(string_times, dtype=np.datetime64)

            # Read latitudes
            self.latitude = np.array(
                f.variables['XLAT'][:],
                dtype=np.float64,
            )[0,:].flatten()

            # Read longitudes
            self.longitude = np.array(
                f.variables['XLONG'][:],
                dtype=np.float64,
            )[0,:].flatten()

            # Read the water vapour for each time step
            vapour = np.array(
                f.variables['QVAPOR'][:],
                dtype=np.float64,
            )
            # Multiply by 1000 to return values in g/Kg
            vapour *= 1E03

            # move the levels on the last position (so the time is the first)
            vapour = np.rollaxis(vapour, 2, 1)
            vapour = np.rollaxis(vapour, 3, 2)

            self._vapour = vapour

            # Read the pressure levels
            p = f.variables['P'][:]
            pb = f.variables['PB'][:]
            pressure = np.array(p + pb, dtype=np.float64)
            pressure = np.rollaxis(pressure, 2, 1)
            pressure = np.rollaxis(pressure, 3, 2)
            self._pressure = pressure

            # Read temperature
            temperature = np.array(
                f.variables['T'][:],
                dtype=np.float64,
            )
            temperature += 300

            temperature = np.rollaxis(temperature, 2, 1)
            temperature = np.rollaxis(temperature, 3, 2)
            temperature *= np.power(self._pressure * 1E-5, 2 / 7.)
            self._temperature = temperature

        self.name = basename(input_file_path)

        # Prepare profiles
        self.profiles = []
        n_of_obs = len(self.latitude)
        n_of_levels = self._pressure.shape[-1]
        self._rh = np.empty( (2,n_of_obs, n_of_levels ) )
        for i, t in enumerate(self.times):
            self._rh[i,:]= mr2rh(self._pressure[i,:].reshape(n_of_obs, n_of_levels)*1e-02,
                                 self._temperature[i,:].reshape(n_of_obs, n_of_levels),
                                 self._vapour[i,:].reshape(n_of_obs, n_of_levels)
                                 )[0].transpose() 

            p = Profile(
                str(t),
                self._pressure[i, :].reshape(n_of_obs, n_of_levels),
                self._temperature[i, :].reshape(n_of_obs, n_of_levels),
                self._vapour[i, :].reshape(n_of_obs, n_of_levels),
                self._rh[i, :]
            )
            self.profiles.append(p)


class MIRTOInputData(InputData):

    def __init__(self, input_file_path):
        from atmos.mirto_atmos_tools import mr2rh

        FIX_TEMP = 273.15

        input_file_name = basename(input_file_path)
        input_file_dir  = dirname(input_file_path)
        
        file_match = re.match(MIRTO_MASK, input_file_name)
        file_type  = file_match.group('type')

        results_file_name = input_file_name.replace(file_type, 'results')
        results_file_path = join(input_file_dir, results_file_name)

        fov_file_name = input_file_name.replace(file_type, 'fov')
        fov_file_path = join(input_file_dir, fov_file_name)

        fg_file_name = input_file_name.replace(file_type, 'fg')
        fg_file_path = join(input_file_dir, fg_file_name)

        print(f"results_file_path: {results_file_path}")
        if not isfile(results_file_path):
            raise IOError('results.nc file {} not found!'.format(results_file_path))

        if not isfile(fov_file_path):
            raise IOError('fov.nc file {} not found!'.format(fov_file_path))

        if not isfile(fg_file_path):
            raise IOError('fg.nc file {} not found!'.format(fg_file_path))

        with Dataset(results_file_path, 'r') as results_file:
            results_vars = results_file.variables

            pressure         = np.ma.masked_invalid(results_vars['p'][:])
            temperature      = np.ma.masked_invalid(results_vars['t'][:])
            water_vapour     = np.ma.masked_invalid(results_vars['q'][:])
            rh               = np.ma.masked_invalid( mr2rh(pressure,
                                                           temperature,
                                                           water_vapour,
                                                           FIX_TEMP)[0] )

            self.d2   =  np.array(results_vars['d2'][:],dtype=np.float64 )
            
        with Dataset(fov_file_path, 'r') as fov_file:
                self.longitude = np.array(fov_file.variables['Longitude'])
                self.latitude = np.array(fov_file.variables['Latitude'])


        with Dataset(fg_file_path, 'r') as fg_file:
            fg_vars = fg_file.variables
            fg_atmos_components = fg_file.groups['atmospheric_components']

            try:
                p_units = fg_atmos_components['p'].units
            except:
                raise ValueError('Pressure units not '
                                 'specified!')
            try:
                T_units = fg_atmos_components['T'].units
            except:
                raise ValueError('Temperature units not '
                                 'specified!')
            try:
                q_units = fg_atmos_components['q'].units
            except:
                raise ValueError('Specific humidity units not '
                                 'specified!')
                
            fg_vars = fg_file.groups['atmospheric_components'].variables
            prior_temperature  = np.ma.masked_invalid(fg_vars['T'][:])
            prior_pressure     = np.ma.masked_invalid(fg_vars['p'][:])
            prior_water_vapour = np.ma.masked_where(fg_vars['q'][:]<0,fg_vars['q'][:])
            prior_rh           = np.ma.masked_invalid( mr2rh(prior_pressure,
                                                             prior_temperature,
                                                             prior_water_vapour,
                                                             FIX_TEMP)[0]  )
        # Prepare profiles
        profile1 = Profile(
                'MIRTO prior',
                prior_pressure,
                prior_temperature,
                prior_water_vapour,
                prior_rh
            )
        profile2 = Profile(
                'MIRTO posterior',
                pressure,
                temperature,
                water_vapour,
                rh
            )
        
        self.profiles = [profile1,profile2]
        
        self.units = [p_units,T_units,q_units]
        
        # Build a string that identifies this input

        print(f"results_file_path: {results_file_path}")
        if results_file_path.split('/')[-2] == 'mirto':
            name = results_file_path.split('/')[-3]
        else:
            name = results_file_path.split('/')[-2]
        if 'mirto' not in name:
                self.name = 'MIRTO_' + name
        else:
                self.name = name
        
        
