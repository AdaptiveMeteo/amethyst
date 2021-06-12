from dobjects.ObservationErrors.Iasi import IasiObservationError
from dobjects.ObservationErrors.Cris import CrisObservationError
from dobjects.ObservationErrors.Iasi import IasiObservationErrorSps
from dobjects.ObservationErrors.Cris import CrisObservationErrorSps

class UnknownObservationType(Exception):
    pass


#PaoloA 25052021
def create_obs_err(data_filename, obs_err_type='iasi'):
#def create_obs_err(data_filename, obs_err_type='cris'):
    if obs_err_type == 'cris':
        return CrisObservationError(data_filename)
    elif obs_err_type == 'iasi':
        return IasiObservationError(data_filename)
    elif obs_err_type == 'iasi_d':
        #PaoloA 12082020
        #return CrisObservationError(data_filename)
        return IasiObservationError(data_filename)
    #elif obs_err_type == 'iasi_sps':
        #return IasiObservationErrorSps(data_filename)
    #elif obs_err_type == 'cris_sps':
        #return CrisObservationErrorSps(data_filename)
    else:
        raise UnknownObservationType

#PaoloA 25052021
def create_obs_err_sps(data_filename, obs_err_type='iasi_d'):
#def create_obs_err_sps(data_filename, obs_err_type='cris'):
    if obs_err_type == 'cris':
        return CrisObservationErrorSps(data_filename)
    elif obs_err_type == 'iasi':
        return IasiObservationErrorSps(data_filename)
    elif obs_err_type == 'iasi_d':
        #PaoloA 12082020
        #return CrisObservationErrorSps(data_filename)
        return IasiObservationErrorSps(data_filename)
    else:
        raise UnknownObservationType

