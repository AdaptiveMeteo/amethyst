from dobjects.ObservationErrors.Iasi import IasiObservationError
from dobjects.ObservationErrors.Cris import CrisObservationError
from dobjects.ObservationErrors.Iasi import IasiObservationErrorSps
from dobjects.ObservationErrors.Cris import CrisObservationErrorSps

class UnknownObservationType(Exception):
    pass


def create_obs_err(data_filename, obs_err_type='cris'):
    if obs_err_type == 'cris':
        return CrisObservationError(data_filename)
    elif obs_err_type == 'iasi':
        return IasiObservationError(data_filename)
    elif obs_err_type == 'iasi_sps':
        return IasiObservationErrorSps(data_filename)
    elif obs_err_type == 'cris_sps':
        return CrisObservationErrorSps(data_filename)
    else:
        raise UnknownObservationType

def create_obs_err_sps(data_filename, obs_err_type='cris'):
    if obs_err_type == 'cris':
        return CrisObservationErrorSps(data_filename)
    elif obs_err_type == 'iasi':
        return IasiObservationErrorSps(data_filename)
    else:
        raise UnknownObservationType

