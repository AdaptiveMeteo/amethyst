from dobjects.ObservationErrors.Iasi import IasiObservationError
from dobjects.ObservationErrors.Cris import CrisObservationError
from dobjects.ObservationErrors.Iasi import IasiObservationErrorSps
from dobjects.ObservationErrors.Cris import CrisObservationErrorSps

class UnknownObservationType(Exception):
    pass


def create_obs_err(data_filename, obs_err_type='cris'):
    if 'cris' in obs_err_type:
        return CrisObservationError(data_filename)
    elif 'iasi' in obs_err_type:
        return IasiObservationError(data_filename)
    else:
        raise UnknownObservationType

def create_obs_err_sps(data_filename, obs_err_type='cris'):
    if 'cris' in obs_err_type:
        return CrisObservationErrorSps(data_filename)
    elif 'iasi' in obs_err_type:
        return IasiObservationErrorSps(data_filename)
    else:
        raise UnknownObservationType

