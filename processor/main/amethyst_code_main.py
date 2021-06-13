#!/usr/bin/env python

"""
Main class providing the inverter method for the UW PHYS retrieval
Implemented 1DVar algorithm using the ForwardModel class
"""

from __future__ import print_function
import numpy as np
from scipy.linalg import lu, solve

from main.ForwardModel import ForwardModel
from main.ForwardModel import NotConvergentIteration
import sys
if "amethyst_config" not in sys.modules:
    import amethyst_config


class amethyst_state(object):
    """
    This class stores the inversion results along with intermediate results
    """
    def __init__(self):
        """ Initialize all object attributes to None type """
        self.Sa_ret = None
        self.SaInv_ret = None
        self.xa = None
        self.xhat = None
        self.xhat_pre = None
        self.xhat_new = None
        self.d2 = None
        self.fm = None
        self.yobs_minus_yhat = None
        #PaoloA 12112018
        #self.jacobian - None
        #self.residuals = None
        #self.fgresiduals = None

class amethyst_core_config(object):
    """ Utility class to keep in one place all needed by Forward Model"""
    def __init__(self, model):
        """ Initialize all class attributes """
        self.model = model
        self.surfacePressure_mb = amethyst_config.processor_vars["constant_pressure"]
        self.observationAltitude_km = amethyst_config.processor_vars["observation_altitude"]["value"]
        self.observationPressure_mb = amethyst_config.processor_vars["observation_pressure"]["value"]
        self.EstimateK = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=int)
        self.Jvar = np.array(amethyst_config.processor_vars["selected_state_vector_variables"])
        
        SKTflag = np.size(np.where(self.Jvar == -1)) > 0
        Tflag = np.size(np.where(self.Jvar == 0)) > 0
        WVflag = np.size(np.where(self.Jvar == 1)) > 0
        CO2flag = np.size(np.where(self.Jvar == 2)) > 0
        O3flag = np.size(np.where(self.Jvar == 3)) > 0
        
        self.variable_selection = [ SKTflag, Tflag, WVflag, CO2flag, O3flag ]

        self.Iteration_limit = sum(self.EstimateK)
        self.retrievalMLGamma      = amethyst_config.processor_vars["ml_gamma"]
        self.retrievalFixGammaZero = amethyst_config.processor_vars["fix_gamma_zero"]
        self.MFRC                  = amethyst_config.processor_vars["minimum_fraction_rate_change"]
        self.MLGIF                 = amethyst_config.processor_vars["ml_gamma_increase_factor"]
        self.MLGDF                 = amethyst_config.processor_vars["ml_gamma_decrease_factor"]
        self.xdim = None
        self.indx = None
        self.state_var_indx = None
        self.fov_latitude = None
        self.fov_longitude = None
        self.fov_time = None
        self.FOVangle = None
        self.wnR = None
        self.R = None
        # Hardcode this, should be function of latitude an time of day....
        self.SunAngle = 90.0
        self.pressure_grid = None
        self.p = None
        self.emiss = None
        self.sfgrd = None
        self.emrf = None


class amethyst_apriori(object):
    """ Utility class to keep apriori estimate and model error """
    def __init__(self):
        """ Initialize all class attributes """
        self.x0 = None
        self.xa = None
        self.Sa = None
        self.SaInv = None

class core(object):
    """
    This class is the core of the inversion system
    """
    def __init__(self, forward_model, obserr):
        """ Load configuration parameters and Input data for retrieval """
        self.cx = amethyst_core_config(forward_model)
        self.obs_err = obserr.obs_err
        self.dot_inv_obs_err = obserr.dot_inv_obs_err
        self.inv_obs_err_dot = obserr.inv_obs_err_dot
        self.apriori = amethyst_apriori()
        self.state = amethyst_state()
        self.yobs_minus_yhat = None
        self.mspo = np.NAN
        self.fm = None
        # Default Marquardt-Levemberg parameter
        self.gamma = 0.0

    def compute_chi_square(self, profile):
        """ Compute X^2 from retrieval residuals """
        if profile.prlinalg is not None:
            profile.prlinalg.enable()
        self.mspo = np.dot(self.yobs_minus_yhat.T,
                           self.inv_obs_err_dot(self.yobs_minus_yhat))
        self.mspo /= len(self.yobs_minus_yhat)

        if profile.prlinalg is not None:
            profile.prlinalg.disable()

    def update_solution(self, fm, profile):
        """ Update state vector at every iteration """
        # From Eq. 3.35 on page 93 of C.Rodgers "Inverse Method for Atmospheric Sounding" 
        if profile.prlinalg is not None:
            profile.prlinalg.enable()
        
        KtSeInv = self.dot_inv_obs_err(fm.K.T)
        KtSeInvK = np.dot(KtSeInv, fm.K)
        A = (KtSeInvK + (1.0+self.gamma)*self.state.SaInv_ret)
        dx = (self.state.xhat - self.state.xa)
        d = (np.dot(KtSeInv, self.yobs_minus_yhat) -
             np.dot(self.state.SaInv_ret, dx))
        
        # Use iterative LU decomposition to determine the solution
        # First iteration
        try:
            [L, U] = lu(A, permute_l=True)
        except ValueError:
            raise NotConvergentIteration("Found an inf or a " 
                                         "NaN in the A matrix!")
        y = solve(L, d)
        x = solve(U, y)
        # Second iteration
        r = d - np.dot(A, x)
        dz = solve(L, r)
        ddx = solve(U, dz)
        # Solution
        totx = x + ddx
        self.state.xhat_new = self.state.xhat + totx
        # Calculate rate of change in the state vector
        # d2 is used to evaluate if convergence criterium is met
        self.state.d2 = np.dot(totx.T, d)
        if profile.prlinalg is not None:
            profile.prlinalg.disable()

    def invert(self, fov, fg, apriori, emiss, obs, log, profile=None):
        """ Invert the measurement to get physical sounding profile """
        L = log
        if self.cx.retrievalFixGammaZero == 0:
            # Marquardt-Levemberg parameter
            self.gamma = self.cx.retrievalMLGamma
        self.mspo = np.NAN
        mspo_hist = []

        self.cx.xdim = fg.xdim
        # Index in file 1 based. Python zero based.
        self.cx.indx = fov.indx
        self.cx.state_var_indx = apriori.varindx
        converge = self.cx.MFRC * len(self.cx.state_var_indx)

        [self.cx.fov_latitude, self.cx.fov_longitude, self.cx.fov_time,
         self.cx.FOVangle, self.cx.wnR, self.cx.R] = fov.fov(obs)

        [self.cx.p, self.apriori.x0, self.apriori.xa] = fg.state_vector(obs)
        self.cx.pressure_grid = self.cx.p[0:self.cx.xdim[0]]
        self.cx.surfacePressure_mb = self.cx.pressure_grid[0]
        self.cx.observationPressure_mb = self.cx.pressure_grid[-1]
        L.log('OBS ' + str(obs) + ': Surface   Pressure = '+
              repr(self.cx.surfacePressure_mb), 4, False)
        L.log('OBS ' + str(obs) + ': Satellite Pressure = '+
              repr(self.cx.observationPressure_mb), 4, False)

        # Load apriori
        [self.apriori.Sa, self.apriori.SaInv] = apriori.covariance_matrix(obs)

        # Assign values for the state vector
        xhat = np.copy(self.apriori.x0)

        # Iteration of the Newton-Gauss method to find the zero of the first
        # derivative of the Gaussian PDF
        Iteration = 0

        # Initialize state vector per previous iteration to apriori.xa
        # (same as apriori.X0)
        xhat_pre = np.copy(self.apriori.xa)
        self.cx.emiss = emiss

        fm = ForwardModel(self.cx)
        self.fm = fm

        ems = fg.xdim[0]+fg.xdim[1]+fg.xdim[2]+fg.xdim[3]+fg.xdim[4]
        eme = fg.xdim[0]+fg.xdim[1]+fg.xdim[2]+fg.xdim[3]+fg.xdim[4]+fg.xdim[5]

        # Selection indexes for variables
        jj = self.cx.state_var_indx

        # Self.state.xxxxx has only the varaibles which will be updated
        # also referred to as reduced or selected variables
        # Starting from October 2014 the input files
        # (apriori.nc, fov.nc, and observationError.nc)
        # store only the selected variables.
        # This was done to reduce the size of the input files.
        self.state.Sa_ret = self.apriori.Sa
        self.state.SaInv_ret = self.apriori.SaInv
        self.state.xa = self.apriori.xa[jj]



        while (Iteration < self.cx.Iteration_limit):
            #
            # Compute F (forward model)
            #
            try:
                [self.cx.sfgrd, self.cx.emrf] = emiss.get(obs, xhat[ems:eme])
            except FloatingPointError:
                raise NotConvergentIteration("Can not get emissivity: overflow!")
            fm.compute_forward(xhat)

            self.state.fm = fm

            self.yobs_minus_yhat = fm.compute_residuals()

            #PaoloA 12112018
            if Iteration == 0:
                self.state.fgresiduals = self.yobs_minus_yhat 
            # Subselect from forward model output (Jacobians)
            # the variables actually retrieved.
            # Radiance and Jacobian calculated with the forward model
            # are already reduced in the wavenumber space (only
            # channels used in the inveersion are received from the
            # forward model)
            fm.K = fm.K[:, jj]
            # self.state.xxxx has only the components of the state vector
            # which are retrieved
            self.state.xhat = xhat[jj]
            self.state.xhat_pre = xhat_pre[jj]

            self.compute_chi_square(profile)
            self.update_solution(fm, profile)

            # By setting update_xhat to True the state vector is updated
            # regardless of the increase/decrease of the cost function
            # By setting update_xhat to False the state vector is updated
            # only if the cost function decreases
            update_xhat = True
            # update_xhat = False
            if Iteration > 0:
                ref_norm = min(mspo_hist)
                L.log('OBS ' + str(obs) + ': Actual value of Norm   : '+repr(self.mspo), 4, False)
                L.log('OBS ' + str(obs) + ': Last low value of Norm : '+repr(ref_norm), 4, False)

                # This part represents the implementation of the
                # Levemberg-Marquardt algorithm.
                # However the true algorithm should use of the ratio
                # between
                # chi_sq_optimal = (yobs_minus_yhat'*
                #        (M*yobs_minus_yhat))/length(yobs_minus_yhat);
                # and
                # chi_sq_linear = (yobs_minus_yhat_lin'*
                #      (M*yobs_minus_yhat_lin))/length(yobs_minus_yhat_lin);
                # where:
                # M = KtSeInv'*Sa_ret*KtSeInv + SeInv;
                # yobs_minus_yhat_lin =
                #    yobs_minus_yhat(:,control.fmiter-1) - K*(xhat-xhat_pre)
                # This implementation instead, uses only
                # self.norm.mspo to decide if and how to
                # change gamma.

                if (self.mspo <= ref_norm):
                    if self.cx.retrievalFixGammaZero == 0:
                        self.gamma = self.gamma / self.cx.MLGDF
                    xxdel = (100.0 * (ref_norm - self.mspo) / self.mspo)
                    L.log('OBS ' + str(obs) + ': MIRTO residuals decreased by ' +
                          repr(xxdel)+'%', 4, False)
                    update_xhat = True
                else:
                    if self.cx.retrievalFixGammaZero == 0:
                        self.gamma = self.gamma * self.cx.MLGIF
                    xxdel = (100.0 * (self.mspo - ref_norm) / self.mspo)
                    L.log('OBS ' + str(obs) + ': MIRTO residuals increased by ' +
                          repr(xxdel)+'%', 4, False)

            L.log('OBS ' + str(obs) + ': Iteration = ' + repr(Iteration), 4, False)
            L.log('OBS ' + str(obs) + ': New Gamma = ' + repr(self.gamma), 4, False)
            L.log('OBS ' + str(obs) + ': Distance  = ' + repr(abs(self.state.d2)), 4, False)
            L.log('OBS ' + str(obs) + ': Wanted    = ' + repr(converge), 4, False)

            if (abs(self.state.d2) < converge):
                # If convergence criterium is met, set gamma to 0 and
                # run one more iteration
                if self.cx.retrievalFixGammaZero == 0:
                    self.gamma = 0.0
                if update_xhat:
                    xhat_pre[jj] = self.state.xhat
                    xhat[jj] = self.state.xhat_new
                [self.cx.sfgrd, self.cx.emrf] = emiss.get(obs,
                                                          xhat[ems:eme])
                fm.compute_forward(xhat)
                self.state.fm = fm
                self.yobs_minus_yhat = fm.compute_residuals()
                fm.K = fm.K[:, jj]
                self.state.xhat = xhat[jj]
                self.state.xhat_pre = xhat_pre[jj]
                self.compute_chi_square(profile)
                self.update_solution(fm, profile)
                xhat[jj] = self.state.xhat
                [self.cx.sfgrd, self.cx.emrf] = emiss.get(obs,
                                                          xhat[ems:eme])
                fm.compute_forward(xhat)

                # Subselect only retrieved state vector variables
                # The below for test transform
                #np.save('fmk',fm.K[:,jj])
                #np.save('fmxhat',xhat[jj])
                #sys.exit(0)
                # The below for test transform
                fm.K = fm.K[:, jj]
                self.yobs_minus_yhat = fm.compute_residuals()
                self.fm.K = fm.K
                # Added to enable one more option in computing
                # transformed retrievals
                self.state.yobs_minus_yhat = self.yobs_minus_yhat
                L.log('**** OBS ' + str(obs) + ' CONVERGED!!! ****', 4, False)
                break
            else:
                if (Iteration < self.cx.Iteration_limit -1):
                    # Assign new value to the solution and continue the
                    # iterative solution
                    if update_xhat:
                        xhat_pre[jj] = self.state.xhat
                        xhat[jj] = self.state.xhat_new
                else:
                    raise NotConvergentIteration("Limit of Iterations reached for "
                                                 " OBS " + str(obs)) 

            mspo_hist.append(self.mspo)
            Iteration = Iteration + 1

        #PaoloA 14112018
        #self.state.jacobian = self.fm.K
        self.state.jacobian = fm.K
        target_length = 249
        axis = 1 
        pad_size = target_length - fm.K.shape[axis]
        axis_nb = len(fm.K.shape)

        if pad_size >=  0:
            npad = [(0, 0) for x in range(axis_nb)]
            npad[axis] = (0, pad_size)
            self.state.jacobian = np.pad(fm.K, pad_width=npad, mode='constant', constant_values=0)
        else:
            self.state.jacobian = self.fm.K 

        #print("fm.K shape {}".format(self.state.jacobian.shape))

        #PaoloA 12112018
        self.state.residuals = self.state.yobs_minus_yhat

        return(self.state)
