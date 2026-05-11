#!/usr/bin/env python

"""
AMETHYST Configuration File — IASI instrument with RTTOV forward model

Working directory: /dev/shm/amethyst_test_iasi/
Test data:         /mnt/miwa/ann/amethyst/test_data/iasi/20220808_134754/mirto/
"""

__author__ = 'Paolo Antonelli'
__copyright__ = "Copyright 2024, Adaptive Meteo S.r.l."

import os

AMETHYST_PATH = os.path.dirname(os.path.realpath(__file__))
RTTOV_PATH    = "/home/paoloa/projects/aurora/ForwardModels/RTTOV"
RTTOV_COEF    = (RTTOV_PATH +
                 "/rtcoef_rttov14/rttov13pred101L/rtcoef_metop_2_iasi_7gas.nc")

### COMMON
common_vars = {
                  "instrument" : "iasi",
                  "fm_version" : 3,      # ignored for RTTOV, kept for compatibility
                  "levels"     : 101,    # RTTOV 101-level predictor
                  "wrfdir"     : "/mnt/miwa/ann/amethyst/test_data/wrf",
                  "wrkdir"     : "/dev/shm/amethyst_test_iasi",

                  # Arctic
                  "geobox"     : { "lonmin" : -180.0,
                                   "lonmax" :  180.0,
                                   "latmin" :  70.0,
                                   "latmax" :  90.0},
             }

#### PRE-PROCESSOR  (unchanged from OSS config)
preprocessor_vars = {
                  "basedir"          : "/mnt/miwa/ann/amethyst/test_data",
                  "iasidir"          : "/mnt/miwa/ann/amethyst/test_data/iasi/20220808_134754",
                  "gcrso"            : "/mnt/miwa/ann/amethyst/test_data/cris/20220820_104410/GCRSO_j01_d20220820_t1046239_e1058057_b24631_c20221011034338465526_cspp_dev.h5",
                  "scris"            : "/mnt/miwa/ann/amethyst/test_data/cris/20220820_104410/SCRIF_j01_d20220820_t1046239_e1058057_b24631_c20221011034338607949_cspp_dev.h5",
                  "l1dir"            : "/mnt/miwa/ann/amethyst/test_data/iasi/",
                  "static_apriori"   : AMETHYST_PATH + "/ancillary/atmosphere/arctic_apriori.nc",
                  "iasi_cmt"         : 5,
                  "cris_cmt"         : .95,
                  "first_guess"      : {
                                          "top_pressure"             : 0.005,
                                          "bottom_pressure"          : 1000,
                                          "surface_pressure"         : 1013,
                                          "min_water_vapor"          : 0.003,
                                          "p_min_water_vapor"        : 45,
                                          "water_vapor_smooth_after" : 60
                                        },
                  "climatology"      : {
                                          "h2o"         : AMETHYST_PATH + "/ancillary/atmosphere/h2o_climatology.nc",
                                          "temperature" : AMETHYST_PATH + "/ancillary/atmosphere/temperature_climatology.nc",
                                          "o3"          : AMETHYST_PATH + "/ancillary/atmosphere/o3_climatology.nc"
                                       }
                  }
preprocessor_vars["bindir"] = preprocessor_vars["basedir"] + "/bin",
preprocessor_vars["logdir"] = preprocessor_vars["basedir"] + "/log/mirto"
preprocessor_vars["rundir"] = {
                                "cris" : preprocessor_vars["basedir"] + "/run_cris",
                                "iasi" : preprocessor_vars["basedir"] + "/run_iasi"
                                }

#### PROCESSOR
processor_vars = {
                # Instrument
                "instrument"      : common_vars["instrument"],
                "noise_file"      : AMETHYST_PATH + "/ancillary/instrument/iasi_obserr.nc",
                "instr_chan_list" : AMETHYST_PATH + "/ancillary/instrument/iasi_chList.dat",

                # Forward model — RTTOV
                "fm_engine"       : "rttov",
                "rttov_path"      : RTTOV_PATH,
                "rttov_coef_file" : RTTOV_COEF,
                "rttov_surftype"  : 2,    # 0=land, 1=sea, 2=sea ice (Arctic)
                # lib_path: prepended to LD_LIBRARY_PATH so pyrttov finds the
                # correct libnetcdf (needed when compiled in a different env).
                "rttov_lib_path"  : "/home/paoloa/anaconda3/envs/rttov/lib",

                # Not used by RTTOV but kept for API compatibility
                "od_file"         : AMETHYST_PATH + "/ancillary/forward_model/iasi_coefficients.nc",

                "fov_file"        : common_vars['wrkdir'] + "/fov.nc",
                "fg_file"         : common_vars['wrkdir'] + "/fg.nc",
                "apriori_file"    : common_vars['wrkdir'] + "/apriori.nc",

                # Geometry
                "observation_altitude" : {"units":"Km",   "value":100.0},
                "observation_pressure" : {"units":"hPa",  "value":0.005},

                # State Vector
                "constant_pressure"                    : {"units":"mb",   "value":1013},
                "constant_co2"                         : {"units":"ppmv", "value":405.0},
                "constant_co2_std"                     : {"units":"ppmv", "value":16.0},
                "constant_skt_std"                     : {"units":"K",    "value":3.0},
                "constant_inflation_surface_emiss_cov" : {"units":1,"value":1},
                "constant_solar_irradiance_file"       : AMETHYST_PATH + "/ancillary/atmosphere/solar_irradiances.nc",
                "eigenforland"         : 4,
                "eigenforsea"          : 3,
                "eigenfortcov"         : 9,
                "clear_outside_diag"   : False,

                # Retrieval
                "selected_state_vector_variables" : [ -2, -1, 0, 1, 3 ],
                "fix_gamma_zero"                  : 0,
                "noise_scaling_factor"            : 1.0,
                "chi_square_limit"                : 0.01,
                "ml_gamma"                        : 3.0,
                "ml_gamma_increase_factor"        : 5.0,
                "ml_gamma_decrease_factor"        : 2.0,
                "minimum_fraction_rate_change"    : 0.03,
                "max_iterations"                  : 30,

                # Transformed Retrieval
                "tr_chan_list"   : AMETHYST_PATH + "/ancillary/instrument/iasi_chList_tr.dat",

                # Amethyst Output
                "output_file"  : "amethyst_output_iasi_rttov.nc",
                "output_vars"  : {
                                     "pressure"                         : True,
                                     "temperature"                      : True,
                                     "water_vapor"                      : True,
                                     "ozone"                            : True,
                                     "surface_temperature"              : True,
                                     "surface_emissivity_coefficients"  : True,
                                     "d2"                               : True,
                                     "da_r"                             : True,
                                     "transformed_retrievals"           : True,
                                     "sn_eigenvalues"                   : True,
                                     "assimilation_operator"            : True,
                                     "residuals"                        : False,
                                     "fg_residuals"                     : False,
                                     "jacobian"                         : False,
                                     "sa"                               : False,
                                     "indices"                          : False},

                # Not used by RTTOV
                "oss_obslevel"   : 0
                }

### POST-PROCESSOR
postprocessor_vars = {}
