#!/usr/bin/env pythonchL

"""
AMETHYST Configuration File 

This Python file contains all default variables used within the Amethyst software.
These are gathered in four dictionaries based on the main sections of the software: 
            - common_vars
            - preprocessor_vars
            - processor_vars
            - postprocessor_vars
            
:copyright: 2021, Adaptive Meteo S.r.l.
"""

__author__ = 'Paolo Scaccia'
__copyright__ = "Copyright 2021, Adaptive Meteo S.r.l."
__credits__ = ["Paolo Scaccia", "Paolo Antonelli"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Paolo Scaccia <paolo.scaccia@adaptivemeteo.com>"
__email__      = "paolo.scaccia@adaptivemeteo.com"

import os

AMETHYST_PATH = os.path.dirname(__file__)

common_vars = {
                  "instrument" : "cris",
                  "fm_version" : 3,
                  "levels"     : 81,
                  "wrfdir"     : "/mnt/satellite/amethyst_test_data/wrf",
                  "wrkdir"     : "/dev/shm/amethyst_test/",
                  
                  # Arctic
                  "geobox"     : { "lonmin" : -180.0,
                                   "lonmax" :  180.0,
                                   "latmin" :  70.0,
                                   "latmax" :  90.0},

                  # Pacific
                  # "geobox"     : { "lonmin" : -166.0,
                  #                  "lonmax" : -148.0,
                  #                  "latmin" :  11.0,
                  #                  "latmax" : 29.0},

             }


preprocessor_vars = {
                  "basedir"   : "/mnt/satellite/amethyst_test_data",
                  "iasidir"   : "/mnt/satellite/amethyst_test_data/iasi/20200820_062454",
                  "gcrso"     : "/mnt/satellite/amethyst_test_data/cris/20220820_104410/GCRSO_j01_d20220820_t1046239_e1058057_b24631_c20221011034338465526_cspp_dev.h5",
                  "scris"     : "/mnt/satellite/amethyst_test_data/cris/20220820_104410/SCRIF_j01_d20220820_t1046239_e1058057_b24631_c20221011034338607949_cspp_dev.h5",
                  "l1dir"     : "/mnt/satellite/amethyst_test_data/cris",
                  "apriori"   : AMETHYST_PATH + "/ancillary/atmosphere/arctic_apriori.nc",
                  "iasi_cmt"  : 5,
                  "cris_cmt"  : .95
            }
preprocessor_vars["bindir"] = preprocessor_vars["basedir"] + "/bin"
preprocessor_vars["logdir"] = preprocessor_vars["basedir"] + "/log/mirto"
preprocessor_vars["rundir"] = { "cris" : preprocessor_vars["basedir"] + "/run_cris",
                                "iasi" : preprocessor_vars["basedir"] + "/run_iasi"
                               }

preprocessor_vars["first_guess"] = { "top_pressure"             : 0.005,
                                     "bottom_pressure"          : 1013,
                                     "min_water_vapor"          : 0.003,
                                     "p_min_water_vapor"        : 45,
                                     "water_vapor_smooth_after" : 60
                                    }
preprocessor_vars["climatology"] = { "h2o"         : AMETHYST_PATH + "/ancillary/atmosphere/h2o_climatology.nc",
                                     "temperature" : AMETHYST_PATH + "/ancillary/atmosphere/temperature_climatology.nc",
                                     "o3"          : AMETHYST_PATH + "/ancillary/atmosphere/o3_climatology.nc"
                                    }



processor_vars = {
                # Instrument
                "instrument"      : common_vars["instrument"],
                "noise_file"      : AMETHYST_PATH + "/ancillary/instrument/cris_obserr.nc",
                "instr_chan_list" : AMETHYST_PATH + "/ancillary/instrument/cris_FSR_chList_862.dat",
                "od_file"         : AMETHYST_PATH + "/ancillary/forward_model/cris_coefficients_oss_v3.nc",
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
                "selected_state_vector_variables" : [ -2,
                                                      -1,
                                                       0,
                                                       1,
                                                       3 ],
                "fix_gamma_zero"                  : 0,
                "noise_scaling_factor"            : 1.0,
                "chi_square_limit"                : 0.01,
                "ml_gamma"                        : 3.0,
                "ml_gamma_increase_factor"        : 5.0,
                "ml_gamma_decrease_factor"        : 2.0,
                "minimum_fraction_rate_change"    : 0.03,
                
                # Transformed Retrieval
                "tr_chan_list"   : AMETHYST_PATH + "/ancillary/instrument/cris_FSR_chList_862_tr.dat",
                
                # Amethyst Output
                "output_file"  : "amethyst_output.nc",
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
                
                # Forward Model
                "oss_obslevel" : 0
                }

postprocessor_vars = {
                """
                      INSERT POSTPROCESSOR VARS  
                """
            }
