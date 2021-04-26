#!/usr/bin/env python
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

common_vars = {
                  "levels"  : 81,
                  "wrfdir"  : "/mnt/satellite/amethyst_test_data/wrf",
                  "geobox"  : { "lonmin" : -166.0,
                                "lonmax" : -148.0,
                                "latmin" :  11.0,
                                "latmax" : 29.0},
             }

preprocessor_vars = {
                  "basedir"   : "/mnt/satellite/amethyst_test_data",
                  "iasidir"   : "/mnt/satellite/amethyst_test_data/iasi/20201120_062454",
                  "gcrso"     : "/mnt/satellite/amethyst_test_data/cris/20201120_002622/GCRSO_npp_d20201120_t0028319_e0040137_b46969_c20210315153613150425_cspp_dev.h5",
                  "scris"     : "/mnt/satellite/amethyst_test_data/cris/20201120_002622/SCRIS_npp_d20201120_t0028319_e0040137_b46969_c20210315153613291913_cspp_dev.h5",
                  "l1dir"     : "/mnt/satellite/amethyst_test_data/cris",
                  "iasi_cmt"  : 5,
                  "cris_cmt"  : .95
            }
preprocessor_vars["bindir"] = preprocessor_vars["basedir"] + "/bin/"
preprocessor_vars["logdir"] = preprocessor_vars["basedir"] + "/log/mirto/"
preprocessor_vars["rundir"] = { "cris" : preprocessor_vars["basedir"] + "/run_cris/",
                                "iasi" : preprocessor_vars["basedir"] + "/run_iasi/"
                               }

processor_vars = {
                """
                      INSERT PROCESSOR VARS  
                """
            }

postprocessor_vars = {
                """
                      INSERT POSTPROCESSOR VARS  
                """
            }
