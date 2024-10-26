#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  1 16:03:49 2018

@author: Paolo Antonelli, Paolo Scaccia
"""
from wrf import rh
import numpy as np

#*************************************************************************
def mr2rh(PRES,TEMP,VAPOUR,*TEMP_CONV):
    """
    
    Parameters
    ----------
    PRES :    np.ndarray    
        Pressure levels in hPa
        NOTE: wrr-py.rh read pressure in Pa

    TEMP :    np.ndarray  
        Temperature profile in K
    VAPOUR :  np.ndarray  
        Water vapour mixing ratio in g/Kg
        NOTE: wrr-py.rh read vapour in Kg/kg

    Returns
    -------
    [ rh , rh ]            Duplicated output to be coherent with MIRTO routines
    rh :      np.ndarray
        Relative humidity
    
    """
    hum = rh(VAPOUR/1000,PRES*100,TEMP,meta=False)
    return [hum, hum]
"""
#             OLD VERSION OF MR2RH
def mr2rh(PRES,TEMP,VAPOUR,*TEMP_CONV):   
    #
    # function [rh1,rh2] = mr2rh(p,t,w,Tconvert);
    #
    # determine relative humidity (#) given
    # reference pressure (mbar), temperature (t,K), and
    # water vapor mass mixing ratio (w,g/kg)
    #
    # Two RHs are returned: rh1 is as the ratio 
    # of water vapor partial pressure to saturation vapor pressure and
    # rh2 is defined as the ratio of water vapor mixing ratio to 
    # saturation mixing ratio.
    #
    # if input, Tconvert is used as the temperature point to switch
    # from using saturation vapor pressure over water to over ice.
    #
    # DCT 3/5/00
    #

    TEMP = np.ma.masked_invalid(TEMP)
    PRES = np.ma.masked_invalid(PRES)
    VAPOUR = np.ma.masked_invalid(VAPOUR)

    

    if len(TEMP_CONV) == 0:
        esat = satvap(TEMP) 
    else:
        esat = satvap(TEMP,TEMP_CONV[0]) 
    
    
    wsat = satmix(PRES,esat)

    # H20 partial pressure
    e = mr2e(PRES,VAPOUR)
    
    # rh using ratios of gas pressure
    rh1 = 100*e/esat
    
    # rh using WMO definition of rel. humidity
    rh2 = 100*VAPOUR/wsat
    
    return rh1,rh2
"""
#*************************************************************************    
def satvap(TEMP,*TEMP_CONV):
    """
    compute saturation vapor pressure [mbar] given temperature, T [K].
    If Tconvert is input, the calculation uses the saturation vapor 
    pressure over ice (opposed to over water) for temperatures less than 
    Tconvert [K].
    
    DCT, updated 3/5/00
    """  

    if len(TEMP_CONV) == 0:
        esat = eswat_goffgratch(TEMP) 
    else:
        esat = esice_goffgratch(np.ma.masked_array(TEMP,TEMP<=TEMP_CONV[0]))

    return esat

#**************************************************************************    
    
def satmix(PRES,ESAT):
    """
        compute saturation mixing ratio [g/kg] given reference pressure, 
     p [mbar] and temperature, T [K].  If Tconvert input, the calculation uses 
     the saturation vapor pressure over ice (opposed to over water) 
     for temperatures less than Tconvert [K].
    
     DCT, updated 3/5/00
    """
    # saturation mixing ratio
    wsat = e2mr(PRES,ESAT)

    return wsat

#*************************************************************************
def e2mr(PRES,e):
    """
    compute H2O mass mixing ratio (mr,g/kg) given
    pressure (p,mbar) and H2O partial pressure (e,mbar)
    """
    eps = 0.621970585
    return eps*1000*e/(PRES-e)
#*************************************************************************
def mr2e(PRES,MR):
    """
    Compute H20 partial pressure (mbar) give pressure (PRES,mbar)
    and H20 mass mixing ratio (MR,g/kg)
    """
    eps = 0.621970585
    return PRES*MR/(1000*eps+MR)

#*************************************************************************
def esice_goffgratch(TEMP):
    """
    Compute water vapour saturation pressure over ice
    using Goff-Gratch formulation. Adopted from PvD's svp_ice.pro
    """
    ewi = 6.1071
    c1 =  9.09718
    c2 = 3.56654
    c3 = 0.876793
    ratio = 273.15 / TEMP
    tmp = ( -c1*( ratio - 1.0 ) )- \
          (  c2*np.log10(ratio) )      + \
          (  c3*( 1.0 - ( 1.0 / ratio ) ) ) +\
          np.log10( ewi )

    return 10.0**tmp
#*************************************************************************
def eswat_goffgratch(TEMP):
    """
    Compute water vapour saturation pressure over water
    using Goff-Gratch formulation. Adopted from PvD's svp_water.iso
    """
    t_sat = 373.16
    t_ratio = t_sat/TEMP
    rt_ratio = 1.0/t_ratio
    sl_pressure = 1013.246
    
    c1 = 7.90298
    c2 = 5.02808
    c3 = 1.3816e-7
    c4 = 11.344
    c5 = 8.1328e-3
    c6 = 3.49149

    tmp = ( -1.0*c1*( t_ratio - 1.0 ) ) + \
          ( c2 * np.log10(t_ratio) )    - \
          ( c3 * ( 10.0**( c4 * ( 1.0 - rt_ratio ) ) - 1.0 ) ) + \
          ( c5 * ( 10.0**( -c6 * ( t_ratio - 1.0 ) ) - 1.0 ) ) + \
          np.log10( sl_pressure )
    return 10.0**tmp
#*************************************************************************
def mr2dp(PRES,TEMP,MR,*T_CONV):
    """
    Compute dew point temperature given watervapour mass mixing ratio (g/kg).
    Uses Goff-Gratch formulation.
    If input, the air temperature (t,K) is used to determine if
    saturation over water (for t > Tconvert) or saturation over ice
    (for t <= Tconvert) is used.
    
    Two dew points are returned: dp1 is with RH defined as the ratio 
    of water vapor partial pressure to saturation vapor pressure and
    dp2 is with RH defined as the ratio of water vapor mixing ratio to 
    saturation mixing ratio.
    
    Notes: results not valid for dew points >= 370 K and  <= 160 K.
    
    Reference: Goff-Gratch formulation from sixth revised 
               edition of Smithsonian Meteorology Tables.

    """
    # compute H20 partial pressure
    e = mr2e(PRES,MR)

    # interpolate (saturation pressure vs T) to desired pressure
    T = np.arange(100,400.2,0.2)
    esat_water = eswat_goffgratch(T)
    dp = np.interp(e,esat_water,T)
    
    if len(T_CONV)!=0 : 
        esat_ice = esice_goffgratch(T);
        dp[TEMP<=T_CONV] = np.interp(e[TEMP<T_CONV],esat_ice,T)

    return dp
#************************************************************************#*****************************************<*******************************

def interpolate_pressure_grid(new_pressure,pressure,observable,*args):

    # Check order
    flip = 1 if pressure[0]>pressure[1] else 0

    xp = np.array(new_pressure) if not flip else np.flip(new_pressure,0)
    x = np.array(pressure) if not flip else np.flip(pressure,0)
    y = np.array(observable) if not flip else np.flip(observable,0)

    if len(args)==0:
        args = 'default'
    1
    # Choose method
    if args == 'default':
        yp = np.interp(np.log(xp),np.log(x),y,left=x[0],right=x[-1])           
    elif args == 'special1':
        raise TypeError('Not yet implemented...')
    elif args == 'special2':
        raise TypeError('Not yet implemented...')
    else:
        raise NameError('Interpolation Function Unknown')
    
    # if flipped return correct array 
    if flip:
        yp = np.flip(yp,0)
        
    return yp

#**********************************************************************************
def hyposometric(pressure,temperature,mixing_ratio,surface_alt):

	# NAME:
	#       geopotential_altitude   - Paul van Delst's name of the routine
	#                       (DT renamed it to this).
	#
	# PURPOSE:
	#       This function calculates geopotential altitudes using the hypsometric equation.
	#    It is like the "hypsometric" function I have implemented, but this one
	#    also accounts for water vapor in the atmosphere.
	#
	#
	# CATEGORY:
	#       Meteorology
	#
	# CALLING SEQUENCE:
	#       altitude = hypsometric( pressure, $          ; Input
	#                               temperature, $       ; Input
	#                               mixing_ratio, $      ; Input
	#                               surface_altitude, $  ; Input
	#
	# INPUTS:
	#       pressure:        Array of pressure in mb
	#       temperature:     Array of temperature in K
	#       mixing_ratio:    Array of water vapor mixing ratio in g/kg
	#       surface_height:  Height of surface in km. 
	#
	# KEYWORD PARAMETERS:
	#       None.
	#
	# OUTPUTS:
	#       altitude:        Geopotential altitudes in km in the same order as the
	#                          input data.
	#
	# CALLS:
	#       None.
	#
	#
	# COMMON BLOCKS:
	#       None
	#
	# SIDE EFFECTS:
	#       None
	#
	# RESTRICTIONS:
	#       - Input pressure, temperature and mixing_ratio MUST be arrays with at
	#         LEAST two elements where the element with the highest pressure corresponds
	#         to the surface altitude passed.
	#       - Input surface height must be > or = to 0.0km and a SCALAR.
	#       - No allowance is made, yet, for the change in the acceleration due to 
	#         gravity with altitude.
	#
	# PROCEDURE:
	#       Geopotential heights are calculated using the hypsometric equation:
	#
	#                         -
	#                    Rd * Tv    [  p1  ]
	#         z2 - z1 = --------- ln[ ---- ]
	#                       g       [  p2  ]
	#
	#       where Rd    = gas constant for dry air (286.9968933 J/K/kg),
	#             g     = acceleration due to gravity (9.80616m/s^2),
	#             Tv    = mean virtual temperature for an atmospheric layer,
	#             p1,p2 = layer boundary pressures, and
	#             z1,z2 = layer boundary heights.
	#
	#       The virtual temperature, the temperature that dry air must have in
	#       order to have the same density as moist air at the same pressure, is
	#       calculated using:
	#
	#                  [      1 - eps      ]
	#         Tv = T * [ 1 + --------- * w ]
	#                  [        eps        ]
	#
	#       where T   = temperature,
	#             w   = water vapor mixing ratio, and
	#             eps = ratio of the molecular weights of water and dry air (0.621970585).
	#
	# EXAMPLE:
	#       Given arrays of pressure, temperature, and mixing ratio:
	#
	#         IDL> PRINT, p, t, mr
	#               1015.42      958.240
	#               297.180      291.060
	#               7.83735      5.71762
	#
	#       the geopotential altitudes can be found by typing:
	#
	#         IDL> result = geopotential_altitude( p, t, mr, 0.0, alt )
	#         IDL> PRINT, result, alt
	#                1      0.00000     0.500970
	#
	# MODIFICATION HISTORY:
	#       Written by:     Paul van Delst, CIMSS/SSEC, 08-Dec-1997
	#
	# Log: hypsometric.pro,m  
	# converted to matlab code - leslie moy 02-02-02
	#
	# $Log: hypsometric2.pro,v $
	# Revision 1.1  2002/01/16 13:28:15  dturner
	# Initial revision
	#
	# Revision 1.1  1999/03/03 16:24:55  paulv
	# Adapted from pressure_height.pro. Improved input argument checking.
	#

	# Be sure to take input as np.array
	pressure = np.array(pressure)
	mixing_ratio = np.array(mixing_ratio)
	temperature = np.array(temperature)
	# Declare constants
	Rd = 286.9968933
	g   = 9.80616;
	eps = 0.621970585;

	# Calculate average dP
	n_levels=len(pressure) ;
	dp_average = np.sum(pressure[:n_levels-1] - pressure[1:]) / n_levels

	# Sort arrays based on  average pressure differential
	sort = 0
	if ( dp_average < 0.0 ):
	# Data being sorted in ascending(descending) altitude(pressure) order!
		pressure = np.flip(pressure,0)
		temperature = np.flip(temperature,0)
		mixing_ratio = np.flip(mixing_ratio,0)
		sort  = 1

    #                   -- Calculate average temperature --
	t_average = 0.5*(temperature[:n_levels-1]+temperature[1:])

    #         -- Calculate average mixing ratio (in kg/kg - hence --
    #         -- the divisor of 2000.0 instead of 2.0 )           --
	mr_average = 0.0005*(mixing_ratio[:n_levels-1]+mixing_ratio[1:])

#                   -- Calculate virtual temperature --
	ratio = ( 1.0 - eps ) / eps
	t_virtual = t_average * (1.0+( ratio*mr_average ) )

#           -- Calculate altitudes (divide by 1000.0 to get km) --
#           -- Make sure the data is returned in an order       --
#           -- consistent with the input                        --
# -- Calculate layer thicknesses

	altitude = np.array([ surface_alt if i == 0 else \
                      0.001*(Rd/g)*t_virtual[i-1]*np.log( pressure[i-1]/pressure[i]) \
                      for i in range(0,n_levels)]
    )

	# -- Add up layer thicknesses
	# stupid
	for i in range(1,n_levels):
		altitude[i] = altitude[i-1]+altitude[i]
#	altitude = altitude+np.concatenate(([0],altitude[:n_levels-1]))

	# -- Reverse array order if required
	if sort == 1:
		altitude = np.flip(altitude,0)

	return altitude
