import numpy as np



def check_units(TEMP,PRES,*units):


    TEMP_UNITS = ['kelvin','K','C','Celsius','C','celsius']
    PRES_UNITS = ['hPa','mb','Pa']
    if len(units) != 0:
        units = units[0]

    if len(units) == 0:
        # default
        units = ['kelvin','hPa']
    elif units[1] not in TEMP_UNITS or \
         units[0] not in PRES_UNITS:
        raise ValueError('Invalid units: %s' %units)

    if TEMP.shape != PRES.shape :
        raise ValueError('Input arrays must have the same shape: \nlen(TEMP):%f len(PRES):%f'%(len(TEMP),len(PRES)))

    # change units if necessary
    if units[1] in TEMP_UNITS[2:]:   # Celsius
        TEMP = TEMP + 273.15
    if units[0] in PRES_UNITS[2:]:   # Pa
        PRES = PRES * 100

    return TEMP,PRES
#******************************************************************************

def compute_temp(THETA,PRES,*units):

    # Be sure to take input as np.array
    THETA,PRES = np.array(THETA),np.array(PRES)

#    THETA,PRES = check_units(THETA,PRES,units)
    THETA,PRES = check_units(THETA,PRES,units) if len(units)!=0 \
            else check_units(THETA,PRES)

    K = 0.2854

    return THETA * (PRES / 1000)**K
#******************************************************************************

def compute_rh(TEMP,PRES,QVAPOR,*units):


    # Be sure to take input as np.array
    TEMP,PRES = np.array(TEMP),np.array(PRES)

    TEMP,PRES = check_units(TEMP,PRES,units) if len(units)!=0 \
            else check_units(TEMP,PRES)

    e_0 = 6.1173 ## mb
    t_0 = 273.16 ## K
    Rv = 461.50 ## J K-1 Kg-1
    Lv_0 = 2.501 * 10**6 ## J Kg-1
    K1 = Lv_0 / Rv ## K 
    K2 = 1 / t_0 ## K-1
    K3 = 1 / TEMP ## K-1
    ## Clausius Clapeyron Equation
    e_s = e_0 * np.exp( K1 * ( K2 - K3 ) )
    w_s = ( 0.622 * e_s ) / ( PRES - e_s )

    # rechange units...to be added...

    return (QVAPOR / w_s ) *100

#*******************************************************************************

def compute_dewp(TEMP,PRES,QVAPOR,*units):


    # Be sure to take input as np.array
    TEMP,PRES = np.array(TEMP),np.array(PRES)

    TEMP,PRES = check_units(TEMP,PRES,units) if len(units)!=0 \
            else check_units(TEMP,PRES)

    # constants for the Clausius Clapeyron Equation
    e_0 = 6.1173 ## mb
    t_0 = 273.16 ## K
    Rv = 461.50 ## J K-1 Kg-1
    Lv_0 = 2.501 * 10**6 ## J Kg-1
    # Compute portions of the equation
    K1 = Lv_0 / Rv ## K 
    K2 = 1 / t_0 ## K-1
    K3 = 1 / TEMP ## K-1
    # Clausius Clapeyron Equation
    e_s = e_0 * np.exp( K1 * ( K2 - K3 ) ) ## mb
    # get saturation mixing ratio for RH
    w_s = ( 0.622 * e_s ) / ( PRES - e_s )
    rh = ( QVAPOR / w_s ) * 100
    # back out the vapor pressure
    e = ( rh / 100 ) * e_s ## mb

    # compute individual terms when solving the equation for Td
    K1_inv = Rv / Lv_0
    print(e/e_0)
    K4 = np.log( e / e_0 )
    term1 = K2 - (K1_inv * K4)

    # rechange units...to be added...

    return 1 / term1

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
#    altitude = altitude+np.concatenate(([0],altitude[:n_levels-1]))

    # -- Reverse array order if required
    if sort == 1:
        altitude = np.flip(altitude,0)

    return altitude

#************************************************************************#*****************************************<*******************************

def interpolate_pressure_grid(input_grid,pressure,observable,*args):

    # Check order
    flip = 1 if pressure[0]>pressure[1] else 0
    
    if type(input_grid) == np.ma.core.MaskedArray :
        mask = 1
        new_pressure = input_grid[~input_grid.mask]
    else:
        mask = 0
        new_pressure = input_grid

    xp = np.array(new_pressure) if not flip else np.flip(new_pressure,0)
    x = np.array(pressure) if not flip else np.flip(pressure,0)
    y = np.array(observable) if not flip else np.flip(observable,0)
    if len(args) == 0 :
        args = 'default'
    
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
    

    if mask:
        out = np.ma.masked_invalid([np.nan for i in range(0,len(input_grid))])
        out[~input_grid.mask] = yp
        return out
    else:
        return yp

#***************************************************************************************************************<*******************************

def cris_viirs_cloudmask(cris_dataset,viirs_dataset,selected_cris_fovs):
    """
    Bla bla bla
    
    INPUT:      - cris_dataset:   CRIS_Wrapper Class defined in data_reader.cris_wrapper
                - viirs_dataset:   VIIRS_Wrapper Class defined in data_reader.viirs_wrapper
                - selected_cris_fovs: CrIS FOV within bounding box 
    OUTPUT:     - cloud_stats:   Matrix od four dimensional vector withs cloud statistics
                                 (N_CRIS_SCAN_LINES x N_CRIS_FOV_4SCAN_LINE x 4) 
    """
    from preprocessor.fov_generator.cris.fov_generator_cris_utilities.collocation import get_collocation
    from preprocessor.fov_generator.cris.fov_generator_cris_utilities.array_reshapers import transform_index
    indices = get_collocation(cris_dataset, viirs_dataset, selected_cris_fovs)

    cris_shape = cris_dataset.longs.shape
    
    # Define matrix with cloudmask statistics. Each element is a four components vector,
    # each component describe the freq. of the relative cloudmask value.
    cloud_stats = np.zeros( shape = cris_shape + (4,) ,dtype=float)

    for n_row, viirs_indices in enumerate(indices):
        
        # Prepare counter for the n_row-th CRIS FOV
        counter = np.array([0,0,0,0] ,dtype=float)
        
        for index in viirs_indices:
            counter[ viirs_dataset.cloudmask[index] ] += 1
       
        #Paolo Antonelli 18.11.2019 
        if len(viirs_indices)!=0:        
           # Normalize for the number of VIIRS FOVs inside current CRIS FOV
           cloud_stats[transform_index(n_row, cris_shape)  ] = counter / len(viirs_indices)
           #print('transform_index({}, {}): {}, counter: {}'.format(n_row, cris_shape,transform_index(n_row, cris_shape),counter/len(viirs_indices)))
           #print('viirs_indices[:10]: {}'.format(viirs_indices[:10]))
        else:
           cloud_stats[transform_index(n_row, cris_shape)  ] = -1 

    # REMOVE INDICES FROM RETURN (FOR TEST ONLY)
    return [ cloud_stats, indices ] 
