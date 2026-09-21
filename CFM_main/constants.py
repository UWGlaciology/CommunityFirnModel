#!/usr/bin/env python
''' 
Constants used in the CFM.
Units are generally mks.
'''

# gas constant used to calculate Arrhenius term
R           = 8.314                          

# number of seconds in a year
S_PER_YEAR  = 31557600.0                     

# cut off density for the first zone densification (kg/m^3)
RHO_1       = 550.0                          

# cut off density for the second zone densification (kg/m^3)
RHO_2       = 815.0                          

# density of ice (kg/m^3)
RHO_I       = 917.0                          

# density of ice (g/m^3)
RHO_I_MGM   = 0.917                          

# cut off density for the first zone densification (g/m^3)
RHO_1_MGM   = 0.550                          

# acceleration due to gravity on Earth
GRAVITY     = 9.8                            

# conversion from Kelvin to Celsius
K_TO_C      = 273.15                         

# melting temperature
T_MELT      = 273.15

# conversion for accumulation rate
BDOT_TO_A   = S_PER_YEAR * RHO_I_MGM         

# density of water
RHO_W_KGM   = 1000.                         

# specific heat of ice at 0C, kJ kg^-1 K^-1
CP_I_kJ     = 2.097                          

# specific heat of ice at 0C, J kg^-1 K^-1 (Cuffey and Patterson, p.400)
CP_I        = 2097.0                        

# latent heat of ice, kJ kg^-1, (Cuffey and Patterson, p.400)
LF_I_kJ     = 333.5                        

# latent heat of ice, J kg^-1
LF_I        = 333500.0                        

# kg/mol
M_AIR       = 28.97e-3                      

# Standard Amtmospheric Pressure, Pa
P_0         = 1.01325e5

# specific heat of water, J kg^-1 K^-1
CP_W		= 4180.0                

VON_KARMAN = 0.4       # von Karman constant [-]
                       # commonly cited range 0.35-0.42; 0.4 is the
                       # standard rounded value used in most glacier/
                       # atmospheric boundary layer literature

CP_AIR = 1005.0        # specific heat of dry air at constant pressure [J/kg/K]

R_DRY = 287.05         # specific gas constant for dry air [J/kg/K]
                       # (R_DRY = R_universal / M_dry_air)

LV_LIQUID = 2.501e6    # latent heat of vaporization, liquid water -> vapor,
                       # at 0degC [J/kg]
                       # VERIFY: slightly temperature-dependent; 2.501e6
                       # is the standard 0degC reference value. If Ts
                       # varies much above 273.15 in your melt cases,
                       # consider whether a temperature-dependent form
                       # is warranted, though the difference is small
                       # over the relevant range.

LS_SUBLIMATION = 2.834e6  # latent heat of sublimation, ice -> vapor,
                          # at 0degC [J/kg]
                          # VERIFY: this is essentially LV_LIQUID + LF_I
                          # (heat of fusion) - worth checking this is
                          # consistent with whatever LF_I value you
                          # already have defined, rather than treating
                          # these as fully independent constants.