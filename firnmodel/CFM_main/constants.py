#!/usr/bin/env python
''' Constants used in the CFM'''

R           = 8.314                          # gas constant used to calculate Arrhenius term
S_PER_YEAR  = 31557600.0                     # number of seconds in a year
RHO_1       = 550.0                          # cut off density for the first zone densification (kg/m^3)
RHO_2       = 815.0                          # cut off density for the second zone densification (kg/m^3)
RHO_I       = 917.0                          # density of ice (kg/m^3)
RHO_I_MGM   = 0.917                          # density of ice (g/m^3)
RHO_1_MGM   = 0.550                          # cut off density for the first zone densification (g/m^3)
GRAVITY     = 9.8                            # acceleration due to gravity on Earth
K_TO_C      = 273.15                         # conversion from Kelvin to Celsius
T_MELT      = 272.65
BDOT_TO_A   = S_PER_YEAR * RHO_I_MGM         # conversion for accumulation rate
RHO_W_KGM   = 1000.                         # density of water
CP_I_kJ     = 2.097                          # specific heat of ice at 0C, kJ kg^-1 K^-1
CP_I        = 2097.0                        # specific heat of ice at 0C, J kg^-1 K^-1 (Cuffey and Patterson, p.400)
LF_I_kJ     = 333.5                        # latent heat of ice, kJ kg^-1, (Cuffey and Patterson, p.400)
LF_I        = 333500.0                        # latent heat of ice, J kg^-1
M_AIR       = 28.97e-3                      # kg/mol
P_0         = 1.01325e5                     # Standard Amtmospheric Pressure, Pa
