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
