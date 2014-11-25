'''
Created on Aug 21, 2013

@author: Max
'''

import numpy as np

def sites(sitechoice):
    
    r_e=6.371e6 #earth radius, m
    g_s=9.807 #standard gravity
    p_0 = 1.01325e5
    R = 8.314 #J/mol/K
    M_air = 28.97e-3 #kg/mol
    
    
    if sitechoice == 'NEEM':     

        elev = 2500. #m
        p_a = 7.45e4  #Pa: kg m^-1 s^-2
        T = -23.5 +273.16 # K
        Accu_0 = 0.25 # m a^-1, ice eq.
        czd = 4. # m
        z_co = 78. # m
        LIZ = 63. # m
        rho0= 360. # kg m^-3 
        hemisphere='NH'
        
    elif sitechoice == 'WAIS':
        elev = 1766. 
        p_a = 7.99e4 
        T = -31 + 273.16 
        Accu_0 = 0.22 
        czd = 4. 
        z_co = 79.5
        LIZ = 67.25
        rho0= 360.
        hemisphere='SH'
        
    elif sitechoice == 'SCENARIO':   

        p_a = 7.45e4
        elev= 2500. 
        T = -24 + 273.15 
        Accu_0 = 0.22
        czd = 1.0 
        z_co = 54.
        LIZ = 51.5
        rho0= 350.
        hemisphere='SCENARIO'
        
    else:  ## something else
    
        pass
    
    g=g_s*(r_e/(r_e+elev))**2
#     p_a=p_0*np.exp(-g_s*M_air*elev/(R*288.15)) #This may not hold at polar sites
    
    return g, p_a, T, Accu_0, czd, z_co, LIZ, rho0, hemisphere