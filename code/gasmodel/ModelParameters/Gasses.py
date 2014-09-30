'''
Created on Aug 21, 2013

@author: Max
'''
from numpy import loadtxt
import os

def gasses(gaschoice, sitechoice, T, p_a, DataPath,hemisphere,loadgas):
    

    #d_0 = 5.e2 # Free air diffusivity, CO2, m**2/yr Schwander, 1988 reports 7.24 mm**2/s =379 m**2/yr
    d_0 = 1.6e-5 # m^2/s :wikipedia value. changed 9/27/13  Schwander, 1988 reports 7.24 mm**2/s = 7.24e-6 m**2/yr
    M_air = 28.97e-3 #kg/mol
    D_ref_CO2 = 5.75E-10*T**1.81*(101325/p_a) #where is this from?
    
    
    if gaschoice == 'CO2':
        D_gas = 1. #free-air diffusivity relative to CO2. Unitless.
        M = 44.01e-3 # molecular mass, kg/mol
        decay = 0.
        omega = 0.0
        
        #if hemisphere == 'SOUTH':
        #    conc1=loadtxt(os.path.join(DataPath,'CO2_NH_history.txt'),skiprows=2) #load data: atmospheric CO2 history.
        #    
        #    firn_meas=loadtxt(os.path.join(DataPath,'CO2_samples_NEEM.txt'),skiprows=2)
        #
        #elif hemisphere == 'NORTH':
        #    conc1=loadtxt(os.path.join(DataPath,'CO2_SH_history.txt'),skiprows=2) #load data: atmospheric CO2 history.
        #    firn_meas=loadtxt(os.path.join(DataPath,'CO2_samples_WAIS.txt'),skiprows=2)   
        #
        #elif hemisphere == 'SCENARIO':
        #    conc1=loadtxt(os.path.join(DataPath,'RampUp2.txt'),skiprows=2) #load data: atmospheric CO2 history.
        #    firn_meas=loadtxt(os.path.join(DataPath,'CO2samples_WAIS.txt'),skiprows=2)   
        #    #conc1=conc1[0:1996,:] # May or may not need this to get time to work...
        
            
    elif gaschoice == 'CH4':
        D_gas = 1.367
        M = 16.04e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'd15N2':
        D_gas = 1.275*0.9912227 # not sure of the origin here... Christo's model?
        M = 1.E-3 + M_air
        decay = 0.
        omega = 0.015/1000

    elif gaschoice == 'SF6':
        D_gas = 0.554
        M = 146.06e-3
        decay = 0.
        omega = 0.
        
    elif gaschoice == 'C14':
        D_gas = 0.991368
        M = 46.01e-3
        decay = 1./8267.
        omega = 0.
        
    elif gaschoice == 'C13':
        D_gas = 0.9955648
        M = 45.01e-3
        decay = 0.
        omega = 0.
        
    elif gaschoice == 'CFC11':
        D_gas = 0.525
        M = 137.37e-3
        decay = 0.

    elif gaschoice == 'CFC12':
        D_gas = 0.596
        M = 120.91e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'C13_CFC12':
        D_gas = 0.59552
        M = 121.91e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'CC14':
        D_gas = 0.470
        M = 153.82e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'CFC113':
        D_gas = 0.453
        M = 187.38e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'CFC115':
        D_gas = 0.532
        M = 154.47e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'R134a':
        D_gas = 0.630
        M = 102.03e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'CH3CCl3':
        D_gas = 0.485
        M = 133.40e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'HCFC22':
        D_gas = 0.710
        M = 86.47e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'C13_CH4':
        D_gas = 1.340806
        M = 17.04e-3
        decay = 0.
        omega = 0.
        
    elif gaschoice == 'd40Ar':
        D_gas = 1.21
        M = 4.e-3 + M_air
        decay = 0.
        omega = 0.0364/1000.

    elif gaschoice == 'FOG':
        D_gas = 1.0
        M = 44e-3
        decay = 1./100.
        omega = 0.
        
        
    
                
    
    ### Load gas history. The file must be located in the correct folder, and have the correct naming convention.
    ### If you want to compare to measured samples, make sure that measurements is on.    
    
    if loadgas:
        gas_string=gaschoice+'_history_'+hemisphere+'.txt'
        meas_string=gaschoice+'_samples_'+sitechoice+'.txt'
    
        conc1=loadtxt(os.path.join(DataPath,gas_string),skiprows=2) #load data: atmospheric CO2 history.
    
#     if measurements=='on':
#         firn_meas=loadtxt(os.path.join(DataPath,meas_string),skiprows=2)
#     else:
#         firn_meas='None'
    
    else:
        conc1=-9999
            
    deltaM = (M-M_air) #delta molecular mass from CO2.
    D_x = D_gas #* D_ref_CO2
    d_0=D_ref_CO2
                   
    return D_x, M, deltaM, conc1, d_0, omega
