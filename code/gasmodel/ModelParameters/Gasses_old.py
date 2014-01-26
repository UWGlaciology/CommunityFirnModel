'''
Created on Aug 21, 2013

@author: Max
'''
from numpy import loadtxt
import os

def gasses(gaschoice, sitechoice, T, p_a, DataPath,hemisphere,measurements):
    

    #d_0 = 5.e2 # Free air diffusivity, CO2, m**2/yr Schwander, 1988 reports 7.24 mm**2/s =379 m**2/yr
    d_0 = 1.6e-5 # m^2/s :wikipedia value. changed 9/27/13  Schwander, 1988 reports 7.24 mm**2/s = 7.24e-6 m**2/yr
    M_air = 28.97e-3 #kg/mol
    D_ref_CO2 = 5.75E-10*T**1.81*(101325/p_a) #where is this from?
    
    
    if gaschoice == 'CO2':
        D_gas = 1. #free-air diffusivity relative to CO2. Unitless.
        M = 44.01e-3 # molecular mass, kg/mol 

        
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
        #    conc1=conc1[0:1996,:]
        
            
    elif gaschoice == 'CH4':
        D_gas = 1.367
        M = 16.04e-3

    elif gaschoice == 'd15N2':
        D_gas = 1.275*0.9912227
        M = 1.E-3 + M_air        
                
    
    ### Load gas history. The file must be located in the correct folder, and have the correct naming convention.
    ### If you want to compare to measured samples, make sure that measurements is on.    
    gas_string=gaschoice+'_history_'+hemisphere+'.txt'
    meas_string=gaschoice+'_samples_'+sitechoice+'.txt'
    
    conc1=loadtxt(os.path.join(DataPath,gas_string),skiprows=2) #load data: atmospheric CO2 history.
    
    if measurements=='on':
        firn_meas=loadtxt(os.path.join(DataPath,meas_string),skiprows=2)
    
    
            
    deltaM = (M-M_air) #delta molecular mass from CO2.
    D_x = D_gas #* D_ref_CO2
    d_0=D_ref_CO2
                   
    return D_x, M, deltaM, conc1, firn_meas, d_0