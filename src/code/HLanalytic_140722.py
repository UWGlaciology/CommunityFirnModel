#7/22/14: Max coded this up for Huong to be able to run simulations for Greenland firn age/depth.

import sys
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def rhoHLAnalytic(T,Accu,rho_surf,z_grid): 
        
    R=8.314        
    rho_i = 0.917
    rho_c = 0.55
    h=z_grid
    #rho_bco = 0.815
    
    rho_bco = 1/( 1/(917.) + T*6.95E-7 - 4.3e-5)/1000.

    k0 = 11  * np.exp(-10160/(R*T))
    k1 = 575 * np.exp(-21400/(R*T))
    h0_55 = 1/(rho_i*k0) * (np.log(rho_c/(rho_i-rho_c))-np.log(rho_surf/(rho_i-rho_surf)))
    Z0 = np.exp(rho_i*k0*h + np.log(rho_surf/(rho_i-rho_surf)))
    t0_55 = 1/(k0*Accu) * np.log((rho_i-rho_surf)/(rho_i-rho_c ))
    rho_h0 = (rho_i* Z0)/(1+Z0)
    t0 =   1/(k0*Accu)*np.log((rho_i-rho_surf)/(rho_i-rho_h0))
    Z1 = np.exp(rho_i*k1*(h-h0_55)/np.sqrt(Accu) +  np.log(rho_c/(rho_i-rho_c)))
    Z = np.concatenate((Z0[h<h0_55], Z1[h>h0_55]))

    rho_h = (rho_i * Z)/(1+Z)
    
    tp = 1/(k1*np.sqrt(Accu)) * np.log((rho_i-rho_c)/(rho_i-rho_h))+ t0_55    
    age = np.concatenate((t0[h<h0_55], tp[h>h0_55]))
    bco_age = min(age[rho_h>=rho_bco])
    bco_dep = min(h[rho_h>=rho_bco])
        
    rhoHL=rho_h
    
    return bco_dep,bco_age,rhoHL,age
    
if __name__ == "__main__":
    #
    dz=0.5
    bottom=1000
    z_grid=np.arange(0,bottom+dz,dz)
    
    rho_surf=0.360; #surface density
    Accu=0.0179; #water equivalent accumulation
    T=213.0; #Temperature in K
    
    bco_dep,bco_age,rhoHL,age=rhoHLAnalytic(T,Accu,rho_surf,z_grid)  
    print 'bco depth = ', bco_dep 