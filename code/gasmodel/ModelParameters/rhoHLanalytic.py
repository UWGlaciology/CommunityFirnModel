import numpy as np


def rhoHLAnalytic(R,T,rho_i,rho0,z_nodes,Accu_0,Accu = None): 
    
    if Accu is None:
        Accu=Accu_0
        
    
    rho_c = 550.0
    h=z_nodes

    k0 = 11  * np.exp(-10160/(R*T))
    k1 = 575 * np.exp(-21400/(R*T))
    
    h0_55 = 1/(rho_i*k0) * (np.log(rho_c/(rho_i-rho_c))-np.log(rho0/(rho_i-rho0)))

    Z0 = np.exp(rho_i*k0*h + np.log(rho0/(rho_i-rho0)))
    
    t0_55 = 1/(k0*Accu) * np.log((rho_i-rho0)/(rho_i-rho_c ))
    
    rho_h0 = (rho_i* Z0)/(1+Z0)

    t0 =   1/(k0*Accu)*np.log((rho_i-rho0)/(rho_i-rho_h0))

    Z1 = np.exp(rho_i*k1*(h-h0_55)/np.sqrt(Accu) +  np.log(rho_c/(rho_i-rho_c)))
    
    Z = np.concatenate((Z0[h<h0_55], Z1[h>h0_55]))

    rho_h = (rho_i * Z)/(1+Z)
    
    tp = 1/(k1*np.sqrt(Accu)) * np.log((rho_i-rho_c)/(rho_i-rho_h))+ t0_55    
    age = np.concatenate((t0[h<h0_55], tp[h>h0_55]))
#     bco = min(age[rho_h>=rho_bco])
    
    rhoHL=rho_h
    
    return rhoHL