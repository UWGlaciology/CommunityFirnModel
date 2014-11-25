import numpy as np

#Herron and Langway requires water eq. accumulation

def rhoHLAnalytic(R,T,rho_i,rho0,rho_bco,z_nodes,Accu_m,Accu = None):
     
    
    if Accu is None:
        Accu=Accu_m
        
    
    rho_c = 550.0
    h=z_nodes
    
    rho_c_HL = rho_c/1000
    rho_i_HL = rho_i/1000
    rho0_HL = rho0/1000
    rho_bco_HL = rho_bco/1000

    k0 = 11  * np.exp(-10160/(R*T))
    k1 = 575 * np.exp(-21400/(R*T))
    
    h0_55 = 1/(rho_i_HL*k0) * (np.log(rho_c_HL/(rho_i_HL-rho_c_HL))-np.log(rho0_HL/(rho_i_HL-rho0_HL)))

    Z0 = np.exp(rho_i_HL*k0*h + np.log(rho0_HL/(rho_i_HL-rho0_HL)))
    
    t0_55 = 1/(k0*Accu) * np.log((rho_i_HL-rho0_HL)/(rho_i_HL-rho_c_HL ))
    
    rho_h0 = (rho_i_HL* Z0)/(1+Z0)

    t0 =   1/(k0*Accu)*np.log((rho_i_HL-rho0_HL)/(rho_i_HL-rho_h0))

    Z1 = np.exp(rho_i_HL*k1*(h-h0_55)/np.sqrt(Accu) +  np.log(rho_c_HL/(rho_i_HL-rho_c_HL)))
    
    Z = np.concatenate((Z0[h<h0_55], Z1[h>h0_55]))

    rho_h = (rho_i_HL * Z)/(1+Z)
    
    tp = 1/(k1*np.sqrt(Accu)) * np.log((rho_i_HL-rho_c_HL)/(rho_i_HL-rho_h))+ t0_55    
    age = np.concatenate((t0[h<h0_55], tp[h>h0_55]))
    #bco = min(age[rho_h>=rho_bco_HL])
    
    rhoHL=rho_h*1000 #return in kg/m^3
    
    return rhoHL