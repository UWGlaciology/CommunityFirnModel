
import numpy as np
import os



def diffusivity(rho_co, por_co, por_tot, por_cl, por_op, d_0, rhoprof): #rhoprof is  density profile
        
    #if rhoprof is None:
    #    rhoprof=rhoHL
    
    ## Constants
    d_eddy_sc=d_0 #Eddy diffusivity in the convective zone
    h=z_nodes
     
    ## Use Severinghaus relationship from Cuffey and Paterson
    d_0_sev=d_0*1.49
    diffu_full_Sev = D_x*d_0_sev*((p_0/p_a)*(T/T_0)**1.85*(2.00*(1-(rhoprof/rho_i))-0.167)) 
    diffu_full_Sev[diffu_full_Sev<=0] = 1e-9
    
    ## Use Schwander 1988, Eq. 2 Diffusivity (does not work very well) use 4e2
    ## for d_0
    k_sch = p_0/p_a*(T/253.16)**1.85 # Constant given in Schwander
    diffu_full_sch =3.72*0.5*k_sch*(23.7*por_tot-2.84)*31.5 # Schwander' diffusivity relationship (for CO2). 31.5 is unit conversion. Added extra 3.72* 9/12/13
    ind = np.nonzero(h>LIZ)
    diffu_full_sch[ind] = 0.001
    diffu_full_sch[diffu_full_sch<0] = 1e-15
    
    ## Use Freitag, 2002, Eq 15 Diffusivity use 9e2 for d_0
    d_0_fre=d_0*4.9
    diffu_full_fre = D_x*d_0_fre*por_op**2.1
    diffu_full_fre[diffu_full_fre<=0] = 1e-15
    
    ## Use Christo's diffusivity data from NEEM-EU
    
    diffu_data=np.loadtxt(os.path.join(DataPath,'c_diffu.txt'))
    h=diffu_data[:,0]
    diffu_full=D_x*d_0*diffu_data[:,1]
      
    ## try a random profile to test if the model is working.
    #h=1:100 
    #diffu_full=d_0*23*ones(length(h))
    diffu_full = np.interp(z_nodes,h,diffu_full)
    diffu_full_Christo=diffu_full
    #d_eddy=1. #hack if not using eddy diffusivity (wrap entirely into D)
    ## Add in high diffusivity in convective zone and low diffusivity below LIZ
    
    diffu_full=diffu_full_Christo #change this line to change your choice of diffusivity
    
    d_eddy=np.zeros(np.size(diffu_full))
    ind = np.nonzero(z_nodes<czd)
    d_eddy[ind] = diffu_full[ind]*10
    ind = np.nonzero(z_nodes>LIZ)
    d_eddy[ind] = diffu_full[ind]
    diffu_full[ind]=1e-15
    
    #d_eddy=np.subtract(d_eddy,diffu)
    diffu=diffu_full
    ## Interpolate diffusivity profile to the finite volume nodes (model space)
    deepnodes = z_nodes>LIZ #leftover line from matlab?
    
    return diffu, deepnodes, d_eddy, diffu_full_fre, diffu_full_sch, diffu_full_Sev, diffu_full_Christo