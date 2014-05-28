# UW Community Firn-Air Transport Model
# Max Stevens
# version 0.20, 22 May 2014
# Working on gravity and commenting, after conversation with Christo at OSU

### In order to run transient, you must put the files from firnmodel.py output into DataImport folder.
### All units should be kilograms, meters, and seconds (MKS system)

import sys
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.sparse.linalg as splin
from scipy.sparse import spdiags
from scipy.sparse.linalg import lsqr
from scipy.integrate import cumtrapz
import math
import ModelParameters.Gasses as MPG
import ModelParameters.Sites as MPS
import ModelParameters.Plotting as plots
import ModelParameters.Diffusivity as MPD
import ModelParameters.density as MPRHO
import csv

# Set path to find all files to import and set output location
spot = os.path.dirname(sys.argv[0]) #Add Folder
print spot 
# os.chdir(spot) #change pwd to location of firnmodel.py
sys.path.insert(1,spot)
ResultsPlace=os.path.join(spot,'Results')
sys.path.insert(3,os.path.join(spot,'DataImport'))
DataPath = os.path.join(spot,'DataImport')

np.set_printoptions(linewidth=300) #set output reading to be wider. A bit easier to read :)

# Downward Advection (of air)        
def w(z_edges_vec,Accu,rho_interface,por_op,T,p_a,por_tot,por_cl): # Function for downward advection of air. 
    
    por_tot_interface=np.interp(z_edges_vec,z_P_vec,por_tot)
    por_cl_interface=np.interp(z_edges_vec,z_P_vec,por_cl)
    por_op_interface=np.interp(z_edges_vec,z_P_vec,por_op)
    teller_co=np.argmax(por_cl_interface)
    w_ice=Accu*rho_i/rho_interface
    
    if ad_method=='ice_vel':
        w_ad=w_ice
    
    elif ad_method=='Christo':
    ### Christo's Method from his thesis (chapter 5). This (maybe) could be vectorized to speed it up.
    
        bubble_pres = np.zeros_like(z_edges_vec)
        dscl = np.append(0, np.diff(por_cl)/dz)    
        C=np.exp(M_air*g*z_edges_vec/(R*T))
        strain = np.gradient(np.log(w_ice),dz )
        s=por_op_interface+por_cl_interface
        
        for teller1 in range (0,teller_co+1): 
            integral = np.zeros(teller1+1)
            integral2 = np.zeros(teller1+1)
            
            for teller2 in range(0,teller1+1):
                integral[teller2] = dscl[teller2]*C[teller2]*(s[teller2]/s[teller1])/(1+np.trapz(strain[teller2:teller1+1],dx=dz)) #need to get this indexing correct
                if dscl[teller2]==0:
                    dscl[teller2]=1e-14
                integral2[teller2] = dscl[teller2]
                
            bubble_pres[teller1] = (dz*np.sum(integral))/(dz*np.sum(integral2))
    
    
        #for teller1 in range (teller_co+1,np.size(z_edges_vec)):
        #    bubble_pres[teller1] = bubble_pres[teller_co+1]*(s[teller_co+1]/s[teller1])/(w_ice[teller1]/w_ice[teller_co+1])
        
        bubble_pres[teller_co+1:] = bubble_pres[teller_co]*(s[teller_co]/s[teller_co+1:])/(w_ice[teller_co+1:]/w_ice[teller_co])
        bubble_pres[0] = 1
    
        flux= w_ice[teller_co]*bubble_pres[teller_co]*por_cl[teller_co]
        velocity = np.minimum(w_ice ,((flux+(1e-10)-w_ice*bubble_pres*por_cl_interface)/((por_op_interface+1e-10)*C)))
        
        w_ad=velocity
    
    ###

    return w_ad

def A(P): # Power-law scheme, Patankar eq. 5.34
    A = np.maximum( (1 - 0.1 * np.abs( P ) )**5, np.zeros(np.size(P) ) )
    return A    
    
def F_upwind(F): # Upwinding scheme
    F_upwind = np.maximum( F, 0 )
    return F_upwind

def solver(a_U,a_D,a_P,b):
    nz=np.size(b)
    
    Diags = (np.append([a_U,-a_P],[a_D],axis=0))
    cols=np.array([1, 0, -1])
       
    big_A=spdiags(Diags,cols,nz,nz,format='csc')
    big_A=big_A.T
    
    #denseA=big_A.todense() #Dense bit for code testing
    
    #print ('a_P='+str(a_P[0:6])) #code testing. Is output what it should be?
    #print ('a_U='+str(a_U[0:6]))
    #print ('a_D='+str(a_D[0:6]))
    #print ('a_Pend='+str(a_P[-5:]))
    #print ('a_Uend='+str(a_U[-5:]))
    #print ('a_Dend='+str(a_D[-5:]))
    #print denseA[0:5,0:6]
    #print denseA[-5:,-6:]    
    
    rhs=-b
    
    phi_t=splin.spsolve(big_A,rhs)
    
    return phi_t

def FirnAir_SS(z_edges_vec,z_P_vec,nt,dt,Gamma_P,bc_u,bc_d,phi_0,rhoHL,deepnodes,R,nz_P,nz_fv,por_op,gas,Accu,T,p_a,diffu,d_eddy):
    
    ConcPath = os.path.join(ResultsPlace, 'conc_out.csv') #Save location.
    
    #Gamma_P=(diffu*por_op+d_eddy*por_op)
    diffu_P=diffu*por_op
    d_eddy_P=d_eddy*por_op
    ##Gamma_P=(diffu*por_op)
    #Gamma_P[Gamma_P<=0]=0#1e-9 # a bit of a hack here - how to deal with diffusivity after close-off?
    
    #nz_P = np.size(z_P_vec)
    #nz_fv = nz_P - 2
    Z_P = z_P_vec # Some things like this remain in case I want to make 2-D.
     
    dZ = np.concatenate(([1],np.diff(z_edges_vec),[1]))
    
    dZ_u = np.diff(Z_P)
    dZ_u = np.append(dZ_u[0], dZ_u)

    dZ_d = np.diff(Z_P)
    dZ_d = np.append(dZ_d,dZ_d[-1])

    f_u = np.append(0, (1 -(z_P_vec[1:] - z_edges_vec)/dZ_u[1:]))
    f_d = np.append(1 - (z_edges_vec - z_P_vec[0:-1])/dZ_d[0:-1], 0)
     
    diffu_U = np.append(diffu_P[0], diffu_P[0:-1] )
    diffu_D = np.append(diffu_P[1:], diffu_P[-1])
    
    diffu_u =  1/ ( (1 - f_u)/diffu_P + f_u/diffu_U )
    diffu_d =  1/ ( (1 - f_d)/diffu_P + f_d/diffu_D )
    
    d_eddy_U = np.append(d_eddy_P[0], d_eddy_P[0:-1] )
    d_eddy_D = np.append(d_eddy_P[1:], d_eddy_P[-1])
    
    d_eddy_u =  1/ ( (1 - f_u)/d_eddy_P + f_u/d_eddy_U )
    d_eddy_d =  1/ ( (1 - f_d)/d_eddy_P + f_d/d_eddy_D )
    
    #Gamma_m = (Gamma_u + Gamma_d)/2
    
    #S_C=0 #use these lines for ignoring gravity
    #S_C=S_C*np.ones(nz_P)
    
    S_C_0=(-diffu_d+diffu_u)*(deltaM*g/(R*T)) #S_C is independent source term in Patankar
    #S_C_01=-Gamma_d*(deltaM*g/(R*T))
    #S_C_02=Gamma_u*(deltaM*g/(R*T))
    
    #S_C_0=(-d_eddy_d+d_eddy_u)*(deltaM*g/(R*T)) #S_C is independent source term in Patankar
    #S_C_01=-d_eddy_d*(deltaM*g/(R*T))*0
    #S_C_02=d_eddy_u*(deltaM*g/(R*T))*0
    
    #S_C_0=(Gamma_m)*(deltaM*g/(R*T)) #S_C is independent source term in Patankar
    #S_C_0=Gamma_P*(deltaM*g/(R*T))/10
    #S_C_0=0
    S_C=S_C_0*phi_0
    #S_C=S_C_01*phi_0+S_C_02*phi_0
    
    S_P=(-diffu_d+diffu_u)*(deltaM*g/(R*T)) #gravity term, S_P is phi-dependent source
    #S_P=(-d_eddy_d+d_eddy_u)*(deltaM*g/(R*T)) #gravity term, S_P is phi-dependent source
    #S_P=Gamma_m*(deltaM*g/(R*T))
    S_P=1.*S_P
    
    b_0 = S_C*dZ
    
    rho_interface=np.interp(z_edges_vec,z_P_vec,rhoHL)
    
    w_edges = w(z_edges_vec,Accu,rho_interface,por_op,T,p_a,por_tot,por_cl)
    
    w_u = np.append(w_edges[0],  w_edges )
    w_d = np.append(w_edges, w_edges[-1])
    
    D_u = ((diffu_u+d_eddy_u) / dZ_u) #check signs
    D_d = ((diffu_d+d_eddy_d) / dZ_d)
    
    F_u =  w_u*por_op #Is this correct?
    F_d =  w_d*por_op 
    
    P_u = F_u/ D_u
    P_d = F_d/ D_d

    a_U = D_u * A( P_u ) + F_upwind(  F_u ) #+ (Gamma_u)*(deltaM*g/(R*T))
    a_D = D_d * A( P_d ) + F_upwind( -F_d ) #- (Gamma_d)*(deltaM*g/(R*T))#need a negative?

    #a_U = D_u + F_u + (Gamma_u)*(deltaM*g/(R*T))/por_op
    #a_D = D_d - F_d - (Gamma_d)*(deltaM*g/(R*T))/por_op

    a_P_0 = por_op*dZ/dt
    #a_P_0[a_P_0<1e-15]=1.e-15
    
    s=(nz_P,nt)
    phi=np.zeros(s)
    a_P_out=np.zeros(s)
    
    phi_t=phi_0
    
    #S_C=S_C_0*phi_t
    #b_0 = S_C*dZ
    
    with open(ConcPath, "w") as f:
        writer = csv.writer(f)       
        writer.writerow(phi_t)    
        
    a_P =  a_U + a_D + a_P_0 - S_P*dZ
    
    for i_time in range(0,nt): 
        
        a_P_out
        
        bc_u_0=gas[i_time]
        bc_type = 1.
        bc_u   = np.array([ bc_u_0, bc_type])
        
        b = b_0 + a_P_0*phi_t       
        
        #Up boundary
        a_P[0] = 1 
        a_U[0] = 0
        a_D[0] = 0
        b[0]=bc_u[0]
        
        #Down boundary
        a_P[-1] = 1 
        #a_U[-1] = 0
        a_D[-1] = 0
        a_U[-1] = 1
        #b[-1]=dZ_u[-1]*bc_d[0] #old
        b[-1]=-dZ_d[-1]*bc_d[0] # probably does not matter as long as it is zero flux.
        
        phi_t = solver(a_U,a_D,a_P,b)
        
        phi[:,i_time]=phi_t
        a_P_out[:,i_time] = a_P
              
        a_P=a_U+a_D+a_P_0 - S_P*dZ
        
        S_C=S_C_0*phi_t
        #S_C=S_C_01*phi_t+S_C_02*phi_t
        b_0 = S_C*dZ
        
        with open(ConcPath, "a") as f:
        
            writer = csv.writer(f)       
            writer.writerow(phi_t)
        
    return phi, a_P_out

    
def FirnAir_TR(z_edges_vec,z_P_vec,nt,dt,Gamma_P,bc_u,bc_d,phi_0,rhoHL,deepnodes,R,nz_P,nz_fv,por_op,gas,p_a,por_tot,por_cl):
    
    ConcPath = os.path.join(ResultsPlace, 'conc_out.csv') #Save location.
    RhoPath = os.path.join(ResultsPlace, 'rho_out.csv') #Save location.
    DiffuPath = os.path.join(ResultsPlace, 'diffu_out.csv') #Save location.
    ZPath = os.path.join(ResultsPlace, 'Znodes_out.csv') #Save location.
    
    f_depth=np.loadtxt(os.path.join(DataPath,'depthS2.csv'),delimiter=',',skiprows=1)
    #f_depth=f_depth[:,1:]
    f_depth=f_depth[0:1996,1:] #what is 1996?

    f_density=np.loadtxt(os.path.join(DataPath,'densityS2.csv'),delimiter=',',skiprows=1)
    #f_density=f_density[:,1:]
    f_density=f_density[0:1996,1:]
        
    for i_time in range(0,nt):
        
        Accu=Accu_m#+0.01*np.sin(i_time)
        
        #rho_prof=rhoHLAnalytic(Accu)
        
        rho_prof=np.interp(z_nodes,f_depth[i_time,:],f_density[i_time,:]) #density profile
        
        rho_co, por_co, por_tot, por_cl, por_op = porosity(rho_prof) #get porosity 
        diffu, deepnodes, d_eddy, diffu_full_fre, diffu_full_sch, diffu_full_Sev, diffu_full_Christo = diffusivity(rho_co, por_co, por_tot, por_cl, por_op, rhoprof=rho_prof) #get diffusivity profile
        
        Gamma_P=(diffu*por_op+d_eddy*por_op)
        #Gamma_P=(diffu*por_op)
        Gamma_P[Gamma_P<=0]=1e-15 # a bit of a hack here - how to deal with diffusivity after close-off?
    
        diffu_P=diffu*por_op
        d_eddy_P=d_eddy*por_op
    
        #nz_P = np.size(z_P_vec)
        #nz_fv = nz_P - 2
        Z_P = z_P_vec

        dZ = np.concatenate(([1],np.diff(z_edges_vec),[1]))
        
        dZ_u = np.diff(Z_P)
        dZ_u = np.append(dZ_u[0], dZ_u)
    
        dZ_d = np.diff(Z_P)
        dZ_d = np.append(dZ_d,dZ_d[-1])
    
        f_u = np.append(0, (1 -(z_P_vec[1:] - z_edges_vec)/dZ_u[1:]))
        f_d = np.append(1 - (z_edges_vec - z_P_vec[0:-1])/dZ_d[0:-1], 0)
        
        diffu_U = np.append(diffu_P[0], diffu_P[0:-1] )
        diffu_D = np.append(diffu_P[1:], diffu_P[-1])
    
        diffu_u =  1/ ( (1 - f_u)/diffu_P + f_u/diffu_U )
        diffu_d =  1/ ( (1 - f_d)/diffu_P + f_d/diffu_D )
    
        d_eddy_U = np.append(d_eddy_P[0], d_eddy_P[0:-1] )
        d_eddy_D = np.append(d_eddy_P[1:], d_eddy_P[-1])
    
        d_eddy_u =  1/ ( (1 - f_u)/d_eddy_P + f_u/d_eddy_U )
        d_eddy_d =  1/ ( (1 - f_d)/d_eddy_P + f_d/d_eddy_D )
        
        #Gamma_U = np.append(Gamma_P[0], Gamma_P[0:-1] )
        #Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])
        #
        #Gamma_u =  1/ ( (1 - f_u)/Gamma_P + f_u/Gamma_U )
        #Gamma_d =  1/ ( (1 - f_d)/Gamma_P + f_d/Gamma_D )
        
        S_C_0=(-diffu_d+diffu_u)*(deltaM*g/(R*T)) #S_C is independent source term in Patankar
        S_C=S_C_0*phi_0
        #S_C=0
        #S_C=S_C*np.ones(nz_P)
        
        #S_C=(-Gamma_d+Gamma_u)*(deltaM*g/(R*T))
        
        S_P=(-diffu_d+diffu_u)*(deltaM*g/(R*T)) #gravity term, S_P is phi-dependent source
        
        b_0 = S_C*dZ
        
        rho_interface=np.interp(z_edges_vec,z_P_vec,rhoHL)
        w_edges = w(z_edges_vec,Accu,rho_interface,por_op,T,p_a,por_tot,por_cl)
        
        w_u = np.append(w_edges[0],  w_edges )
        w_d = np.append(w_edges, w_edges[-1])
        
        D_u = ((diffu_u+d_eddy_u) / dZ_u) #check signs
        D_d = ((diffu_d+d_eddy_d) / dZ_d)
        
        F_u =  w_u*por_op #Is this correct?
        F_d =  w_d*por_op 
        
        P_u = F_u/ D_u
        P_d = F_d/ D_d
        
        a_U = D_u * A( P_u ) + F_upwind(  F_u )
        a_D = D_d * A( P_d ) + F_upwind( -F_d )
    
        a_P_0 = por_op*dZ/dt

        
        if i_time==0:
            s=(nz_P,nt)
            phi=np.zeros(s)
            diffu_hold=np.zeros(s)
            rho_hold=np.zeros(s)
            phi_t=phi_0
            with open(ConcPath, "w") as f:
                writer = csv.writer(f)       
                writer.writerow(phi_t)
            with open(RhoPath, "w") as f:
                writer = csv.writer(f)       
                writer.writerow(rho_prof)
            with open(DiffuPath, "w") as f:
                writer = csv.writer(f)       
                writer.writerow(diffu)
            with open(ZPath, "w") as f:
                writer = csv.writer(f)       
                writer.writerow(z_nodes)
                                
        a_P =  a_U + a_D + a_P_0
           
        bc_u_0=gas[i_time]
        bc_type = 1.
        bc_u   = np.array([ bc_u_0, bc_type])
        
        b = b_0 + a_P_0*phi_t       
        
        #Upper boundary
        a_P[0] = 1 
        a_U[0] = 0
        a_D[0] = 0
        b[0]=bc_u[0]
        
        #Down boundary
        a_P[-1] = 1 
        #a_U[-1] = 0
        a_D[-1] = 0
        a_U[-1] = 1
        b[-1]=dZ_u[-1]*bc_d[0]

        phi_t = solver(a_U,a_D,a_P,b)
        
        phi[:,i_time] = phi_t
        diffu_hold[:,i_time]=diffu
        rho_hold[:,i_time]=rho_prof
    
        
        with open(ConcPath, "a") as f:
            writer = csv.writer(f)       
            writer.writerow(phi_t)
        with open(RhoPath, "a") as f:
            writer = csv.writer(f)       
            writer.writerow(rho_prof)
        with open(DiffuPath, "a") as f:
            writer = csv.writer(f)       
            writer.writerow(diffu)
        #with open(ZPath, "a") as f:
        #    writer = csv.writer(f)       
        #    writer.writerow(z_nodes)
        
        #plt.figure(13)
        #plt.clf
        #if np.mod(i_time,50)==0:    
        #    
        #    plt.clf
        #    plt.pause(2)
        #    plt.plot(z_nodes,phi_t,'b')
        #    plt.show()    
            
        a_P=a_U+a_D+a_P_0
        
    return phi, diffu_hold, rho_hold
        
        
def space(depth,z_res):
    Lz=depth
    nz_fv=np.around(depth/z_res)
    nz_P=nz_fv+2
    
    dz=Lz/nz_fv
    z_edges_vec=dz*np.arange(0,nz_fv+1)
        
    z_P_vec=np.concatenate(([z_edges_vec[0]],z_edges_vec[0:-1]+np.diff(z_edges_vec)/2,[z_edges_vec[-1]]))
    
    z_nodes = z_P_vec
    nodes = np.size(z_nodes)
    
    return dz, z_edges_vec, z_P_vec, z_nodes, nodes, nz_P, nz_fv
    
#def rhoHLAnalytic(Accu = None): 
#    
#    if Accu is None:
#        Accu=Accu_0
#        
#    
#    rho_c = 550.0
#    h=z_nodes
#
#    k0 = 11  * np.exp(-10160/(R*T))
#    k1 = 575 * np.exp(-21400/(R*T))
#    
#    h0_55 = 1/(rho_i*k0) * (np.log(rho_c/(rho_i-rho_c))-np.log(rho0/(rho_i-rho0)))
#
#    Z0 = np.exp(rho_i*k0*h + np.log(rho0/(rho_i-rho0)))
#    
#    t0_55 = 1/(k0*Accu) * np.log((rho_i-rho0)/(rho_i-rho_c ))
#    
#    rho_h0 = (rho_i* Z0)/(1+Z0)
#
#    t0 =   1/(k0*Accu)*np.log((rho_i-rho0)/(rho_i-rho_h0))
#
#    Z1 = np.exp(rho_i*k1*(h-h0_55)/np.sqrt(Accu) +  np.log(rho_c/(rho_i-rho_c)))
#    
#    Z = np.concatenate((Z0[h<h0_55], Z1[h>h0_55]))
#
#    rho_h = (rho_i * Z)/(1+Z)
#    
#    tp = 1/(k1*np.sqrt(Accu)) * np.log((rho_i-rho_c)/(rho_i-rho_h))+ t0_55    
#    age = np.concatenate((t0[h<h0_55], tp[h>h0_55]))
#    bco = min(age[rho_h>=rho_bco])
#    
#    rhoHL=rho_h
#    
#    return rhoHL

def porosity(rho_prof):
    
    ## Porosity, from Goujon et al., 2003, equations 9 and 10
    por_tot = 1-rho_prof/rho_i # Total porosity
    rho_co = rho_bco #use Martinerie close-off criteria
    #rho_co = 0.815 # User chosen close off-density (put in site-specific in sites?)
    por_co = 1 - rho_co/rho_i # Porosity at close-off
    alpha = 0.37 # constant determined in Goujon
    por_cl = alpha*por_tot*(por_tot/por_co)**(-7.6)
    ind=por_cl>por_tot
    por_cl[ind]=por_tot[ind]
    
    #por_cl[por_cl>1]=1
    #por_cl = por_cl*por_tot # Final closed porosity
    por_op = por_tot - por_cl # Open Porosity
    por_op[por_op<=0] = 1.e-10
    
    return rho_co, por_co, por_tot, por_cl, por_op
      
def diffusivity(rho_co, por_co, por_tot, por_cl, por_op, rhoprof = None): #rhoprof is density profile
        
    if rhoprof is None:
        rhoprof=rhoHL
    
    ## Constants
    d_eddy_sc=d_0 #Eddy diffusivity in the convective zone
    h=z_nodes
     
    ## Use Severinghaus relationship from Cuffey and Paterson
    #d_0_sev=d_0*1.49
    d_0_sev=d_0*0.1
    diffu_full_Sev = D_x*d_0_sev*((p_0/p_a)*(T/T_0)**1.85*(2.00*(1-(rhoprof/rho_i))-0.167)) 
    diffu_full_Sev[diffu_full_Sev<=0] = 1e-9
    
    ## Use Schwander 1988, Eq. 2 Diffusivity (does not work very well) use 4e2
    ## for d_0
    k_sch = p_0/p_a*(T/253.16)**1.85 # Constant given in Schwander
    diffu_full_sch =3.72*0.5*k_sch*(23.7*por_tot-2.84)*31.5 # Schwander' diffusivity relationship (for CO2). 31.5 is unit conversion. Added extra 3.72* 9/12/13
    diffu_full_sch =k_sch*(23.7*por_tot-2.84)/(1000**2) # Schwander' diffusivity relationship (for CO2). 31.5 is unit conversion. Added extra 3.72* 9/12/13
    #ind = np.nonzero(h>LIZ)
    #diffu_full_sch[ind] = 0.001
    diffu_full_sch[diffu_full_sch<0] = 1.e-15
    
    ## Use Freitag, 2002, Eq 15 Diffusivity use 9e2 for d_0
    d_0_fre=d_0*4.9
    diffu_full_fre = D_x*d_0_fre*por_op**2.1
    diffu_full_fre[diffu_full_fre<=0] = 1e-15
    
    ## Use Christo's diffusivity data from NEEM-EU
    if sitechoice == 'NEEM':
        diffu_data=np.loadtxt(os.path.join(DataPath,'c_diffu_NEEM.txt'))
        h=diffu_data[:,0]
        diffu_full=D_x*d_0*diffu_data[:,1]

    elif sitechoice == 'WAIS':
        diffu_data=np.loadtxt(os.path.join(DataPath,'c_diffu_WAIS.txt'))
        h=diffu_data[:,0]
        diffu_full=D_x*d_0*diffu_data[:,1]
        
    elif sitechoice == 'SCENARIO':
        diffu_data=np.loadtxt(os.path.join(DataPath,'c_diffu_NEEM.txt'))
        h=diffu_data[:,0]
        diffu_full=D_x*d_0*diffu_data[:,1]        
        
      
    ## try a random profile to test if the model is working.
    #h=1:100 
    #diffu_full=d_0*23*ones(length(h))
    diffu_full = np.interp(z_nodes,h,diffu_full)
    diffu_full_Christo=diffu_full
    #d_eddy=1. #hack if not using eddy diffusivity (wrap entirely into D)
    
    ## Add in high diffusivity in convective zone and low diffusivity below LIZ
    
    diffu_full=diffu_full_Sev #change this line to change your choice of diffusivity
    
    #Add eddy diffusivity terms: convective zone and non-diffusive zone
    d_eddy=np.zeros(np.size(diffu_full))
    ind = np.nonzero(z_nodes<czd)
    d_eddy_surf=2.426405E-5 #Kawamura, 2006
    H_scale=4.4
    d_eddy_up=d_eddy_surf*np.exp(-1*z_nodes/H_scale)
    #d_eddy[ind] = diffu_full[ind]*10
    ind = np.flatnonzero(z_nodes>LIZ)
    ind2 = np.flatnonzero(z_nodes<z_co)
    ind3 = np.intersect1d(ind,ind2)
    d_eddy[ind3] = diffu_full[ind]/10
    d_eddy=d_eddy+d_eddy_up
    diffu_full[ind]=1e-15 #set molecular diffusivity equal to zero after LIZ - eddy diffusivity term drives diffusion below
    
    #diffu_full[:]=1.e-10 #Use these lines for advection testing.
    #d_eddy[:]=1.e-10
    
    #d_eddy=np.subtract(d_eddy,diffu)
    diffu=diffu_full
    ## Interpolate diffusivity profile to the finite volume nodes (model space)
    deepnodes = z_nodes>LIZ #leftover line from matlab?
    
    return diffu, deepnodes, d_eddy, diffu_full_fre, diffu_full_sch, diffu_full_Sev, diffu_full_Christo
    
def boundaries(gas_org):
    bc_u_0 = gas_org[0] #this is concentration on upper boundary (i.e. atmosphere) It is not actually necessary here, because it gets looped over.
    bc_type = 1
    bc_u   = np.concatenate( ([bc_u_0], [bc_type]))

    bc_d_0 = 0 
    bc_type = 2
    bc_d   = np.concatenate(([ bc_d_0 ], [ bc_type ]))
     
    return bc_u,bc_d, bc_u_0

def gasinterp(time_yr,gas_org,model_time):
    
    gas=np.interp(model_time,time_yr,gas_org)
    
    return gas
    
        
if __name__ == "__main__":
    #globals
    
    reload(MPG)
    reload(MPS)
    reload(MPD)
    reload(plots)
    reload(MPRHO)
    import time
    
    tic=time.time()
    
    rho_i = 917.0 #kg/m^3
    R = 8.314 #J/mol/K
    M_air = 28.97e-3 #kg/mol
    rho_bco = 815. #kg/m^3
    p_0 = 1.01325e5 # Standard Amtmospheric Pressure, Pa
    T_0 = 273.15 # Standard Temp, K
    sPerYear = 365.25*24*3600 #seconds per year
    depth = 120. # m
    
    ad_method="Christo" #advection method 
    measurements = 'on'
    
    
    # Set up parameters for different sites.
    sitechoice = 'NEEM'
    g, p_a, T, Accu_0, czd, z_co, LIZ, rho0, hemisphere = MPS.sites(sitechoice)   
    Accu_m=Accu_0 #Accumulation in m/year
    Accu_0=Accu_0/sPerYear #accumulation in m/second

    
    # Set up Gas Properties. Just CO2 for now.
    #gaschoice='CO2'
    #gaschoice='CH4'
    #gaschoice='SF6'    
    gaschoice='d15N2'
    D_x, M, deltaM, conc1, firn_meas, d_0 = MPG.gasses(gaschoice, sitechoice,T,p_a,DataPath,hemisphere,measurements)

    time_yr=conc1[:,0] # Atmospheric measurements times
    
    #timestop=400.0
    #iin=np.flatnonzero(time_yr==timestop)
    #time_yr=time_yr[0:iin+1] 
       
    time_yr_s=time_yr*sPerYear
    gas_org=conc1[:,1] # Atmospheric measurement concentrations 
    #gas_org=conc1[0:iin+1,1] #Added 4/2/14 for gala runs. Makes gas_org the same length as the time.
    if firn_meas != 'None':
        meas_depth=firn_meas[:,0]
        meas_conc=firn_meas[:,1] # measured gas concentrations in firn
        #meas_conc=firn_meas[:,2] # gravity corrected
        meas_uncert=firn_meas[:,3]
    
    #Space and time. Time is in seconds!
    z_res=0.5 #resolution of grid, m
    #dt=0.01
    
    yrs=np.around(time_yr[-1]-time_yr[0]) #should I/can I round here? 9/10

    time_total=yrs*sPerYear
    stpsperyear=5 #If this is for transient this number must (for now) be the same time steps as the input density/depth files.s
    t_steps=np.int(yrs*stpsperyear)
    #dt=0.2 #time step size.
    dt=time_total/t_steps #time step size. 
    model_time=np.arange(np.around(time_yr[0]*sPerYear),np.around(time_yr[-1]*sPerYear),dt) #set model time steps
    nt=np.size(model_time) #number of time steps
    
    dz, z_edges_vec, z_P_vec, z_nodes, nodes, nz_P, nz_fv = space(depth,z_res) #call on space function to set up spatial grid
    
    rhoHL = MPRHO.rhoHLAnalytic(R,T,rho_i,rho0,rho_bco,z_nodes,Accu_m) # Get density profile from H&L analytic
    rho_co, por_co, por_tot, por_cl, por_op = porosity(rhoHL)
    diffu, deepnodes, d_eddy, diffu_full_fre, diffu_full_sch, diffu_full_Sev, diffu_full_Christo = diffusivity(rho_co, por_co, por_tot, por_cl, por_op) #get diffusivity profiles
  
    ind_co=np.argmax(z_edges_vec>=z_co) #index of close-off depth in z_edges vec. A bit of a hack for now...
    ind_co=np.argmax(rhoHL>=rho_co) #index of close-off depth in z_edges vec. A bit of a hack for now...    
    #gas = gasinterp(time_yr_s,gas_org,model_time) #interpolate atmospheric gas history to model time.
    gas=np.interp(model_time,time_yr_s,gas_org) #interpolate atmospheric gas history to model time.
    bc_u, bc_d, bc_u_0 = boundaries(gas_org) #set boundary and initial conditions: bc_u is atmospheric condition, bc_d is zero gradient.
    
    Z_P = z_P_vec
     
    dZ = np.concatenate(([1],np.diff(z_edges_vec),[1]))
    
    dZ_u = np.diff(Z_P)
    dZ_u = np.append(dZ_u[0], dZ_u)

    dZ_d = np.diff(Z_P)
    dZ_d = np.append(dZ_d,dZ_d[-1])
    
    phi_0 = np.zeros(nz_P)
    phi_0[:]=bc_u_0
    
    Gamma_P=(diffu*por_op+d_eddy*por_op)
    #Gamma_P=(diffu*por_op)
    Gamma_P[Gamma_P<=0]=0#1e-9 # a bit of a hack here - how to deal with diffusivity after close-off?
    
    rho_interface=np.interp(z_edges_vec,z_P_vec,rhoHL)
    
    transdiffu = 'off'
    
    if transdiffu == 'on':
        phi, diffu_hold, rho_hold = FirnAir_TR(z_edges_vec,z_P_vec,nt,dt,Gamma_P,bc_u,bc_d,phi_0,rhoHL,deepnodes, R,nz_P,nz_fv,por_op,gas,p_a,por_tot,por_cl)
    
    elif transdiffu == 'off':
        phi, a_P_out =FirnAir_SS(z_edges_vec,z_P_vec,nt,dt,Gamma_P,bc_u,bc_d,phi_0,rhoHL,deepnodes, R,nz_P,nz_fv,por_op,gas,Accu_0,T,p_a,diffu,d_eddy)

    else:
        sys.exit("Oops. You have a typo at transdiffu. Try again.")
        
    #phi_final=np.append([phi_0],[phi],axis=0)
    #phi=phi.T
    
    phi_toplot=phi[:,1:-1:20]
    aa=np.shape(phi_toplot)
    aa=aa[1]
    
    if firn_meas != 'None':
        gas_meas=np.interp(z_nodes, meas_depth, meas_conc);
    
    plotting = 'on'
    
    plots.makeplots(plotting,Z_P,phi,gas_meas,meas_depth,meas_conc,ResultsPlace,
                    diffu_full_Sev,diffu_full_fre,diffu_full_sch,diffu_full_Christo, meas_uncert=meas_uncert)
        
    elapsed=time.time()-tic
    elapsed_min=elapsed/60.
    mins=np.floor(elapsed_min)
    secs=(elapsed_min-mins)*60
    print mins, 'min', secs, 'sec elapsed'
        