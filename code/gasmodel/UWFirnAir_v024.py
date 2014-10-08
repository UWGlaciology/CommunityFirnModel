# UW Community Firn-Air Transport Model
# Max Stevens
# version 0.24, 20 August 2014
# Working on getting transient to interact smoothly with input from firnmodel.py 

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
import json

# Set path to find all files to import and set output location
spot = os.path.dirname(sys.argv[0]) #Add Folder
print spot 
# os.chdir(spot) #change pwd to location of firnmodel.py
sys.path.insert(1,spot)
ResultsPlace=os.path.join(spot,'Results')
sys.path.insert(3,os.path.join(spot,'DataImport'))
DataPath = os.path.join(spot,'DataImport')
DataPathUser = os.path.join(spot,'DataImport/user_input')

np.set_printoptions(linewidth=300) #set output reading to be wider. A bit easier to read

#Set Globals
rho_i = 917.0 #kg/m^3
R = 8.314 #J/mol/K
M_air = 28.97e-3 #kg/mol
rho_bco = 815. #kg/m^3
p_0 = 1.01325e5 # Standard Amtmospheric Pressure, Pa
T_0 = 273.15 # Standard Temp, K
sPerYear = 365.25*24*3600 #seconds per year
g=9.8

# Downward Advection (of air)        
def w(z_edges,Accu,rho_interface,por_op,T,p_a,por_tot,por_cl,z_nodes,ad_method,dz): # Function for downward advection of air and also calculates total air content. 
    
    global teller_co, por_cl_interface
    
    por_tot_interface=np.interp(z_edges,z_nodes,por_tot)
    por_cl_interface=np.interp(z_edges,z_nodes,por_cl)
    por_op_interface=np.interp(z_edges,z_nodes,por_op)
    teller_co=np.argmax(por_cl_interface)
    w_ice=Accu*rho_i/rho_interface #Check this - is there a better way?
    
    if ad_method=='ice_vel':
        w_ad=w_ice
        trapped = 0.0
        bubble_pres = np.zeros_like(z_edges)
    
    elif ad_method=='Christo':
    ### Christo's Method from his thesis (chapter 5). This (maybe) could be vectorized to speed it up.
    
        bubble_pres = np.zeros_like(z_edges)
        dscl = np.append(0, np.diff(por_cl)/dz)    
        C=np.exp(M_air*g*z_edges/(R*T))
        strain = np.gradient(np.log(w_ice),dz )
        s=por_op_interface+por_cl_interface
        
        for teller1 in range (0,teller_co+1): 
            integral = np.zeros(teller1+1)
            integral2 = np.zeros(teller1+1)
            
            for teller2 in range(0,teller1+1):
                integral[teller2] = dscl[teller2]*C[teller2]*(s[teller2]/s[teller1])/(1+np.trapz(strain[teller2:teller1+1],dx=dz)) #need to get this indexing correct 6/19/14: I think it is fine.
                if dscl[teller2]==0:
                    dscl[teller2]=1e-14
                integral2[teller2] = dscl[teller2]
                
            bubble_pres[teller1] = (dz*np.sum(integral))/(dz*np.sum(integral2))
        
        bubble_pres[teller_co+1:] = bubble_pres[teller_co]*(s[teller_co]/s[teller_co+1:])/(w_ice[teller_co+1:]/w_ice[teller_co])
        
        bubble_pres[0] = 1
        #print 'bubble pressure = %s' % bubble_pres
        
        flux= w_ice[teller_co]*bubble_pres[teller_co]*por_cl[teller_co]
        velocity = np.minimum(w_ice ,((flux+(1e-10)-w_ice*bubble_pres*por_cl_interface)/((por_op_interface+1e-10)*C)))
            
        w_ad=velocity

    return w_ad, bubble_pres

def A(P): # Power-law scheme, Patankar eq. 5.34
    A = np.maximum( (1 - 0.1 * np.abs( P ) )**5, np.zeros(np.size(P) ) )
    return A    
    
def F_upwind(F): # Upwinding scheme
    F_upwind = np.maximum( F, 0 )
    return F_upwind

def solver(a_U,a_D,a_P,b): #routine to solve Ax=b
    nz=np.size(b)
    Diags = (np.append([a_U,-a_P],[a_D],axis=0))
    cols=np.array([1, 0, -1])      
    big_A=spdiags(Diags,cols,nz,nz,format='csc')
    big_A=big_A.T        
    rhs=-b    
    phi_t=splin.spsolve(big_A,rhs)    
    return phi_t

def FirnAir_SS(cc,gaschoice):
        
    depth = cc["depth"] # m
    ad_method = cc["ad_method"] #advection method 
    
    # Set up parameters for different sites.
    sitechoice = cc["sitechoice"]
    g, p_a, T, Accu_0, czd, z_co, LIZ, rho0, hemisphere = MPS.sites(sitechoice)   
    Accu_m=Accu_0 #Accumulation in m/year 
    Accu_0=Accu_0/sPerYear #accumulation in m/second
    
    loadgas = True        
    D_x, M, deltaM, conc1, d_0, omega = MPG.gasses(gaschoice, sitechoice,T,p_a,DataPath,hemisphere,loadgas)
    
    time_yr=conc1[:,0] # Atmospheric measurements times
                
    time_yr_s=time_yr*sPerYear
    
    gas_org=conc1[:,1] # Atmospheric measurement concentrations
      
    #Space and time. Time is in seconds!
    z_res=cc['z_resolution'] #resolution of grid, m
    
    yrs=np.around(time_yr[-1]-time_yr[0])

    time_total=yrs*sPerYear #total model run time in seconds
    stpsperyear=cc['StepsPerYear'] #If this is for transient this number must (for now) be the same time steps as the input density/depth files. Make sure that it is a float.
    t_steps=np.int(yrs*stpsperyear)
    dt=time_total/t_steps #time step size. 
    model_time=np.arange(np.around(time_yr[0]*sPerYear),np.around(time_yr[-1]*sPerYear),dt) #set model time steps
    model_time_years=model_time/sPerYear
    nt=np.size(model_time) #number of time steps
    
    dz, z_edges, z_nodes, nodes, nz_P, nz_fv = space(depth,z_res) #call on space function to set up spatial grid
    
    rhoHL = MPRHO.rhoHLAnalytic(R,T,rho_i,rho0,rho_bco,z_nodes,Accu_m) # Get density profile from H&L analytic
    rho_co, por_co, por_tot, por_cl, por_op, bcoRho, LIDRho = porosity(rhoHL,T)
    
    if sitechoice=='SCENARIO': 
        z_co = min(z_nodes[rhoHL>=(bcoRho)]) #close-off depth; bcoRho is close off density
        LIZ = min(z_nodes[rhoHL>=(LIDRho)]) #lock in depth; LIDRho is lock-in density
    
    diffu,  d_eddy, diffu_full_fre, diffu_full_sch, diffu_full_sev, diffu_full_data = diffusivity(cc,rho_co, por_co, por_tot, por_cl, por_op, z_co, czd, LIZ,d_0,D_x,p_a,z_nodes,T,sitechoice, rhoHL) #get diffusivity profiles
    
    dcon=1.0 
    diffu=diffu*dcon   
    gas=np.interp(model_time,time_yr_s,gas_org) #interpolate atmospheric gas history to model time.
    bc_u, bc_d, bc_u_0 = boundaries(gas_org) #set boundary and initial conditions: bc_u is atmospheric condition, bc_d is zero gradient.
    
    phi_0 = np.zeros(nz_P) #phi is mixing ratio of gas.
    phi_0[:]=bc_u_0
    
    ConcPath = os.path.join(ResultsPlace, 'concentration.csv') #Save location.
    DepthPath = os.path.join(ResultsPlace, 'depth.csv') #Save location.
    
    diffu_P=diffu*por_op
    d_eddy_P=d_eddy*por_op
     
    dZ = np.concatenate(([1],np.diff(z_edges),[1])) #Not sure why there is the 1 on either side... from Ed's Patankar code.
    
    dZ_u = np.diff(z_nodes)
    dZ_u = np.append(dZ_u[0], dZ_u)

    dZ_d = np.diff(z_nodes)
    dZ_d = np.append(dZ_d,dZ_d[-1])

    f_u = np.append(0, (1 -(z_nodes[1:] - z_edges)/dZ_u[1:]))
    f_d = np.append(1 - (z_edges - z_nodes[0:-1])/dZ_d[0:-1], 0)
     
    diffu_U = np.append(diffu_P[0], diffu_P[0:-1] )
    diffu_D = np.append(diffu_P[1:], diffu_P[-1])
    
    diffu_u =  1/ ( (1 - f_u)/diffu_P + f_u/diffu_U )
    diffu_d =  1/ ( (1 - f_d)/diffu_P + f_d/diffu_D )
    
    d_eddy_U = np.append(d_eddy_P[0], d_eddy_P[0:-1] )
    d_eddy_D = np.append(d_eddy_P[1:], d_eddy_P[-1])
    
    d_eddy_u =  1/ ( (1 - f_u)/d_eddy_P + f_u/d_eddy_U )
    d_eddy_d =  1/ ( (1 - f_d)/d_eddy_P + f_d/d_eddy_D )
    
    S_P=0.0 # Source dependent term. Should always be zero, but leaving it in in case there becomes a need.
    #S_P=(-diffu_d+diffu_u)*(deltaM*g/(R*T)) #OLD gravity term, S_P is phi-dependent source
    S_P=1.*S_P
    
    if cc['gravity']=="off" and cc['thermal']=="off":
        print 'gravity and thermal are off'
        S_C_0=0.0
    
    elif cc['gravity']=='on' and cc['thermal']=='off':
        S_C_0=(-diffu_d+diffu_u)*(deltaM*g/(R*T))/dz #S_C is independent source term in Patankar
    
    elif cc['gravity']=='on' and cc['thermal']=='on':
        print 'thermal on'
        dTdz=np.zeros(np.size(diffu_d))
        dTdz[0:100]=-0.6 #K/m. Negative gradient here means that it is colder deep (z is positive down)
        S_C_0=(diffu_d-diffu_u)*((-deltaM*g/(R*T))+(omega*dTdz))/dz #S_C is independent source term in Patankar
    
    S_C=S_C_0*phi_0
    
    b_0 = S_C*dZ
    
    rho_interface=np.interp(z_edges,z_nodes,rhoHL) #Density at finite-volume interfaces
    
    w_edges, bubble_pres = w(z_edges,Accu_0,rho_interface,por_op,T,p_a,por_tot,por_cl,z_nodes,ad_method,dz)
    
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
    
    s=(nz_P,nt)
    phi=np.zeros(s)
    a_P_out=np.zeros(s)
    
    phi_t=phi_0
    
    with open(ConcPath, "w") as f:
        writer = csv.writer(f)       
        writer.writerow(np.append(0,phi_t))    
        
    with open(DepthPath, "w") as f:
        writer = csv.writer(f)       
        writer.writerow(np.append(0,z_nodes)) 
            
    a_P =  a_U + a_D + a_P_0 - S_P*dZ
    
    for i_time in range(0,nt): 
        
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
        a_D[-1] = 0
        a_U[-1] = 1
        b[-1]=-dZ_d[-1]*bc_d[0] # probably does not matter as long as it is zero flux.
        
        phi_t = solver(a_U,a_D,a_P,b)
        
        phi[:,i_time]=phi_t
        a_P_out[:,i_time] = a_P
              
        a_P=a_U+a_D+a_P_0 - S_P*dZ
        
        S_C=S_C_0*phi_t
        b_0 = S_C*dZ
        
        with open(ConcPath, "a") as f:        
            writer = csv.writer(f)       
            writer.writerow(np.append(model_time_years[i_time],phi_t))
        
    return phi, a_P_out, bubble_pres, z_nodes

    
def FirnAir_TR(cc,gaschoice,jj):
    
    g=9.8
    FID='conc_out_%s.csv' %gaschoice
    
    ConcPath = os.path.join(ResultsPlace, FID) #Save location.
    RhoPath = os.path.join(ResultsPlace, 'rho_out.csv') #Save location.
    DiffuPath = os.path.join(ResultsPlace, 'diffu_out.csv') #Save location.
    ZPath = os.path.join(ResultsPlace, 'Znodes_out.csv') #Save location.
    
    if (cc['UserData']):
    
        f_depth=np.loadtxt(os.path.join(DataPathUser,'depth.csv'),delimiter=',',skiprows=0) #load data from firnmodel.py. This is what should be streamlined.    
        f_density=np.loadtxt(os.path.join(DataPathUser,'density.csv'),delimiter=',',skiprows=0)
        f_dcon=np.loadtxt(os.path.join(DataPathUser,'Dcon.csv'),delimiter=',',skiprows=0)
        #f_dcon[f_dcon==0.9]=0.2        
        f_temp=np.loadtxt(os.path.join(DataPathUser,'temp.csv'),delimiter=',',skiprows=0)        
        f_clim=np.loadtxt(os.path.join(DataPathUser,'Clim.csv'),delimiter=',',skiprows=0)
        f_BCO=np.loadtxt(os.path.join(DataPathUser,'BCO.csv'),delimiter=',',skiprows=0)
        f_LID=np.loadtxt(os.path.join(DataPathUser,'LID.csv'),delimiter=',',skiprows=0)
        f_gas=np.loadtxt(os.path.join(DataPathUser,'GasHistory.csv'),delimiter=',',skiprows=0)
        
        print "Gasses Loaded"
        
        Accu_vec=f_clim[:,1]
        
        T_vec=f_clim[:,2]
        #Accu_vec=np.ones(len(f_density[:,0])) #nned to be aware of accumulation per second or year. Accu_m is per year. Accu_0 is per s.
       
        depth = cc["depth"] # m
        ad_method = cc["ad_method"] #advection method
        p_a=cc['pressure'] #pressure at site. Could be vectorized to change through time.
        czd=cc['ConZoneDepth'] #could also be a vector
        z_co_vec=f_BCO[2,:] #2 for mart, 4 for 815
        LIZ_vec = f_LID[2,:]
        rho0=cc['rho0']
        hemisphere='SCENARIO'
        loadgas = False
        sitechoice=cc['sitechoice']
        T_DX=T_vec[0]
        D_x, M, deltaM, conc1, d_0, omega = MPG.gasses(gaschoice, sitechoice,T_DX,p_a,DataPath,hemisphere,loadgas)
        
        time_yr=f_clim[:,0]
        time_yr_s=time_yr*sPerYear
        gas_org=f_gas[jj+1,:] #this needs to match up with the order of gasses specified in the config.
        z_res=cc['z_resolution']
        yrs=time_yr[-1]-time_yr[0]
        time_total=yrs*sPerYear #total model run time in seconds
        stpsperyear=1/(time_yr[2]-time_yr[1])
        #t_steps=np.int(yrs*stpsperyear)
        t_steps=yrs*stpsperyear
        dt=time_total/t_steps #time step size.
        print dt
        #model_time=np.arange(np.around(time_yr[0]*sPerYear),np.around(time_yr[-1]*sPerYear),dt) #set model time steps
        model_time=np.arange(time_yr[0]*sPerYear,time_yr[-1]*sPerYear+dt,dt) #set model time steps
        
        model_time_years=model_time/sPerYear
        nt=np.size(model_time) #number of time steps

        if nt>len(time_yr): #this is because occasionally np.arrange will include end points.
            model_time=np.arange(time_yr[0]*sPerYear,time_yr[-1]*sPerYear,dt)
            model_time_years=model_time/sPerYear
            nt=np.size(model_time) #number of time steps        
                        
        dz, z_edges, z_nodes, nodes, nz_P, nz_fv = space(depth,z_res)
        
        bc_u, bc_d, bc_u_0 = boundaries(gas_org) #set boundary and initial conditions: bc_u is atmospheric condition, bc_d is zero gradient.        
        phi_0 = np.zeros(nz_P) #phi is mixing ratio of gas.
        phi_0[:]=bc_u_0
        gas=np.interp(model_time,time_yr_s,gas_org) #interpolate atmospheric gas history to model time.
        Accu_vec=np.interp(model_time,time_yr,Accu_vec)
        #Accu_vec_s=Accu_vec/sPerYear
        T_vec=np.interp(model_time,time_yr,T_vec)
        #rho_init=f_
        #rho_co, por_co, por_tot, por_cl, por_op, bcoRho, LIDRho = porosity(rhoHL,T)
    
    #######
    else:    
        depth = cc["depth"] # m
        ad_method = cc["ad_method"] #advection method 
    
    # Set up parameters for different sites.
        sitechoice = cc["sitechoice"]
        g, p_a, T, Accu_0, czd, z_co, LIZ, rho0, hemisphere = MPS.sites(sitechoice)   
        Accu_m=Accu_0 #Accumulation in m/year 
        Accu_0=Accu_0/sPerYear #accumulation in m/second
        
        loadgas = True
        
        D_x, M, deltaM, conc1, d_0, omega = MPG.gasses(gaschoice, sitechoice,T,p_a,DataPath,hemisphere,loadgas)
    
        time_yr=conc1[:,0] # Atmospheric measurements times
                
        time_yr_s=time_yr*sPerYear
    
        gas_org=conc1[:,1] # Atmospheric measurement concentrations
        
        #Space and time. Time is in seconds!
        z_res=cc['z_resolution'] #resolution of grid, m
        
        yrs=np.around(time_yr[-1]-time_yr[0])
    
        time_total=yrs*sPerYear #total model run time in seconds
        stpsperyear=cc['StepsPerYear'] #If this is for transient this number must (for now) be the same time steps as the input density/depth files. Make sure that it is a float.
        t_steps=np.int(yrs*stpsperyear)
        dt=time_total/t_steps #time step size. 
        model_time=np.arange(np.around(time_yr[0]*sPerYear),np.around(time_yr[-1]*sPerYear),dt) #set model time steps
        model_time_years=model_time/sPerYear
        nt=np.size(model_time) #number of time steps
        
        dz, z_edges, z_nodes, nodes, nz_P, nz_fv = space(depth,z_res) #call on space function to set up spatial grid
        
        rhoHL = MPRHO.rhoHLAnalytic(R,T,rho_i,rho0,rho_bco,z_nodes,Accu_m) # Get density profile from H&L analytic
        rho_co, por_co, por_tot, por_cl, por_op, bcoRho, LIDRho = porosity(rhoHL,T)
        
        if sitechoice=='SCENARIO':
            z_co = min(z_nodes[rhoHL>=(bcoRho)]) #close-off depth; bcoRho is close off density
            LIZ = min(z_nodes[rhoHL>=(LIDRho)]) #lock in depth; LIDRho is lock-in density
        
        gas=np.interp(model_time,time_yr_s,gas_org) #interpolate atmospheric gas history to model time.
        bc_u, bc_d, bc_u_0 = boundaries(gas_org) #set boundary and initial conditions: bc_u is atmospheric condition, bc_d is zero gradient.        
        phi_0 = np.zeros(nz_P) #phi is mixing ratio of gas.
        phi_0[:]=bc_u_0                    
                                                               
    for i_time in range(0,nt): #6/18/14: need to fix nt so that it is consistent with output from firnmodel.py

        #Accu=Accu_m#+0.01*np.sin(i_time) #6/18: what is this used for? Should probably use output from firnmodel.py
        if cc['UserData']: #the 1: is because the first value in the csv rows is the time.
            rho_prof=np.interp(z_nodes,f_depth[i_time,1:],f_density[i_time,1:]) #density profile, interpolating onto consistent grid (z_nodes); 6/18/14: need to make sure that z_nodes go deep enough to track.       
            dconint=np.interp(z_nodes,f_depth[i_time,1:],f_dcon[i_time,1:]) #this is the diffusivity constant
            Tprof=np.interp(z_nodes,f_depth[i_time,1:],f_temp[i_time,1:])
            T=Tprof[0]
            Accu_m=Accu_vec[i_time] #Accumulation in m/year 
            Accu_0=Accu_m/sPerYear #accumulation in m/second
        
        else:
            rho_prof=rhoHL
            
        rho_co, por_co, por_tot, por_cl, por_op, bcoRho, LIDRho = porosity(rho_prof,T) #get porosity
        
        z_co = min(z_nodes[rho_prof>=(bcoRho)]) #close-off depth; bcoRho is close off density
        LIZ = min(z_nodes[rho_prof>=(LIDRho)]) #lock in depth; LIDRho is lock-in density
              
        diffu,  d_eddy, diffu_full_fre, diffu_full_sch, diffu_full_sev, diffu_full_data = diffusivity(cc,rho_co, por_co, por_tot, por_cl, por_op, z_co, czd, LIZ,d_0,D_x,p_a,z_nodes,T,sitechoice, rhoprof=rho_prof) #get diffusivity profile
        
        #diffu=diffu*dconint
        
        if i_time==0:
            s=(nt,nz_P)
            phi=np.zeros(s)
            diffu_hold=np.zeros(s)
            rho_hold=np.zeros(s)
            phi_t=phi_0
            with open(ConcPath, "w") as f:
                writer = csv.writer(f)       
                writer.writerow(np.append(model_time_years[i_time],phi_t))
            with open(RhoPath, "w") as f:
                writer = csv.writer(f)       
                writer.writerow(np.append(model_time_years[i_time],rho_prof))
            with open(DiffuPath, "w") as f:
                writer = csv.writer(f)       
                writer.writerow(np.append(model_time_years[i_time],diffu))
            with open(ZPath, "w") as f:
                writer = csv.writer(f)       
                writer.writerow(np.append(model_time_years[i_time],z_nodes))
    
        diffu_P=diffu*por_op
        d_eddy_P=d_eddy*por_op

        dZ = np.concatenate(([1],np.diff(z_edges),[1]))
        
        dZ_u = np.diff(z_nodes)
        dZ_u = np.append(dZ_u[0], dZ_u)
    
        dZ_d = np.diff(z_nodes)
        dZ_d = np.append(dZ_d,dZ_d[-1])
    
        f_u = np.append(0, (1 -(z_nodes[1:] - z_edges)/dZ_u[1:]))
        f_d = np.append(1 - (z_edges - z_nodes[0:-1])/dZ_d[0:-1], 0)
        
        diffu_U = np.append(diffu_P[0], diffu_P[0:-1] )
        diffu_D = np.append(diffu_P[1:], diffu_P[-1])
    
        diffu_u =  1/ ( (1 - f_u)/diffu_P + f_u/diffu_U )
        diffu_d =  1/ ( (1 - f_d)/diffu_P + f_d/diffu_D )
    
        d_eddy_U = np.append(d_eddy_P[0], d_eddy_P[0:-1] )
        d_eddy_D = np.append(d_eddy_P[1:], d_eddy_P[-1])
    
        d_eddy_u =  1/ ( (1 - f_u)/d_eddy_P + f_u/d_eddy_U )
        d_eddy_d =  1/ ( (1 - f_d)/d_eddy_P + f_d/d_eddy_D )
        
        if cc['gravity']=="off" and cc['thermal']=="off":
            #print 'gravity and thermal are off'
            S_C_0=0.0
    
        elif cc['gravity']=='on' and cc['thermal']=='off':
            S_C_0=(-diffu_d+diffu_u)*(deltaM*g/(R*T))/dz #S_C is independent source term in Patankar
    
        elif cc['gravity']=='on' and cc['thermal']=='on':
            #print "thermal on"
            if cc["UserData"]:
                dTdz=np.gradient(Tprof)/dz
            else:
                dTdz=np.ones(np.size(diffu_d))
                dTdz=cc["Tgrad"]
                
            S_C_0=(diffu_d-diffu_u)*((-deltaM*g/(R*T))+(omega*dTdz))/dz #S_C is independent source term in Patankar. 8/28: Check thermal - should it still work in LIZ? if so use d_eddy+diffu
        
        #S_C_0=(-diffu_d+diffu_u)*(deltaM*g/(R*T)) #S_C is independent source term in Patankar
        
        S_C=S_C_0*phi_t #this line might be the troublesome one! Should it be phi_0 instead?
        
        #S_P=(-diffu_d+diffu_u)*(deltaM*g/(R*T)) #gravity term, S_P is phi-dependent source
        S_P=0.0
        
        b_0 = S_C*dZ
        
        rho_interface=np.interp(z_edges,z_nodes,rho_prof)
        
        w_edges, bubble_pres = w(z_edges,Accu_0,rho_interface,por_op,T,p_a,por_tot,por_cl,z_nodes,ad_method,dz)

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
                                
        a_P =  a_U + a_D + a_P_0 - S_P*dZ
           
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
        
        phi[i_time,:] = phi_t
        diffu_hold[i_time,:]=diffu
        rho_hold[i_time,:]=rho_prof
        
        if i_time>0:        
            with open(ConcPath, "a") as f:
                writer = csv.writer(f)       
                writer.writerow(np.append(model_time_years[i_time],phi_t))
            with open(RhoPath, "a") as f:
                writer = csv.writer(f)       
                writer.writerow(np.append(model_time_years[i_time],rho_prof))
            with open(DiffuPath, "a") as f:
                writer = csv.writer(f)       
                writer.writerow(np.append(model_time_years[i_time],diffu))
        
    return phi, diffu_hold, rho_hold, z_nodes
        
        
def space(depth,z_res):
    Lz=depth
    nz_fv=np.around(depth/z_res)
    nz_P=nz_fv+2
    
    dz=Lz/nz_fv
    z_edges=dz*np.arange(0,nz_fv+1)
        
    z_nodes=np.concatenate(([z_edges[0]],z_edges[0:-1]+np.diff(z_edges)/2,[z_edges[-1]]))
    nodes = np.size(z_nodes)
    
    return dz, z_edges, z_nodes, nodes, nz_P, nz_fv
    
def porosity(rho_prof,T):
    
    bcoRho = 1/( 1/(rho_i) + T*6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
    LIDRho = bcoRho - 14.0 #LIZ depth (Blunier and Schwander, 2000)
    
    ## Porosity, from Goujon et al., 2003, equations 9 and 10
    por_tot = 1-rho_prof/rho_i # Total porosity
    rho_co = bcoRho #use Martinerie close-off criteria
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
    
    return rho_co, por_co, por_tot, por_cl, por_op, bcoRho, LIDRho
      
def diffusivity(cc, rho_co, por_co, por_tot, por_cl, por_op, z_co, czd, LIZ,d_0,D_x,p_a,z_nodes,T,sitechoice, rhoprof = None): #rhoprof is density profile
        
    if rhoprof is None:
        rhoprof=rhoHL
    
    ## Constants
    d_eddy_sc=d_0 #Eddy diffusivity in the convective zone
    h=z_nodes
    dind=np.min(np.where(z_nodes>LIZ))
     
    ## Use Severinghaus relationship from Cuffey and Paterson
    d_0_sev=d_0*1.7
    #d_0_sev=d_0
    diffu_full_sev = D_x*d_0_sev*((p_0/p_a)*(T/T_0)**1.85*(2.00*(1-(rhoprof/rho_i))-0.167))
    diffu_full_sev = diffu_full_sev-diffu_full_sev[dind]
    diffu_full_sev[diffu_full_sev<=0] = 1e-15
    
    ## Use Schwander 1988, Eq. 2 Diffusivity (does not work very well) use 4e2
    ## for d_0
    k_sch = p_0/p_a*(T/253.16)**1.85 # Constant given in Schwander
    #diffu_full_sch =3.72*0.5*k_sch*(23.7*por_tot-2.84)*31.5 # Schwander' diffusivity relationship (for CO2). 31.5 is unit conversion. Added extra 3.72* 9/12/13
    diffu_full_sch =k_sch*(23.7*por_tot-2.84)/(1000**2) # Schwander' diffusivity relationship (for CO2). 1/1000**2 is unit conversion. Added extra 3.72* 9/12/13
    #ind = np.nonzero(h>LIZ)
    #diffu_full_sch[ind] = 0.001
    diffu_full_sch = diffu_full_sch-diffu_full_sch[dind]
    diffu_full_sch[diffu_full_sch<0] = 1.e-15
    
    ## Use Freitag, 2002, Eq 15 Diffusivity use 9e2 for d_0
    d_0_fre=d_0*4.9
    diffu_full_fre = D_x*d_0_fre*por_op**2.1
    diffu_full_fre = diffu_full_fre - diffu_full_fre[dind]
    diffu_full_fre[diffu_full_fre<=0] = 1e-15
    
    ## Use Christo's diffusivity data from NEEM-EU
    if sitechoice == 'NEEM':
        diffu_data=np.loadtxt(os.path.join(DataPath,'c_diffu_NEEM.txt'))
        h=diffu_data[:,0]
        diffu_full_data=D_x*d_0*diffu_data[:,1]

    elif sitechoice == 'WAIS':
        diffu_data=np.loadtxt(os.path.join(DataPath,'c_diffu_WAIS.txt'))
        h=diffu_data[:,0]
        diffu_full_data=D_x*d_0*diffu_data[:,1]
        
    elif sitechoice == 'SCENARIO':
        diffu_data=np.loadtxt(os.path.join(DataPath,'c_diffu_NEEM.txt'))
        h=diffu_data[:,0]
        diffu_full_data=D_x*d_0*diffu_data[:,1]        
        
      
    ## try a random profile to test if the model is working.
    #h=1:100 
    #diffu_full=d_0*23*ones(length(h))
    diffu_full_data = np.interp(z_nodes,h,diffu_full_data)
    #diffu_full_data=diffu_full
    #d_eddy=1. #hack if not using eddy diffusivity (wrap entirely into D)
    
    ## Add in high diffusivity in convective zone and low diffusivity below LIZ
    if cc["diffu"]=="data":
        diffu_full=diffu_full_data #change this line to change your choice of diffusivity
    
    if cc["diffu"]=="Sev":
        diffu_full=diffu_full_sev #change this line to change your choice of diffusivity
        
    if cc["diffu"]=="Fre":
        diffu_full=diffu_full_fre #change this line to change your choice of diffusivity
    
    if cc["diffu"]=="Sch":
        diffu_full=diffu_full_sch #change this line to change your choice of diffusivity
    
    #Add eddy diffusivity terms: convective zone and non-diffusive zone
    d_eddy=np.zeros(np.size(diffu_full))
    ind = np.nonzero(z_nodes<czd)
    d_eddy_surf=2.426405E-5 #Kawamura, 2006
    H_scale=czd
    d_eddy_up=d_eddy_surf*np.exp(-1*z_nodes/H_scale)
    #d_eddy[ind] = diffu_full[ind]*10
    ind = np.flatnonzero(z_nodes>LIZ)
    ind2 = np.flatnonzero(z_nodes<z_co)
    ind3 = np.intersect1d(ind,ind2)
    d_eddy[ind3] = diffu_full[ind]/10
    d_eddy=d_eddy+d_eddy_up
    diffu_full[ind]=1e-15 #set molecular diffusivity equal to zero after LIZ - eddy diffusivity term drives diffusion below
    
    diffu=diffu_full
    ## Interpolate diffusivity profile to the finite volume nodes (model space)
    #deepnodes = z_nodes>LIZ #leftover line from matlab?
    
    return diffu,  d_eddy, diffu_full_fre, diffu_full_sch, diffu_full_sev, diffu_full_data  
            
def boundaries(gas_org):
    bc_u_0 = gas_org[0] #this is concentration on upper boundary (i.e. atmosphere) It is not actually necessary here, because it gets looped over.
    bc_type = 1
    bc_u   = np.concatenate( ([bc_u_0], [bc_type]))

    bc_d_0 = 0 
    bc_type = 2
    bc_d   = np.concatenate(([ bc_d_0 ], [ bc_type ]))
     
    return bc_u,bc_d, bc_u_0
    
def firnairmodel(airconfig):
    
    with open(airconfig, "r") as f:
        json_data=f.read()
        cc = json.loads(json_data)
    
    #global depth, T, Accu_0, g, d_0, z_nodes, D_x, p_a, sitechoice, deltaM, por_tot, por_cl, ad_method,options,dz,omega,z_co,z_edges,rho_co, por_co, por_tot, por_cl, por_op, bcoRho, LIDRho
    gaschoice_all=cc["gaschoice"]
    nogas=len(gaschoice_all)
    d={}
    zz={}
    
    for jj in xrange(nogas):
        print jj
        
        gaschoice=gaschoice_all[jj]
        
        print jj, gaschoice
        
        runtype = cc["runtype"] #this line chooses transient or steady-state
        
        if runtype == 'transient':
            phi, diffu_hold, rho_hold, z_nodes = FirnAir_TR(cc,gaschoice,jj)
            print 'maximum = %s' % np.max(phi)
            
        elif runtype == 'steady':
            phi, a_P_out, bubble_pres,z_nodes =FirnAir_SS(cc,gaschoice)
            print 'maximum = %s' % np.max(phi)
    
        d[gaschoice]=phi
        print '%s done' % gaschoice
        d['nodes']=z_nodes
            
    return d
    
#    depth = cc["depth"] # m
#    ad_method = cc["ad_method"]
#    measurements = cc["measurements"] 
#    options=cc["options"]    
#    # Set up parameters for different sites.
#    sitechoice = cc["sitechoice"]
#    g, p_a, T, Accu_0, czd, z_co, LIZ, rho0, hemisphere = MPS.sites(sitechoice)   
#    Accu_m=Accu_0 #Accumulation in m/year 
#    Accu_0=Accu_0/sPerYear #accumulation in m/second
#        
#    D_x, M, deltaM, conc1, firn_meas, d_0, omega = MPG.gasses(gaschoice, sitechoice,T,p_a,DataPath,hemisphere,measurements)
#    
#    time_yr=conc1[:,0] # Atmospheric measurements times
#                
#    time_yr_s=time_yr*sPerYear
#    gas_org=conc1[:,1] # Atmospheric measurement concentrations 
#    if firn_meas != 'None': #incomplete
#        meas_depth=firn_meas[:,0]
#        meas_conc=firn_meas[:,1] # measured gas concentrations in firn
#        #meas_conc=firn_meas[:,2] # gravity corrected
#        meas_uncert=firn_meas[:,3]
#        
#    #Space and time. Time is in seconds!
#    z_res=cc['z_resolution'] #resolution of grid, m
#    
#    yrs=np.around(time_yr[-1]-time_yr[0])
#
#    time_total=yrs*sPerYear #total model run time in seconds
#    stpsperyear=cc['StepsPerYear'] #If this is for transient this number must (for now) be the same time steps as the input density/depth files. Make sure that it is a float.
#    t_steps=np.int(yrs*stpsperyear)
#    dt=time_total/t_steps #time step size. 
#    model_time=np.arange(np.around(time_yr[0]*sPerYear),np.around(time_yr[-1]*sPerYear),dt) #set model time steps
#    nt=np.size(model_time) #number of time steps
#    
#    dz, z_edges, z_nodes, nodes, nz_P, nz_fv = space(depth,z_res) #call on space function to set up spatial grid
#    
#    rhoHL = MPRHO.rhoHLAnalytic(R,T,rho_i,rho0,rho_bco,z_nodes,Accu_m) # Get density profile from H&L analytic
#    rho_co, por_co, por_tot, por_cl, por_op, bcoRho, LIDRho = porosity(rhoHL)
#    
#    if sitechoice=='SCENARIO':
#        z_co = min(z_nodes[rhoHL>=(bcoRho)]) #close-off depth; bcoRho is close off density
#        LIZ = min(z_nodes[rhoHL>=(LIDRho)]) #lock in depth; LIDRho is lock-in density
#    
#    diffu,  d_eddy, diffu_full_fre, diffu_full_sch, diffu_full_Sev, diffu_full_data = diffusivity(rho_co, por_co, por_tot, por_cl, por_op, z_co, czd, LIZ, rhoHL) #get diffusivity profiles
#    
#    dcon=1.0
#    diffu=diffu*dcon   
#    gas=np.interp(model_time,time_yr_s,gas_org) #interpolate atmospheric gas history to model time.
#    bc_u, bc_d, bc_u_0 = boundaries(gas_org) #set boundary and initial conditions: bc_u is atmospheric condition, bc_d is zero gradient.
        

        
    #phi_0 = np.zeros(nz_P) #phi is mixing ratio of gas.
    #phi_0[:]=bc_u_0
        
    #trapped = 1000*(bubble_pres*por_cl_interface)*(p_a/101325)*(273.15/T)/(rho_interface/1000) #amount of trapped air in ml STP/cm^3 (from Christo)
        #print trapped    
                    
        #phi_final=np.append([phi_0],[phi],axis=0)
        #phi=phi.T
        
        #phi_toplot=phi[:,1:-1:20]
        #aa=np.shape(phi_toplot)
        #aa=aa[1]
        #d[gaschoice]=phi
        #d['bubblePres']=bubble_pres
        #d['TAC']=trapped
        

    
    #d15max=np.max(phi[:,-1])
    #d15LID=np.log(d15max)*R*T/(g*deltaM)
    
    #print 'dcon =', dcon
    #print 'd15N2 LID =', d15LID
    #print 'LID from density/temperature (Martinerie) =', LIZ    
    #print 'Close off depth =', z_co
    #print 'Difference (close-off - d15N2 LID)=', z_co-d15LID
    

    
if __name__ == "__main__":
    
    reload(MPG)
    reload(MPS)
    reload(MPD)
    reload(plots)
    reload(MPRHO)
    import time
    
    tic=time.time()
    
    airconfig = os.path.join(os.path.dirname(__file__),'FirnAirConfig.json')

    d = firnairmodel(airconfig)

    #if firn_meas != 'None':
    #    gas_meas=np.interp(z_nodes, meas_depth, meas_conc);
    #
    #plotting = 'on'
    #
    #if plotting != 'off':
    #    plots.makeplots(plotting,z_nodes,phi,gas_meas,meas_depth,meas_conc,ResultsPlace,
    #                    diffu_full_Sev,diffu_full_fre,diffu_full_sch,diffu_full_data, meas_uncert=meas_uncert)
    
    #d40Ar=d['d40Ar']
    #d15N2=d['d15N2']
    #N2_p=d15N2[:,-1]
    #Ar_p=(d40Ar[:,-1]-1)/4+1
    #plt.figure(1)
    #plt.plot(z_nodes,N2_p,'b')
    #plt.plot(z_nodes,Ar_p,'r')
    #plt.show()
    
    
        
                
    elapsed=time.time()-tic
    elapsed_min=elapsed/60.
    mins=np.floor(elapsed_min)
    secs=(elapsed_min-mins)*60
    print mins, 'min', secs, 'sec elapsed'
        