# UW Community Firn-Air Transport Model
# Max Stevens
# version 0.1, 6/3/13

import sys
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.sparse.linalg as splin
from scipy.sparse import spdiags
from scipy.sparse.linalg import lsqr
import math

sys.path.insert(0, '~/Documents/Grad_School/PIRE/FMbeta/code/gasmodel')
os.chdir('/Users/Max/Documents/Grad_School/PIRE/FMbeta/code/gasmodel')
np.set_printoptions(linewidth=300) #set output reading to be wider. A bit easier to read :)

def w(z_edges_vec,Accu,rho_interface,por_op,T,p_a,por_tot,por_cl):
    M_air = 29e-3
    rho_i=0.917
    g=9.81
    R = 8.314
    rho_bco = 0.815
    por_op_star=por_op*np.exp(M_air*g*z_nodes/(R*T)) #Christo Thesis 5.3 Do I want nodes or edges for z?
    
    w_ice=Accu*rho_i/rho_interface
    w_ice_co=w_ice[ind_co]
    
    zeta=(por_co/por_tot)/(1+np.log(w_ice/w_ice_co))
    
    p_ratio=np.zeros(nz_P)
    
    p_ratio[0:ind_co+1]=
    #P_z=p_a*np.exp(M_air*g*z_edges_vec/(R*T)) #Barometric Equation
    
    #need s_cl, 
    
    w = 1.7*Accu/rho_interface #probably need the interfaces!
    #w = A*rho_i/(por_op_star)*(por_co/rho_bco*
    
    return w

def A(P):
    A = np.maximum( (1 - 0.1 * np.abs( P ) )**5, np.zeros(np.size(P) ) )
    return A    
    
def F_upwind(F):
    F_upwind = np.maximum( F, 0 )
    return F_upwind

def solver(a_U,a_D,a_P,b):
    nz=np.size(b)
    
    Diags = (np.append([a_U,-a_P],[a_D],axis=0))
    cols=np.array([1, 0, -1])
    
    
    big_A=spdiags(Diags,cols,nz,nz,format='csc')
    big_A=big_A.T
    denseA=big_A.todense()
    #print ('a_P='+str(a_P[0:6]))
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
    
def transient_solve(z_edges_vec,z_P_vec,nt,dt,Gamma_P,bc_u,bc_d,phi_0,rhoHL,deepnodes,R,nz_P,nz_fv,por_op,gas,Accu,T,p_a,por_tot,por_cl):
    
    #nz_P = np.size(z_P_vec)
    #nz_fv = nz_P - 2
    Z_P = z_P_vec
     
    dZ = np.concatenate(([1],np.diff(z_edges_vec),[1]))
    
    S_C=0
    S_C=S_C*np.ones(nz_P)
    
    b_0 = S_C*dZ
    
    dZ_u = np.diff(Z_P)
    dZ_u = np.append(dZ_u[0], dZ_u)

    dZ_d = np.diff(Z_P)
    dZ_d = np.append(dZ_d,dZ_d[-1])


    f_u = np.append(0, (1 -(z_P_vec[1:] - z_edges_vec)/dZ_u[1:]))
    f_d = np.append(1 - (z_edges_vec - z_P_vec[0:-1])/dZ_d[0:-1], 0)
     
    Gamma_U = np.append(Gamma_P[0], Gamma_P[0:-1] )
    Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])
    
    Gamma_u =  1/ ( (1 - f_u)/Gamma_P + f_u/Gamma_U )
    Gamma_d =  1/ ( (1 - f_d)/Gamma_P + f_d/Gamma_D )
    
    #Z_w_velo = z_edges_vec #relic of old code?
    
    rho_interface=np.interp(z_edges_vec,z_P_vec,rhoHL)
    w_edges = w(z_edges_vec,Accu,rho_interface,por_op,T,p_a,por_tot,por_cl)
    
    w_u = np.append(w_edges[0],  w_edges )
    w_d = np.append(w_edges, w_edges[-1])
    
    D_u = (Gamma_u / dZ_u)
    D_d = (Gamma_d / dZ_d)
    
    F_u =  w_u 
    F_d =  w_d 
    
    P_u = F_u/ D_u
    P_d = F_d/ D_d
    print P_u
    
    a_U = D_u * A( P_u ) + F_upwind(  F_u )
    a_D = D_d * A( P_d ) + F_upwind( F_d )

    a_P_0 = por_op*dZ/dt

    s=(nz_P,nt)
    phi=np.zeros(s)
    
    phi_t=phi_0
    
    #S_P = (M*g/(R*T))*phi_0
    
    a_P =  a_U + a_D + a_P_0
    
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
        #a_U[-1] = 0
        a_D[-1] = 0
        a_U[-1] = 1
        b[-1]=dZ_u[-1]*bc_d[0]
        
        #print ("a_U=" + str(a_U[-1]))
        phi_t = solver(a_U,a_D,a_P,b)
        
        phi[:,i_time]=phi_t
        
        a_P=a_U+a_D+a_P_0
        
    return phi
        
        
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
    
def rhoHL(T,Accu,rho0,z_nodes,rho_i,R): 
    
    rho_c = 0.550
    h=z_nodes
    rho_bco = 0.815
    
    
    k0 = 11  * np.exp(-10160/(R*T ))
    k1 = 575 * np.exp(-21400/(R*T ))
    
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
    bco = min(age[rho_h>=rho_bco])
    
    rhoHL=rho_h
    
    return rhoHL
    
#def gasses():
#    conc1=np.loadtxt('SCENARIO_NEEM08_CO2.txt',skiprows=1)
#    firn_meas=np.loadtxt('CO2samples.txt',skiprows=1)
#    meas_conc=firn_meas.data(:,2)
#    
def diffusivity(rho_i,rhoHL,z_nodes):
    
    ## Constants
    p_0 = 1.01325e5 # Standard Amtmospheric Pressure, Pa
    T_0 = 273.15 # Standard Temp
    d_0 = 5e2 # Free air diffusivity, CO2, m**2/yr Schwander, 1988 reports 7.24 mm**2/s =379 m**2/yr
    d_eddy_sc=d_0 #Eddy diffusivity in the convective zone
    h=z_nodes
    
    ## Porosity, from Goujon et al., 2003, equations 9 and 10
    por_tot = 1-rhoHL/rho_i # Total porosity
    # rho_co = rho_bco #use Martinerie close-off criteria
    rho_co = 0.815 # User chosen close off-density (put in site-specific in sites?)
    por_co = 1 - rho_co/rho_i # Porosity at close-off
    alpha = 0.37 # constant determined in Goujon
    por_cl = alpha*por_tot*(por_tot/por_co)**(-7.6)
    por_cl[por_cl>1]=1
    por_cl = por_cl*por_tot # Final closed porosity
    por_op = por_tot - por_cl # Open Porosity
    
    
    ## Use Severinghaus relationship from Cuffey and Paterson
    d_0_sev=d_0*1.5
    diffu_full_Sev = D_x*d_0_sev*((p_0/p_a)*(T/T_0)**1.85*(2.00*(1-(rhoHL/rho_i))-0.167)) 
    diffu_full_Sev[diffu_full_Sev<=0] = 1e-9
    
    ## Use Schwander 1988, Eq. 2 Diffusivity (does not work very well) use 4e2
    ## for d_0
    k_sch = p_0/p_a*(T/253.16)**1.85 # Constant given in Schwander
    diffu_full_sch =4*0.5*k_sch*(23.7*por_tot-2.84)*31.5 # Schwander' diffusivity relationship (for CO2). 31.5 is unit conversion. Added extra 4* 6/10/13
    ind = np.nonzero(h>LIZ)
    diffu_full_sch[ind] =0.001
    diffu_full_sch[diffu_full_sch<0] = 1e-15
    
    ## Use Freitag, 2002, Eq 15 Diffusivity use 9e2 for d_0
    d_0_fre=d_0*4
    diffu_full_fre = D_x*d_0_fre*por_op**2.1
    diffu_full_fre[diffu_full_fre<=0] = 1e-15
    
    ## Use Christo's diffusivity data from NEEM-EU
    
    diffu_data=np.loadtxt('c_diffu.txt')
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
    
    return diffu,deepnodes,por_op,d_eddy,diffu_full_fre,diffu_full_sch,diffu_full_Sev,diffu_full_Christo,por_cl,por_tot,por_co
    
def boundaries(gas_org):
    bc_u_0 = gas_org[0] #this is concentration on upper boundary (i.e. atmosphere) It is not actually necessary here, because it gets looped over.
    bc_type = 1
    bc_u   = np.concatenate( ([bc_u_0], [bc_type]))

    bc_d_0 = 0 #how do I handle this?
    bc_type = 2
    bc_d   = np.concatenate(([ bc_d_0 ], [ bc_type ]))
     
    return bc_u,bc_d, bc_u_0

def gasinterp(time_yr,gas_org,model_time):
    gas=np.interp(model_time,time_yr,gas_org)
    
    return gas
    
        
if __name__ == "__main__":
    #globals
    rho_i = 0.917
    R = 8.314
    
    #these parameters are for NEEM - 6/9/13 need to set up a "sites" script
    g=9.81        
    depth = 80. 
    p_a = 7.45e4 
    T = -28.9 +273.16 
    Accu = 0.20 
    czd = 4 
    z_co = 78
    LIZ = 63
    rho0=0.36
    
    #Gas Properties. Just CO2 for now.Should make a definition for this.
    D_x = 1. #free-air diffusivity relative to CO2
    M = 44.01e-3 # molecular mass
    conc1=np.loadtxt('SCENARIO_NEEM08_CO2.txt',skiprows=1) #load data: atmospheric CO2 history.
    firn_meas=np.loadtxt('CO2samples.txt',skiprows=2) #load data: firn air samples.
    meas_conc=firn_meas[:,2] # measured CO2 concentration in firn
    time_yr=conc1[:,0] 
    gas_org=conc1[:,1]
    meas_depth=firn_meas[:,0]
    meas_uncert=firn_meas[:,3]
    
    #Space and time
    z_res=0.5 #resolution of grid
    #t_res=0.01
    t_res=0.08 #time step size. 
    model_time=np.arange(time_yr[0],time_yr[-1],t_res) #set model time steps
    dt=t_res
    nt=np.size(model_time) #number of time steps
    
    dz, z_edges_vec, z_P_vec, z_nodes, nodes, nz_P, nz_fv = space(depth,z_res) #call on space function to set up spatial grid
    
    rhoHL = rhoHL(T,Accu,rho0,z_nodes,rho_i,R) # Get density profile from H&L analytic
    diffu, deepnodes, por_op, d_eddy, diffu_full_fre,diffu_full_sch,diffu_full_Sev, diffu_full_Christo, por_tot, por_cl, por_co = diffusivity(rho_i,rhoHL,z_nodes) #get diffusivity profiles
    ind_co=np.argmax(z_edges_vec>=z_co) #index of close-off depth in z_edges vec. A bit of a hack for now...
    gas = gasinterp(time_yr,gas_org,model_time) #interpolate atmospheric gas history to model time.
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
    
    phi=transient_solve(z_edges_vec,z_P_vec,nt,dt,Gamma_P,bc_u,bc_d,phi_0,rhoHL,deepnodes, R,nz_P,nz_fv,por_op,gas,Accu,T,p_a,por_tot,por_cl)
    #phi_final=np.append([phi_0],[phi],axis=0)
    #phi=phi.T
    
    phi_toplot=phi[:,1:-1:20]
    aa=np.shape(phi_toplot)
    aa=aa[1]
    
    neem_meas=np.interp(z_nodes, meas_depth, meas_conc);
    
    ## Plotting
    #
    #fig1=plt.figure(1)
    #plt.clf()
    #plt.plot(Z_P,phi[:,-3],'b')
    #plt.plot(Z_P,neem_meas,'r')
    #plt.plot(meas_depth,meas_conc,'k+')
    #plt.xlabel('Depth (m)')
    #plt.ylabel('$CO_{2}$ (ppm)')
    ##plt.title('Concentration of $CO_{2}$ in firn at NEEM')
    #plt.legend(('Max\'s model','Measurements'))
    #plt.grid()
    #fname='ESS524_fig1.eps'
    #plt.savefig(fname,dpi=100)
    #
    #fig2=plt.figure(2)
    #plt.clf()
    #plt.plot(Z_P,diffu_full_Sev,'r')
    #plt.plot(Z_P,diffu_full_fre,'b')    
    #plt.plot(Z_P,diffu_full_sch,'g')    
    #plt.plot(Z_P,diffu_full_Christo,'m')
    #plt.xlabel('Depth (m)')
    #plt.ylabel('Diffusivity ($m^{2} yr^{-1})$')
    ##plt.title('Diffusivity with depth for different parameterizations')
    #plt.legend(('Severinghaus et al., 2001','Freitag et al., 2002','Schwander, 1988','Buizert,2011 (tuned from data)'))  
    #plt.grid()
    #fname='ESS524_fig2.eps'
    #plt.savefig(fname,dpi=100)    
    ##plt.show()  
    #
    #fig3=plt.figure(3)
    #plt.clf()
    #plt.plot(time_yr,gas_org,'b')
    #plt.xlabel('Time (yrs)')
    #plt.ylabel('CO$_{2}$ (ppm)')
    ##plt.title('Northern Hemisphere Atmospheric CO$_{2}$ Concentration')
    #plt.grid()
    #fig3.set_size_inches(10,5)
    #fname='ESS524_fig3.eps'
    #
    #plt.savefig(fname,dpi=100)    
    #plt.show()  
    #files = []
    #fig2 = plt.figure()
    #ax = fig2.add_subplot(111)
#    for i in range(aa):  # nt frames
#        plt.plot(Z_P,phi_toplot[:,i])
#        plt.axis((Z_P[0],Z_P[-1],270,400))
#        plt.xlabel('Depth (m)')
#        plt.ylabel('$CO_{2}$ (ppm)')
#        plt.title('Concentration of $CO_{2}$ in firn at NEEM')
#
#        fname = 'saved_figs/'+str('%03d' %i)+'.png'
#        plt.savefig(fname,dpi=100)
#        plt.clf()
        

        