# -*- coding: utf-8 -*-
"""
Functions required for the Darcy-type liquid water flow scheme

@author: Vincent
"""

import numpy as np


def hydrconducsat_Calonne(rad,rho):
    '''Saturated hydraulic conductivity of Calonne et al. Eq. (6)'''
    mu = 0.001792 #dynamic viscosity of water at 273.15K [kg m-1 s-1]
    bigksat = 3*(rad)**2 * 1000*9.81/mu * np.exp(-0.013*rho) #[m s-1]
    return(bigksat)

def vG_Yama(rad,rho,thetaeff):
    '''Pressure head and relative hydraulic conductivity computations
    from the van Genuchten (1980) model with the Yamaguchi et al. (2012)
    parameterisation'''
    alpha = 4.4e6*(rho/(2*rad))**(-0.98) #.alpha parameter, Yamaguchi 2012 Eq.(6)
    n     = 1+2.7e-3*(rho/(2*rad))**0.61 #n parameter, Yamaguchi 2012 Eq. (7)
    m     = 1-(1/n) #m parameter, Yamaguchi 2012 p.7
    head  = 1/alpha * (thetaeff**(-1/m)-1)**(1/n) #head pressure, Hirashima 2010 (9)
    bigkrel = thetaeff**(1/2) * (1-(1-thetaeff**(1/m))**m)**2 # Hirashima 2010 (10)
    return(head,bigkrel)

def thetae_update(absfl,th_i,th_s,LWC,dz):
    '''Updates effective saturation for a given total water flux absfl [m]'''
    lw_in  = np.append(0,absfl)
    lw_out = np.append(absfl,0)
    th_w   = (LWC+lw_in-lw_out)/dz
    th_e   = (th_w-th_i)/(th_s-th_i) #effective water saturatin, Hirashima 2010 (5)
    stab_e    = 1e-9 #stabilisation theta_e
    th_e   = np.maximum(stab_e,th_e) #avoid negative effective saturation
    th_e   = np.minimum(1-stab_e,th_e) #avoid effective saturation equal to 1
    return(th_e)

def thetaeff_equaliser(th_i2,th_s2,LWC2,dz2):
    '''
    Computes the total flow needed to equalise saturation between neighbouring nodes
    All input arrays must be of two elements
    '''
    th_w   = LWC2/dz2
    th_e0  = (th_w-th_i2)/(th_s2-th_i2) #effective water saturatin, Hirashima 2010 (5)
    # lwflux from index[0] to index[1] ensures equal saturation between both volumes #
    lwflux = ((dz2[0]*(th_s2[0]-th_i2[0]))**(-1)+(dz2[1]*(th_s2[1]-th_i2[1]))**(-1))**(-1) * (th_e0[0]-th_e0[1])
    return(lwflux)

def vG_Yama_params(rad,rho):
    '''Computes the van Genuchten parameters following the parameterisation
    of Yamaguchi et al. (2012)'''
    alpha = 4.4e6*(rho/(2*rad))**(-0.98) #.alpha parameter, Yamaguchi 2012 Eq.(6)
    n     = 1+2.7e-3*(rho/(2*rad))**0.61 #n parameter, Yamaguchi 2012 Eq. (7)
    m     = 1-(1/n) #m parameter, Yamaguchi 2012 p.7
    return(alpha,n,m)

def phead_vG(alpha,n,m,thetaeff):
    '''Computes pressure head according to the van Genuchten model'''
    head  = 1/alpha * (thetaeff**(-1/m)-1)**(1/n) #head pressure, Hirashima 2014 (3)
    return(head) 

def krel_vG(m,thetaeff):
    '''Computes relative hydraulic conductivity according to the van Genuchten model'''
    bigkrel = thetaeff**(1/2) * (1-(1-thetaeff**(1/m))**m)**2 # Hirashima 2010 (10)
    return(bigkrel) 

def dfdg_derivative(th_sfull,th_ifull,th_efull,alphafull,nfull,mfull,dzfull):
    '''
    Computes the derivative of the equilibrium variable (f_eq, defined by setting all terms
    of Eq.(20) of Hirashima et al. (2010) on the right-hand-side) with respect to the water
    flux at the interface of two neighbouring nodes
    '''
    th_s,th_i,th_e    = th_sfull[0:-1],th_ifull[0:-1],th_efull[0:-1]
    alpha,n,m,dz      =    alphafull[0:-1],nfull[0:-1],mfull[0:-1],dzfull[0:-1]
    th_sd,th_id,th_ed = th_sfull[1:],th_ifull[1:],th_efull[1:]
    alphad,nd,md,dzd  =    alphafull[1:],nfull[1:],mfull[1:],dzfull[1:] 
    dfdg = 1/((th_s-th_i)*alpha*n*m*dz) * th_e**(-1*(1+1/m)) * (th_e**(-1/m)-1)**((1-n)/n) + \
        1/((th_sd-th_id)*alphad*nd*md*dzd) * th_ed**(-1*(1+1/md)) * (th_ed**(-1/md)-1)**((1-nd)/nd)
    return(dfdg)



def flux_bisection(gc,LWCav,glwcacm,th_i,th_s,lwc,dz,avG,nvG,mvG,eps_cvg):
    '''
    Bisection algorithm to find guess of water flux (gc) between two neighbouring nodes that brings
    head pressures close to equilibrium (f_eq, defined by setting all terms of Eq.(20) of Hirashima
    et al. (2010) on the right-hand-side)
    '''
    bisitmax = 100 #maximum number of iteration for bisection algorithm
    bisit    = 0 #iteration number for bisection algorithm
    cvg_bis = False #convergence criterion for Bisection
    dltz    = 1/2*sum(dz) #distance between the centres of the two nodes
    gth_e = thetae_update(gc,th_i,th_s,lwc,dz)
    ghd   = phead_vG(avG,nvG,mvG,gth_e) #pressure head [m]
    f_eq = ghd[0]-ghd[1]-dltz #Hirashima 2010 Eq.(20) evaluated at interface
    # Initialise bisection bounds #
    g0   = 0. #lower bisection bound
    g1   = min(LWCav[0],glwcacm[1]) #upper bisection bound (limited by lwc available and pore space available)
    while (cvg_bis==False and bisit<bisitmax): #start Bisection (if cvg_bis is False)
        gprev0    = np.copy(gc) #flux guess at previous iteration
        if f_eq<0: #hd[0] too low -> increase outgoing flux
            g0 = np.copy(gc) #updated lower bound
            gc = (g1+g0)/2 #updated guess
        elif f_eq>0: #hd[j1] too high -> decrease outgoing flux
            g1 = np.copy(gc) #updated lower bound
            gc = (g1+g0)/2 #updated guess
        gth_e = thetae_update(gc,th_i,th_s,lwc,dz)
        ghd   = phead_vG(avG,nvG,mvG,gth_e) #pressure head [m]
        # Evaluate equilibrium #
        f_eq = ghd[0]-ghd[1]-dltz #Hirashima 2010 Eq.(20) evaluated at interface
        if (f_eq<0 and gth_e[0]<1e-8): 
            # Equilibrium requires higher outflow but max outflow already prescribed (CV[0] dried) #
            cvg_bis = True
        elif (f_eq<0 and gth_e[1]>0.95):
            # Equilibrium requires higher outflow but underlying node is aturated (CV[j1+1] saturated) #
            cvg_bis = True
        elif (f_eq>0 and gc<=1e-6):
            # Equilibrium requires less outflow but min outflow already prescribed
            gc = 0. #set outlfow to 0
            cvg_bis = True
        if abs(f_eq)<eps_cvg or abs(gc-gprev0)<1e-6:
            # Flux estimate has converged
            cvg_bis = True
        bisit += 1 #increase iteration number
        if bisit==bisitmax:
            print('Maximum iteration number reached in bisection algorithm')  
    return(gc)


def flux_newtonraphson(gc,LWCav,glwcacm,th_i,th_s,lwc,dz,avG,nvG,mvG,eps_cvg):
    '''
    Newton-Raphson algorithm to find guess of water flux (gc) between two neighbouring nodes that brings
    head pressures close to equilibrium (f_eq, defined by setting all terms of Eq.(20) of Hirashima
    et al. (2010) on the right-hand-side)
    If the algorithm diverges: call for the bisection algorithm
    '''
    nritmax = 20 #maximum number of iteration for Newton-Raphson algorithm
    nrit    = 0 #iteration number for Newton-Raphson algorithm
    cvg_nr = False #convergence criterion for Newton-Raphson
    dltz    = 1/2*sum(dz) #distance between the centres of the two nodes
    gth_e = thetae_update(gc,th_i,th_s,lwc,dz)
    ghd   = phead_vG(avG,nvG,mvG,gth_e) #pressure head [m]
    f_eq = ghd[0]-ghd[1]-dltz #Hirashima 2010 Eq.(20) evaluated at interface
    while ((cvg_nr==False) and nrit<nritmax):
        # Values of previous iteration #
        gprev0    = np.copy(gc)
        fprev0    = np.copy(f_eq)
        
        if (f_eq<0 and gth_e[0]<1e-8): 
            # Equilibrium requires higher outflow but max outflow already prescribed (CV[0] dried) #
            cvg_nr = True
        elif (f_eq<0 and gth_e[1]>0.95):
            # Equilibrium requires higher outflow but underlying node is aturated (CV[j1+1] saturated) #
            cvg_nr = True
        elif (f_eq>0 and gc<=1e-6):
            # Equilibrium requires less outflow but min outflow already prescribed
            gc = 0. #set outlfow to 0
            cvg_nr = True

        else: #Newton-Raphson to improve current guess of glw
            #Computation for derivative d(f_eq)/d(glw[j1])
            dfdg = dfdg_derivative(th_s,th_i,gth_e,avG,nvG,mvG,dz)
            deltaglw = -1*fprev0/dfdg #step in glw guess
            gc = gprev0+deltaglw #adjust glw guess
            gth_e = thetae_update(gc,th_i,th_s,lwc,dz)
            ghd   = phead_vG(avG,nvG,mvG,gth_e) #pressure head [m]
            f_eq = ghd[0]-ghd[1]-dltz #Hirashima 2010 Eq.(20) evaluated at interface
            if ((abs(f_eq)>abs(fprev0)) or (abs(dfdg)>1e6)):
                #Newton-Raphson diverges: use bisection algorithm #
                gc = flux_bisection(gprev0,LWCav,glwcacm,th_i,th_s,lwc,dz,avG,nvG,mvG,eps_cvg)
                cvg_nr = True
        
        nrit += 1
        if abs(f_eq)<eps_cvg or abs(gc-gprev0)<1e-6:
            # Flux estimate has converged
            cvg_nr = True
        if nrit==nritmax:
            print('Maximum iteration number reached in Newton-Raphson algorithm')  
    return(gc)


def runoffZuoOerlemans(dt,slope,lwcexcess,inds):
    '''
    Computes runoff according to the Zuo and Oerlemans (1996) parameterisation
    lwcexcess: liquid water content in excess of irreducible water content [m]
    inds:      indices of the firn column domain where runoff should be computed
    '''
    c1zuo = 1.5*24*3600 #constant from Zuo and Oerlemans (1996), converted in [s]
    c2zuo = 25.*24*3600 #constant from Zuo and Oerlemans (1996), converted in [s]
    c3zuo = 140. #constant from Zuo and Oerlemans (1996) [/]
    tstar = c1zuo + c2zuo*np.exp(-1*c3zuo*slope) # Eq.(22) Zuo and Oerlemans 1996 [s]
    rfout = np.zeros(len(lwcexcess)) #initialise runoff
    rfout[inds] = dt*lwcexcess[inds]/tstar #from Eq.(21) Zuo and Oerlemans 1996 [m]
    return(rfout)

def runoffDarcy(dt,slope,bigk,inds):
    '''
    Computes lateral runoff according to Darcy's law assuming a certain slope
    and that horizontal firn hydraulic properties are homogeneous.
    bigk: hydraulic conductivity [m s-1]
    inds: indices of the firn column domain where runoff should be computed
    '''
    rfout = np.zeros(len(bigk)) #initialise runoff
    rfout[inds] = dt*bigk[inds]*slope #total outgoing Darcy flux
    return(rfout)







