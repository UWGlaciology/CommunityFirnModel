#!/usr/bin/env python
'''
Functions to solve the diffusion equation
'''

import numpy as np
from scipy import interpolate
import scipy.integrate
from scipy.sparse import spdiags
import scipy.sparse.linalg as splin
from constants import *
import sys


def solver(a_U, a_D, a_P, b):
    '''
    function for solving matrix problem

    :param a_U:
    :param a_D:
    :param a_P:
    :param b:

    :return phi_t:
    '''

    nz = np.size(b)

    diags = (np.append([a_U, -a_P], [a_D], axis = 0))
    cols = np.array([1, 0, -1])

    big_A = spdiags(diags, cols, nz, nz, format = 'csc')
    big_A = big_A.T

    rhs = -b
    phi_t = splin.spsolve(big_A, rhs)

    return phi_t

####!!!!

def transient_solve_TR(z_edges, Z_P, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s, tot_rho, c_vol, airdict=None):
    '''
    transient 1-d diffusion finite volume method
    :param z_edges:
    :param Z_P:
    :param nt:
    :param dt:
    :param Gamma_P:
    :param phi_0:
    :param nz_P:
    :param nz_fv:
    :param phi_s:
    :return phi_t:
    '''

    phi_t = phi_0

    for i_time in range(nt):

        dZ = np.diff(z_edges) #width of nodes

        deltaZ_u = np.diff(Z_P)
        deltaZ_u = np.append(deltaZ_u[0], deltaZ_u)
        
        deltaZ_d = np.diff(Z_P)
        deltaZ_d = np.append(deltaZ_d, deltaZ_d[-1])

        f_u = 1 - (Z_P[:] - z_edges[0:-1]) / deltaZ_u[:]
        f_d = 1 - (z_edges[1:] - Z_P[:]) / deltaZ_d[:]

        #######################################
        # this part is for gas diffusion, which takes a bit more physics
        if airdict!=None: 
            Gamma_Po    = Gamma_P * airdict['por_op'] #This is the diffusivity times the open porosity.

            Gamma_U     = np.append(Gamma_Po[0], Gamma_Po[0: -1] )
            Gamma_D     = np.append(Gamma_Po[1:], Gamma_Po[-1])
            Gamma_u     =  1 / ((1 - f_u) / Gamma_Po + f_u / Gamma_U) #Patankar Eq. 4.11
            Gamma_d     =  1 / ((1 - f_d) / Gamma_Po + f_d / Gamma_D)

            d_eddy_P    = airdict['d_eddy'] * airdict['por_op']
            d_eddy_U    = np.append(d_eddy_P[0], d_eddy_P[0:-1] )
            d_eddy_D    = np.append(d_eddy_P[1:], d_eddy_P[-1])
            d_eddy_u    =  1/ ( (1 - f_u)/d_eddy_P + f_u/d_eddy_U )
            d_eddy_d    =  1/ ( (1 - f_d)/d_eddy_P + f_d/d_eddy_D )
            
            if airdict['gravity']=="off" and airdict['thermal']=="off":
                S_C_0   = 0.0

            elif airdict['gravity']=='on' and airdict['thermal']=='off':
                S_C_0   = (-Gamma_d + Gamma_u) * (airdict['deltaM'] * GRAVITY / (R * airdict['Tz'])) / airdict['dz'] #S_C is independent source term in Patankar

            elif airdict['gravity']=='on' and airdict['thermal']=='on':
                # dTdz    = np.gradient(airdict['Tz'])/airdict['dz']
                dTdz    = np.gradient(airdict['Tz'], Z_P)
                Gamma_del = (Gamma_d-Gamma_u)
                # Gamma_del[Gamma_del<0]=1e-65
                S_C_0   = Gamma_del * (-(airdict['deltaM'] * GRAVITY / (R * airdict['Tz'])) + (airdict['omega'] * dTdz)) / airdict['dz'] # should thermal still work in LIZ? if so use d_eddy+diffu
                # S_C_0[S_C_0<0]=1.e-40
            else:
                print('Error at in solver.py at 119')
                sys.exit()
            
            S_C         = S_C_0 * phi_t
            b_0         = S_C * dZ 

            rho_edges = np.interp(z_edges,Z_P,airdict['rho'])
            
            w_edges = w(airdict, z_edges, rho_edges, Z_P, dZ) # advection term (upward relative motion due to porosity changing)

            w_p = np.interp(Z_P,z_edges,w_edges) # Units m/s
            w_edges[z_edges>airdict['z_co']] = 0.0          
            w_u = w_edges[0:-1]
            w_d = w_edges[1:]

            D_u = ((Gamma_u+d_eddy_u) / deltaZ_u) # Units m/s
            D_d = ((Gamma_d+d_eddy_d) / deltaZ_d)

            F_u =  w_u * airdict['por_op'] # Units m/s
            F_d =  w_d * airdict['por_op']
            
            P_u = F_u / D_u
            P_d = F_d / D_d

            op_ind              = np.where(z_edges<=airdict['z_co'])[0] #indices of all nodes wiht open porosity (shallower than CO)
            op_ind2             = np.where(z_edges<=airdict['z_co']+20)[0] # a bit deeper
            co_ind              = op_ind[-1]
            
            a_U = D_u * A( P_u ) + F_upwind(  F_u )
            a_D = D_d * A( P_d ) + F_upwind( -F_d )
        
            a_P_0 = airdict['por_op'] * dZ / dt
        #######################################
        ### end gas physics portion ###########
        #######################################

        #######################################        
        else: # just for heat, enthalpy, isotope diffusion 
            Gamma_U = np.append(Gamma_P[0], Gamma_P[0: -1] )
            Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])

            Gamma_u =  1 / ((1 - f_u) / Gamma_P + f_u / Gamma_U) # Patankar eq. 4.9
            Gamma_d =  1 / ((1 - f_d) / Gamma_P + f_d / Gamma_D)

            S_C = 0
            S_C = S_C * np.ones(nz_P)

            D_u = (Gamma_u / deltaZ_u)
            D_d = (Gamma_d / deltaZ_d)

            b_0 = S_C * dZ # first term of Patankar eq. 4.41d

            a_U = D_u # Patankar eq. 4.41a,b
            a_D = D_d # Patankar eq. 4.41a,b

            # a_P_0 = dZ / dt
            # a_P_0 = tot_rho * dZ / dt #  (old)
            a_P_0 = c_vol * dZ / dt # (new) Patankar eq. 4.41c
            # a_P_0 = RHO_I * c_firn * dZ / dt
        #######################################

        S_P     = 0.0
        a_P     = a_U + a_D + a_P_0 - S_P*dZ

        bc_u_0  = phi_s # need to pay attention to surface boundary for gas
        bc_type = 1
        bc_u    = np.concatenate(([ bc_u_0], [bc_type]))

        bc_d_0  = 0
        bc_type = 2
        bc_d    = np.concatenate(([ bc_d_0 ], [ bc_type ]))

        b       = b_0 + a_P_0 * phi_t #Patankar 4.41d

        #Upper boundary
        a_P[0]  = 1
        a_U[0]  = 0
        a_D[0]  = 0
        b[0]    = bc_u[0]

        #Down boundary
        a_P[-1] = 1
        a_D[-1] = 0
        a_U[-1] = 1
        b[-1]   = deltaZ_u[-1] * bc_d[0]

        phi_t = solver(a_U, a_D, a_P, b)
        a_P = a_U + a_D + a_P_0

    if airdict!=None:
        return phi_t, w_p
    else:
        return phi_t

###################################
### end transient_solve_TR ########
###################################


def transient_solve_EN(z_edges, Z_P, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s, tot_rho, c_vol, g_liq, deltaH):
    '''
    transient 1-d diffusion finite volume method for enthalpy

    :param z_edges:
    :param Z_P:
    :param nt:
    :param dt:
    :param Gamma_P:
    :param phi_0:
    :param nz_P:
    :param nz_fv:
    :param phi_s:
    :param g_liq

    :return phi_t:

    The source terms S_P and S_C come from the linearization described
    in Voller and Swaminathan, 1991, equations 31 and 32.
    and Voller, Swaminathan, and Thomas, 1990, equation 61
    '''

    phi_t = phi_0.copy()
    g_liq_old = g_liq.copy()
    g_lig_iter = g_liq.copy() 
    # g_ice_old = g_ice.copy() 

    for i_time in range(nt):
        phi_t = phi_0.copy()
        g_lig_iter = g_liq.copy()
        dZ = np.diff(z_edges) #width of nodes

        deltaZ_u = np.diff(Z_P)
        deltaZ_u = np.append(deltaZ_u[0], deltaZ_u)
        
        deltaZ_d = np.diff(Z_P)
        deltaZ_d = np.append(deltaZ_d, deltaZ_d[-1])

        f_u = 1 - (Z_P[:] - z_edges[0:-1]) / deltaZ_u[:]
        f_d = 1 - (z_edges[1:] - Z_P[:]) / deltaZ_d[:]

        Gamma_U = np.append(Gamma_P[0], Gamma_P[0: -1] )
        Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])

        Gamma_u =  1 / ((1 - f_u) / Gamma_P + f_u / Gamma_U) # Patankar eq. 4.9
        Gamma_d =  1 / ((1 - f_d) / Gamma_P + f_d / Gamma_D)

        # beta_P = np.zeros_like(Gamma_P)
        # beta_P[(g_liq>=0) & (g_liq<=1)]= -1.0e9
        # S_P     = beta_P * deltaH * g_liq_old
        # S_C = 1 * S_P * 273.15 + deltaH * g_liq_old - deltaH * g_liq 

        beta_P = np.zeros_like(Gamma_P)
        beta_P[(g_liq>=0) & (g_liq<=1)]= 1.0e12

        S_P = -1 * beta_P * deltaH * g_liq_old
        S_C = 1 * S_P * 273.15 + deltaH * g_liq_old - deltaH * g_liq
        # Vincent Verjans pointed out that Voller et al 1990 use -1*S_P; but it does not
        # converge as quickly and gives the same answer. Not sure how I figured that out.

        D_u = (Gamma_u / deltaZ_u)
        D_d = (Gamma_d / deltaZ_d)

        b_0 = S_C * dZ * dt

        a_U = D_u * dt
        a_D = D_d * dt

        # a_P_0 = dZ / dt
        # a_P_0 = tot_rho * dZ / dt #  (old)
        a_P_0 = c_vol * dZ #/ dt # (new) Patankar eq. 4.41c, this is b_p in Voller (1990; Eq. 30)
        # a_P_0 = RHO_I * c_firn * dZ / dt
            
        a_P     = a_U + a_D + a_P_0 - S_P*dZ*dt

        bc_u_0  = phi_s # need to pay attention for gas
        bc_type = 1
        bc_u    = np.concatenate(([ bc_u_0], [bc_type]))

        bc_d_0  = 0
        bc_type = 2
        bc_d    = np.concatenate(([ bc_d_0 ], [ bc_type ]))

        b       = b_0 + a_P_0 * phi_t

        #Upper boundary
        a_P[0]  = 1
        a_U[0]  = 0
        a_D[0]  = 0
        b[0]    = bc_u[0]

        #Down boundary
        a_P[-1] = 1
        a_D[-1] = 0
        a_U[-1] = 1
        b[-1]   = deltaZ_u[-1] * bc_d[0]

        #####
        phi_t = solver(a_U, a_D, a_P, b)
        #####

        # a_P = a_U + a_D + a_P_0

        g_liq_old = g_liq.copy()
        g_liq = g_liq + a_P * ((phi_t - 273.15) / (dZ * deltaH))
        g_liq[g_liq>1] = 1
        g_liq[g_liq<0] = 0

    return phi_t, g_liq

###################################
### end transient_solve_EN ########
###################################

'''
Functions below are for firn air
Works, but consider to be in beta
'''
def w(airdict, z_edges, rho_edges, Z_P, dZ): 
    '''
    Function for downward advection of air and also calculates total air content. 
    '''
    if airdict['advection_type']=='Darcy':
        por_op_edges=np.interp(z_edges,airdict['z'],airdict['por_op'])
        T_edges = np.interp(z_edges,airdict['z'],airdict['Tz'])
        p_star = por_op_edges * np.exp(M_AIR *GRAVITY*z_edges/(R*T_edges))
        dPdz = np.gradient(airdict['air_pressure'],airdict['z'])
        dPdz_edges=np.interp(z_edges,airdict['z'],dPdz)

        # perm = 10.0**(-7.29) * por_op_edges**3.71 # Adolph and Albert, 2014, eq. 5, units m^2
        # perm = 10.0**(-7.7) * por_op_edges**3.4 #Freitag, 2002 
        perm = 10.0**(-7.7) * p_star**3.4 #Freitag, 2002 
        visc = 1.5e-5 #kg m^-1 s^-1, dynamic viscosity, source?
        flux = -1.0 * perm / visc * dPdz_edges # units m/s
        # w_ad = flux / airdict['dt']  / por_op_edges # where did I get this?
        w_ad = flux / p_star / airdict['dt']
        # w_ad = flux / por_op_edges / airdict['dt']

    elif airdict['advection_type']=='Christo':
        por_tot_edges       = np.interp(z_edges,Z_P,airdict['por_tot'])
        por_cl_edges        = np.interp(z_edges,Z_P,airdict['por_cl'])
        por_op_edges        = np.interp(z_edges,Z_P,airdict['por_op'])
        w_firn_edges        = np.interp(z_edges,Z_P,airdict['w_firn']) # units m/s
        T_edges             = np.interp(z_edges,Z_P,airdict['Tz'])
        p_star              = por_op_edges * np.exp(M_AIR *GRAVITY*z_edges/(R*T_edges))
        dscl                = np.gradient(por_cl_edges,z_edges)
        C                   = np.exp(M_AIR*GRAVITY*z_edges/(R*T_edges))

        op_ind              = np.where(z_edges<=airdict['z_co'])[0] #indices of all nodes wiht open porosity (shallower than CO)
        op_ind2             = np.where(z_edges<=airdict['z_co']+20)[0] # a bit deeper
        co_ind              = op_ind[-1] 
        cl_ind1             = np.where(z_edges>airdict['z_co'])[0] #closed indices
        cl_ind              = np.intersect1d(cl_ind1,op_ind2)

        # print('depth co_ind',z_edges[co_ind])

        Xi                  = np.zeros((len(op_ind2),len(op_ind2)))
        Xi_up               = por_op_edges[op_ind2]/np.reshape(por_op_edges[op_ind2], (-1,1))
        Xi_down             = (1 + np.log( np.reshape(w_firn_edges[op_ind2], (-1,1))/ w_firn_edges[op_ind2] ))
        Xi                  = Xi_up / Xi_down # Equation 5.10 in Christo's thesis; Xi[i,j] is the pressure increase (ratio) for bubbles at depth[i] that were trapped at depth[j]

        integral_matrix     = (Xi.T*dscl[op_ind2]*C[op_ind2]).T 
        integral_matrix_sum = integral_matrix.sum(axis=1)

        p_ratio_t           = np.zeros_like(op_ind2)
        p_ratio             = np.zeros_like(z_edges)
        p_ratio[op_ind]         = integral_matrix_sum[op_ind]   #5.11
        p_ratio[cl_ind]         = p_ratio[co_ind]*Xi[cl_ind, co_ind] # 5.12
        p_ratio[cl_ind[-1]+1:]  = p_ratio[cl_ind[-1]]

        flux                = w_firn_edges[co_ind-1] * p_ratio[co_ind-1] * por_cl_edges[co_ind-1]

        velocity            = np.minimum(w_firn_edges ,((flux + 1e-10 - w_firn_edges * p_ratio * por_cl_edges) / ((por_op_edges + 1e-10 * C))))

        # velocity            = (flux + 1e-10 - w_firn_edges * p_ratio * por_cl_edges) / ((por_op_edges + 1e-10 * C))
        # velocity = flux / p_star# / airdict['dt']

        # w_ad = velocity
        w_ad              = (velocity - w_firn_edges)
        
        # w_ad[w_ad>0] = 0

        # w_ad[co_ind:+1] = 0

        # veldiff = velocity
        # print('dt',airdict['dt'])
        # print('flux',flux)
        # print('vel',velocity[co_ind-5:co_ind+5])
        # print('w_firn',w_firn_edges[co_ind-5:co_ind+5])
        # print('w_ad',w_ad[co_ind-5:co_ind+5])
        # print('w_ad',w_ad)
        # input() 

    elif airdict['advection_type']=='zero':
        w_ad = np.zeros_like(rho_edges)

    return w_ad


def A(P):
    '''Power-law scheme, Patankar eq. 5.34'''
    A = np.maximum( (1 - 0.1 * np.abs( P ) )**5, np.zeros(np.size(P) ) )
    return A    

def F_upwind(F): 
    ''' Upwinding scheme '''
    F_upwind = np.maximum( F, 0 )
    return F_upwind


# def w(z_edges,rho_edges,por_op,T,p_a,por_tot,por_cl,Z_P,dz, w_firn): # Function for downward advection of air and also calculates total air content. 
    
#     por_tot_edges=np.interp(z_edges,Z_P,por_tot)
#     por_cl_edges=np.interp(z_edges,Z_P,por_cl)
#     por_op_edges=np.interp(z_edges,Z_P,por_op)
#     teller_co=np.argmax(por_cl_edges)
#     # w_firn_edges=Accu*rho_i/rho_edges #Check this - is there a better way?
#     w_firn_edges=np.interp(z_edges,Z_P,w_firn)
    
#     # if ad_method=='ice_vel':
#     #     w_ad=w_firn_edges
#     #     trapped = 0.0
#     #     bubble_pres = np.zeros_like(z_edges)
    

#     ### Christo's Method from his thesis (chapter 5). This (maybe) could be vectorized to speed it up.
    
#     bubble_pres = np.zeros_like(z_edges)
#     # print(len(np.diff(por_cl)))
#     # print(len(dz))
#     # dscl = np.append(0, np.diff(por_cl)/dz)
#     dscl = np.append(0, np.gradient(por_cl,dz))
#     T_edges = np.interp(z_edges,Z_P,T)
#     C=np.exp(M_AIR*GRAVITY*z_edges/(R*T_edges))
#     strain = np.gradient(np.log(w_firn),dz)
#     s=por_op_edges+por_cl_edges
    
#     for teller1 in range (0,teller_co+1): 
#         integral = np.zeros(teller1+1)
#         integral2 = np.zeros(teller1+1)
        
#         for teller2 in range(0,teller1+1):
#             # integral[teller2] = dscl[teller2]*C[teller2]*(s[teller2]/s[teller1])/(1+scipy.integrate.trapz(strain[teller2:teller1+1],dz)) #need to get this indexing correct 6/19/14: I think it is fine.
#             integral[teller2] = dscl[teller2]*C[teller2]*(s[teller2]/s[teller1])/(1+scipy.integrate.trapz(strain[teller2:teller1+1],z_edges[teller2:teller1+1])) #need to get this indexing correct 6/19/14: I think it is fine.
#             if dscl[teller2]==0:
#                 dscl[teller2]=1e-14
#             integral2[teller2] = dscl[teller2]
            
#         bubble_pres[teller1] = (np.mean(dz)*np.sum(integral))/(np.mean(dz)*np.sum(integral2))
    
#     bubble_pres[teller_co+1:] = bubble_pres[teller_co]*(s[teller_co]/s[teller_co+1:])/(w_firn_edges[teller_co+1:]/w_firn_edges[teller_co])
    
#     bubble_pres[0] = 1
#     #print 'bubble pressure = %s' % bubble_pres
    
#     flux= w_firn_edges[teller_co]*bubble_pres[teller_co]*por_cl[teller_co]

#     velocity = np.minimum(w_firn_edges ,((flux+(1e-10)-w_firn_edges*bubble_pres*por_cl_edges)/((por_op_edges+1e-10)*C)))
#     #velocity = velocity * 2   
#     w_ad=velocity


#     return w_ad #, bubble_pres