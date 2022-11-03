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

        #######################################
        ### Boundary conditions:
        ### type 1 is a specified value, type 2 is a specified gradient
        ### (units for gradient are degrees/meter)
        bc_u_0  = phi_s # need to pay attention to surface boundary for gas
        bc_type_u = 1
        bc_u    = np.concatenate(([ bc_u_0], [bc_type_u]))

        bc_d_0  = 0
        bc_type_d = 2
        # bc_d_0  = 273.149
        # bc_type_d = 1
        bc_d    = np.concatenate(([ bc_d_0 ], [ bc_type_d ]))
        #########################################

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
        # print(phi_t[-1])
        # input('waiting')
        a_P = a_U + a_D + a_P_0

    if airdict!=None:
        return phi_t, w_p
    else:
        return phi_t

###################################
### end transient_solve_TR ########
###################################

def transient_solve_EN_h(z_edges, Z_P, nt, dt, Gamma_P, Tc, nz_P, nz_fv, phi_s, mix_rho, c_vol, LWC, mass_sol, dz, ICT,iii=0):
    '''
    transient 1-d diffusion finite volume method for enthalpy

    :param z_edges:
    :param Z_P:
    :param nt: number of iterations; depricated (now uses while loop)
    :param dt: time step size
    :param Gamma_P:
    :param phi_0:
    :param nz_P:
    :param nz_fv:
    :param phi_s:
    :param g_liq
    :param c_vol: [J/m3/K] 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)

    :return phi_t:

    The source terms S_P and S_C come from the linearization described
    in Voller and Swaminathan, 1991, equations 31 and 32.
    and Voller, Swaminathan, and Thomas, 1990, equation 61
    '''


    # h = 
    # phi_t = Tc.copy()
    # phi_t_old = phi_t.copy()
    LWC_in = LWC.copy()

    vol_S       = mass_sol / RHO_I     # volume_Solid, i.e. volume of the ice (solid) portion of each control volume
    vol_SL      = vol_S + LWC    # volume of solid and liquid in each control volume
    mass_liq    = LWC * RHO_W_KGM  # mass of liquid water in each control
    mass_tot    = mass_liq + mass_sol # total mass (solid +liquid) of each control
    rho_liq_eff = mass_liq / dz      # effective density of the liquid portion
    mix_rho     = (mass_sol + mass_liq) / dz # mixture, or total density of volume (solid plus liquid), see Aschwanden (2012)
    g_liq       = LWC / dz    #  use liquid volume fraction of total volume of the control, which will net us the enthalpy/volume
    # g_liq       = LWC / vol_SL  # alternatively, liquid volume fraction (of the material portion, porosity ignored), (I don't think this one is correct. See code at end of while loop to swap if you want to use this)
   
    rho_firn = mass_sol/dz

    h = ((rho_firn * Tc * CP_I) + LF_I * rho_liq_eff) * dz #[J]

    phi_t = h.copy()
    phi_t_old = phi_t.copy()

    g_liq_old = g_liq.copy()

    itercheck = 0.9
    count = 0

    # while itercheck>ICT:
    phi_iter = phi_t.copy()
    g_liq_iter = g_liq.copy() #unitless
    # deltaH_l = rho_liq_eff * LF_I # [kg/m3 * J/kg = J/m3]
    # deltaH = deltaH_l # deltaH is zero anywhere that has LWC = 0
    # phi_t = phi_0.copy() # T is in C
    dZ = np.diff(z_edges) #width of nodes

    deltaZ_u = np.diff(Z_P) # [m]
    deltaZ_u = np.append(deltaZ_u[0], deltaZ_u)
    
    deltaZ_d = np.diff(Z_P)
    deltaZ_d = np.append(deltaZ_d, deltaZ_d[-1])

    f_u = 1 - (Z_P[:] - z_edges[0:-1]) / deltaZ_u[:] # unitless
    f_d = 1 - (z_edges[1:] - Z_P[:]) / deltaZ_d[:]

    ### Gamma has units J/s/m/K (W/m/K)
    Gamma_U = np.append(Gamma_P[0], Gamma_P[0: -1])
    Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])

    Gamma_u =  1 / ((1 - f_u) / Gamma_P + f_u / Gamma_U) # Patankar eq. 4.9
    Gamma_d =  1 / ((1 - f_d) / Gamma_P + f_d / Gamma_D)

    # #### eqs. 31, 32 from Voller 1991
    # dFdT = np.zeros_like(Gamma_P) # [1/K] dFdT is slope of the liquid fraction/temperature curve. Large at T=0C, i.e. fraction increases rapidly at T=0
    # dFdT[(g_liq>=0) & (g_liq<=1)]= 1.0e10
    # # ### Finv has units of temperature: T=F^-1(g_liq). When g_liq>0, T=0C, and Finv=0. For g_liq<=0, T<=0 - but
    # # ### the equation for S_C contains deltaH multiplying Finv. deltaH is zero for nodes with T<0, so we can 
    # # ### just set Finv=0 at those nodes.
    # Finv = np.zeros_like(dFdT)


    S_P = np.zeros_like(Gamma_P)
    S_C = np.zeros_like(Gamma_P)
    # S_C = LF_I * (rho_liq_eff) 

    D_u = (Gamma_u / deltaZ_u) # [W/m2/K]
    D_d = (Gamma_d / deltaZ_d)
    
    a_U   = D_u        #* dt # [W/m2/K]
    a_D   = D_d        #* dt # [W/m2/K]
    a_P_0 = c_vol * dZ / dt # [W/m2/K] (new) Patankar eq. 4.41c, this is b_p in Voller (1990; Eq. 30)           
    a_P   = a_U + a_D + a_P_0 - S_P * dZ #* dt # check the multiply on the S_P
    
    b_0   = S_C  * dZ #* dt # [W/m2]
    b     = b_0 + a_P_0 * phi_t_old # is this one right? change 22/10/27

    ### Boundary conditions:
    ### type 1 is a specified value, type 2 is a specified gradient
    ### (units for gradient are degrees/meter)
    # bc_u_0  = phi_s 
    bc_u_0  = h[0] 
    bc_type_u = 1
    bc_u    = np.concatenate(([ bc_u_0], [bc_type_u]))

    bc_d_0  = 0
    bc_type_d = 1
    # bc_d_0  = 0
    # bc_type_d = 2
    bc_d    = np.concatenate(([ bc_d_0 ], [ bc_type_d ]))

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

    # h = (rho_firn * Tc * CP_I)

    # Delta_gl = a_P * (phi_t - )

    # deltaH_new = deltaH + 0.95 * a_P/a_P_0*c_vol*phi_iter
    # deltaH_new[deltaH_new<0] = 0
    # deltaH_new[deltaH_new>1] = 1
    # rho_liq_eff = deltaH_new / LF_I
    # mass_liq = rho_liq_eff*dz
    # LWC = mass_liq/RHO_W_KGM
    # g_liq = LWC/dz

    # LWC = g_liq*dz
    # mass_liq = LWC * RHO_W_KGM
    # rho_liq_eff = mass_liq / dz #dz or dZ? Only a very small difference.


    T_new = np.zeros_like(LWC)
    T_new[LWC==0] = phi_t[LWC==0] / (rho_firn[LWC==0] * CP_I *dz[LWC==0])


    delH = np.zeros_like(LWC)
    g_liq_out = np.zeros_like(LWC)
    delH[LWC>0] = phi_t_old[LWC>0] - phi_t_old[LWC>0]
    E_to_freeze = LF_I * LWC * RHO_W_KGM #energy to freeze all water
    iLWC = np.where(LWC>0)[0]
    for jj in iLWC:
        # print('jj',jj)
        if delH[jj]>E_to_freeze[jj]: #all refrozen
            # g_liq_out[jj] = 0
            T_new[jj] = -1*(delH[jj]-E_to_freeze[jj])/(mass_tot[jj]*CP_I)
            # print('option1: ',T_new[jj])
            mass_liq[jj] = 0
        else:
            T_new[jj] = 0.0
            refrozen = delH[jj]/LF_I
            mass_liq[jj] = mass_liq[jj]-refrozen
            # LWC[]
    
    LWC = mass_liq/RHO_W_KGM
    g_liq_out = LWC/dz

    iterdiff=0
    print('iii',iii)
    print(T_new[0:10])
    print(LWC[0:10])
    print(g_liq_out[0:10])
    print(g_liq_old[0:10])
    input('paused.')
    return T_new, g_liq_out, count, iterdiff

###################################
### end transient_solve_EN ########
###################################


def transient_solve_EN_save(z_edges, Z_P, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s, mix_rho, c_vol, LWC, mass_sol, dz, ICT,iii=0):
    '''
    transient 1-d diffusion finite volume method for enthalpy

    :param z_edges:
    :param Z_P:
    :param nt: number of iterations; depricated (now uses while loop)
    :param dt: time step size
    :param Gamma_P:
    :param phi_0:
    :param nz_P:
    :param nz_fv:
    :param phi_s:
    :param g_liq
    :param c_vol: [J/m3/K] 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)

    :return phi_t:

    The source terms S_P and S_C come from the linearization described
    in Voller and Swaminathan, 1991, equations 31 and 32.
    and Voller, Swaminathan, and Thomas, 1990, equation 61
    '''

    # phi_t = phi_0.copy()
    LWC_old = LWC.copy()

    vol_S       = mass_sol / RHO_I     # volume_Solid, i.e. volume of the ice (solid) portion of each control volume
    vol_SL      = vol_S + LWC    # volume of solid and liquid in each control volume
    mass_liq    = LWC * RHO_W_KGM  # mass of liquid water in each control
    mass_tot    = mass_liq + mass_sol # total mass (solid +liquid) of each control
    # rho_liq_eff = mass_liq / dz      # effective density of the liquid portion
    rho_liq_eff = RHO_W_KGM*dz #new 10/31
    mix_rho     = (mass_sol + mass_liq) / dz # mixture, or total density of volume (solid plus liquid), see Aschwanden (2012)
    g_liq       = LWC / dz    #  use liquid volume fraction of total volume of the control, which will net us the enthalpy/volume
    g_sol       = vol_S / dz
    # g_liq       = LWC / vol_SL  # alternatively, liquid volume fraction (of the material portion, porosity ignored), (I don't think this one is correct. See code at end of while loop to swap if you want to use this)
 
    
    H_L_liq = RHO_W_KGM*LF_I #volumetric latent enthalpy
    phi_t = phi_0*CP_I*RHO_I*g_sol #sensible enthalpy
    bigH = phi_t + H_L_liq*g_liq #total enthalpy, voller 1990b eq 4a

    phi_t_old = phi_t.copy()
    g_liq_old = g_liq.copy()

    itercheck = 0.9
    count = 0

    # deltaH = deltaH_l = RHO_W_KGM*LF_I # [J/m3]

    # deltaH = -g_sol*RHO_I*CP_I*phi_t + LF_I*RHO_W_KGM*g_liq
    
    # g_liq_iter = g_liq

    # print('iii:',iii)
    # print('max T in:', np.max(phi_t))

    H_tot = phi_t + H_L_liq*g_liq

    while itercheck>ICT:
        H_tot_iter = H_tot.copy()
        phi_iter = phi_t.copy()
        g_liq_iter = g_liq.copy() #unitless
        # deltaH_l = rho_liq_eff * LF_I # [kg/m3 * J/kg = J/m3]
        # deltaH_l = rho_liq_eff * LF_I # [kg/m2 * J/kg = J/m2]
        # deltaH = deltaH_l # deltaH is zero anywhere that has LWC = 0
        # phi_t = phi_0.copy() # T is in C


        dZ = np.diff(z_edges) #width of nodes

        deltaZ_u = np.diff(Z_P) # [m]
        deltaZ_u = np.append(deltaZ_u[0], deltaZ_u)
        
        deltaZ_d = np.diff(Z_P)
        deltaZ_d = np.append(deltaZ_d, deltaZ_d[-1])

        f_u = 1 - (Z_P[:] - z_edges[0:-1]) / deltaZ_u[:] # unitless
        f_d = 1 - (z_edges[1:] - Z_P[:]) / deltaZ_d[:]

        ### Gamma has units J/s/m/K (W/m/K)
        Gamma_U = np.append(Gamma_P[0], Gamma_P[0: -1])
        Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])

        Gamma_u =  1 / ((1 - f_u) / Gamma_P + f_u / Gamma_U) # Patankar eq. 4.9
        Gamma_d =  1 / ((1 - f_d) / Gamma_P + f_d / Gamma_D)

        # #### eqs. 31, 32 from Voller 1991
        dFdT = np.zeros_like(Gamma_P) # [1/K] dFdT is slope of the liquid fraction/temperature curve. Large at T=0C, i.e. fraction increases rapidly at T=0
        dFdT[(g_liq>=0) & (g_liq<=1)]= 1.0e10
        # ### Finv has units of temperature: T=F^-1(g_liq). When g_liq>0, T=0C, and Finv=0. For g_liq<=0, T<=0 - but
        # ### the equation for S_C contains deltaH multiplying Finv. deltaH is zero for nodes with T<0, so we can 
        # ### just set Finv=0 at those nodes.
        Finv = np.zeros_like(dFdT)

        # method='Voller'
        # if method=='Voller':
        #     # deltaH is [J/m3]
        # S_P = (-1 * deltaH * dFdT) / dt #[J/m3/K/s = W/m3/K]
        # S_P = (deltaH * dFdT) / dt #[J/m3/K/s = W/m3/K]
        # S_C = -1*(deltaH * (g_liq_old - g_liq_iter) + deltaH*dFdT*Finv) / dt # [J/m3/s = W/m3]
        # ### end 1991

        ### S_C is independent
        S_P = np.zeros_like(Gamma_P)
        S_C = H_L_liq  * (g_liq_old - g_liq_iter) / dt # J/m3/s = W/m3 

        D_u = (Gamma_u / deltaZ_u) # [W/m2/K]
        D_d = (Gamma_d / deltaZ_d)
        
        a_U   = D_u        #* dt # [W/m2/K]
        a_D   = D_d        #* dt # [W/m2/K]
        a_P_0 = c_vol * dZ / dt # [W/m2/K] (new) Patankar eq. 4.41c, this is b_p in Voller (1990; Eq. 30)           
        a_P   = a_U + a_D + a_P_0 - S_P * dZ #* dt # check the multiply on the S_P

        # a_P[g_liq_old>0] = 1e15
        
        b_0   = S_C  * dZ #* dt # [W/m2]
        b     = b_0 + a_P_0 * phi_t_old # is this one right? change 22/10/27
        # b       = b_0 + a_P_0 * phi_t

        # print(iii,count)
        # print('S_C:',np.max(b_0))
        # print('a_P_0 * phi_t_old:',np.max(a_P_0 * phi_t_old))
        # print('a_U', a_U)
        # input('paused.')

        ### Boundary conditions:
        ### type 1 is a specified value, type 2 is a specified gradient
        ### (units for gradient are degrees/meter)
        bc_u_0  = phi_s
        bc_type_u = 1
        bc_u    = np.concatenate(([ bc_u_0], [bc_type_u]))

        bc_d_0  = 0
        bc_type_d = 1
        # bc_d_0  = 0
        # bc_type_d = 2
        bc_d    = np.concatenate(([ bc_d_0 ], [ bc_type_d ]))

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
        # print(count)
        # print(np.max(phi_t))
        # input('waiting.')
        #####

        # deltaH_new = -g_sol*RHO_I*CP_I*phi_t + LF_I*RHO_W_KGM*g_liq

        phi_iter = phi_t

        UCF = 1
        # dh_sens = CP_I * RHO_I * phi_t

        # Delta_g_liq=a_P*phi_t/(dz*deltaH)*dt
        # Delta_g_liq = np.zeros_like(LWC)
        # Delta_g_liq[g_liq_old>0] = (a_P_0 * phi_t[g_liq_old>0]) / (dz[g_liq_old>0] * deltaH)

        Delta_g_liq = a_P_0 * phi_t / (H_L_liq) 
     
        g_liq = g_liq + UCF * Delta_g_liq
        # g_liq = g_liq_iter + (a_P * (phi_t) / (RHO_W_KGM*LF_I*dz)) 
        g_liq[g_liq<0] = 0
        g_liq[g_liq>1] = 1

        # Delta

        H_tot = phi_t + H_L_liq*g_liq

        # deltaLWC = (g_liq - g_liq_iter)/dz
        # delta_vol_S = -1*deltaLWC/0.917
        # g_sol = (vol_S + delta_vol_S) / dz

        # deltaH = deltaH_new
        # deltaH = -g_sol*RHO_I*CP_I*phi_t + LF_I*RHO_W_KGM*g_liq

        # deltaH_new = deltaH + 0.95 * a_P/a_P_0*c_vol*phi_iter
        # deltaH_new[deltaH_new<0] = 0
        # deltaH_new[deltaH_new>1] = 1
        # rho_liq_eff = deltaH_new / LF_I
        # mass_liq = rho_liq_eff*dz
        # LWC = mass_liq/RHO_W_KGM
        # g_liq = LWC/dz

        # LWC = g_liq*dz
        # mass_liq = LWC * RHO_W_KGM
        # rho_liq_eff = mass_liq / dz #dz or dZ? Only a very small difference.

        # iterdiff = (np.sum(g_liq_iter) - np.sum(g_liq))
        iterdiff = np.abs(np.sum(H_tot_iter/(CP_I*RHO_I*g_sol)) - np.sum(H_tot/(CP_I*RHO_I*g_sol)))
        if iterdiff<1e-4:
            itercheck = 0
        else:
            itercheck = iterdiff

        # if iterdiff==0:
        #     itercheck = 0
        #     # print(f'g_liq_iter: {np.sum(g_liq_iter)}')
        #     # print(f'g_liq: {np.sum(g_liq)}')
        #     # print('break!')
        #     break
        # else:
        #     itercheck = np.abs( iterdiff/np.sum(g_liq_iter))
        count += 1
        if count==100:
            if ICT == 0:
                ICT = 1e-12
            else:
                pass
        if ((count==200) and (ICT==1e-12)):
            # if ICT > 1e-12:
            #     pass
            # else:
            ICT = ICT * 10
        if count>1000:
            print('breaking loop')
            print('iii', iii)
            print('itercheck',itercheck)
            print('iterdiff', iterdiff)
            print('np.sum(g_liq_iter)', np.sum(g_liq_iter))
            break
    ### END ITERATION

    # print(count)
    # print('###############')
    # input('waiting')
    # print('count',count)

    phi_t_out = H_tot/(CP_I*RHO_I*g_sol)
    phi_t_out[g_liq>0] = 0
    # print('max out:', np.max(phi_t_out))
    return phi_t_out, g_liq, count, iterdiff

###################################
### end transient_solve_EN ########
###################################

def transient_solve_EN(z_edges, Z_P, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s, mix_rho, c_vol, LWC, mass_sol, dz, ICT,iii=0):
    '''
    transient 1-d diffusion finite volume method for enthalpy
    Diffuse sensible enthalpy version.

    :param z_edges:
    :param Z_P:
    :param nt: number of iterations; depricated (now uses while loop)
    :param dt: time step size
    :param Gamma_P:
    :param phi_0: temperature profile [C]
    :param nz_P:
    :param nz_fv:
    :param phi_s: surface temperature [C]
    :param g_liq
    :param c_vol: [J/m3/K] 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)

    :return phi_t:

    The source terms S_P and S_C come from the linearization described
    in Voller and Swaminathan, 1991, equations 31 and 32.
    and Voller, Swaminathan, and Thomas, 1990, equation 61
    '''

    # phi_t = phi_0.copy()
    LWC_old = LWC.copy()
    phi_in = phi_0.copy()

    vol_S       = mass_sol / RHO_I     # volume_Solid, i.e. volume of the ice (solid) portion of each control volume
    vol_SL      = vol_S + LWC    # volume of solid and liquid in each control volume
    mass_liq    = LWC * RHO_W_KGM  # mass of liquid water in each control
    mass_tot    = mass_liq + mass_sol # total mass (solid +liquid) of each control
    # rho_liq_eff = mass_liq / dz      # effective density of the liquid portion
    rho_liq_eff = RHO_W_KGM*dz #new 10/31
    # mix_rho     = (mass_sol + mass_liq) / dz # mixture, or total density of volume (solid plus liquid), see Aschwanden (2012)
    g_liq       = LWC / dz    #  use liquid volume fraction of total volume of the control, which will net us the enthalpy/volume
    g_sol       = vol_S / dz

    g_liq_in = g_liq.copy()
      
    H_L_liq = RHO_W_KGM*LF_I #volumetric latent enthalpy
    phi_t = phi_0*CP_I*RHO_I*g_sol #sensible enthalpy
    
    phi_t_old = phi_t.copy() #initial enthalpy
    g_liq_old = g_liq.copy() # initial liquid fraction

    itercheck = 0.9
    count = 0

    H_latent = H_L_liq*g_liq
    H_latent_old = H_latent.copy()
    H_tot = phi_t + H_latent  #total enthalpy, voller 1990b eq 4a

    while itercheck>ICT:

        H_tot_iter = H_tot.copy()
        phi_iter = phi_t.copy()
        g_liq_iter = g_liq.copy() #unitless
        H_latent_iter = H_latent.copy()

        dZ = np.diff(z_edges) #width of nodes

        deltaZ_u = np.diff(Z_P) # [m]
        deltaZ_u = np.append(deltaZ_u[0], deltaZ_u)
        
        deltaZ_d = np.diff(Z_P)
        deltaZ_d = np.append(deltaZ_d, deltaZ_d[-1])

        f_u = 1 - (Z_P[:] - z_edges[0:-1]) / deltaZ_u[:] # unitless
        f_d = 1 - (z_edges[1:] - Z_P[:]) / deltaZ_d[:]

        ### Gamma has units J/s/m/K (W/m/K)
        Gamma_U = np.append(Gamma_P[0], Gamma_P[0: -1])
        Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])

        Gamma_u =  1 / ((1 - f_u) / Gamma_P + f_u / Gamma_U) # Patankar eq. 4.9
        Gamma_d =  1 / ((1 - f_d) / Gamma_P + f_d / Gamma_D)

        ### S_C is independent
        S_P = np.zeros_like(Gamma_P)
        S_C = H_L_liq  * (g_liq_old - g_liq_iter) / dt # J/m3/s = W/m3 
        # S_C = H_L_liq  * g_liq_old / dt # J/m3/s = W/m3 

        D_u = (Gamma_u / deltaZ_u) # [W/m2/K]
        D_d = (Gamma_d / deltaZ_d)
        
        a_U   = D_u        #* dt # [W/m2/K]
        a_D   = D_d        #* dt # [W/m2/K]
        a_P_0 = c_vol * dZ / dt # [W/m2/K] (new) Patankar eq. 4.41c, this is b_p in Voller (1990; Eq. 30)           
        a_P   = a_U + a_D + a_P_0 - S_P * dZ #* dt # check the multiply on the S_P

        b_0   = S_C  * dZ #* dt # [W/m2]
        b     = b_0 + a_P_0 * phi_t_old # is this one right? change 22/10/27

        ### Boundary conditions:
        ### type 1 is a specified value, type 2 is a specified gradient
        ### (units for gradient are degrees/meter)
        # bc_u_0  = phi_s
        bc_u_0 = phi_t_old[0]
        bc_type_u = 1
        bc_u    = np.concatenate(([ bc_u_0], [bc_type_u]))

        bc_d_0  = 0
        bc_type_d = 1
        # bc_d_0  = 0
        # bc_type_d = 2
        bc_d    = np.concatenate(([ bc_d_0 ], [ bc_type_d ]))

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

        delta_h = phi_t - phi_iter #this will always be negative for layers with LWC

        ndh = -1*delta_h 

        ### everything refreezes
        cond1 = ((ndh>=H_latent) & (g_liq>0)) #layers where the energy change is larger than the energy required to refreeze, and where there is water
        phi_t[cond1] = delta_h[cond1] + H_latent[cond1] #sensible enthalpy will be the difference between the latent heat and the calculated change in energy
        g_liq[cond1] = 0 #no liquid left

        ### partial refreezing
        cond2 = ((ndh<H_latent) & (g_liq>0))
        phi_t[cond2] = 0
        g_liq[cond2] = (H_latent[cond2] + delta_h[cond2])/H_L_liq

        delta_g_liq = g_liq - g_liq_old # change in liquid fraction. Should always be negative.

        g_sol = g_sol + -1*delta_g_liq/0.917
        H_tot = phi_t + H_L_liq*g_liq

        iterdiff = (np.sum(g_liq_iter) - np.sum(g_liq))

        if iterdiff==0:
            itercheck = 0
            break
        else:
            itercheck = np.abs( iterdiff/np.sum(g_liq_iter))


        count += 1
        if count==100:
            if ICT == 0:
                ICT = 1e-12
            else:
                pass
        if ((count==200) and (ICT==1e-12)):
            # if ICT > 1e-12:
            #     pass
            # else:
            ICT = ICT * 10
        if count>1000:
            print('breaking loop')
            print('iii', iii)
            print('itercheck',itercheck)
            print('iterdiff', iterdiff)
            print('np.sum(g_liq_iter)', np.sum(g_liq_iter))
            break
        ### END ITERATION

    # iterdiff=0
    phi_t_out = H_tot/(CP_I*RHO_I*g_sol)
    phi_t_out[g_liq>0] = 0
    phi_t_out[(phi_t_out>0)] = 0 

    bcw = ((phi_t_out<0)&(g_liq>0))
    # bcw = ((g_liq - g_liq_in)>0)
    if np.any(bcw):
        icw = np.where(bcw)[0]
        print('##########')
        print('solver returning cold and wet firn.')
        print('count', count)
        print('itercheck',itercheck)
        print('phi_in',phi_in[icw])
        print('g_liq_in',g_liq_in[icw])
        print('phi_t_out',phi_t_out[icw])
        print('g_liq',g_liq[icw])
        print('g_liq diff: ',(g_liq - g_liq_in)[icw])
        print('depths', Z_P[icw])
        input("waiting. (solver.py)")

    return phi_t_out, g_liq, count, iterdiff

###################################
### end transient_solve_EN ########
###################################


def transient_solve_EN_noloop(z_edges, Z_P, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s, mix_rho, c_vol, LWC, mass_sol, dz, ICT,iii=0):
    '''
    transient 1-d diffusion finite volume method for enthalpy
    Diffuse sensible enthalpy version.

    :param z_edges:
    :param Z_P:
    :param nt: number of iterations; depricated (now uses while loop)
    :param dt: time step size
    :param Gamma_P:
    :param phi_0:
    :param nz_P:
    :param nz_fv:
    :param phi_s:
    :param g_liq
    :param c_vol: [J/m3/K] 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)

    :return phi_t:

    The source terms S_P and S_C come from the linearization described
    in Voller and Swaminathan, 1991, equations 31 and 32.
    and Voller, Swaminathan, and Thomas, 1990, equation 61
    '''

    # phi_t = phi_0.copy()
    LWC_old = LWC.copy()

    vol_S       = mass_sol / RHO_I     # volume_Solid, i.e. volume of the ice (solid) portion of each control volume
    vol_SL      = vol_S + LWC    # volume of solid and liquid in each control volume
    mass_liq    = LWC * RHO_W_KGM  # mass of liquid water in each control
    mass_tot    = mass_liq + mass_sol # total mass (solid +liquid) of each control
    # rho_liq_eff = mass_liq / dz      # effective density of the liquid portion
    rho_liq_eff = RHO_W_KGM*dz #new 10/31
    # mix_rho     = (mass_sol + mass_liq) / dz # mixture, or total density of volume (solid plus liquid), see Aschwanden (2012)
    g_liq       = LWC / dz    #  use liquid volume fraction of total volume of the control, which will net us the enthalpy/volume
    g_sol       = vol_S / dz
      
    H_L_liq = RHO_W_KGM*LF_I #volumetric latent enthalpy
    phi_t = phi_0*CP_I*RHO_I*g_sol #sensible enthalpy
    
    phi_t_old = phi_t.copy() #initial enthalpy
    g_liq_old = g_liq.copy() # initial liquid fraction

    itercheck = 0.9
    count = 0

    H_latent = H_L_liq*g_liq
    H_latent_old = H_latent.copy()
    H_tot = phi_t + H_latent  #total enthalpy, voller 1990b eq 4a

    # while itercheck>ICT:

    H_tot_iter = H_tot.copy()
    phi_iter = phi_t.copy()
    g_liq_iter = g_liq.copy() #unitless
    H_latent_iter = H_latent.copy()

    dZ = np.diff(z_edges) #width of nodes

    deltaZ_u = np.diff(Z_P) # [m]
    deltaZ_u = np.append(deltaZ_u[0], deltaZ_u)
    
    deltaZ_d = np.diff(Z_P)
    deltaZ_d = np.append(deltaZ_d, deltaZ_d[-1])

    f_u = 1 - (Z_P[:] - z_edges[0:-1]) / deltaZ_u[:] # unitless
    f_d = 1 - (z_edges[1:] - Z_P[:]) / deltaZ_d[:]

    ### Gamma has units J/s/m/K (W/m/K)
    Gamma_U = np.append(Gamma_P[0], Gamma_P[0: -1])
    Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])

    Gamma_u =  1 / ((1 - f_u) / Gamma_P + f_u / Gamma_U) # Patankar eq. 4.9
    Gamma_d =  1 / ((1 - f_d) / Gamma_P + f_d / Gamma_D)

    ### S_C is independent
    S_P = np.zeros_like(Gamma_P)
    # S_C = H_L_liq  * (g_liq_old - g_liq_iter) / dt # J/m3/s = W/m3 
    S_C = H_L_liq  * g_liq_old / dt # J/m3/s = W/m3 

    D_u = (Gamma_u / deltaZ_u) # [W/m2/K]
    D_d = (Gamma_d / deltaZ_d)
    
    a_U   = D_u        #* dt # [W/m2/K]
    a_D   = D_d        #* dt # [W/m2/K]
    a_P_0 = c_vol * dZ / dt # [W/m2/K] (new) Patankar eq. 4.41c, this is b_p in Voller (1990; Eq. 30)           
    a_P   = a_U + a_D + a_P_0 - S_P * dZ #* dt # check the multiply on the S_P

    # a_P[g_liq_old>0] = 1e10
    
    b_0   = S_C  * dZ #* dt # [W/m2]
    b     = b_0 + a_P_0 * phi_t_old # is this one right? change 22/10/27

    ### Boundary conditions:
    ### type 1 is a specified value, type 2 is a specified gradient
    ### (units for gradient are degrees/meter)
    # bc_u_0  = phi_s
    bc_u_0 = phi_t_old[0]
    bc_type_u = 1
    bc_u    = np.concatenate(([ bc_u_0], [bc_type_u]))

    bc_d_0  = 0
    bc_type_d = 1
    # bc_d_0  = 0
    # bc_type_d = 2
    bc_d    = np.concatenate(([ bc_d_0 ], [ bc_type_d ]))

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

    delta_h = phi_t - phi_t_old #this will always be negative for layers with LWC

    ndh = -1*delta_h 

    cond1 = ((ndh>=H_latent) & (g_liq>0)) #layers where the energy change is larger than the energy required to refreeze, and where there is water
    phi_t[cond1] = delta_h[cond1] + H_latent[cond1] #sensible enthalpy will be the difference be
    g_liq[cond1] = 0 #no liquid left

    cond2 = ((ndh<H_latent) & (g_liq>0))
    phi_t[cond2] = 0
    g_liq[cond2] = (H_latent[cond2] + delta_h[cond2])/H_L_liq

    delta_g_liq = g_liq - g_liq_old # change in liquid fraction. Should always be negative.

    g_sol = g_sol + -1*delta_g_liq/0.917

    # g_liq = g_liq_iter + delta_h/(RHO_W_KGM*LF_I)
    # g_liq[g_liq<0] = 0
    # g_liq[g_liq>1] = 1
    # delta_gl = g_liq - g_liq_iter
    
    # H_excess = H_latent_old + delta_h
    # phi_t[H_excess<0] = phi_t[H_excess<0] + H_excess[H_excess<0]

    # UCF = 0.9
    # Delta_g_liq = a_P/a_P_0 * phi_t / (H_L_liq) 
 
    # g_liq = g_liq + UCF * Delta_g_liq
    # g_liq[g_liq<0] = 0
    # g_liq[g_liq>1] = 1

    H_tot = phi_t + H_L_liq*g_liq
    
    # iterdiff = np.abs(np.sum(H_tot_iter/(CP_I*RHO_I*g_sol)) - np.sum(H_tot/(CP_I*RHO_I*g_sol)))
    # if iterdiff<1e-4:
    #     itercheck = 0
    # else:
    #     itercheck = iterdiff

    # iterdiff = (np.sum(g_liq_iter) - np.sum(g_liq))
    # if iterdiff==0:
    #     itercheck = 0
    #     # print(f'g_liq_iter: {np.sum(g_liq_iter)}')
    #     # print(f'g_liq: {np.sum(g_liq)}')
    #     # print('break!')
    #     break
    # else:
    #     itercheck = np.abs( iterdiff/np.sum(g_liq_iter))
    # count += 1
    # if count==100:
    #     if ICT == 0:
    #         ICT = 1e-12
    #     else:
    #         pass
    # if ((count==200) and (ICT==1e-12)):
    #     # if ICT > 1e-12:
    #     #     pass
    #     # else:
    #     ICT = ICT * 10
    # if count>1000:
    #     print('breaking loop')
    #     print('iii', iii)
    #     print('itercheck',itercheck)
    #     print('iterdiff', iterdiff)
    #     print('np.sum(g_liq_iter)', np.sum(g_liq_iter))
    #     break
    ### END ITERATION

    iterdiff=0
    phi_t_out = H_tot/(CP_I*RHO_I*g_sol)
    phi_t_out[g_liq>0] = 0
    phi_t_out[(phi_t_out>0)] = 0 
    # print('max out:', np.max(phi_t_out))
    # print('iii',iii)
    # print('count',count)
    return phi_t_out, g_liq, count, iterdiff

###################################
### end transient_solve_EN ########
###################################


def transient_solve_EN_VOL(z_edges, Z_P, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s, mix_rho, c_vol, LWC, mass_sol, dz, ICT,iii=0):
    '''
    transient 1-d diffusion finite volume method for enthalpy

    :param z_edges:
    :param Z_P:
    :param nt: number of iterations; depricated (now uses while loop)
    :param dt: time step size
    :param Gamma_P:
    :param phi_0:
    :param nz_P:
    :param nz_fv:
    :param phi_s:
    :param g_liq
    :param c_vol: [J/m3/K] 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)

    :return phi_t:

    The source terms S_P and S_C come from the linearization described
    in Voller and Swaminathan, 1991, equations 31 and 32.
    and Voller, Swaminathan, and Thomas, 1990, equation 61
    '''

    phi_t = phi_0.copy()
    phi_t_old = phi_t.copy()

    vol_S       = mass_sol / RHO_I     # volume_Solid, i.e. volume of the ice (solid) portion of each control volume
    vol_SL      = vol_S + LWC    # volume of solid and liquid in each control volume
    mass_liq    = LWC * RHO_W_KGM  # mass of liquid water in each control
    mass_tot    = mass_liq + mass_sol # total mass (solid +liquid) of each control
    rho_liq_eff = mass_liq / dz      # effective density of the liquid portion
    mix_rho     = (mass_sol + mass_liq) / dz # mixture, or total density of volume (solid plus liquid), see Aschwanden (2012)
    g_liq       = LWC / dz    #  use liquid volume fraction of total volume of the control, which will net us the enthalpy/volume
    # g_liq_alt   = LWC / vol_SL  # alternatively, liquid volume fraction (of the material portion, porosity ignored), (I don't think this one is correct. See code at end of while loop to swap if you want to use this)
   
    g_liq_old = g_liq.copy()

    itercheck = 0.9
    count = 0

    # ICT = ICT # itercheck threshold 
    # print('###############')
    # print(f'iii = {iii}')
    # print(f'g_liq_old: {np.sum(g_liq_old)}')

    while itercheck>ICT:
        g_liq_iter = g_liq.copy() #unitless
        deltaH_l = rho_liq_eff * LF_I # [kg/m3 * J/kg = J/m3]
        deltaH = deltaH_l # deltaH is zero anywhere that has LWC = 0
        phi_t = phi_0.copy() # T in C
        dZ = np.diff(z_edges) #width of nodes

        deltaZ_u = np.diff(Z_P) # [m]
        deltaZ_u = np.append(deltaZ_u[0], deltaZ_u)
        
        deltaZ_d = np.diff(Z_P)
        deltaZ_d = np.append(deltaZ_d, deltaZ_d[-1])

        f_u = 1 - (Z_P[:] - z_edges[0:-1]) / deltaZ_u[:] # unitless
        f_d = 1 - (z_edges[1:] - Z_P[:]) / deltaZ_d[:]

        ### Gamma has units J/s/m/K (W/m/K)
        Gamma_U = np.append(Gamma_P[0], Gamma_P[0: -1])
        Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])

        Gamma_u =  1 / ((1 - f_u) / Gamma_P + f_u / Gamma_U) # Patankar eq. 4.9
        Gamma_d =  1 / ((1 - f_d) / Gamma_P + f_d / Gamma_D)

        ### Voller 1990 # Should be same result as 31, 32 from Voller 1991.
        # beta_P = np.zeros_like(Gamma_P)
        # beta_P[(g_liq>=0) & (g_liq<=1)]= -1.0e16
        # S_P = beta_P * deltaH * g_liq_old #* dZ * dt
        # S_C = -S_P*phi_t + deltaH*g_liq_old - deltaH*g_liq #+ deltaH*beta_P*
        ### end 1990 

        #### eqs. 31, 32 from Voller 1991
        dFdT = np.zeros_like(Gamma_P) # [1/K] dFdT is slope of the liquid fraction/temperature curve. Large at T=0C, i.e. fraction increases rapidly at T=0
        dFdT[(g_liq>=0) & (g_liq<=1)]= 1.0e10
        ### Finv has units of temperature: T=F^-1(g_liq). When g_liq>0, T=0C, and Finv=0. For g_liq<=0, T<=0 - but
        ### the equation for S_C contains deltaH multiplying Finv. deltaH is zero for nodes with T<0, so we can 
        ### just set Finv=0 at those nodes.
        Finv = np.zeros_like(dFdT)

        method='Voller'
        if method=='Voller':
            # deltaH is [J/m3]
            S_P = (-1 * deltaH * dFdT) / dt #[J/m3/K/s = W/m3/K]
            S_C = (deltaH * (g_liq_old - g_liq_iter) + deltaH*dFdT*Finv) / dt # [J/m3/s = W/m3]
        ### end 1991

        ### Brent method
        if method=='Brent':
            S_P = np.zeros_like(Gamma_P)
            S_C = deltaH * (g_liq_old - g_liq_iter) / dt


        D_u = (Gamma_u / deltaZ_u) # [W/m2/K]
        D_d = (Gamma_d / deltaZ_d)

        b_0   = S_C   * dZ #* dt # [W/m2]

        a_U   = D_u        #* dt # [W/m2/K]
        a_D   = D_d        #* dt # [W/m2/K]
        a_P_0 = c_vol * dZ / dt # [W/m2/K] (new) Patankar eq. 4.41c, this is b_p in Voller (1990; Eq. 30)           
        
        a_P   = a_U + a_D + a_P_0 - S_P * dZ #* dt # check the multiply on the S_P

        # b       = b_0 + a_P_0 * phi_t
        b       = b_0 + a_P_0 * phi_t_old # is this one right? change 22/10/27

        ### Boundary conditions:
        ### type 1 is a specified value, type 2 is a specified gradient
        ### (units for gradient are degrees/meter)
        bc_u_0  = phi_s
        bc_type_u = 1
        bc_u    = np.concatenate(([ bc_u_0], [bc_type_u]))

        bc_d_0  = 0
        bc_type_d = 1
        # bc_d_0  = 0
        # bc_type_d = 2
        bc_d    = np.concatenate(([ bc_d_0 ], [ bc_type_d ]))

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

        ### Use underrelaxation
        ### Brent
        if method=='Brent':
            g_delta = np.zeros_like(g_liq_old)
            g_delta[g_liq>0] = (a_P[g_liq>0] * ((phi_t[g_liq>0]) / (dZ[g_liq>0] * deltaH[g_liq>0]))) # dz or dZ?
            # print('g_delta',max(g_delta),min(g_delta))
            UF = 0.9 # underrelaxation parameter
            g_liq[g_liq>0] = g_liq[g_liq>0] + UF * g_delta[g_liq>0]

            ### Do not use underrelaxation (Voller says should not need to)
            # g_liq[g_liq>0] = g_liq[g_liq>0] + (a_P[g_liq>0] * ((phi_t[g_liq>0]) / (dZ[g_liq>0] * deltaH[g_liq>0]))) # dz or dZ?
        

        ### taken from v1.1.2; bit in v1.1.5 does not work!
        ### Voller 91
        elif method=='Voller':
            g_liq[g_liq>0] = (g_liq[g_liq>0] + a_P[g_liq>0] * ((phi_t[g_liq>0]) / (dZ[g_liq>0] * deltaH[g_liq>0]))) # dz or dZ?
            g_liq[g_liq<=0] = 0
            g_liq[g_liq>1] = 1
       
        LWC = g_liq*dz
        mass_liq = LWC * RHO_W_KGM
        rho_liq_eff = mass_liq / dz #dz or dZ? Only a very small difference.

        ### Use the following if you use g_liq_alt above
        # new_vol = mass_tot/(g_liq*RHO_W_KGM + RHO_I - g_liq*RHO_I)
        # LWC = g_liq * new_vol
        # mass_liq = LWC * RHO_W_KGM
        # rho_liq_eff = mass_liq / dZ #dz or dZ?

        iterdiff = (np.sum(g_liq_iter) - np.sum(g_liq))
        if iterdiff==0:
            itercheck = 0
            # print(f'g_liq_iter: {np.sum(g_liq_iter)}')
            # print(f'g_liq: {np.sum(g_liq)}')
            # print('break!')
            break
        else:
            itercheck = np.abs( iterdiff/np.sum(g_liq_iter))
        # print(f'count: {count}')
        # print(f'itercheck: {itercheck}')
        # print(f'g_liq_old: {np.sum(g_liq_old)}')
        # print(f'g_liq_iter: {np.sum(g_liq_iter)}')
        # print(f'g_liq: {np.sum(g_liq)}')
        # print(f'ICT: {ICT}')
        count += 1
        if count==100:
            if ICT == 0:
                ICT = 1e-12
            else:
                pass
        if ((count==200) and (ICT==1e-12)):
            # if ICT > 1e-12:
            #     pass
            # else:
            ICT = ICT * 10
        if count>1000:
            print('breaking loop')
            print(iii)
            print('itercheck',itercheck)
            break

    # print(count)
    # print('###############')
    # input('waiting')
    return phi_t, g_liq, count, iterdiff

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