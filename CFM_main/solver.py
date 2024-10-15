#!/usr/bin/env python
'''
Functions to solve the diffusion equation
'''

import numpy as np
np.set_printoptions(precision=4)
from scipy import interpolate
import scipy.integrate
from scipy.sparse import spdiags
import scipy.sparse.linalg as splin
from constants import *
import sys
from scipy.linalg import lapack


def solver(a_U, a_D, a_P, b):
    '''
    function for solving matrix problem

    :param a_U:
    :param a_D:
    :param a_P:
    :param b:

    :return phi_t:
    '''

    #NOTE: lapack contains a dedicated solver that handles tridagonal
    # solutions.  It's quite a bit faster than using the general spsolve
    # function.  I've included two boolean variables that let you choose
    # which function to solve the equation-- if both are enabled, the
    # normalized difference between the results will be reported.

    #use_splin=True
    #use_dgtsv=False

    use_splin=False
    use_dgtsv=True
    if use_splin:
        nz = np.size(b)

        diags = (np.append([a_U, -a_P], [a_D], axis = 0))
        cols = np.array([1, 0, -1])

        big_A = spdiags(diags, cols, nz, nz, format = 'csc')
        big_A = big_A.T

        rhs = -b
        phi_t = splin.spsolve(big_A, rhs)
        if use_dgtsv:
            phi_t_spsolve=phi_t.copy()
    if use_dgtsv:
        dl=np.ascontiguousarray(a_U[1:])
        d=np.ascontiguousarray(-a_P)
        du=np.ascontiguousarray(a_D[:-1])
        rhs = np.ascontiguousarray(-b)
        _, _, _, phi_t, _ = lapack.dgtsv(dl, d, du, rhs)
    if use_splin and use_dgtsv:
        print(np.max(np.abs(phi_t-phi_t_spsolve))/np.median(np.abs(phi_t)))
    return phi_t

####!!!!

def transient_solve_TR(z_edges, Z_P, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s, tot_rho, c_vol, airdict=None):
    '''
    transient 1-d diffusion finite volume method

    This is for standard heat (no liquid water), isotope, and air diffusion.
    If there is liquid water is should use the enthalpy solver.

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
    phi_t_old = phi_t.copy()

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
        bc_type_d = 2 # 2 is gradient
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
        if bc_type_d==2:
            a_U[-1] = 1
            b[-1]   = deltaZ_u[-1] * bc_d[0]
        elif bc_type_d==1:
            a_U[-1] = 0
            b[-1]   = bc_d[0]

        phi_t = solver(a_U, a_D, a_P, b)

        a_P = a_U + a_D + a_P_0

    if airdict!=None:
        return phi_t, w_p
    else:
        return phi_t

###################################
### end transient_solve_TR ########
###################################

def transient_solve_EN(z_edges, Z_P, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s, mix_rho, c_vol, LWC, mass_sol, dz, ICT, rho_firn, iii=0):
    '''
    transient 1-d diffusion finite volume method for enthalpy

    This is for heat diffusion when there is liquid water present.
    It uses a source term for the latent heat associated with the liquid water.

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
    rho_liq_eff = RHO_W_KGM*dz #new 10/31
    g_liq       = LWC / dz    #  use liquid volume fraction of total volume of the control, which will net us the enthalpy/volume
    g_sol       = vol_S / dz #unitless

    dZ = np.diff(z_edges) #width of nodes

    H_L_liq = RHO_W_KGM*LF_I #volumetric latent enthalpy [J/m3]

    phi_t = phi_0 # phi_t is just the temperature

    phi_t_old = phi_t.copy() # initial temperature
    g_liq_old = g_liq.copy() # initial liquid fraction
    g_sol_old = g_sol.copy() # initial solid fraction

    itercheck = 0.9
    count = 0


    ### Big H stands for latent enthalpy, little h is sensible enthalpy

    ### H_tot is the sum of latent and sensible enthalpy
    H_lat      = H_L_liq*g_liq # Latent enthalpy for each layer
    H_lat_old  = H_lat.copy()
    H_tot         = phi_t * g_sol * RHO_I * CP_I + H_lat  #total enthalpy, voller 1990b eq 4a
    H_tot_old     = H_tot.copy()
    h_old         = phi_t * g_sol * RHO_I * CP_I
    h_updated     = h_old.copy()

    update_gsol = True

    for i_time in range(100): # Testing indicates that this should never need this many iterations

        H_tot_iter  = H_tot.copy()
        phi_iter    = phi_t.copy()
        g_liq_iter  = g_liq.copy() #unitless
        g_sol_iter  = g_sol.copy()
        H_lat_iter  = H_lat.copy()
        h_iter      = h_updated.copy()

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

        #### version with working dt ####
        ### S_C is independent
        S_P = np.zeros_like(Gamma_P)
        S_C = H_L_liq  * (g_liq_old - g_liq_iter) # J/m3/s = W/m3 Latent heat as source term.
        # S_C = RHO_I * CP_I * phi_iter * (g_sol_old - g_sol_iter) + (RHO_W_KGM * CP_W * phi_iter + H_L_liq)  * (g_liq_old - g_liq_iter)

        D_u = (Gamma_u / deltaZ_u) # [W/m2/K]
        D_d = (Gamma_d / deltaZ_d)

        a_U   = D_u        #* dt # [W/m2/K]
        a_D   = D_d        #* dt # [W/m2/K]

        c_vol1 = RHO_I * CP_I # gets multiplied by g_sol below

        a_P_0 = c_vol1 * dZ / dt # [W/m2/K] (new) Patankar eq. 4.41c, this is b_p in Voller (1990; Eq. 30)

        if update_gsol:
            a_P   = a_U + a_D + a_P_0 * g_sol_iter - S_P * dZ #* dt # check the multiply on the S_P
        else:
            a_P   = a_U + a_D + a_P_0 * g_sol_old - S_P * dZ #* dt # check the multiply on the S_P

        b_0   = S_C * dZ/dt #* dt # [W/m2]
        b     = b_0 + a_P_0 * g_sol_old * phi_t_old # By this phi_t_old has to be in K

        ###############

        ### Boundary conditions:
        ### type 1 is a specified value, type 2 is a specified gradient
        ### (units for gradient are degrees/meter)
        # bc_u_0  = phi_s
        bc_u_0 = phi_t_old[0]
        # bc_u_0 = phi_iter[0]
        bc_type_u = 1
        bc_u      = np.concatenate(([ bc_u_0], [bc_type_u]))

        bc_d_0    = 0 #Value
        bc_type_d = 2 # Type, 1 for fixed temperature, 2 for gradient
        bc_d      = np.concatenate(([ bc_d_0 ], [ bc_type_d ]))

        #Upper boundary
        a_P[0]  = 1
        a_U[0]  = 0
        a_D[0]  = 0
        b[0]    = bc_u[0]

        #Down boundary
        a_P[-1] = 1
        a_D[-1] = 0
        if bc_type_d==2: # Gradient (flux)
            a_U[-1] = 1
            b[-1]   = deltaZ_u[-1] * bc_d[0]
        elif bc_type_d==1: # Fixed value
            a_U[-1] = 0
            b[-1]   = bc_d[0]

        #####
        phi_t = solver(a_U, a_D, a_P, b) #sensible enthalpy. 0 for layers at freezing (have LWC), negative for dry layers
        #####

        '''
        ####
        The crux is to adjust liquid fraction and temperture field based on solution
        Calculations are (partially) based on the fact that the freezing temp is 0,
        which means that if H_tot<0 there is no liquid and you can calculate temperature, and if H_tot>0 there is liquid and T is 0.

        Note previous ways of solving in dev branch and previous releases.
        ###
        '''

        ### The best way to solve: calculate new g_liq, and apply overshoot correction on g_liq
        ### Break loop if iteration (prior to overshoot) is the same as previous solution

        h_updated = phi_t * CP_I * RHO_I * g_sol_old # updated sensible enthalpy after solver. g_sol is volume_solid/dz
        delta_h = h_updated - h_old # change in sensible enthalpy, relative to initial (not iteration)

        ### Figure out what delta_h and g_liq should be based on different conditions
        cond0 = ((delta_h>0) & (LWC_old>0)) # Layers where there was sensible enthalpy increased and there is liquid water
        delta_h[cond0] = 0 # sensible enthalpy should not increase if LWC present (should either stay at 0C or cool down)

        ndh = -1*delta_h #negative delta_h (which makes it positive), makes corrections below easy

        ### everything refreezes if the calculated change in enthalpy is greater than the latent enthalpy
        cond1 = ((ndh>=H_lat_old) & (g_liq_old>0)) #layers where energy change is larger than the needed to refreeze, and where there is water
        H_tot[cond1] = (delta_h[cond1] + H_tot_old[cond1]) # total enthalpy in those layers is change+total, should be net negative
        g_liq[cond1] = 0 #no liquid left

        ### partial refreezing if the delta_h is less than the latent enthalpy
        cond2 = ((ndh<H_lat) & (g_liq_old>0))
        H_tot[cond2] = 0
        g_liq[cond2] = (H_lat_old[cond2] + delta_h[cond2])/H_L_liq #remaining liquid

        ### Make sure that there is no liquid in layers that did not have liquid at start
        cond3 = (g_liq_old<=0)
        g_liq[cond3] = 0
        H_tot[cond3] = h_updated[cond3]

        ### if this iteration gave the same solution as the last iteration, break
        iterdiff = (np.sum(g_liq_iter) - np.sum(g_liq)) # Deprecated? Used to use to calulate time to break loop
        if ((np.allclose(phi_iter,phi_t,rtol=1e-4,atol=1e-3))):
            break
        elif ((np.allclose(g_liq_iter,g_liq,rtol=1e-4,atol=1e-4))):
            break

        ### otherwise apply overshoot correction to g_liq and iterate again.
        delta_g_liq = g_liq - g_liq_iter # change in liquid fraction. Should always be negative.
        g_liq = g_liq_iter + 0.6*delta_g_liq
        g_liq[g_liq<0]=0
        g_liq[g_liq>1]=1
        ################

        g_sol = g_sol + -1 * delta_g_liq/0.917 # Update solid fraction

        ## Now update temperatures after liquid corrections
        phi_t[H_tot>=0] = 0 # H_tot>0 means liquid present, T=0
        phi_t[H_tot<0] = H_tot[H_tot<0] / (CP_I * RHO_I * g_sol_old[H_tot<0]) # H_tot<0 means no liquid; all enthalpy is sensible

        phi_t[g_liq>0] = 0
        #############

        count += 1

        ### END ITERATION LOOP
        ######################

    phi_t_out = phi_t
    phi_t_out[g_liq>0] = 0
    phi_t_out[(phi_t_out>0)] = 0

    return phi_t_out, g_liq, count, iterdiff,g_sol

###################################
### end transient_solve_EN ########
###################################

def apparent_heat(z_edges, Z_P, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s, mix_rho, c_vol, LWC, mass_sol, dz, ICT, rho_firn, iii=0):
    '''
    transient 1-d diffusion finite volume method

    This is for standard heat (no liquid water), isotope, and air diffusion.
    If there is liquid water is should use the enthalpy solver.

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

    LWC_old = LWC.copy()
    phi_in = phi_0.copy()

    vol_S       = mass_sol / RHO_I     # volume_Solid, i.e. volume of the ice (solid) portion of each control volume
    vol_SL      = vol_S + LWC    # volume of solid and liquid in each control volume
    mass_liq    = LWC * RHO_W_KGM  # mass of liquid water in each control
    mass_tot    = mass_liq + mass_sol # total mass (solid +liquid) of each control
    rho_tot = mass_tot/dz
    rho_liq_eff = RHO_W_KGM*dz #new 10/31
    g_liq       = LWC / dz    #  use liquid volume fraction of total volume of the control, which will net us the enthalpy/volume
    g_sol       = vol_S / dz #unitless

    dZ = np.diff(z_edges) #width of nodes

    H_L_liq = RHO_W_KGM*LF_I #volumetric latent enthalpy [J/m3]

    phi_t = phi_0 # phi_t is just the temperature, C
    phi_t_old = phi_t.copy()

    phi_t_old = phi_t.copy() # initial temperature
    g_liq_old = g_liq.copy() # initial liquid fraction
    g_sol_old = g_sol.copy() # initial solid fraction

    H_lat = g_liq * H_L_liq
    h_sens = phi_t * CP_I * mass_sol
    H_tot = H_lat + h_sens

    H_tot_old = H_tot.copy()
    h_old = h_sens.copy()
    H_lat_old = H_lat.copy()

    count = 0

    for i_time in range(20):
        phi_t_iter = phi_t.copy()
        H_tot_iter = H_tot.copy()
        g_liq_iter = g_liq.copy()

        delH = np.gradient(H_tot,Z_P)
        delT = np.gradient(phi_t,Z_P)

        # C_app = (delH**2/delT**2)**0.5
        # if i_time == 0:

        # C_app = np.zeros_like(phi_t)
        C_app_0 = g_liq_iter * H_L_liq * 50

        c_vol = c_vol + C_app_0

        # else:
        #     C_app = (H_tot - H_tot_old)/(phi_t_iter - phi_t_old)

        # lwc_layers  = np.where(g_liq>0)[0]
        # c_vol[lwc_layers] = C_app[lwc_layers]

        dZ = np.diff(z_edges) #width of nodes

        deltaZ_u = np.diff(Z_P)
        deltaZ_u = np.append(deltaZ_u[0], deltaZ_u)

        deltaZ_d = np.diff(Z_P)
        deltaZ_d = np.append(deltaZ_d, deltaZ_d[-1])

        f_u = 1 - (Z_P[:] - z_edges[0:-1]) / deltaZ_u[:]
        f_d = 1 - (z_edges[1:] - Z_P[:]) / deltaZ_d[:]

        #######################################
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
        if bc_type_d==2:
            a_U[-1] = 1
            b[-1]   = deltaZ_u[-1] * bc_d[0]
        elif bc_type_d==1:
            a_U[-1] = 0
            b[-1]   = bc_d[0]

        phi_t = solver(a_U, a_D, a_P, b)

        h_updated = phi_t * CP_I * RHO_I * g_sol_old # updated sensible enthalpy after solver

        delta_h = (h_updated - h_old) #* 0.6 # overshoot correction. delta_h should always be negative for layers with LWC
        h_updated_trun = h_old + delta_h # sensible enthalpy after overshoot

        cond0 = ((delta_h>0) & (LWC_old>0)) # Layers where there was sensible enthalpy increased and there is liquid water
        delta_h[cond0] = 0 # sensible enthalpy should not increase if LWC present (should either stay at 0C or cool down)

        ndh = -1*delta_h #negative delta_h (which makes it positive), makes corrections below easy

        ### everything refreezes if the calculated change in enthalpy is greater than the latent enthalpy
        cond1 = ((ndh>=H_lat_old) & (g_liq_old>0)) #layers where energy change is larger than the needed to refreeze, and where there is water
        H_tot[cond1] = (delta_h[cond1] + H_tot_old[cond1]) # total enthalpy in those layers is change+total, should be net negative
        g_liq[cond1] = 0 #no liquid left

        ### partial refreezing if the delta_h is less than the latent enthalpy
        cond2 = ((ndh<H_lat) & (g_liq_old>0))
        H_tot[cond2] = 0
        g_liq[cond2] = (H_lat_old[cond2] + delta_h[cond2])/H_L_liq #remaining liquid

        ### Make sure that there is no liquid in layers that did not have liquid at start
        cond3 = (g_liq_old<=0)
        g_liq[cond3] = 0
        H_tot[cond3] = h_updated_trun[cond3]

        g_liq[g_liq<0]=0
        g_liq[g_liq>1]=1

        delta_g_liq = g_liq - g_liq_old # change in liquid fraction. Should always be negative.
        #######

        g_sol = g_sol + -1 * delta_g_liq/0.917 # Update solid fraction

        ### Now update temperatures after liquid corrections
        phi_t[H_tot>=0] = 0 # H_tot>0 means liquid present, T=0
        phi_t[H_tot<0] = H_tot[H_tot<0] / (CP_I * RHO_I * g_sol_old[H_tot<0]) # H_tot<0 means no liquid; all enthalpy is sensible

        iterdiff = (np.sum(g_liq_iter) - np.sum(g_liq)) # Deprecated? Used to use to calulate time to break loop

        if ((np.allclose(phi_t_iter,phi_t,rtol=1e-4,atol=1e-3))):
            break

        elif ((np.allclose(g_liq_iter,g_liq,rtol=1e-4,atol=1e-4))):
            break

        count+=1

    return phi_t, g_liq, count, iterdiff,g_sol

###################################
### end apparent_heat ########
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
