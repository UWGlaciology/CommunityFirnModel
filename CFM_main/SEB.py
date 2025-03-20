#!/usr/bin/env python

import numpy as np 
# import scipy.interpolate as interpolate
from scipy import optimize
from constants import *
# import os
import sys
import math, cmath
# try:
# from numba import njit
    # usejit = True
# except:

import pandas as pd

usejit = False
### working on jit implementation. not working 24/10/22

class SurfaceEnergyBudget:
    '''
    Class to handle surface energy balance in the CFM

    The energy balance equation:

    E_net = SW_d + SW_u + LW_d + LW_u + G + QH + QL + EP

    Parameters
    ------------
    E_net: value
        Net surface energy [W/m2]
    SW_d: value
        Downward shortwave radiation [W/m2]
    SW_u: value
        Upward shortwave radiation [W/m2]
    LW_d: value
        Downward longwave radiation [W/m2]
    LW_u: value
        Upward longwave radiation [W/m2]
    G: value
        Subsurface energy flux
    QH: value
        sensible heat flux
    QL: value
        latent heat flux
    EP: value
        heat flux from precipitation
    ALBEDO: value
        surface ALBEDO
    T2m: value
        2-m temperature (K)

    '''

    def __init__(self,config,climateTS,start_ind,SEBfluxes, start_ind_EF):
        '''
        intialalize seb
        consider to be in beta

        in any case, climateTS needs to have T2m in it, 2 m temperature at main model run time resolution

        '''

        self.c = config
        
        if SEBfluxes is not None:
            self.time_in = SEBfluxes['time'][start_ind_EF:]
            self.SW_d    = SEBfluxes['SW_d'][start_ind_EF:]
            self.LW_d    = SEBfluxes['LW_d'][start_ind_EF:]
            self.ALBEDO  = SEBfluxes['ALBEDO'][start_ind_EF:]
            self.ALBEDO[np.isnan(self.ALBEDO)] = np.nanmean(self.ALBEDO)
            try:
                self.albedo_factor = self.c['albedo_factor']
                self.ALBEDO = self.ALBEDO * self.albedo_factor
                print(f'ALBEDO is multipled by {self.albedo_factor}')
            except:
                print('"albedo_factor" not defined in .json. Using 1')
            self.T2m     = SEBfluxes['T2m'][start_ind_EF:]
            self.TSKIN   = SEBfluxes['TSKIN'][start_ind_EF:]
            self.QH      = SEBfluxes['QH'][start_ind_EF:] #MERRA fluxes are upward positive, so multiply by -1
            self.QL      = SEBfluxes['QL'][start_ind_EF:]
            self.RAIN    = SEBfluxes['RAIN'][start_ind_EF:] # [m i.e./year]
            if 'LW_u' in SEBfluxes:
                self.LW_u_input = True
                self.LW_u = SEBfluxes['LW_u'][start_ind_EF:]
            self.dtRATIO = SEBfluxes['dtRATIO'] # This is the ratio of dt in SEBfluxes to dt in climateTS (how many sub steps to do, e.g it is 6 with daily res climateTS and 4-hourly SEB fields)
            self.T2m_daily = climateTS['T2m'][start_ind:]

            self.df_CLIM = pd.DataFrame({'SW_d':self.SW_d,'LW_d':self.LW_d,'ALBEDO':self.ALBEDO, 'T2m':self.T2m, 
                'TSKIN':self.TSKIN, 'QH':self.QH, 'QL':self.QL, 'RAIN':self.RAIN},index=self.time_in)
            self.df_CLIM['flux'] = ((self.df_CLIM['SW_d'] * (1 - self.df_CLIM['ALBEDO'])) + self.df_CLIM['LW_d'] + self.df_CLIM['QH'] + self.df_CLIM['QL'])

        else:       
            self.time_in = climateTS['time'][start_ind:]
            self.SW_d    = climateTS['SW_d'][start_ind:]
            self.LW_d    = climateTS['LW_d'][start_ind:]
            self.ALBEDO  = climateTS['ALBEDO'][start_ind:]
            self.ALBEDO[np.isnan(self.ALBEDO)] = np.nanmean(self.ALBEDO)
            try:
                self.albedo_factor = self.c['albedo_factor']
                self.ALBEDO = self.ALBEDO * self.albedo_factor
                print(f'ALBEDO is multipled by {self.albedo_factor}')
            except:
                print('"albedo_factor" not defined in .json. Using 1')
            self.T2m     = climateTS['T2m'][start_ind:]
            self.TSKIN   = climateTS['TSKIN'][start_ind:]
            if (('EVAP' in climateTS.keys()) and ('SUBLIM' in climateTS.keys())):
                self.EVAP    = climateTS['EVAP'][start_ind:]
                self.SUBLIM    = climateTS['SUBLIM'][start_ind:]
            elif 'SUBLIM' in climateTS.keys():
                self.SUBLIM    = climateTS['SUBLIM'][start_ind:]
                self.EVAP    = np.zeros_like(self.SUBLIM)
            elif 'EVAP' in climateTS.keys():
                self.EVAP    = climateTS['EVAP'][start_ind:]
                self.SUBLIM    = np.zeros_like(self.EVAP)
            self.QH      = -1*climateTS['QH'][start_ind:] #MERRA fluxes are upward positive, so multiply by -1
            self.QL      = -1*climateTS['QL'][start_ind:]
            self.RAIN    = climateTS['RAIN'][start_ind:] # [m i.e./year]
            if 'LW_u' in climateTS:
                self.LW_u_input = True
                self.LW_u = climateTS['LW_u'][start_ind:]
            self.dtRATIO = 1
        # need to account for cold snow falling on warmer surface

        # self.EP      = np.zeros_like(self.SW_d)
        # self.G       = np.zeros_like(self.SW_d) # For now we are not considering any flux in/out of the upper model node from below
        try:
            self.TL_thick = self.c['SEB_TL_thick']
        except:
            print('Top layer thickness ("SEB_TL_thick") not defined in .json')
            print('Setting to 5 cm')
            self.TL_thick = 0.05

        
        self.SBC = 5.67e-8 # Stefan-Boltzmann constant [W K^-4 m^-2
        self.emissivity_air = 1
        self.emissivity_snow = 0.98
        
        # self.D_sh = 15 # Sensible heat flux coefficient, Born et al. 2019 [W m^-2 K^-1]

    def SEB_fqs(self,PhysParams,iii,T_old):
        '''
        Calculate surface energy using fqs solver
        Positive fluxes are into the surface, negative are out
        SEBparams: mass, Tz, dt
        '''

        Tz   = PhysParams['Tz']
        mass = PhysParams['mass']
        dt   = PhysParams['dt'] # [s]
        Tguess = self.T2m[iii]
        dz  = PhysParams['dz']
        z = PhysParams['z']
        mtime = PhysParams['mtime']

        T_rain = np.max((self.T2m[iii],T_MELT))
        # Qrain_i = 0

        rain_mass = self.RAIN[iii] * RHO_I / S_PER_YEAR * dt #[kg] of rain at this timestep
        Qrain_i =  CP_W * rain_mass * (T_rain - T_MELT) # Assume rain temperature is air temp, Hock 2005, eq 19
        #latent heat for rain falling on top of cold snow should be handled in melt.py

        Q_SW_net = self.SW_d[iii] * (1-self.ALBEDO[iii])
        # Q_LW_d = self.SBC * (self.emissivity_air * self.T2m[iii]**4)
        Q_LW_d = self.emissivity_air * self.LW_d[iii]

        i_GL = np.where(z>=1)[0][0]
        z_GL = z[i_GL]
        m_GL = np.cumsum(mass)[i_GL]
        T_GL = np.cumsum(mass*Tz)[i_GL]/m_GL
        rho_GL = m_GL/z_GL
        K_ice   = 9.828 * np.exp(-0.0057 * T_GL) #[W/m/K]

        K_GL  = K_ice * (rho_GL/RHO_I) ** (2 - 0.5 * (rho_GL/RHO_I))
        
        G = (K_GL * (Tz[i_GL] - Tz[0])/z_GL) # estimated temperature flux in firn due to temperature gradient
        # G = 0

        # TL_thick = 0.1 # thickness of snow/firn "Top Layer" that energy goes into. Reducing results in higher melt.
        iTL = np.where(z>=self.TL_thick)[0][0]

        for kk in range(10): # this loop is to make sure that the toplayer is thick enough        
            dTL = np.cumsum(dz)[iTL]

            m = np.cumsum(mass)[iTL] #mass of the TL
            
            TTL = np.cumsum(mass*Tz)[iTL]/m # mean temperature of top X cm (weighted mean)
            cold_content_TL = CP_I * m * (T_MELT - TTL) # cold content [J], positive quantity if T<T_melt

            Qnet = Q_SW_net + Q_LW_d + self.QH[iii] + self.QL[iii] + Qrain_i + G
            # fqs = FQS()

            pmat = np.zeros(5)

            a = self.emissivity_snow * self.SBC * dt/(CP_I*m)
            b = 0
            c = 0
            d = 1
            e = -1 * (Qnet*dt/(CP_I*m)+TTL)

            pmat[0] = a
            pmat[3] = d
            pmat[4] = e
            pmat[np.isnan(pmat)] = 0

            r = quartic_roots(pmat)
            Tsurface = (r[((np.isreal(r)) & (r>0))].real)
            
            if Tsurface>=273.15:
                Tsurface = 273.15
                meltmass = (Qnet - self.SBC*273.15**4) * dt / LF_I #multiply by dt to put in units per time step
                # do not need to subtract cold content to calculate cold content b/c Q_melt = sum(energies), Q_melt=0 if energy can be balanced, i.e. sum(energies)=0
                # melt_mass = (Qnet - self.SBC*273.15**4) * dt / LF_I
            ### meltmass has units [kg/m2/s]
            else:
                meltmass = 0

            if meltmass<=m: #if the melt mass is greater than the mass of the TL layer, we need a thicker TL because the next layer could be below freezing, and it needs to warm before melting
                break
            else:
                iTL = np.where(np.cumsum(mass)>=meltmass)[0][0]

        Tz[0:iTL+1] = Tsurface

        return Tsurface, Tz, meltmass, self.TSKIN[iii]
    ############################
    ### end SEB_fqs
    ############################

    def SEB_fqs_subdt(self,PhysParams,iii,T_old):
        '''
        Same as above, but loops to find Ts, melt when SEB fluxes have finer time resolution
        Calculate surface energy using fqs solver
        Positive fluxes are into the surface, negative are out
        SEBparams: mass, Tz, dt

        the indices in SEBfluxes should start at iii*dtRATIO
        '''

        Tz   = PhysParams['Tz']
        mass = PhysParams['mass']
        dtMAIN = PhysParams['dt']
        dt   = PhysParams['dt']/self.dtRATIO # [s]
        Tguess = self.T2m[iii]
        dz  = PhysParams['dz']
        z = PhysParams['z']
        mtime = PhysParams['mtime']

        iTL = np.where(z>=self.TL_thick)[0][0]

        i_GL = np.where(z>=1)[0][0]
        z_GL = z[i_GL]
        m_GL = np.cumsum(mass)[i_GL]
        T_GL = np.cumsum(mass*Tz)[i_GL]/m_GL
        rho_GL = m_GL/z_GL
        K_ice   = 9.828 * np.exp(-0.0057 * T_GL) #[W/m/K]

        K_GL  = K_ice * (rho_GL/RHO_I) ** (2 - 0.5 * (rho_GL/RHO_I))
        
        G = (K_GL * (Tz[i_GL] - Tz[0])/z_GL) # estimated temperature flux in firn due to temperature gradient
        # G = 0

        iiisub = (iii)*self.dtRATIO
        
        flux_df1 = self.df_CLIM.iloc[iiisub:iiisub+self.dtRATIO]
        flux_df1_r = flux_df1.flux.values + G

        Tcalc = np.zeros_like(flux_df1_r)
        meltmass = np.zeros_like(flux_df1_r)

        dts = flux_df1.TSKIN.values

        tindex = flux_df1.index

        dTL = np.cumsum(dz)[iTL]
        m = np.cumsum(mass)[iTL] #mass of the TL
        
        TTL = np.cumsum(mass*Tz)[iTL]/m # Temp Top Layer, mean temperature of top X cm (weighted mean)
        
        # if usejit:
        #     fqs = FQS_nojit()
        # else:
        #     pass

        for kk,mdate in enumerate(tindex): #loop through all time steps, calculate melt for all lat/lon pairs at that time
            pmat = np.zeros(5) # p matrix to put into FQS solver
            if kk==0:
                T_0 = TTL
            else:
                T_0 = Tcalc[kk-1]

            a = np.array(self.SBC * dt / (CP_I*m))
            b = np.array(0.0)
            c = np.array(0.0)
            d = np.array(1.0)
            e = np.array(-1 * (flux_df1_r[kk]*dt/(CP_I*m)+T_0))

            pmat[0] = a
            pmat[3] = d
            pmat[4] = e
            pmat[np.isnan(pmat)] = np.array(0.0)
            
            # if usejit:
            #     r = fqs.quartic_roots(pmat)
            # else:
            r = quartic_roots(pmat)
            Tnew = (r[((np.isreal(r)) & (r>0))].real)
            Tnew[np.isnan(e)] = np.nan
            try:
                if Tnew>=273.15:
                    Tcalc[kk] = 273.15000000000000        
                    meltmass[kk] = (flux_df1_r[kk] - self.SBC*273.15**4) / LF_I * dt #multiply by dt to put in units per time step

                else:
                    Tcalc[kk] = Tnew
                    meltmass[kk] = 0
            except:
                print(f'r: {r}')
                print(f'Tnew:{Tnew}')
                print(f'T_0:{Tcalc[kk-10:kk]}')
                
                print(iii)
                sys.exit()

        Tsurface_out = np.mean(Tcalc)
        meltmass_out = np.sum(meltmass) #if you sum the mass, there will be a too much melt because model resolution does not also calculate the temperature profile during the colder hours

        # if meltmass_out<=m: #if the melt mass is greater than the mass of the TL layer, we need a thicker TL because the next layer could be below freezing, and it needs to warm before melting
        #     break
        # else:
        #     iTL = np.where(np.cumsum(mass)>=meltmass_out)[0][0]

        if meltmass_out>0:
            Tsurface_out = 273.15

        Tz[0:iTL+1] = Tsurface_out

        return Tsurface_out, Tz, meltmass_out, self.TSKIN[iii]
    ############################
    ### end SEB_fqs_subdt
    ############################

    def SEB_loop(self,PhysParams,iii,T_old):
        '''
        This is developement code - do not use for science

        Calculate surface energy using fqs solver
        Positive fluxes are into the surface, negative are out
        SEBparams: mass, Tz, dt
        '''

        Tz   = PhysParams['Tz']
        mass = PhysParams['mass']
        dt   = PhysParams['dt'] # [s]
        Tguess = self.T2m[iii]
        dz  = PhysParams['dz']
        z = PhysParams['z']

        T_rain = np.max((self.T2m[iii],T_MELT))
        # Qrain_i = 0

        rain_mass = self.RAIN[iii] * RHO_I / S_PER_YEAR * dt #[kg] of rain at this timestep
        Qrain_i =  CP_W * rain_mass * (T_rain - T_MELT) # Assume rain temperature is air temp, Hock 2005, eq 19
        #latent heat for rain falling on top of cold snow should be handled in melt.py

        Q_SW_net = self.SW_d[iii] * (1-self.ALBEDO[iii])
        # Q_LW_d = self.SBC * (self.emissivity_air * self.T2m[iii]**4)
        Q_LW_d = self.emissivity_air * self.LW_d[iii]

        # TL_thick = 0.1 # thickness of snow/firn "Top Layer" that energy goes into. Reducing results in higher melt.
        iTL = np.where(z>=self.TL_thick)[0][0] 
        dTL = z[iTL]

        i_GL = np.where(z>=1)[0][0]
        z_GL = z[i_GL]
        m_GL = np.cumsum(mass)[i_GL]
        T_GL = np.cumsum(mass*Tz)[i_GL]/m_GL
        rho_GL = m_GL/z_GL
        K_ice   = 9.828 * np.exp(-0.0057 * T_GL) #[W/m/K]

        K_GL  = K_ice * (rho_GL/RHO_I) ** (2 - 0.5 * (rho_GL/RHO_I))
        G = (K_GL * (Tz[i_GL] - Tz[0])/z_GL) # estimated temperature flux in firn due to temperature gradient

        m = np.cumsum(mass)[iTL] #mass of the top layer
        TTL = np.cumsum(mass*Tz)[iTL]/m # mean temperature of top layer (weighted mean)
        cold_content_TL = CP_I * m * (T_MELT - TTL) # cold content [J], positive quantity if T<T_melt

        Tnew = TTL.copy()

        Q_sum = Q_SW_net + Q_LW_d + self.QH[iii] + self.QL[iii] + Qrain_i + G #sum of all flux terms that are not Temperature dependent

        def Qnet(Ts,Qsum):
            Qout = np.abs(-1*self.SBC*Ts**4 + Qsum)
            return Qout

        sol = optimize.minimize(Qnet,TTL,args=Q_sum,method='Nelder-Mead')

        Tsurface = sol.x[0]

       
        if Tsurface>=273.15:
            Tsurface = 273.15
            meltmass = (Q_sum - self.SBC*273.15**4) * dt / LF_I #*dt #multiply by dt to put in units per day
            # do not need to subtract cold content to calculate cold content b/c Q_melt = sum(energies), Q_melt=0 if energy can be balanced, i.e. sum(energies)=0
            # melt_mass = (Qnet - self.SBC*273.15**4) * dt / LF_I
        ### meltmass has units [kg/m2/s]
        else:
            meltmass = 0

        Tz[0:iTL+1] = Tsurface            

        return Tsurface, Tz, meltmass
    ############################
    ### end SEB_loop
    ############################

    def SEB(self, PhysParams,iii,T_old):
        '''
        Development-level code

        Calculate the surface energy budget
        Positive fluxes are into the surface, negative are out
        SEBparams: mass, Tz, dt
        '''

        Tz   = PhysParams['Tz']
        mass = PhysParams['mass']
        dt   = PhysParams['dt'] # [s]
        Tguess = self.T2m[iii]
        dz  = PhysParams['dz']
        z = PhysParams['z']

        T_rain = np.max((self.T2m[iii],T_MELT))
        # Qrain_i = 0

        rain_mass = self.RAIN[iii] * RHO_I / S_PER_YEAR * dt #[kg] of rain at this timestep
        Qrain_i =  CP_W * rain_mass * (T_rain - T_MELT) # Assume rain temperature is air temp, Hock 2005, eq 19
        #latent heat for rain falling on top of cold snow should be handled in melt.py

        Q_SW_net = self.SW_d[iii] * (1-self.ALBEDO[iii])
        Q_LW_d = self.emissivity_air * self.LW_d[iii]

        # TL_thick = 0.1 # thickness of snow/firn "Top Layer" that energy goes into. Reducing results in higher melt.
        iTL = np.where(z>=self.TL_thick)[0][0] 
        dTL = z[iTL]

        i_GL = np.where(z>=1)[0][0]
        z_GL = z[i_GL]
        m_GL = np.cumsum(mass)[i_GL]
        T_GL = np.cumsum(mass*Tz)[i_GL]/m_GL
        rho_GL = m_GL/z_GL
        K_ice   = 9.828 * np.exp(-0.0057 * T_GL) #[W/m/K]

        K_GL  = K_ice * (rho_GL/RHO_I) ** (2 - 0.5 * (rho_GL/RHO_I))
        
        G = (K_GL * (Tz[i_GL] - Tz[0])/z_GL) # estimated temperature flux in firn due to temperature gradient

        iTL = np.where(z>=0.1)[0][0]
        dTL = z[iTL]

        m = np.cumsum(mass)[iTL] # this was on staging

        # G=0
 
        Qnet = Q_SW_net + Q_LW_d + self.QH[iii] + self.QL[iii] + Qrain_i + G

        a = self.SBC*dt/(CP_I*m)
        b = 0
        c = 0
        d = 1
        e = -1 * (Qnet*dt/(CP_I*m)+T_old)
        p = np.poly1d([a,b,c,d,e])
        r = np.roots(p)
        Tnew = (r[((np.isreal(r)) & (r>0))].real)[0]
        if Tnew<273.15:
            Tsurface = Tnew
            melt_mass = 0
        else:
            Tsurface = 273.15
            melt_mass = (Qnet - self.SBC*273.15**4) * dt / LF_I 

        Tz[0] = Tsurface
        if melt_mass<0:
            melt_mass = 0

        return Tsurface, Tz, melt_mass

### FQS below ###########
#########################

# class FQS_jit:
#     '''
#     Fast Quartic Solver: analytically solves quartic equations (needed to calculate melt)
#     Takes methods from fqs package (@author: NKrvavica)
#     full documentation: https://github.com/NKrvavica/fqs/blob/master/fqs.py
#     '''
#     def __init__(self):
#         pass

#     @njit
#     def single_quadratic(self, a0, b0, c0):
#         ''' 
#         Analytical solver for a single quadratic equation
#         '''
#         a, b = b0 / a0, c0 / a0

#         # Some repating variables
#         a0 = -0.5*a
#         delta = a0*a0 - b
#         sqrt_delta = cmath.sqrt(delta)

#         # Roots
#         r1 = a0 - sqrt_delta
#         r2 = a0 + sqrt_delta

#         return r1, r2


#     @njit
#     def single_cubic(self, a0, b0, c0, d0):
#         ''' 
#         Analytical closed-form solver for a single cubic equation
#         '''
#         a, b, c = b0 / a0, c0 / a0, d0 / a0

#         # Some repeating constants and variables
#         third = 1./3.
#         a13 = a*third
#         a2 = a13*a13
#         sqr3 = math.sqrt(3)

#         # Additional intermediate variables
#         f = third*b - a2
#         g = a13 * (2*a2 - b) + c
#         h = 0.25*g*g + f*f*f

#         def cubic_root(x):
#             ''' Compute cubic root of a number while maintaining its sign'''
#             if x.real >= 0:
#                 return x**third
#             else:
#                 return -(-x)**third

#         if f == g == h == 0:
#             r1 = -cubic_root(c)
#             return r1, r1, r1

#         elif h <= 0:
#             j = math.sqrt(-f)
#             k = math.acos(-0.5*g / (j*j*j))
#             m = math.cos(third*k)
#             n = sqr3 * math.sin(third*k)
#             r1 = 2*j*m - a13
#             r2 = -j * (m + n) - a13
#             r3 = -j * (m - n) - a13
#             return r1, r2, r3

#         else:
#             sqrt_h = cmath.sqrt(h)
#             S = cubic_root(-0.5*g + sqrt_h)
#             U = cubic_root(-0.5*g - sqrt_h)
#             S_plus_U = S + U
#             S_minus_U = S - U
#             r1 = S_plus_U - a13
#             r2 = -0.5*S_plus_U - a13 + S_minus_U*sqr3*0.5j
#             r3 = -0.5*S_plus_U - a13 - S_minus_U*sqr3*0.5j
#             return r1, r2, r3


#     @njit
#     def single_cubic_one(self, a0, b0, c0, d0):
#         ''' 
#         Analytical closed-form solver for a single cubic equation
#         '''
#         a, b, c = b0 / a0, c0 / a0, d0 / a0

#         # Some repeating constants and variables
#         third = 1./3.
#         a13 = a*third
#         a2 = a13*a13

#         # Additional intermediate variables
#         f = third*b - a2
#         g = a13 * (2*a2 - b) + c
#         h = 0.25*g*g + f*f*f

#         def cubic_root(x):
#             ''' Compute cubic root of a number while maintaining its sign
#             '''
#             if x.real >= 0:
#                 return x**third
#             else:
#                 return -(-x)**third

#         if f == g == h == 0:
#             return -cubic_root(c)

#         elif h <= 0:
#             j = math.sqrt(-f)
#             k = math.acos(-0.5*g / (j*j*j))
#             m = math.cos(third*k)
#             return 2*j*m - a13

#         else:
#             sqrt_h = cmath.sqrt(h)
#             S = cubic_root(-0.5*g + sqrt_h)
#             U = cubic_root(-0.5*g - sqrt_h)
#             S_plus_U = S + U
#             return S_plus_U - a13


#     @njit
#     def single_quartic(self, a0, b0, c0, d0, e0):
#         '''
#         Analytical closed-form solver for a single quartic equation
#         '''
#         a, b, c, d = b0/a0, c0/a0, d0/a0, e0/a0

#         # Some repeating variables
#         a0 = 0.25*a
#         a02 = a0*a0

#         # Coefficients of subsidiary cubic euqtion
#         p = 3*a02 - 0.5*b
#         q = a*a02 - b*a0 + 0.5*c
#         r = 3*a02*a02 - b*a02 + c*a0 - d

#         # One root of the cubic equation
#         z0 = self.single_cubic_one(1, p, r, p*r - 0.5*q*q)

#         # Additional variables
#         s = cmath.sqrt(2*p + 2*z0.real + 0j)
#         if s == 0:
#             t = z0*z0 + r
#         else:
#             t = -q / s

#         # Compute roots by quadratic equations
#         r0, r1 = self.single_quadratic(1, s, z0 + t)
#         r2, r3 = self.single_quadratic(1, -s, z0 - t)

#         return r0 - a0, r1 - a0, r2 - a0, r3 - a0


#     def multi_quadratic(self, a0, b0, c0):
#         ''' 
#         Analytical solver for multiple quadratic equations
#         '''
#         a, b = b0 / a0, c0 / a0

#         # Some repating variables
#         a0 = -0.5*a
#         delta = a0*a0 - b
#         sqrt_delta = np.sqrt(delta + 0j)

#         # Roots
#         r1 = a0 - sqrt_delta
#         r2 = a0 + sqrt_delta

#         return r1, r2


#     def multi_cubic(self, a0, b0, c0, d0, all_roots=True):
#         '''
#         Analytical closed-form solver for multiple cubic equations
#         '''
#         a, b, c = b0 / a0, c0 / a0, d0 / a0

#         # Some repeating constants and variables
#         third = 1./3.
#         a13 = a*third
#         a2 = a13*a13
#         sqr3 = math.sqrt(3)

#         # Additional intermediate variables
#         f = third*b - a2
#         g = a13 * (2*a2 - b) + c
#         h = 0.25*g*g + f*f*f

#         # Masks for different combinations of roots
#         m1 = (f == 0) & (g == 0) & (h == 0)     # roots are real and equal
#         m2 = (~m1) & (h <= 0)                   # roots are real and distinct
#         m3 = (~m1) & (~m2)                      # one real root and two complex

#         def cubic_root(x):
#             ''' Compute cubic root of a number while maintaining its sign
#             '''
#             root = np.zeros_like(x)
#             positive = (x >= 0)
#             negative = ~positive
#             root[positive] = x[positive]**third
#             root[negative] = -(-x[negative])**third
#             return root

#         def roots_all_real_equal(c):
#             ''' Compute cubic roots if all roots are real and equal
#             '''
#             r1 = -cubic_root(c)
#             if all_roots:
#                 return r1, r1, r1
#             else:
#                 return r1

#         def roots_all_real_distinct(a13, f, g, h):
#             ''' Compute cubic roots if all roots are real and distinct
#             '''
#             j = np.sqrt(-f)
#             k = np.arccos(-0.5*g / (j*j*j))
#             m = np.cos(third*k)
#             r1 = 2*j*m - a13
#             if all_roots:
#                 n = sqr3 * np.sin(third*k)
#                 r2 = -j * (m + n) - a13
#                 r3 = -j * (m - n) - a13
#                 return r1, r2, r3
#             else:
#                 return r1

#         def roots_one_real(a13, g, h):
#             ''' Compute cubic roots if one root is real and other two are complex
#             '''
#             sqrt_h = np.sqrt(h)
#             S = cubic_root(-0.5*g + sqrt_h)
#             U = cubic_root(-0.5*g - sqrt_h)
#             S_plus_U = S + U
#             r1 = S_plus_U - a13
#             if all_roots:
#                 S_minus_U = S - U
#                 r2 = -0.5*S_plus_U - a13 + S_minus_U*sqr3*0.5j
#                 r3 = -0.5*S_plus_U - a13 - S_minus_U*sqr3*0.5j
#                 return r1, r2, r3
#             else:
#                 return r1

#         # Compute roots
#         if all_roots:
#             roots = np.zeros((3, len(a))).astype(complex)
#             roots[:, m1] = roots_all_real_equal(c[m1])
#             roots[:, m2] = roots_all_real_distinct(a13[m2], f[m2], g[m2], h[m2])
#             roots[:, m3] = roots_one_real(a13[m3], g[m3], h[m3])
#         else:
#             roots = np.zeros(len(a))  # .astype(complex)
#             roots[m1] = roots_all_real_equal(c[m1])
#             roots[m2] = roots_all_real_distinct(a13[m2], f[m2], g[m2], h[m2])
#             roots[m3] = roots_one_real(a13[m3], g[m3], h[m3])

#         return roots


#     def multi_quartic(self, a0, b0, c0, d0, e0):
#         ''' 
#         Analytical closed-form solver for multiple quartic equations
#         '''
#         a, b, c, d = b0/a0, c0/a0, d0/a0, e0/a0

#         # Some repeating variables
#         a0 = 0.25*a
#         a02 = a0*a0

#         # Coefficients of subsidiary cubic euqtion
#         p = 3*a02 - 0.5*b
#         q = a*a02 - b*a0 + 0.5*c
#         r = 3*a02*a02 - b*a02 + c*a0 - d

#         # One root of the cubic equation
#         z0 = self.multi_cubic(1, p, r, p*r - 0.5*q*q, all_roots=False)

#         # Additional variables
#         s = np.sqrt(2*p + 2*z0.real + 0j)
#         t = np.zeros_like(s)
#         mask = (s == 0)
#         t[mask] = z0[mask]*z0[mask] + r[mask]
#         t[~mask] = -q[~mask] / s[~mask]

#         # Compute roots by quadratic equations
#         r0, r1 = self.multi_quadratic(1, s, z0 + t) - a0
#         r2, r3 = self.multi_quadratic(1, -s, z0 - t) - a0

#         return r0, r1, r2, r3


#     def cubic_roots(self, p):
#         '''
#         A caller function for a fast cubic root solver (3rd order polynomial).
#         '''
#         # Convert input to array (if input is a list or tuple)
#         p = np.asarray(p)

#         # If only one set of coefficients is given, add axis
#         if p.ndim < 2:
#             p = p[np.newaxis, :]

#         # Check if four coefficients are given
#         if p.shape[1] != 4:
#             raise ValueError('Expected 3rd order polynomial with 4 '
#                              'coefficients, got {:d}.'.format(p.shape[1]))

#         if p.shape[0] < 100:
#             roots = [self.single_cubic(*pi) for pi in p]
#             return np.array(roots)
#         else:
#             roots = self.multi_cubic(*p.T)
#             return np.array(roots).T


#     def quartic_roots(self, p):
#         '''
#         A caller function for a fast quartic root solver (4th order polynomial).
#         '''
#         # Convert input to an array (if input is a list or tuple)
#         p = np.asarray(p)

#         # If only one set of coefficients is given, add axis
#         if p.ndim < 2:
#             p = p[np.newaxis, :]

#         # Check if all five coefficients are given
#         if p.shape[1] != 5:
#             raise ValueError('Expected 4th order polynomial with 5 '
#                              'coefficients, got {:d}.'.format(p.shape[1]))

#         if p.shape[0] < 100:
#             roots = [self.single_quartic(*pi) for pi in p]
#             return np.array(roots)
#         else:
#             roots = self.multi_quartic(*p.T)
#             return np.array(roots).T

#########################
#########################

# class FQS_nojit:
#     '''
#     # Fast Quartic Solver: analytically solves quartic equations (needed to calculate melt)
#     # Takes methods from fqs package (@author: NKrvavica)
#     # full documentation: https://github.com/NKrvavica/fqs/blob/master/fqs.py
#     '''    

#     def __init__(self):
#         pass

def single_quadratic(a0, b0, c0):
    ''' 
    Analytical solver for a single quadratic equation
    '''
    a, b = b0 / a0, c0 / a0

    # Some repating variables
    a0 = -0.5*a
    delta = a0*a0 - b
    sqrt_delta = np.sqrt(delta)

    # Roots
    r1 = a0 - sqrt_delta
    r2 = a0 + sqrt_delta

    return r1, r2

def single_cubic(a0, b0, c0, d0):
    ''' 
    Analytical closed-form solver for a single cubic equation
    '''
    a, b, c = b0 / a0, c0 / a0, d0 / a0

    # Some repeating constants and variables
    third = 1./3.
    a13 = a*third
    a2 = a13*a13
    sqr3 = np.sqrt(3)

    # Additional intermediate variables
    f = third*b - a2
    g = a13 * (2*a2 - b) + c
    h = 0.25*g*g + f*f*f

    def cubic_root(x):
        ''' Compute cubic root of a number while maintaining its sign'''
        if x.real >= 0:
            return x**third
        else:
            return -(-x)**third

    if f == g == h == 0:
        r1 = -cubic_root(c)
        return r1, r1, r1

    elif h <= 0:
        j = np.sqrt(-f)
        k = np.arccos(-0.5*g / (j*j*j))
        m = np.cos(third*k)
        n = sqr3 * np.sin(third*k)
        r1 = 2*j*m - a13
        r2 = -j * (m + n) - a13
        r3 = -j * (m - n) - a13
        return r1, r2, r3

    else:
        sqrt_h = np.sqrt(h)
        S = cubic_root(-0.5*g + sqrt_h)
        U = cubic_root(-0.5*g - sqrt_h)
        S_plus_U = S + U
        S_minus_U = S - U
        r1 = S_plus_U - a13
        r2 = -0.5*S_plus_U - a13 + S_minus_U*sqr3*0.5j
        r3 = -0.5*S_plus_U - a13 - S_minus_U*sqr3*0.5j
        return r1, r2, r3



def single_cubic_one(a0, b0, c0, d0):
    ''' 
    Analytical closed-form solver for a single cubic equation
    '''
    a, b, c = b0 / a0, c0 / a0, d0 / a0

    # Some repeating constants and variables
    third = 1./3.
    a13 = a*third
    a2 = a13*a13

    # Additional intermediate variables
    f = third*b - a2
    g = a13 * (2*a2 - b) + c
    h = 0.25*g*g + f*f*f

    def cubic_root(x):
        ''' Compute cubic root of a number while maintaining its sign
        '''
        if x.real >= 0:
            return x**third
        else:
            return -(-x)**third

    if f == g == h == 0:
        return -cubic_root(c)

    elif h <= 0:
        j = np.sqrt(-f)
        k = np.arccos(-0.5*g / (j*j*j))
        m = np.cos(third*k)
        return 2*j*m - a13

    else:
        sqrt_h = np.sqrt(h)
        S = cubic_root(-0.5*g + sqrt_h)
        U = cubic_root(-0.5*g - sqrt_h)
        S_plus_U = S + U
        return S_plus_U - a13

def single_quartic(a0, b0, c0, d0, e0):
    '''
    Analytical closed-form solver for a single quartic equation
    '''
    a, b, c, d = b0/a0, c0/a0, d0/a0, e0/a0

    # Some repeating variables
    a0 = 0.25*a
    a02 = a0*a0

    # Coefficients of subsidiary cubic euqtion
    p = 3*a02 - 0.5*b
    q = a*a02 - b*a0 + 0.5*c
    r = 3*a02*a02 - b*a02 + c*a0 - d

    # One root of the cubic equation
    z0 = single_cubic_one(1, p, r, p*r - 0.5*q*q)

    # Additional variables
    s = np.sqrt(2*p + 2*z0.real + 0j)
    if s == 0:
        t = z0*z0 + r
    else:
        t = -q / s

    # Compute roots by quadratic equations
    r0, r1 = single_quadratic(1, s, z0 + t)
    r2, r3 = single_quadratic(1, -s, z0 - t)

    return r0 - a0, r1 - a0, r2 - a0, r3 - a0

def multi_quadratic(a0, b0, c0):
    ''' 
    Analytical solver for multiple quadratic equations
    '''
    a, b = b0 / a0, c0 / a0

    # Some repating variables
    a0 = -0.5*a
    delta = a0*a0 - b
    sqrt_delta = np.sqrt(delta + 0j)

    # Roots
    r1 = a0 - sqrt_delta
    r2 = a0 + sqrt_delta

    return r1, r2


def multi_cubic(a0, b0, c0, d0, all_roots=True):
    '''
    Analytical closed-form solver for multiple cubic equations
    '''
    a, b, c = b0 / a0, c0 / a0, d0 / a0

    # Some repeating constants and variables
    third = 1./3.
    a13 = a*third
    a2 = a13*a13
    sqr3 = np.sqrt(3)

    # Additional intermediate variables
    f = third*b - a2
    g = a13 * (2*a2 - b) + c
    h = 0.25*g*g + f*f*f

    # Masks for different combinations of roots
    m1 = (f == 0) & (g == 0) & (h == 0)     # roots are real and equal
    m2 = (~m1) & (h <= 0)                   # roots are real and distinct
    m3 = (~m1) & (~m2)                      # one real root and two complex

    def cubic_root(x):
        ''' Compute cubic root of a number while maintaining its sign
        '''
        root = np.zeros_like(x)
        positive = (x >= 0)
        negative = ~positive
        root[positive] = x[positive]**third
        root[negative] = -(-x[negative])**third
        return root

    def roots_all_real_equal(c):
        ''' Compute cubic roots if all roots are real and equal
        '''
        r1 = -cubic_root(c)
        if all_roots:
            return r1, r1, r1
        else:
            return r1

    def roots_all_real_distinct(a13, f, g, h):
        ''' Compute cubic roots if all roots are real and distinct
        '''
        j = np.sqrt(-f)
        k = np.arccos(-0.5*g / (j*j*j))
        m = np.cos(third*k)
        r1 = 2*j*m - a13
        if all_roots:
            n = sqr3 * np.sin(third*k)
            r2 = -j * (m + n) - a13
            r3 = -j * (m - n) - a13
            return r1, r2, r3
        else:
            return r1

    def roots_one_real(a13, g, h):
        ''' Compute cubic roots if one root is real and other two are complex
        '''
        sqrt_h = np.sqrt(h)
        S = cubic_root(-0.5*g + sqrt_h)
        U = cubic_root(-0.5*g - sqrt_h)
        S_plus_U = S + U
        r1 = S_plus_U - a13
        if all_roots:
            S_minus_U = S - U
            r2 = -0.5*S_plus_U - a13 + S_minus_U*sqr3*0.5j
            r3 = -0.5*S_plus_U - a13 - S_minus_U*sqr3*0.5j
            return r1, r2, r3
        else:
            return r1

    # Compute roots
    if all_roots:
        roots = np.zeros((3, len(a))).astype(complex)
        roots[:, m1] = roots_all_real_equal(c[m1])
        roots[:, m2] = roots_all_real_distinct(a13[m2], f[m2], g[m2], h[m2])
        roots[:, m3] = roots_one_real(a13[m3], g[m3], h[m3])
    else:
        roots = np.zeros(len(a))  # .astype(complex)
        roots[m1] = roots_all_real_equal(c[m1])
        roots[m2] = roots_all_real_distinct(a13[m2], f[m2], g[m2], h[m2])
        roots[m3] = roots_one_real(a13[m3], g[m3], h[m3])

    return roots


def multi_quartic(a0, b0, c0, d0, e0):
    ''' 
    Analytical closed-form solver for multiple quartic equations
    '''
    a, b, c, d = b0/a0, c0/a0, d0/a0, e0/a0

    # Some repeating variables
    a0 = 0.25*a
    a02 = a0*a0

    # Coefficients of subsidiary cubic euqtion
    p = 3*a02 - 0.5*b
    q = a*a02 - b*a0 + 0.5*c
    r = 3*a02*a02 - b*a02 + c*a0 - d

    # One root of the cubic equation
    z0 = multi_cubic(1, p, r, p*r - 0.5*q*q, all_roots=False)

    # Additional variables
    s = np.sqrt(2*p + 2*z0.real + 0j)
    t = np.zeros_like(s)
    mask = (s == 0)
    t[mask] = z0[mask]*z0[mask] + r[mask]
    t[~mask] = -q[~mask] / s[~mask]

    # Compute roots by quadratic equations
    r0, r1 = multi_quadratic(1, s, z0 + t) - a0
    r2, r3 = multi_quadratic(1, -s, z0 - t) - a0

    return r0, r1, r2, r3


def cubic_roots(p):
    '''
    A caller function for a fast cubic root solver (3rd order polynomial).
    '''
    # Convert input to array (if input is a list or tuple)
    p = np.asarray(p)

    # If only one set of coefficients is given, add axis
    if p.ndim < 2:
        p = p[np.newaxis, :]

    # Check if four coefficients are given
    if p.shape[1] != 4:
        raise ValueError('Expected 3rd order polynomial with 4 '
                        'coefficients, got {:d}.'.format(p.shape[1]))

    if p.shape[0] < 100:
        roots = [single_cubic(*pi) for pi in p]
        return np.array(roots)
    else:
        roots = multi_cubic(*p.T)
        return np.array(roots).T


def quartic_roots(p):
    '''
    A caller function for a fast quartic root solver (4th order polynomial).
    p[0]*x^4 + p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4] = 0
    '''
    # Convert input to an array (if input is a list or tuple)
    p = np.asarray(p)

    # If only one set of coefficients is given, add axis
    if p.ndim < 2:
        p = p[np.newaxis, :]

    # Check if all five coefficients are given
    if p.shape[1] != 5:
        raise ValueError('Expected 4th order polynomial with 5 '
                        'coefficients, got {:d}.'.format(p.shape[1]))

    if p.shape[0] < 100:
        roots = [single_quartic(*pi) for pi in p]
        return np.array(roots)
    else:
        roots = multi_quartic(*p.T)
        return np.array(roots).T

def multi_quadratic(a0, b0, c0):
    ''' 
    Analytical solver for multiple quadratic equations
    '''
    a, b = b0 / a0, c0 / a0

    # Some repating variables
    a0 = -0.5*a
    delta = a0*a0 - b
    sqrt_delta = np.sqrt(delta + 0j)

    # Roots
    r1 = a0 - sqrt_delta
    r2 = a0 + sqrt_delta

    return r1, r2


def multi_cubic(a0, b0, c0, d0, all_roots=True):
    '''
    Analytical closed-form solver for multiple cubic equations
    '''
    a, b, c = b0 / a0, c0 / a0, d0 / a0

    # Some repeating constants and variables
    third = 1./3.
    a13 = a*third
    a2 = a13*a13
    sqr3 = math.sqrt(3)

    # Additional intermediate variables
    f = third*b - a2
    g = a13 * (2*a2 - b) + c
    h = 0.25*g*g + f*f*f

    # Masks for different combinations of roots
    m1 = (f == 0) & (g == 0) & (h == 0)     # roots are real and equal
    m2 = (~m1) & (h <= 0)                   # roots are real and distinct
    m3 = (~m1) & (~m2)                      # one real root and two complex

    def cubic_root(x):
        ''' Compute cubic root of a number while maintaining its sign
        '''
        root = np.zeros_like(x)
        positive = (x >= 0)
        negative = ~positive
        root[positive] = x[positive]**third
        root[negative] = -(-x[negative])**third
        return root

    def roots_all_real_equal(c):
        ''' Compute cubic roots if all roots are real and equal
        '''
        r1 = -cubic_root(c)
        if all_roots:
            return r1, r1, r1
        else:
            return r1

    def roots_all_real_distinct(a13, f, g, h):
        ''' Compute cubic roots if all roots are real and distinct
        '''
        j = np.sqrt(-f)
        k = np.arccos(-0.5*g / (j*j*j))
        m = np.cos(third*k)
        r1 = 2*j*m - a13
        if all_roots:
            n = sqr3 * np.sin(third*k)
            r2 = -j * (m + n) - a13
            r3 = -j * (m - n) - a13
            return r1, r2, r3
        else:
            return r1

    def roots_one_real(a13, g, h):
        ''' Compute cubic roots if one root is real and other two are complex
        '''
        sqrt_h = np.sqrt(h)
        S = cubic_root(-0.5*g + sqrt_h)
        U = cubic_root(-0.5*g - sqrt_h)
        S_plus_U = S + U
        r1 = S_plus_U - a13
        if all_roots:
            S_minus_U = S - U
            r2 = -0.5*S_plus_U - a13 + S_minus_U*sqr3*0.5j
            r3 = -0.5*S_plus_U - a13 - S_minus_U*sqr3*0.5j
            return r1, r2, r3
        else:
            return r1

    # Compute roots
    if all_roots:
        roots = np.zeros((3, len(a))).astype(complex)
        roots[:, m1] = roots_all_real_equal(c[m1])
        roots[:, m2] = roots_all_real_distinct(a13[m2], f[m2], g[m2], h[m2])
        roots[:, m3] = roots_one_real(a13[m3], g[m3], h[m3])
    else:
        roots = np.zeros(len(a))  # .astype(complex)
        roots[m1] = roots_all_real_equal(c[m1])
        roots[m2] = roots_all_real_distinct(a13[m2], f[m2], g[m2], h[m2])
        roots[m3] = roots_one_real(a13[m3], g[m3], h[m3])

    return roots


def multi_quartic(a0, b0, c0, d0, e0):
    ''' 
    Analytical closed-form solver for multiple quartic equations
    '''
    a, b, c, d = b0/a0, c0/a0, d0/a0, e0/a0

    # Some repeating variables
    a0 = 0.25*a
    a02 = a0*a0

    # Coefficients of subsidiary cubic euqtion
    p = 3*a02 - 0.5*b
    q = a*a02 - b*a0 + 0.5*c
    r = 3*a02*a02 - b*a02 + c*a0 - d

    # One root of the cubic equation
    z0 = multi_cubic(1, p, r, p*r - 0.5*q*q, all_roots=False)

    # Additional variables
    s = np.sqrt(2*p + 2*z0.real + 0j)
    t = np.zeros_like(s)
    mask = (s == 0)
    t[mask] = z0[mask]*z0[mask] + r[mask]
    t[~mask] = -q[~mask] / s[~mask]

    # Compute roots by quadratic equations
    r0, r1 = multi_quadratic(1, s, z0 + t) - a0
    r2, r3 = multi_quadratic(1, -s, z0 - t) - a0

    return r0, r1, r2, r3


def cubic_roots(p):
    '''
    A caller function for a fast cubic root solver (3rd order polynomial).
    '''
    # Convert input to array (if input is a list or tuple)
    p = np.asarray(p)

    # If only one set of coefficients is given, add axis
    if p.ndim < 2:
        p = p[np.newaxis, :]

    # Check if four coefficients are given
    if p.shape[1] != 4:
        raise ValueError('Expected 3rd order polynomial with 4 '
                        'coefficients, got {:d}.'.format(p.shape[1]))

    if p.shape[0] < 100:
        roots = [single_cubic(*pi) for pi in p]
        return np.array(roots)
    else:
        roots = multi_cubic(*p.T)
        return np.array(roots).T


def quartic_roots(p):
    '''
    A caller function for a fast quartic root solver (4th order polynomial).
    p[0]*x^4 + p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4] = 0
    '''
    # Convert input to an array (if input is a list or tuple)
    p = np.asarray(p)

    # If only one set of coefficients is given, add axis
    if p.ndim < 2:
        p = p[np.newaxis, :]

    # Check if all five coefficients are given
    if p.shape[1] != 5:
        raise ValueError('Expected 4th order polynomial with 5 '
                        'coefficients, got {:d}.'.format(p.shape[1]))

    if p.shape[0] < 100:
        roots = [single_quartic(*pi) for pi in p]
        return np.array(roots)
    else:
        roots = multi_quartic(*p.T)
        return np.array(roots).T


'''
References
Cuffey and Paterson, p. 140-150
Hock 2005
van As, 2005
Van Pelt, 2012
Klok, 2002 
Born 2019
'''