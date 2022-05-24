#!/usr/bin/env python

import numpy as np 
# from solver import solver
from solver import transient_solve_TR
# from Gasses import gasses 
# from Sites import sites 
from reader import read_input, read_init
import json
import scipy.interpolate as interpolate
from scipy import optimize
from constants import *
import os
import sys

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

    def __init__(self,config,climateTS,start_ind):
        '''
        intialalize seb
        '''

        self.c = config
        
        self.time_in = climateTS['time'][start_ind:]
        self.SW_d    = climateTS['SW_d'][start_ind:]
        self.ALBEDO  = climateTS['ALBEDO'][start_ind:]
        self.T2m     = climateTS['T2m'][start_ind:]
        if (('EVAP' in climateTS.keys()) and ('SUBLIM' in climateTS.keys())):
            self.EVAP    = climateTS['EVAP'][start_ind:]
            self.SUBLIM    = climateTS['SUBLIM'][start_ind:]
        elif 'SUBLIM' in climateTS.keys():
            self.SUBLIM    = climateTS['SUBLIM'][start_ind:]
            self.EVAP    = np.zeros_like(self.SUBLIM)
        elif 'EVAP' in climateTS.keys():
            self.EVAP    = climateTS['EVAP'][start_ind:]
            self.SUBLIM    = np.zeros_like(self.EVAP)
        self.QH      = climateTS['QH'][start_ind:]
        self.QL      = climateTS['QL'][start_ind:]
        self.RAIN    = climateTS['RAIN'][start_ind:]
        # need to account for cold snow falling on warmer surface

        # self.EP      = np.zeros_like(self.SW_d)
        # self.G       = np.zeros_like(self.SW_d) # For now we are not considering any flux in/out of the upper model node from below
        
        
        self.SBC = 5.67e-8 # Stefan-Boltzmann constant [W K^-4 m^-2
        self.emissivity_air = 1
        self.emissivity_snow = 1
        
        # self.D_sh = 15 # Sensible heat flux coefficient, Born et al. 2019 [W m^-2 K^-1] 


    def SEB(self, PhysParams,iii):
        '''
        Calculate the surface energy budget
        Positive fluxes are into the surface, negative are out
        SEBparams: mass, Tz, dt
        '''

        # def enet(Tsurf,Qn):
        #     return np.abs(Qn - 5.67e-8 * Tsurf**4)

        def enet(Tsurf,Qn,Tz,dz):
            # e1 = np.abs(Qn - 5.670374419e-8 * Tsurf**4)
            gflux = np.abs(0.3*(Tz - Tsurf)/dz)
            # gflux = 0
            e1 = np.abs(Qn - 5.670374419e-8 * Tsurf**4) + gflux 
            return e1

        Tz   = PhysParams['Tz']
        mass = PhysParams['mass']
        dt   = PhysParams['dt']

        T_rain = np.max((self.T2m[iii],T_MELT))
        Qrain_i = 0

        ### CHECK THAT THE RAIN UNITS ARE CORRECT!
        # Qrain_i = RHO_W_KGM * CP_W * self.RAIN[iii] * (T_rain - T_MELT) # Assume rain temperature is air temp, Hock 2005, eq 19
            # heat for rain falling on top of cold snow should be handled in melt.py
        ### CHECK THAT THE RAIN UNITS ARE CORRECT!

        Q_SW_net = self.SW_d[iii] * (1-self.ALBEDO[iii])
        Q_LW_d = self.SBC * (self.emissivity_air * self.T2m[iii]**4)
        Qnet = Q_SW_net + Q_LW_d - self.QH[iii] - self.QL[iii] + Qrain_i

        Tlast=Tz[1]
        dz=PhysParams['dz'][1]
        cold_content = CP_I * mass * (T_MELT - Tz) # cold content [J]
        # LH_melt = LF_I * mass
        # esum = cold_content + LH_melt
        # ccs = np.cumsum(esum)

        mresults = optimize.minimize(enet,method = 'Nelder-Mead',x0=self.T2m[iii],args=(Qnet,Tlast,dz),tol=1e-6)
        Tsurface = mresults.x[0]

        if Tsurface>T_MELT:
            print('One')
            input('waiting')
            Q_LW_me = (self.SBC * (self.emissivity_air * self.T2m[iii]**4 - self.emissivity_snow * (T_MELT)**4))
            G = 0.3*(Tlast - T_MELT)/dz
            E_melt = (Q_SW_net + Q_LW_me - QH - QL + Qrain_i + G)

            E_iter = E_melt.copy()
            kk = 0
            melt_mass = 0
            
            while E_iter > 0:
                if E_iter <= cold_content[kk]: # all energy warms/cools this firn; no melt in it
                    Tz[kk] = Tz[kk] + E_net / (CP_I * mass[kk]) * dt
                    E_iter = 0

                elif E_net > cold_content[kk]:
                    Tz[kk] = T_MELT
                    E_melt = E_net - cold_content[kk] #energy available for melting
                    melt_mass_pot = E_melt/LF_I # mass that can be melted with the energy available              

                    if melt_mass_pot <= mass[kk]: #less than the entire node melts
                        melt_mass += melt_mass_pot
                        E_iter = 0

                    else: #entire node melts
                        melt_mass += mass[kk]
                        E_iter -= mass[kk]*LF_I

                kk+=1
            Tsurface = T_MELT
        else:
            # print('Two')
            melt_mass = 0
            # Tz[0] = Tsurface

        # print(iii)
        # print('Tz', Tz[0])
        # print('Melt mass:',melt_mass)
        # input('waiting')
        return Tsurface, melt_mass


    def SEB_direct(self, PhysParams,iii):
        '''
        Calculate the surface energy budget
        Positive fluxes are into the surface, negative are out
        SEBparams: mass, Tz, dt
        '''

        # max_iter = 20
        # ijk = 0
        Tz   = PhysParams['Tz']
        mass = PhysParams['mass']
        dt   = PhysParams['dt']
        print('dt',dt)
        
        T_rain = np.max((self.T2m[iii],T_MELT))
        Qrain_i = RHO_W_KGM * CP_W * self.RAIN[iii] * (T_rain - T_MELT) # Assume rain temperature is air temp, Hock 2005, eq 19
        # heat for rain falling on top of cold snow should be handled in melt.py

        Q_SW_i = self.SW_d[iii] * (1-self.ALBEDO[iii])

        emissivity_air = emissivity_snow = 1
        Q_LW_i = self.SBC * (emissivity_air * self.T2m[iii]**4 - emissivity_snow * (Tz[0])**4)

        E_net = Q_SW_i + Q_LW_i + self.QH[iii] + self.QL[iii] + Qrain_i + self.G[iii]
        print('E_net',E_net)
        # time dimension for conversion?
        # dTz[0] = E_net * dt[iii] / CP_I / mass

        cold_content    = CP_I * mass * (T_MELT - Tz)           # cold content [J]
        print('cold_content',cold_content)
        print('mass',mass[0:5])


        E_iter = E_net.copy()
        kk = 0
        melt_mass = 0
        
        while E_iter > 0:
            if E_iter <= cold_content[kk]: # all energy warms/cools the firn; no melt
                Tz[kk] = Tz[kk] + E_net / (CP_I * mass[kk]) * dt
                E_iter = 0

            elif E_net > cold_content[kk]:
                Tz[kk] = T_MELT
                E_melt = E_net - cold_content[kk] #energy available for melting
                melt_mass_pot = E_melt/LF_I # mass that can be melted                

                if melt_mass_pot <= mass[kk]: #less than the entire node melts
                    melt_mass += melt_mass_pot
                    E_iter = 0

                else: #entire node melts
                    melt_mass += mass[kk]
                    E_iter -= mass[kk]*LF_I

            kk+=1

        print(Tz[0:5])
        input('SEB')
        return Tz, melt_mass








'''
References
Cuffey and Paterson, p. 140-150
Hock 2005
van As, 2005
Van Pelt, 2012
Klok, 2002 
Born 2019
'''