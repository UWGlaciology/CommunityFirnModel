#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
physics.py
==========

The Community Firn Model physics module.

equations for each of the densification models.
'''


import numpy as np
from constants import *
from scipy import interpolate
import sys
import numpy.polynomial.polynomial as poly



class FirnPhysics:

    '''
    class for the various firn densification schemes

    The standard parameters that get passed are: iii, steps, gridLen, bdotSec, bdot_mean, bdot_type, Tz, T10m, rho, sigma, dt, Ts, r2, physGrain

    if you want to add physics that require more parameters, you need to change the 'PhysParams' dictionary in both the spin and nospin classes.

    :param steps: # of steps per year
    :param gridLen:
    :param bdotSec:
    :param Tz:
    :param rho:
    :param sigma:

    :return drho_dt:
    :return viscosity:

    '''

    def __init__(self,PhysParams):
        '''
        The init docstring
        bdot_mean units are m I.E./year
        bdotSec units are m I.E./time step
        '''
        for k,v in list(PhysParams.items()):
            setattr(self,k,v)
        self.RD = {} # RD = Return Dictionary, set up this way so that more things can be returned easily if needed.


    def HL_dynamic(self):
        '''
        Accumulation units are m W.E. per year
        viscosity doesn't work
        '''
        Q1  = 10160.0
        Q2  = 21400.0
        k1  = 11.0
        k2  = 575.0
        aHL = 1.0
        bHL = 0.5

        A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM # Accumulation in units m W.E. per year
        A_mean = self.bdot_mean * RHO_I_MGM

        drho_dt = np.zeros(self.gridLen)
        viscosity = np.zeros(self.gridLen)
        
        if A_instant<0:
            A_instant = 0.0

        if self.bdot_type == 'instant':
            if A_instant<0:
                A_instant = 0.0
            if (self.FirnAir and self.AirRunType=='steady'):
                Tcon = self.steady_T * np.ones_like(self.Tz)
                drho_dt[self.rho < RHO_1]     = k1 * np.exp(-Q1 / (R * Tcon[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * A_instant**aHL * 1000 / S_PER_YEAR
                drho_dt[self.rho >= RHO_1]    = k2 * np.exp(-Q2 / (R * Tcon[self.rho >= RHO_1])) * (RHO_I_MGM - self.rho[self.rho >= RHO_1] / 1000) * A_instant**bHL * 1000 / S_PER_YEAR
            else:
                drho_dt[self.rho < RHO_1]     = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * A_instant**aHL * 1000 / S_PER_YEAR
                drho_dt[self.rho >= RHO_1]    = k2 * np.exp(-Q2 / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I_MGM - self.rho[self.rho >= RHO_1] / 1000) * A_instant**bHL * 1000 / S_PER_YEAR
                
                viscosity[self.rho < RHO_1]   = (self.rho[self.rho < RHO_1]* self.sigma[self.rho < RHO_1])/ (2 )/drho_dt[self.rho < RHO_1] 
                viscosity[self.rho >= RHO_1]  = (self.rho[self.rho >= RHO_1] * self.sigma[self.rho >= RHO_1] ) / (2)/drho_dt[self.rho >= RHO_1]

        elif self.bdot_type == 'mean':
            drho_dt[self.rho < RHO_1]     = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * (A_mean[self.rho < RHO_1])**aHL * 1000 / S_PER_YEAR
            drho_dt[self.rho >= RHO_1]    = k2 * np.exp(-Q2 / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I_MGM - self.rho[self.rho >= RHO_1] / 1000) * (A_mean[self.rho >= RHO_1])**bHL * 1000 / S_PER_YEAR

            # viscosity[self.rho < RHO_1]   = (self.rho[self.rho < RHO_1]* self.sigma[self.rho < RHO_1])/ (2 )/drho_dt[self.rho < RHO_1] 
            # viscosity[self.rho >= RHO_1]  = (self.rho[self.rho >= RHO_1] * self.sigma[self.rho >= RHO_1] ) / (2)/drho_dt[self.rho >= RHO_1]

        self.RD['drho_dt']   = drho_dt

        self.RD['viscosity'] = viscosity
        return self.RD
    ### end HL_dynamic ###
    ######################

    def HL_Sigfus(self):
        '''
        Accumulation units are m W.E. per year (zone 1); uses stress for zone 2 
        '''
        Q1  = 10160.0
        Q2  = 21400.0
        k1  = 11.0
        k2  = 575.0
        aHL = 1.0

        A_instant   = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean      = self.bdot_mean * RHO_I_MGM

        drho_dt     = np.zeros(self.gridLen)
        viscosity   = np.zeros(self.gridLen)
        f550        = interpolate.interp1d(self.rho, self.sigma)
        sigma550    = f550(RHO_1)
        rhoDiff     = (RHO_I_MGM - self.rho / 1000)

        z1mask = (self.rho < RHO_1)
        z2mask = ((self.rho >= RHO_1) & (self.rho < RHO_I))
        zImask = (self.rho >= RHO_I)

        kSig        = (k2 * np.exp(-1 * Q2 / (R * self.Tz[z2mask])))**2 / S_PER_YEAR
        sigmaDiff   = (self.sigma[z2mask] - sigma550)

        if self.bdot_type == 'instant':
            drho_dt[z1mask] = k1 * np.exp(-Q1 / (R * self.Tz[z1mask])) * (RHO_I_MGM - self.rho[z1mask] / 1000) * A_instant**aHL * 1000 / S_PER_YEAR
        elif self.bdot_type == 'mean':
            drho_dt[z1mask] = k1 * np.exp(-Q1 / (R * self.Tz[z1mask])) * (RHO_I_MGM - self.rho[z1mask] / 1000) * (A_mean[z1mask])**aHL * 1000 / S_PER_YEAR

        drho_dt[z2mask]  = kSig * (sigmaDiff * rhoDiff[z2mask]) / (GRAVITY * np.log((RHO_I_MGM - RHO_1 / 1000) / (rhoDiff[(self.rho >= RHO_1) & (self.rho < RHO_I)])))
        drho_dt[zImask]  = 0
        
        viscosity[self.rho < RHO_1]   = (self.rho[self.rho < RHO_1]/ (2 * self.sigma[self.rho < RHO_1]))/drho_dt[self.rho < RHO_1] 
        viscosity[self.rho >= RHO_1]  = (self.rho[self.rho >= RHO_1] / (2 * self.sigma[self.rho >= RHO_1] ))/drho_dt[self.rho >= RHO_1]

        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        return self.RD
    ### end HL_Sigfus ###
    #####################

    def Li_2004(self):
        '''
        Accumulation units are m W.E. per year (?)
        Equation from Arthern, 2010 (eq. 2): not sure where Rob got that? 
        (Arthern implies accumulation is m I.E./year for bdot; not doing that here.)
        Paper would lead me to believe that it is m W.E.
        
        Needs to have the vapor flux coded in if we want to use these physics properly.
        '''

        A_instant   = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean      = self.bdot_mean * RHO_I_MGM
        viscosity   = np.zeros(self.gridLen)
        
        if self.bdot_type == 'instant':
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Li and Zwally 2004 physics")
            dr_dt = (RHO_I - self.rho) * A_instant * (139.21 - 0.542 * self.T_mean[self.iii]) * 8.36 * (K_TO_C - self.Tz) ** -2.061
        elif self.bdot_type == 'mean':
            dr_dt = (RHO_I - self.rho) * A_mean * (139.21 - 0.542 * self.T_mean[self.iii]) * 8.36 * (K_TO_C - self.Tz) ** -2.061   

        drho_dt = dr_dt / S_PER_YEAR
        
        viscosity[self.rho < RHO_1]     = (self.rho[self.rho < RHO_1] * self.sigma[self.rho < RHO_1])/ (2)/drho_dt[self.rho < RHO_1] 
        viscosity[self.rho >= RHO_1]    = (self.rho[self.rho >= RHO_1] * self.sigma[self.rho >= RHO_1] ) / (2)/drho_dt[self.rho >= RHO_1]  
    
        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        return self.RD
    ### end Li_2004 ###
    ###################

    def Li_2011(self):
        '''
        Accumulation units are m W.E. per year (email correspondence with J. Li, 12/3/13)
        Temperature in the equation for beta is in C.
        Temperature in the 8.36(273.15-T)**-2.061 is in K. (implied in Arthern)

        beta should be calculated with the long-term mean accumulation rate
        '''

        A_instant   = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean      = self.bdot_mean * RHO_I_MGM

        TmC   = self.T_mean[self.iii] - K_TO_C

        dr_dt       = np.zeros(self.gridLen)
        viscosity   = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Li and Zwally 2011 physics")
            beta1 = -9.788 + 8.996 * A_instant - 0.6165 * TmC
            beta2 = beta1 / (-2.0178 + 8.4043 * A_instant - 0.0932 * TmC)

            dr_dt[self.rho <= RHO_1] = (RHO_I - self.rho[self.rho <= RHO_1]) * A_instant * beta1 * 8.36 * (K_TO_C - self.Tz[self.rho <= RHO_1]) ** -2.061
            dr_dt[self.rho > RHO_1]  = (RHO_I - self.rho[self.rho > RHO_1]) * A_instant * beta2 * 8.36 * (K_TO_C - self.Tz[self.rho > RHO_1]) ** -2.061
        
        elif self.bdot_type == 'mean':

            ### These lines are using a different beta for each node
            # beta1 = -9.788 + 8.996 * A_mean - 0.6165 * TmC
            # beta2 = beta1 / (-2.0178 + 8.4043 * A_mean - 0.0932 * TmC)

            ### These lines are for a single value of beta based on long-term accumulation rate

            beta1a = -9.788 + 8.996 * np.mean(A_mean) - 0.6165 * TmC
            beta2a = beta1a / (-2.0178 + 8.4043 * np.mean(A_mean) - 0.0932 * TmC)
            beta1 = np.ones(len(A_mean))*beta1a
            beta2 = np.ones(len(A_mean))*beta2a

            dr_dt[self.rho <= RHO_1] = (RHO_I - self.rho[self.rho <= RHO_1]) * A_mean[self.rho <= RHO_1] * beta1[self.rho <= RHO_1] * 8.36 * (K_TO_C - self.Tz[self.rho <= RHO_1]) ** -2.061
            dr_dt[self.rho > RHO_1]  = (RHO_I - self.rho[self.rho > RHO_1]) * A_mean[self.rho > RHO_1] * beta2[self.rho > RHO_1] * 8.36 * (K_TO_C - self.Tz[self.rho > RHO_1]) ** -2.061

        drho_dt = dr_dt / S_PER_YEAR
        
        viscosity[self.rho < RHO_1]   = (self.rho[self.rho < RHO_1]* self.sigma[self.rho < RHO_1])/ (2)/drho_dt[self.rho < RHO_1] 
        viscosity[self.rho >= RHO_1]  = (self.rho[self.rho >= RHO_1] * self.sigma[self.rho >= RHO_1] ) / (2)/drho_dt[self.rho >= RHO_1]  
          
        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        return self.RD
    ### end Li_2011 ###
    ###################

    def Li_2015(self):

        
        '''
        Accumulation units are m W.E. per year (email correspondence with J. Li, 12/3/13)
        Temperature in the equation for beta is in C.
        Temperature in the 8.36(273.15-T)**-2.061 is in K. (implied in Arthern)

        beta should be calculated with the long-term mean accumulation rate
        '''

        #Thanks to Jiajen Chen for pointing out that there were updated beta parameters in Li and Zwally, 2015

        A_instant     = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean = self.bdot_mean * RHO_I_MGM
        # A_mean = 0.45*np.ones_like(A_mean)
        # print(A_mean)

        # TmC   = self.T10m - K_TO_C
        TmC   = self.T_mean[self.iii] - K_TO_C 
        # TmC = -17.267
        # A_meanLZ = 0.45

        dr_dt = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Li and Zwally 2011 physics")
#            beta1 = -9.788 + 8.996 * A_instant - 0.6165 * TmC
#            beta2 = beta1 / (-2.0178 + 8.4043 * A_instant - 0.0932 * TmC)
            beta1a = -1.218 - 0.403 * TmC
            beta2a = beta1a * (0.792 - 1.080 * np.mean(A_mean) + 0.00465 * TmC)

            dr_dt[self.rho <= RHO_1] = (RHO_I - self.rho[self.rho <= RHO_1]) * A_instant * beta1 * 8.36 * (K_TO_C - self.Tz[self.rho <= RHO_1]) ** -2.061
            dr_dt[self.rho > RHO_1]  = (RHO_I - self.rho[self.rho > RHO_1]) * A_instant * beta2 * 8.36 * (K_TO_C - self.Tz[self.rho > RHO_1]) ** -2.061
        
        elif self.bdot_type == 'mean':

            ### These lines are using a different beta for each node
            # beta1 = -9.788 + 8.996 * A_mean - 0.6165 * TmC
            # beta2 = beta1 / (-2.0178 + 8.4043 * A_mean - 0.0932 * TmC)

            ### These lines are for a single value of beta based on long-term accumulation rate

#            beta1a = -9.788 + 8.996 * np.mean(A_mean) - 0.6165 * TmC
#            beta2a = beta1a / (-2.0178 + 8.4043 * np.mean(A_mean) - 0.0932 * TmC)
            beta1a = -1.218 - 0.403 * TmC
            beta2a = beta1a * (0.792 - 1.080 * np.mean(A_mean) + 0.00465 * TmC)
            beta1 = np.ones(len(A_mean))*beta1a
            beta2 = np.ones(len(A_mean))*beta2a

            dr_dt[self.rho <= RHO_1] = (RHO_I - self.rho[self.rho <= RHO_1]) * A_mean[self.rho <= RHO_1] * beta1[self.rho <= RHO_1] * 8.36 * (K_TO_C - self.Tz[self.rho <= RHO_1]) ** -2.061
            dr_dt[self.rho > RHO_1]  = (RHO_I - self.rho[self.rho > RHO_1]) * A_mean[self.rho > RHO_1] * beta2[self.rho > RHO_1] * 8.36 * (K_TO_C - self.Tz[self.rho > RHO_1]) ** -2.061

        drho_dt = dr_dt / S_PER_YEAR
        # self.viscosity = np.ones(self.gridLen)
        self.RD['drho_dt'] = drho_dt
        return self.RD
    ### end Li_2015 ###
    ###################

    def Arthern_2010S(self):
        '''
        This is the steady-state solution described in the main text of Arthern et al. (2010)
        Accumulation units are kg/m^2/year

        '''

        ar1 = 0.07
        ar2 = 0.03
        Ec  = 60.0e3
        Eg  = 42.4e3

        A_instant   = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM * 1000
        A_mean_1    = self.bdot_mean[self.rho < RHO_1] * RHO_I_MGM * 1000
        A_mean_2    = self.bdot_mean[self.rho >= RHO_1] * RHO_I_MGM * 1000
        dr_dt       = np.zeros(self.gridLen)
        viscosity   = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Arthern 2010 physics")
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean[self.iii]))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean[self.iii]))
        elif self.bdot_type == 'mean':
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean[self.iii]))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean[self.iii]))

        drho_dt = dr_dt / S_PER_YEAR
        
        viscosity[self.rho < RHO_1]   = (self.rho[self.rho < RHO_1]/ (2 * self.sigma[self.rho < RHO_1]))/drho_dt[self.rho < RHO_1] 
        viscosity[self.rho >= RHO_1]  = (self.rho[self.rho >= RHO_1] / (2 * self.sigma[self.rho >= RHO_1] ))/drho_dt[self.rho >= RHO_1]  
          
        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        return self.RD
    ### end Arthern_2010S ###
    #########################

    def Arthern_2010T(self):
        '''
        This is the transient solution described in the appendix of Arthern et al. (2010)

        Uses stress rather than accumulation rate.

        '''

        kc1 = 9.2e-9
        kc2 = 3.7e-9
        Ec  = 60.0e3       

        if not self.physGrain:
           print("Grain growth should be on for Arthern Transient")
           return

        drho_dt     = np.zeros(self.gridLen)
        viscosity   = np.zeros(self.gridLen)
    
        drho_dt[self.rho < RHO_1]  = kc1 * (RHO_I - self.rho[self.rho < RHO_1]) * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1])) * self.sigma[self.rho < RHO_1] / (self.r2[self.rho < RHO_1]) 
        drho_dt[self.rho >= RHO_1] = kc2 * (RHO_I - self.rho[self.rho >= RHO_1]) * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1])) * self.sigma[self.rho >= RHO_1] / (self.r2[self.rho >= RHO_1])
       
        viscosity[self.rho < RHO_1]   = (self.rho[self.rho < RHO_1]* self.sigma[self.rho < RHO_1])/ (2 )/drho_dt[self.rho < RHO_1] 
        viscosity[self.rho >= RHO_1]  = (self.rho[self.rho >= RHO_1] * self.sigma[self.rho >= RHO_1] )/ (2 )/drho_dt[self.rho >= RHO_1]  
          
        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        return self.RD
    ### end Arthern_2010T ###
    #########################

    def Helsen_2008(self):
        '''
        Accumulation units are m W.E. per year (?)
        Equation is from Arthern et al. 2010 (2)
        (Arthern implies units are m I.E.; not doing that here) 

        '''

        A_instant   = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean      = self.bdot_mean * RHO_I_MGM
        viscosity   = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Helsen 2008 physics")            
            dr_dt = (RHO_I - self.rho) * A_instant * (76.138 - 0.28965 * self.T_mean[self.iii]) * 8.36 * (K_TO_C - self.Tz) ** -2.061
        elif self.bdot_type == 'mean':
            dr_dt = (RHO_I - self.rho) * A_mean * (76.138 - 0.28965 * self.T_mean[self.iii]) * 8.36 * (K_TO_C - self.Tz) ** -2.061

        drho_dt = dr_dt / S_PER_YEAR    #To get into (kg m^3)/seconds
      
        viscosity = ((self.rho* self.sigma )/2)* (1/(drho_dt))  #In Pascal seconds

        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        return self.RD
    ### end Helsen_2008 ###
    #######################

    def Simonsen_2013(self):
        '''
        Accumulation units are kg/m^2/year

        '''
        ar1 = 0.07
        ar2 = 0.03
        Ec  = 60.0e3
        Eg  = 42.4e3
        # F0  = 0.68 #firnmice value?
        # F1  = 1.03 #firnmice value?
        F0=0.8 # Simonsen's recommended (email correspondence)
        F1=1.25 # Simonsen's recommended (email correspondence) (See email from 4/7/15)


        A_instant   = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM * 1000
        A_mean_1    = self.bdot_mean[self.rho < RHO_1] * RHO_I_MGM * 1000
        A_mean_2    = self.bdot_mean[self.rho >= RHO_1]* RHO_I_MGM * 1000
        dr_dt       = np.zeros(self.gridLen)
        viscosity   = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Simonsen physics")
            gamma = 61.7 / (A_instant ** (0.5)) * np.exp(-3800. / (R * self.T_mean[self.iii]))
            dr_dt[self.rho < RHO_1]  = F0 * (RHO_I - self.rho[self.rho < RHO_1]) * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean[self.iii]))
            dr_dt[self.rho >= RHO_1] = F1 * gamma * (RHO_I - self.rho[self.rho >= RHO_1]) * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean[self.iii]))
        elif self.bdot_type == 'mean':
            gamma = 61.7 / (A_mean_2 ** (0.5)) * np.exp(-3800.0 / (R * self.T_mean[self.iii]))
            dr_dt[self.rho < RHO_1]  = F0 * (RHO_I - self.rho[self.rho < RHO_1]) * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean[self.iii]))
            dr_dt[self.rho >= RHO_1] = F1 * gamma * (RHO_I - self.rho[self.rho >= RHO_1]) * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean[self.iii]))

        drho_dt = dr_dt / S_PER_YEAR
        
        viscosity[self.rho < RHO_1]   = ((self.rho[self.rho < RHO_1]*(self.sigma[self.rho < RHO_1]))/2)/drho_dt[self.rho < RHO_1] 
        viscosity[self.rho >= RHO_1]  = ((self.rho[self.rho >= RHO_1]* self.sigma[self.rho >= RHO_1] )/2)/drho_dt[self.rho >= RHO_1] 
        
        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        return self.RD
    ### end Simonsen_2013 ###
    #########################

    def Ligtenberg_2011(self):
        '''
        Units are mm W.E. per year
        b_dot is meant to be accumulation over a reference period (20 years for spin up, 1 year for regular?) (not mean over the lifetime  of a parcel)

        '''
        ar1 = 0.07
        ar2 = 0.03
        Ec  = 60.0e3
        Eg  = 42.4e3

        dr_dt     = np.zeros(self.gridLen)
        viscosity = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Ligtenberg 2011 physics")
            A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM * 1000
            M_0 = 1.435 - 0.151 * np.log(A_instant)
            M_1 = 2.366 - 0.293 * np.log(A_instant)
            M_0 = np.max((0.25,M_0))
            M_1 = np.max((0.25,M_1))
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean[self.iii]))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1])+ Eg / (R * self.T_mean[self.iii]))
        elif self.bdot_type == 'mean':
            A_mean_1 = self.bdot_mean[self.rho < RHO_1] * RHO_I
            A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * RHO_I
            M_0 = 1.435 - 0.151 * np.log(A_mean_1)
            M_1 = 2.366 - 0.293 * np.log(A_mean_2)
            M_0[M_0<0.25]=0.25
            M_1[M_1<0.25]=0.25
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean[self.iii]))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean[self.iii]))

        elif self.bdot_type == 'stress':
            A_mean_1 = self.bdot_mean[self.rho < RHO_1] * RHO_I
            sigma_1 = self.sigma[self.rho < RHO_1]
            A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * RHO_I
            sigma_2 = self.sigma[self.rho>= RHO_1]
            M_0 = 1.435 - 0.151 * np.log(A_mean_1)
            M_1 = 2.366 - 0.293 * np.log(A_mean_2)
            M_0[M_0<0.25]=0.25
            M_1[M_1<0.25]=0.25
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean[self.iii]))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean[self.iii]))


        drho_dt = dr_dt / S_PER_YEAR
        
        viscosity[self.rho < RHO_1]   = (self.rho[self.rho < RHO_1]* self.sigma[self.rho < RHO_1])/ (2) /drho_dt[self.rho < RHO_1] 
        viscosity[self.rho >= RHO_1]  = (self.rho[self.rho >= RHO_1]* self.sigma[self.rho >= RHO_1] )/ (2)/drho_dt[self.rho >= RHO_1] 

        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        return self.RD
    ### end Ligtenberg_2011 ###
    ###########################

    def Barnola_1991(self):
        '''

        uses m W.E. (zone 1) and stress (zone 2)

        '''
        Q1              = 10160.0
        k1              = 11.0
        aHL             = 1.0
        alphaBarnola    = -37.455
        betaBarnola     = 99.743
        deltaBarnola    = -95.027
        gammaBarnola    = 30.673
        A0b             = 2.54e4
        n               = 3.0
        QBarnola        = 60.0e3
        closeOff        = 800.0

        self.rho[self.rho > RHO_I] = RHO_I # The Barnola model will go a fraction over the ice density (oself.RDer 10^-3), so this stops that.
        drho_dt     = np.zeros(self.gridLen)
        viscosity   = np.zeros(self.gridLen)
        D           = self.rho / RHO_I
        nBa         = n * np.ones(self.gridLen)
        A0          = A0b * np.ones(self.gridLen) / 1.e18 #this is for the n=3 region.

        sigmaEff = self.sigma
        # nBa[sigmaEff<1.0e5]=1.0


        ### Zone 1 ###
        A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean_1 = self.bdot_mean[self.rho < RHO_1] * RHO_I_MGM

        if self.bdot_type == 'instant':
            drho_dt[self.rho < RHO_1] = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * A_instant ** aHL * 1000 / S_PER_YEAR
        elif self.bdot_type == 'mean':
            drho_dt[self.rho < RHO_1] = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * A_mean_1 ** aHL * 1000 / S_PER_YEAR

        ### Zone 2 ###
        condition_zn2 = ((self.rho >= RHO_1) & (self.rho <= closeOff))
        fe = 10.0 ** (alphaBarnola * (self.rho[condition_zn2] / 1000) ** 3. + betaBarnola * (self.rho[condition_zn2] / 1000) ** 2. + deltaBarnola * self.rho[condition_zn2] / 1000 + gammaBarnola)
        drho_dt[condition_zn2] = self.rho[condition_zn2] * A0[condition_zn2] * np.exp(-QBarnola / (R * self.Tz[condition_zn2])) * fe * (sigmaEff[condition_zn2] ** nBa[condition_zn2])

        # zone 3
        condition_zn3 = (self.rho > closeOff)
        fs = (3. / 16.) * (1 - self.rho[condition_zn3] / RHO_I) / (1 - (1 - self.rho[condition_zn3] / RHO_I) ** (1. / 3.)) ** 3.
        drho_dt[condition_zn3] = self.rho[condition_zn3] * A0[condition_zn3] * np.exp(-QBarnola / (R * self.Tz[condition_zn3])) * fs * (sigmaEff[condition_zn3] ** nBa[condition_zn3])
        
        viscosity[self.rho < RHO_1]   = (self.rho[self.rho < RHO_1]* self.sigma[self.rho < RHO_1])/ (2)/drho_dt[self.rho < RHO_1] 
        viscosity[self.rho >= RHO_1]  = (self.rho[self.rho >= RHO_1]* self.sigma[self.rho >= RHO_1] )/ (2)/drho_dt[self.rho >= RHO_1] 
    
        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        return self.RD
    ### end Barnola_1991 ###
    ########################
    
    def Morris_HL_2014(self):
        '''

        Uses stress instead of accumulation.

        Need to choose physics for zone 2. Herron and Langway here.
        '''

        QMorris = self.QMorris
        # QMorris = 60.e3

        if self.iii ==0:
            print('CAUTION: MORRIS PHYSICS ARE STILL UNDER CODE DEVELOPMENT!')
            print('see physics.py for more details')          
            print('QMorris is %s' %QMorris)
            print('If you want to change this, you must in physics.py and spin.py')

        ### Figure out what the densification constant, k, should be.
        ### I calculated coefficients for slope and intercept by regressing the 
        ### (E_alpha - E_H) for each site and each activation energy provided by
        ### Liz Morris (the fix of the erroneous Table 2 in original paper)
        slope = -0.0009667915546575245 * self.QMorris /1.e3 + 0.001947860800510695
        intercept = 0.29455063899108685* self.QMorris /1.e3 - 2.652540535829697
        deltaE = (slope*self.T_mean[self.iii] + intercept) * 1000
        kHL = 11.0 # units m water eq.
        Estar = 10.16e3
        kMorris = kHL * np.exp(-1 * (Estar - deltaE) / (R * self.T_mean[self.iii]))
        if self.iii==0:
            print('kMorris=',kMorris)
        ##

        ### calculations for the (1-M_0*m) term
        ind_mor = np.where(self.rho<RHO_1)[0] #indices in zone 1
        rho_mor=self.rho[ind_mor]
        dep_mor=self.z[ind_mor]
        coefs = poly.polyfit(dep_mor, rho_mor, 2)
        rho_fit = poly.polyval(dep_mor, coefs) #rho_0 in Morris
        M0bar = 3.3
        m = (self.rho[ind_mor]-rho_fit)/(RHO_I-rho_fit)
        ###

        drho_dt   = np.zeros(self.gridLen)
        viscosity = np.zeros(self.gridLen)
 
        drho_dt[self.rho < RHO_1] = (kMorris / (RHO_W_KGM * GRAVITY)) * ((RHO_I - self.rho[self.rho < RHO_1])) * (1 / self.Hx[self.rho < RHO_1]) * np.exp(-1*self.QMorris / (R * self.Tz[self.rho < RHO_1])) * self.sigma[self.rho < RHO_1] * (1 - M0bar*m)

        # Use HL Dynamic for zone 2 b/c Morris does not specify zone 2.
        Q2  = 21400.0
        k2  = 575.0
        bHL = 0.5
        A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * RHO_I_MGM
        if self.bdot_type == 'instant':
            drho_dt[self.rho >= RHO_1]   = k2 * np.exp(-Q2 / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I_MGM - self.rho[self.rho >= RHO_1] / 1000) * A_instant ** bHL * 1000 / S_PER_YEAR
        elif self.bdot_type == 'mean':
            drho_dt[self.rho >= RHO_1]   = k2 * np.exp(-Q2 / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I_MGM - self.rho[self.rho >= RHO_1] / 1000) * A_mean_2 ** bHL * 1000 / S_PER_YEAR
        
        viscosity[self.rho < RHO_1]   = -(self.rho[self.rho < RHO_1]* self.sigma[self.rho < RHO_1])/ (2)/drho_dt[self.rho < RHO_1] 
        viscosity[self.rho >= RHO_1]  = (self.rho[self.rho >= RHO_1]* self.sigma[self.rho >= RHO_1]) / (2)/drho_dt[self.rho >= RHO_1] 

        self.Hx = self.Hx + np.exp(-1*self.QMorris / (R * self.Tz)) * self.dt
        Hx_new  = np.exp(-1 * self.QMorris / (R * self.Tz[0])) * self.dt 
        self.Hx = np.concatenate(([Hx_new],self.Hx[:-1]))

        # input('enter to continue')
        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        self.RD['Hx'] = self.Hx
        return self.RD
    ### end Morris_HL_2014 ###
    ##########################

    def KuipersMunneke_2015(self):
        '''

        Units are mm W.E. per year
        b_dot is meant to be accumulation over a reference period (20 years for spin up, 1 year for regular?) (not mean over the lifetime  of a parcel)

        '''
        ar1 = 0.07
        ar2 = 0.03
        Ec  = 60.0e3
        Eg  = 42.4e3

        dr_dt     = np.zeros(self.gridLen)
        viscosity = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM * 1000
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Ligtenberg 2011 physics")
            M_0 = 1.042 - 0.0916 * np.log(A_instant)
            M_1 = 1.734 - 0.2039 * np.log(A_instant)
            M_0 = np.max((0.25,M_0))
            M_1 = np.max((0.25,M_1))
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean[self.iii]))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1])+ Eg / (R * self.T_mean[self.iii]))

        elif self.bdot_type == 'mean':
            A_mean_1 = self.bdot_mean[self.rho < RHO_1] * RHO_I
            A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * RHO_I

            M_0 = 1.042 - 0.0916 * np.log(A_mean_1)
            M_1 = 1.734 - 0.2039 * np.log(A_mean_2)

            M_0[M_0<0.25]=0.25
            M_1[M_1<0.25]=0.25

            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean[self.iii]))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean[self.iii]))

            # dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T50))
            # dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T50))

        elif self.bdot_type == 'stress':

            A_mean_1 = self.mass[self.rho < RHO_1]*10
            A_mean_2 = self.mass[self.rho >= RHO_1]*10

            M_0 = 1.042 - 0.0916 * np.log(A_mean_1)
            M_1 = 1.734 - 0.2039 * np.log(A_mean_2)

            M_0[M_0<0.25] = 0.25
            M_1[M_1<0.25] = 0.25

            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean[self.iii]))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean[self.iii]))

        drho_dt = dr_dt / S_PER_YEAR
        drho_dt[self.rho>=RHO_I] = 0
        
        viscosity[self.rho < RHO_1]   = (self.rho[self.rho < RHO_1]* self.sigma[self.rho < RHO_1])/ (2 )/drho_dt[self.rho < RHO_1] 
        viscosity[self.rho >= RHO_1]  = (self.rho[self.rho >= RHO_1] * self.sigma[self.rho >= RHO_1] )/ (2 )/drho_dt[self.rho >= RHO_1] 
        
        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        return self.RD
    ### end KuipersMunneke_2015 ###
    ###############################

    def Goujon_2003(self):
        '''
        Uses stress
        '''

        atmosP      = 101325.0 # Atmospheric Pressure
        dDdt        = np.zeros(self.gridLen) # Capital D is change in relative density
        # top2m       = np.nonzero(self.z <= 2.)
        
        # self.rho[top2m] = self.rhos0 # top 2 meters of Goujon model are always set to surface density
        
        sigma_MPa   = self.sigma / (1.0e6)
        sigma_bar   = self.sigma / (1.0e5)
        Qgj         = 60.0e3
        n           = 3.0 

        rhoi2cgs    = .9165 * (1.-1.53e-4 * (self.T_mean[self.iii] - 273.15)) # Density of ice, temperature dependent, g/cm^3
        rhoi2       = rhoi2cgs * 1000.0 # density of ice, kg/m^3
        
        D           = self.rho / rhoi2 # Relative density
        Dm23        = 0.9 #transition from zone 2 to 3
        rho23       = Dm23 * rhoi2 #density of 2/3 transition
        
        D0          = 0.00226 * self.T_mean[self.iii] + 0.03 #D0 is the relative density of transition from zone 1 to 2. Here is from Arnaud et al. (2000) eq. 8, not Goujon (Goujon uses C, not K)

        if D0 > 0.59: #Model requires zone 1/2 transition to be less than D=0.6
            D0      = 0.59
        Dms         = D0 + 0.009 #D0 + epsilon factor, maximum is 0.599
        Dmsrho      = Dms * rhoi2 # density of zone 1/2 transition
        
        ind1        = np.argmax(D >= Dms) #first index where the density is greater than or equal to the transition
        Dm          = D[ind1] #actual transition relative density. Use this so that the transition falls at a node
        Dmrho       = Dm * rhoi2 #density of first node, zone 2

        A           = 7.89e3 * np.exp(-Qgj/(R * self.Tz)) * 1.0e-3 # A given in MPa^-3 s^-1, Goujon uses bar as pressure unit. Eq. A5 in Goujon
        ccc         = 15.5 # no units given, equation A7, given as c
        Z0g         = 110.2 * D0 ** 3.-148.594 * D0 ** 2.+87.6166 * D0-17. # from Anais' code           
        lp          = (D/D0) ** (1.0/3.0) # A6
        Zg          = Z0g+ccc * (lp-1.0) # A7
        lpp_n       = (4.0 * Z0g * (lp-1.0) ** 2.0 * (2.0 * lp+1.0) + ccc * (lp-1.0) ** 3.0 * (3.0 * lp + 1.0)) # A8, numerator
        lpp_d       = (12.0 * lp * (4.0 * lp - 2.0 * Z0g * (lp-1.0) - ccc * (lp-1.0) ** 2.0)) # A8, denominator
        lpp         = lp + (lpp_n/lpp_d) # A8 (l double prime)
        a           = (np.pi/(3.0 * Zg * lp ** 2.0)) * (3.0 * (lpp ** 2.0 - 1.0) * Z0g + lpp ** 2.0 * ccc * (2.0 * lpp-3.0)+ccc) # A9
        sigmastar   = (4.0 * np.pi * sigma_bar)/(a * Zg * D) # A4

        # gamma_An    = (5.3*A[ind1] * (Dms**2*D0)**(1.0/3.0) * (a[ind1]/np.pi)**(1.0/2.0) * (sigmastar[ind1]/3.0)**n) / ((sigma_bar[ind1]/(Dms**2))*(1-(5.0/3.0*Dms))) #this is the analytic solution of what gamma should be by combining equations A1 and A3 and solving for gamma (densification rate should be smooth at the zone 1/2 transition). Does not get used.
        # gamma_An    = (5.3*A[ind1+1] * (Dms**2*D0)**(1.0/3.0) * (a[ind1+1]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1]/3.0)**n) / ((sigma_bar[ind1+1]/(Dms**2))*(1-(5.0/3.0*Dms))) #this is the analytic solution of what gamma should be by combining equations A1 and A3 and solving for gamma (densification rate should be smooth at the zone 1/2 transition). Does not get used.

        if self.iii == 0 or ind1 != self.ind1_old:
            self.Gamma_Gou       = 0.5 / S_PER_YEAR
            self.Gamma_old_Gou   = self.Gamma_Gou
        else:
            self.Gamma_Gou       = self.Gamma_old_Gou

        dDdt[0:ind1+1]  = self.Gamma_Gou*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
        dDdt[ind1+1:]   = 5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n         
        gfrac           = 0.03
        gam_div         = 1 + gfrac #change this if want: making it larger will make the code run faster. Must be >=1.
        
        ########## iterate to increase gamma first if not in steady state    
        if self.iii != 0 and dDdt[ind1] <= dDdt[ind1+1] and self.Gamma_Gou!=self.Gamma_old2_Gou:
            cc = 1
            while dDdt[ind1] < dDdt[ind1 + 1]:
                self.Gamma_Gou       = self.Gamma_Gou * (gam_div)
                dDdt[0:ind1+1]  = self.Gamma_Gou*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
                dDdt[ind1+1:]   = 5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n

                cc += 1
                if cc>10000:
                    print('Goujon is not converging. exiting')
                    sys.exit()

        ### then iterate to find the maximum value of gamma that will make a continuous drho/dt

        counter = 1
        while dDdt[ind1] >= dDdt[ind1 + 1]:
            
            self.Gamma_Gou      = self.Gamma_Gou / (1 + gfrac/2.0)
            dDdt[0:ind1+1] = self.Gamma_Gou*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
            dDdt[ind1+1:]  = 5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n
            counter += 1

            if counter>10000:
                print('Goujon is not converging. exiting')
                sys.exit()

        #####
        # dDdt[0:ind1+1] = gamma_An*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
        # dDdt[ind1+1:]  = 5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n
        #####

        # if self.iii<10:
            # print('dDdt',dDdt[ind1:ind1+2])
        self.Gamma_old2_Gou  = self.Gamma_old_Gou
        self.Gamma_old_Gou   = self.Gamma_Gou
        self.ind1_old        = ind1
        #####################
        
        rhoC        = RHO_2 #should be Martinerie density
        frho2       = interpolate.interp1d(self.rho,sigma_bar,bounds_error=False,fill_value='extrapolate')
        sigmarho2   = frho2(rhoC) #pressure at close off

        ind2 = np.argmax(D >= Dm23)
        
        # sigma_b = sigmarho2 * (D*(1-rhoC/rhoi)) / (rhoC/rhoi*(1-D))
        # sigma_b = (sigma_MPa[ind2] * (D*(1-Dm23)) / (Dm23*(1-D)))/10. #works for Ex2
        # sigma_b = ((sigma_bar + atmosP/1.0e5) * (D*(1-Dm23)) / (Dm23*(1-D)))
        sigma_b                 = ((atmosP / 1.0e5) * (D * (1 - Dm23)) / (Dm23 * (1 - D)))
        sigmaEff                = (sigma_bar + atmosP / 1.0e5 - sigma_b)
        sigmaEff[sigmaEff <= 0] = 1.0e-9
        
        ind2 = np.argmax(D >= Dm23)
        ind3 = ind2 + 10

        dDdt[D>Dm23]    = 2. * A[D>Dm23] * ( (D[D>Dm23] * (1 - D[D>Dm23])) / (1 - (1 - D[D>Dm23])**(1/n))**n ) * (2 * sigmaEff[D>Dm23]/n)**3.0
        Ad              = 1.2e3 * np.exp(-Qgj / (R * self.Tz)) * 1.0e-1
        T34             = 0.98
        dDdt[D>T34]     = 9/4 * Ad[D>T34] * (1-D[D>T34]) * sigmaEff[D>T34]

        dDdt_old        = dDdt
        drho_dt         = dDdt*rhoi2
        # drho_dt[top2m]  = 0.0
        
        viscosity = np.zeros(self.gridLen)
        
        self.RD['drho_dt'] = drho_dt
        # global Gamma_Gou, Gamma_old_Gou, Gamma_old2_Gou, ind1_old
        self.RD['Gamma_Gou'] = self.Gamma_Gou
        self.RD['Gamma_old_Gou'] = self.Gamma_old_Gou
        self.RD['Gamma_old2_Gou'] = self.Gamma_old2_Gou
        self.RD['ind1_old'] = self.ind1_old
        return self.RD
    ### end Goujon_2003 ###
    #######################

    def Breant2017(self):
        '''
        Bréant, C., Martinerie, P., Orsi, A., Arnaud, L., and Landais, A.: Modelling firn thickness evolution during the last deglaciation: constraints on sensitivity to temperature and impurities, Clim. Past, 13, 833–853, https://doi.org/10.5194/cp-13-833-2017, 2017.

        '''

        atmosP      = 101325.0 # Atmospheric Pressure
        dDdt        = np.zeros(self.gridLen) # Capital D is change in relative density
        # top2m       = np.nonzero(self.z <= 2.)       
        # self.rho[top2m] = self.rhos0 # top 2 meters of Goujon model are always set to surface density
        
        sigma_MPa   = self.sigma / (1.0e6)
        sigma_bar   = self.sigma / (1.0e5)
        Qgj         = 60.0e3
        n           = 3.0 

        rhoi2cgs    = .9165 * (1.-1.53e-4 * (self.T_mean[self.iii] - 273.15)) # Density of ice, temperature dependent, g/cm^3
        rhoi2       = rhoi2cgs * 1000.0 # density of ice, kg/m^3
        
        D           = self.rho / rhoi2 # Relative density
        Dm23        = 0.9 #transition from zone 2 to 3
        rho23       = Dm23 * rhoi2 #density of 2/3 transition
        
        D0          = 0.56        

        Dmsrho      = D0 * rhoi2 # density of zone 1/2 transition
        
        ind1        = np.argmax(D >= D0) #first index where the density is greater than or equal to the transition
        Dm          = D[ind1] #actual transition relative density. Use this so that the transition falls at a node
        Dmrho       = Dm * rhoi2 #density of first node, zone 2

        # A           = 7.89e3 * np.exp(-Qgj/(R * self.Tz)) * 1.0e-3 # A given in MPa^-3 s^-1, Goujon uses bar as pressure unit. Eq. A5 in Goujon
        A0 = 7.89e3 * 1.0e-3 #[bar^-3 s^-1]
        a1 = 1.05912848e-02
        a2 = 1400
        a3 = 6.0e-15
        Q1 = 110000
        Q2 = 75000
        Q3 = 1500
        A = A0 * (a1*np.exp(-Q1/(R*self.Tz)) + a2*np.exp(-Q2/(R*self.Tz)) + a1*np.exp(-Q3/(R*self.Tz)))

        ccc         = 15.5 # no units given, equation A7, given as c
        Z0g         = 110.2 * D0 ** 3.-148.594 * D0 ** 2.+87.6166 * D0-17. # from Anais' code           
        lp          = (D/D0) ** (1.0/3.0) # A6
        Zg          = Z0g + ccc * (lp-1.0) # A7
        lpp_n       = (4.0 * Z0g * (lp-1.0) ** 2.0 * (2.0 * lp+1.0) + ccc * (lp-1.0) ** 3.0 * (3.0 * lp + 1.0)) # A8, numerator
        lpp_d       = (12.0 * lp * (4.0 * lp - 2.0 * Z0g * (lp-1.0) - ccc * (lp-1.0) ** 2.0)) # A8, denominator
        lpp         = lp + (lpp_n/lpp_d) # A8 (l double prime)
        a           = (np.pi/(3.0 * Zg * lp ** 2.0)) * (3.0 * (lpp ** 2.0 - 1.0) * Z0g + lpp ** 2.0 * ccc * (2.0 * lpp-3.0)+ccc) # A9
        sigmastar   = (4.0 * np.pi * sigma_bar)/(a * Zg * D) # A4

        # gamma_An    = (5.3*A[ind1] * (Dms**2*D0)**(1.0/3.0) * (a[ind1]/np.pi)**(1.0/2.0) * (sigmastar[ind1]/3.0)**n) / ((sigma_bar[ind1]/(Dms**2))*(1-(5.0/3.0*Dms))) #this is the analytic solution of what gamma should be by combining equations A1 and A3 and solving for gamma (densification rate should be smooth at the zone 1/2 transition). Does not get used.
        # gamma_An    = (5.3*A[ind1+1] * (Dms**2*D0)**(1.0/3.0) * (a[ind1+1]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1]/3.0)**n) / ((sigma_bar[ind1+1]/(Dms**2))*(1-(5.0/3.0*Dms))) #this is the analytic solution of what gamma should be by combining equations A1 and A3 and solving for gamma (densification rate should be smooth at the zone 1/2 transition). Does not get used.

        if self.iii == 0 or ind1 != self.ind1_old:
            self.Gamma_Gou       = 0.5 / S_PER_YEAR
            self.Gamma_old_Gou   = self.Gamma_Gou
        else:
            self.Gamma_Gou       = self.Gamma_old_Gou

        dDdt[0:ind1+1]  = self.Gamma_Gou*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
        dDdt[ind1+1:]   = 5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n         
        gfrac           = 0.03
        gam_div         = 1 + gfrac #change this if want: making it larger will make the code run faster. Must be >=1.
        
        ########## iterate to increase gamma first if not in steady state    
        if self.iii != 0 and dDdt[ind1] <= dDdt[ind1+1] and self.Gamma_Gou!=self.Gamma_old2_Gou:
            cc = 1
            while dDdt[ind1] < dDdt[ind1 + 1]:
                self.Gamma_Gou       = self.Gamma_Gou * (gam_div)
                dDdt[0:ind1+1]  = self.Gamma_Gou*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
                dDdt[ind1+1:]   = 5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n

                cc += 1
                if cc>10000:
                    print('Goujon is not converging. exiting')
                    sys.exit()

        ### then iterate to find the maximum value of gamma that will make a continuous drho/dt

        counter = 1
        while dDdt[ind1] >= dDdt[ind1 + 1]:
            
            self.Gamma_Gou      = self.Gamma_Gou / (1 + gfrac/2.0)
            dDdt[0:ind1+1] = self.Gamma_Gou*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
            dDdt[ind1+1:]  = 5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n
            counter += 1

            if counter>10000:
                print('Goujon is not converging. exiting')
                sys.exit()

        #####
        # dDdt[0:ind1+1] = gamma_An*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
        # dDdt[ind1+1:]  = 5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n
        #####

        # if self.iii<10:
            # print('dDdt',dDdt[ind1:ind1+2])
        self.Gamma_old2_Gou  = self.Gamma_old_Gou
        self.Gamma_old_Gou   = self.Gamma_Gou
        self.ind1_old        = ind1
        #####################
        
        rhoC        = RHO_2 #should be Martinerie density
        frho2       = interpolate.interp1d(self.rho,sigma_bar,bounds_error=False,fill_value='extrapolate')
        sigmarho2   = frho2(rhoC) #pressure at close off

        ind2 = np.argmax(D >= Dm23)
        
        # sigma_b = sigmarho2 * (D*(1-rhoC/rhoi)) / (rhoC/rhoi*(1-D))
        # sigma_b = (sigma_MPa[ind2] * (D*(1-Dm23)) / (Dm23*(1-D)))/10. #works for Ex2
        # sigma_b = ((sigma_bar + atmosP/1.0e5) * (D*(1-Dm23)) / (Dm23*(1-D)))
        sigma_b                 = ((atmosP / 1.0e5) * (D * (1 - Dm23)) / (Dm23 * (1 - D)))
        sigmaEff                = (sigma_bar + atmosP / 1.0e5 - sigma_b)
        sigmaEff[sigmaEff <= 0] = 1.0e-9
        
        ind2 = np.argmax(D >= Dm23)
        ind3 = ind2 + 10

        dDdt[D>Dm23]    = 2. * A[D>Dm23] * ( (D[D>Dm23] * (1 - D[D>Dm23])) / (1 - (1 - D[D>Dm23])**(1/n))**n ) * (2 * sigmaEff[D>Dm23]/n)**3.0
        Ad              = 1.2e3 * np.exp(-Qgj / (R * self.Tz)) * 1.0e-1
        T34             = 0.98
        dDdt[D>T34]     = 9/4 * Ad[D>T34] * (1-D[D>T34]) * sigmaEff[D>T34]

        dDdt_old        = dDdt
        drho_dt         = dDdt*rhoi2
        # drho_dt[top2m]  = 0.0
        
        viscosity = np.zeros(self.gridLen)
        
        self.RD['drho_dt'] = drho_dt
        # global Gamma_Gou, Gamma_old_Gou, Gamma_old2_Gou, ind1_old
        self.RD['Gamma_Gou'] = self.Gamma_Gou
        self.RD['Gamma_old_Gou'] = self.Gamma_old_Gou
        self.RD['Gamma_old2_Gou'] = self.Gamma_old2_Gou
        self.RD['ind1_old'] = self.ind1_old
        return self.RD


    def Crocus(self):
        '''
        Uses stress

        '''

        f1 = 1.0 # unitless
        f1 = (1+60*self.LWC/self.dz)**(-1) #From Vincent
        f2 = 4.0 # unitless
        nu_0 = 7.62237e6 # kg s^-1
        a_n = 0.1 # K^-1
        b_n = 0.023 # m^3 kg^-1
        # c_n = 250 # kg m^-3 # original
        c_n = 358. # kg m^-3 #VV eq (8) van Kampenhout et al. (2017)

        viscosity = f1 * f2 * nu_0 * self.rho / c_n * np.exp(a_n * (273.15 - self.Tz) + b_n * self.rho)

        dr_dt = self.rho * self.sigma / viscosity

        dr_dt[self.rho>=RHO_I] = 0

        drho_dt = dr_dt #/ S_PER_YEAR

        ##########
        CHybrid = True
        if CHybrid:
            if self.iii==0:
                print("You are using CROCUS hybrid")
            RHO_CC = 450
            ar1 = 0.07
            ar2 = 0.03
            Ec  = 60.0e3
            Eg  = 42.4e3
            A_mean_1 = self.bdot_mean[self.rho < RHO_CC] * RHO_I
            A_mean_2 = self.bdot_mean[self.rho >= RHO_CC] * RHO_I

            M_0 = 1.042 - 0.0916 * np.log(A_mean_1)
            M_1 = 1.734 - 0.2039 * np.log(A_mean_2)

            M_0[M_0<0.25]=0.25
            M_1[M_1<0.25]=0.25

            # dr_dt[self.rho < RHO_CC]  = (RHO_I - self.rho[self.rho < RHO_CC]) * M_0 * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_CC]) + Eg / (R * self.T_mean[self.iii]))
            dr_dt[self.rho >= RHO_CC] = ((RHO_I - self.rho[self.rho >= RHO_CC]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_CC]) + Eg / (R * self.T_mean[self.iii])))/ S_PER_YEAR
            drho_dt[self.rho>=RHO_I] = 0

        # dr_dt[self.rho >= RHO_CC] = ((RHO_I - self.rho[self.rho >= RHO_CC]) * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_CC]) + Eg / (R * self.T_mean[self.iii])))/ S_PER_YEAR


        #########
        
        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        
        return self.RD
    ### end Crocus ###
    ##################

    def GSFC2020(self):
        '''
        Taken from Medley et al. (2020)
        Accumulation units are kg/m^2/year
        Subscript here are changed to 1 and 2 (correspond to densification stage)
        (in Medley et al. 2020 they are 0 and 1)
        '''

        ar1 = 0.07
        ar2 = 0.03
        Eg  = 42.4e3

        # old ones
        # Ec1  = 60.0e3
        # Ec2 = 56973.0
        # alpha1 = 0.9250
        # alpha2 = 0.6354

        # New:
        Ec1  = 59442.0
        Ec2 = 56827.0
        alpha1 = 0.9021
        alpha2 = 0.6372

        # beta1 = -0.1483 # Changed from Smith 2020, which used subcripts 0 and 1
        # beta2 = -0.3510 # we want subscripts to reflect zone so we use 1 and 2
        # E1 = -731
        # E2 = -2861

        A_mean_1    = self.bdot_mean[self.rho < RHO_1] * RHO_I_MGM * 1000 # Convert to kg/m^2/time
        A_mean_2    = self.bdot_mean[self.rho >= RHO_1] * RHO_I_MGM * 1000
        dr_dt       = np.zeros(self.gridLen)
        viscosity   = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            print('instant bdot_type does not work with GSFC physics. exiting')
            sys.exit()
        elif self.bdot_type == 'mean':
            # R1 = A_mean_1**beta1 * np.exp(-E1 / (R * self.Tz[self.rho < RHO_1]))
            # R2 = A_mean_2**beta2 * np.exp(-E2 / (R * self.Tz[self.rho >= RHO_1]))
            indS1 = np.where(self.rho < RHO_1)[0]
            indS2 = np.where(self.rho >= RHO_1)[0]

            c1 = ar1 * A_mean_1**alpha1 * GRAVITY * np.exp(-Ec1 / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean[self.iii]))
            c2 = ar2 * A_mean_2**alpha2 * GRAVITY * np.exp(-Ec2 / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean[self.iii]))

            dr_dt[self.rho < RHO_1]  = c1 * (RHO_I - self.rho[self.rho < RHO_1])
            dr_dt[self.rho >= RHO_1] = c2 * (RHO_I - self.rho[self.rho >= RHO_1])

        drho_dt = dr_dt / S_PER_YEAR
        
        viscosity[self.rho < RHO_1]   = (self.rho[self.rho < RHO_1]/ (2 * self.sigma[self.rho < RHO_1]))/drho_dt[self.rho < RHO_1] 
        viscosity[self.rho >= RHO_1]  = (self.rho[self.rho >= RHO_1] / (2 * self.sigma[self.rho >= RHO_1] ))/drho_dt[self.rho >= RHO_1]  
          
        self.RD['drho_dt']   = drho_dt
        self.RD['viscosity'] = viscosity
        return self.RD
    ### end GSFC2020 ###
    #########################



    def Max2018b(self):
        '''
        Experimental: based on firn compaction data from FirnCover campaigns
        '''

        # k0 = 0.15 # units Pa^-1 s^-1
        # k1 = 0.05 # units Pa^-1 s^-1
        k0 = 6.0e7
        Q = 60000.0
        dr_dt = np.zeros(self.gridLen)
        Q2  = 21400.0        
        k2  = 575.0
        aHL = 1.0
        bHL = 0.5

        A_mean = self.bdot_mean * RHO_I_MGM
        
        # dr_dt[self.rho < RHO_1] = k0 * np.exp(-1*Q / (R * self.Tz[self.rho < RHO_1])) * (RHO_I - self.rho[self.rho < RHO_1]) * self.sigma[self.rho < RHO_1]
        # dr_dt[self.rho >= RHO_1] = k1 * np.exp(-1*Q / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I - self.rho[self.rho >= RHO_1]) * self.sigma[self.rho >= RHO_1]

        msk = ((self.age>0) & (self.rho < RHO_1))
        msk2 = self.rho>=RHO_1
        dr_dt[msk] = k0 * np.exp(-1*Q / (R * self.Tz[msk])) * (RHO_I - self.rho[msk]) * self.sigma[msk] / self.age[msk]
        # dr_dt[msk] = k0 * np.exp(-1*Q / (R * self.Tz[msk])) * (RHO_I - self.rho[msk]) * self.sigma[msk] / self.age[msk]

        dr_dt[self.rho >= RHO_1]    = k2 * np.exp(-Q2 / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I_MGM - self.rho[self.rho >= RHO_1] / 1000) * (A_mean[self.rho >= RHO_1])**bHL * 1000 / S_PER_YEAR

        self.RD['drho_dt'] = dr_dt # units are (kg m^-3) s^-1
        return self.RD
    ### end Max2018b ###
    ####################

    def Max2018(self):
        '''
        Experimental: based on firn compaction data from FirnCover campaigns
        '''

        # k0 = 0.15 # units Pa^-1 s^-1
        # k1 = 0.05 # units Pa^-1 s^-1
        # k0 = 8.5e9
        if self.iii == 0:
            print('Caution: Max2018 physics are still in development.')
        Q = 70000.0
        dr_dt = np.zeros(self.gridLen)
        Q2  = 21400.0        
        k2  = 575.0
        aHL = 1.0
        bHL = 0.5

        A_mean = self.bdot_mean * RHO_I_MGM

        k0 = -1.387e10 * np.nanmean(self.bdot_mean) + 1.042e10
        # print(np.nanmean(self.bdot_mean))
        
        # dr_dt[self.rho < RHO_1] = k0 * np.exp(-1*Q / (R * self.Tz[self.rho < RHO_1])) * (RHO_I - self.rho[self.rho < RHO_1]) * self.sigma[self.rho < RHO_1]
        # dr_dt[self.rho >= RHO_1] = k1 * np.exp(-1*Q / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I - self.rho[self.rho >= RHO_1]) * self.sigma[self.rho >= RHO_1]

        msk = ((self.age>0) & (self.rho < RHO_1))
        msk2 = self.rho>=RHO_1
        dr_dt[msk] = k0 * np.exp(-1*Q / (R * self.Tz[msk])) * (RHO_I - self.rho[msk]) * self.sigma[msk] / self.age[msk]

        # Helsen
        # dr_dt[self.rho >= RHO_1] = ((RHO_I - self.rho[self.rho >= RHO_1]) * A_mean[self.rho >= RHO_1] * (76.138 - 0.28965 * self.T_mean) * 8.36 * (K_TO_C - self.Tz[self.rho >= RHO_1]) ** -2.061)/S_PER_YEAR

        # Li and Zwally
        # TmC   = self.T_mean - K_TO_C
        # beta1a = -9.788 + 8.996 * np.mean(A_mean) - 0.6165 * TmC
        # beta2a = beta1a / (-2.0178 + 8.4043 * np.mean(A_mean) - 0.0932 * TmC)
        # beta1 = np.ones(len(A_mean))*beta1a
        # beta2 = np.ones(len(A_mean))*beta2a
        # dr_dt[self.rho >= RHO_1]  = ((RHO_I - self.rho[self.rho >= RHO_1]) * A_mean[self.rho >= RHO_1] * beta2[self.rho >= RHO_1] * 8.36 * (K_TO_C - self.Tz[self.rho >= RHO_1]) ** -2.061)/S_PER_YEAR

        #use KM physics below 550
        ar2 = 0.03
        Ec = 60.0e3
        Eg = 42.4e3
        A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * RHO_I
        M_1 = 1.734 - 0.2039 * np.log(A_mean_2)
        M_1[M_1<0.25]=0.25

        dr_dt[self.rho >= RHO_1] = ((RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean[self.iii])))/S_PER_YEAR


        self.RD['drho_dt'] = dr_dt # units are (kg m^-3) s^-1
        return self.RD
    ### end Max2018 ###
    ###################

    def MaxSP(self):

        Q = 60.0e3
        # kz1 = 19.37
        # kz2 = 50.01
        # L = kz2-kz1
        # rho_mid = 550
        # kL = 0.1
        # kM = (L / (1 + np.exp(-kL*(self.rho - rho_mid)))) + kz1

        # def k_log(L,kL,x0,xa,d,os):
        #     denom = 1 + np.exp(-kL*(xa - x0))
        #     return L/denom + d*xa + ic 
        # kz1=np.median(kmed) #0.087 19.37
        # kz2 = visc_df.iloc[-4:-1]['k'].median() #0.16 50.01
        # L=kz2-kz1
        # L=25
        # rho_mid = 550
        # # xa = np.arange(300,801)
        # d = 0.075
        # ic = - 20
        # kL = 0.1
        # kM = k_log(L,kL,rho_mid,self.rho,d,ic)
        # kM[kM<8] = 8

        ###
        # # kV4 (these work pretty well)
        # def k_log_cf(xa,L,kL,os,rhoc):
        #     denom = 1 + np.exp(-kL*(xa - rhoc))
        #     return L/denom + os
        # L  = 44.0 #4.5e+01 
        # kL = 3.3284e-2 #3.12e-02
        # os = 9.65
        # rhoc = 542.8
        # kM = k_log_cf(self.rho,L,kL,os,rhoc)
        # [4.50000000e+01 3.27472804e-02 1.00000000e+01 5.45780400e+02]
        ####



        ## kV6
        def k_log_cf_d(xa,L,kL,os,d,rhoc):
            denom = 1 + np.exp(-kL*(xa - rhoc))
            return L/denom + d*xa + os
        # L  = 3.54937382e+01
        # kL = 4.18361570e-02
        # os = 9.22241253
        # d  = 2.19293203e-29
        # rhoc = 5.22153672e+02

        # L,kL,os,d,rhoc = 2.86331016e+11, 3.32122415e-02, 7.10378487e+10, 9.00002301e+01,5.26031683e+02 # good with 60 kJ/mol and Arthern dual activtion energy. Derived with age max.
        
        # L,kL,os,d,rhoc = 2.9e11, 3.4e-2,9.5e10,100e2,550
        # L,kL,os,d,rhoc = 2.94203531e+11, 3.38197506e-02, 9.52135030e+10, 3.00000000e+02,5.49900000e+02
        # L,kL,os,d,rhoc = 2.20486112e+01, 4.26473629e-02, 6.55373968e-01, 2.00000000e-02, 5.11228713e+02
        # L,kL,os,d,rhoc = 2.52880882e+01, 4.01341086e-02, 4.64676820e+00, 1.00000000e-02,5.13089576e+02 #good for 10m, little toofast for others
        # L,kL,os,d,rhoc = 2.86587915e+01, 3.76819172e-02, 8.54495830e+00, 9.67012244e-30, 5.14586341e+02 # real 
        # L,kL,os,d,rhoc = 3.09885751e+01, 3.78208923e-02, 7.89969584e+00, 7.53593340e-33, 5.14133552e+02
        Q=60.0e3
        # L,kL,os,d,rhoc = 2.20570722e+01, 3.69831602e-02, 1.00000000e+01, 3.46114982e-38,5.15734776e+02

        L,kL,os,d,rhoc = 22.13, 3.95e-02, 1.00000000e+01, 0 ,511.2 #Set cutoff to 60 for K values.
        
        kM = k_log_cf_d(self.rho,L,kL,os,d,rhoc)

# [3.54937382e+01 4.18361570e-02 9.22241253e+00 2.19293203e-29
 # 5.22153672e+02]

        # [4.11166976e+01 3.39248620e-02 0.00000000e+00 2.17863365e-02 5.47146094e+02]
        # 70kJ numbers: [ 9.06842e-02  1.5756e-01 -4.56035e-02  2.90756e-04 5.6125e+02]
        # 50kJ numbeer: [ 4.50626058e+03  1.60444934e-01 -2.79657914e+03  1.55511701e+01 5.60840881e+02]
        # 60kJ numbers: [ 1.97234192e1  1.65074446e-01 -1.29809477e1  7.05912848e-02 5.60986300e2]

        # kM = np.zeros_like(self.rho)
        # kM[self.rho<=550.0] = 165348093000.8031
        # kM[self.rho>550.0] = 363691127138.42395

        ####
        na = 1.0 #n exponent for age
        
        viscosity = (kM * (self.age / S_PER_YEAR)**na) / ((RHO_I - self.rho) * np.exp(-Q/(R*self.Tz)))
        
        # viscosity = (kM * (self.age / S_PER_YEAR)**na) / ((RHO_I - self.rho) * np.exp((-Q/(R*self.Tz)+(42.2e3/(R*222.15)))))
        
        viscosity[0] = viscosity[1]
        # print(viscosity[0])
        # top1m = np.nonzero(self.z <= 1.)
        # viscosity[top1m] = np.nanmean(viscosity[top1m])
        dr_dt = (0.5 * self.rho * self.sigma / viscosity)
        # dr_dt[top1m] = 
        
        self.RD['drho_dt'] = dr_dt # units are (kg m^-3) s^-1
        # if self.iii<10:
        #     print(dr_dt[0:20])
        return self.RD

    # def grainGrowth(self):
    #     '''
    #     :param Tz:
    #     :param Ts:
    #     :param iii:
    #     :param dt:
    #     :return r2: radius squared in m^2
    #     '''

    #     kgr = 1.3e-7 # grain growth rate from Arthern (2010), m^2/s
    #     Eg  = 42.4e3 # kJ/mol

    #     if self.MELT:
    #         porosity = 1 - self.rho / RHO_I 
    #         porespace = porosity * self.dz # meters

    #         # sat = self.LWC / porespace 
    #         sat = np.zeros_like(self.dz) # use 0 sat if our porespace is 0
    #         sat[np.where(porespace>0)[0]] = self.LWC[np.where(porespace>0)[0]] / porespace[np.where(porespace>0)[0]]

    #         if self.GrGrowPhysics == 'Katsushima':
    #             dr2_dt = 1e-9/(4*(self.r2)**0.5)*np.minimum(2/(np.pi)*(1.28e-8+4.22e-10*(sat*((1000*(RHO_I-self.rho)/(self.rho*RHO_I))*100))**3),6.94e-8)
    #         elif self.GrGrowPhysics == 'Arthern':
    #             dr2_dt = kgr * np.exp(-Eg / (R * self.Tz))

    #     else: # no MELT
    #         dr2_dt = kgr * np.exp(-Eg / (R * self.Tz)) #Arthern et al., 2010 grain growth, units are m^2/s

    #     r2 = self.r2 + dr2_dt * self.dt

    #     if ((self.calcGrainSize) and (self.bdotSec[self.iii]>0)): #VV if there is a new layer and we use Linow param
    #     #if self.calcGrainSize: # Apply initial grain size parameterisation from Linow et al., 2012: eqs (11) and (12)
    #         # uses mean annual T in [C] and mean annual bdot in [m w.e. yr-1]

    #         b0Lnw = 0.781
    #         b1Lnw = 0.0085
    #         b2Lnw = -0.279
    #         #r2_surface = ((b0Lnw+b1Lnw*(self.Ts[self.iii]-K_TO_C) + b2Lnw*(self.bdot_mean[0]*RHO_I/1000))*10**(-3))**2
    #         #r2_surface = ((b0Lnw+b1Lnw*(self.Ts[self.iii]-K_TO_C) + b2Lnw*(self.bdot_mean[-1]*RHO_I/1000))*10**(-3))**2 # More accurate to use bdot_mean[-1] as value of mean accumulation ??
    #         r2_surface = ((b0Lnw+b1Lnw*(self.T_mean[self.iii] - K_TO_C) + b2Lnw*(self.bdot_mean[-1]*RHO_I/1000))*10**(-3))**2 #VV
    #         r2 = np.concatenate(([r2_surface], r2[:-1]))

    #         # r2 = np.concatenate(([-2.42e-9 * self.Ts[self.iii] + 9.46e-7], r2[:-1])) # legacy code. Not sure where this equation is from. Gow 1967ish?

    #     elif (self.bdotSec[self.iii]>0): #VV if there is a new layer but we don't use Linow param
    #     #else: # use a fixed surface value, r2s0.

    #         #VV
    #         r2 = np.concatenate(([self.r2s0], r2[:-1])) # Rob Arthern's recommended value, personal communication.
    #         #r2 = np.concatenate(([self.r2s0 ** 2], r2[:-1])) # Rob Arthern's recommended value, personal communication.
    #     else: #VV if no new layer, no need for a new grain size at surface
    #         pass

    #     return r2, dr2_dt
    ### end grainGrowth ###
    #######################

    def surfacegrain(self):

        '''
        function to calculate the surface grain size at each time step
        '''

        if self.calcGrainSize: #VV if there is a new layer and we use Linow param
            #if self.calcGrainSize: # Apply initial grain size parameterisation from Linow et al., 2012: eqs (11) and (12)
                # uses mean annual T in [C] and mean annual bdot in [m w.e. yr-1]
            b0Lnw = 0.781
            b1Lnw = 0.0085
            b2Lnw = -0.279
            #r2_surface = ((b0Lnw+b1Lnw*(self.Ts[self.iii]-K_TO_C) + b2Lnw*(self.bdot_mean[0]*RHO_I/1000))*10**(-3))**2
            #r2_surface = ((b0Lnw+b1Lnw*(self.Ts[self.iii]-K_TO_C) + b2Lnw*(self.bdot_mean[-1]*RHO_I/1000))*10**(-3))**2 # More accurate to use bdot_mean[-1] as value of mean accumulation ??
            #r2_surface = ((b0Lnw+b1Lnw*(self.T_mean[self.iii]-K_TO_C) + b2Lnw*(self.bdot_mean[-1]*RHO_I/1000))*10**(-3))**2 #VV
            r2_surface = ((b0Lnw+b1Lnw*(self.T_mean[self.iii] - K_TO_C) + b2Lnw*(self.bdot_av[self.iii]*RHO_I/1000))*10**(-3))**2 #VV
            
            #r2 = np.concatenate(([r2_surface], r2[:-1])) #VV commented

            # r2 = np.concatenate(([-2.42e-9 * self.Ts[self.iii] + 9.46e-7], r2[:-1])) # legacy code. Not sure where this equation is from. Gow 1967ish?

        else: #VV if there is a new layer but we don't use Linow param
        # use a fixed surface value, r2s0.
            r2_surface = self.r2s0

        return r2_surface
    ### end surfacegrain ###
    ########################
        
    def graincalc(self,iii):
        
        '''
        Evolve the grain size

        This is the same as graingrowth except that we do not calculate for the surface grain, which is done by surfacegrain() function      
        '''

        kgr = 1.3e-7 # grain growth rate from Arthern (2010), m^2/s
        Eg  = 42.4e3 # kJ/mol

        if self.MELT:
            porosity = 1 - self.rho / RHO_I 
            porespace = porosity * self.dz # meters
            
            sat = np.zeros_like(self.dz) #VV use 0 sat if our porespace is 0
            sat[np.where(porespace>0)[0]] = self.LWC[np.where(porespace>0)[0]] / porespace[np.where(porespace>0)[0]] #VV 

            if self.GrGrowPhysics == 'Katsushima':
                dr2_dt = 1e-9/(4*(self.r2)**0.5)*np.minimum(2/(np.pi)*(1.28e-8+4.22e-10*(sat*((1000*(RHO_I-self.rho)/(self.rho*RHO_I))*100))**3),6.94e-8)
            elif self.GrGrowPhysics == 'Arthern':
                dr2_dt = kgr * np.exp(-Eg / (R * self.Tz))

        else: # no MELT
            dr2_dt = kgr * np.exp(-Eg / (R * self.Tz)) #Arthern et al., 2010 grain growth, units are m^2/s

        if np.any(dr2_dt<0):
            print('Negative grain growth at layers:',np.where(dr2_dt<0)[0])

        r2 = self.r2 + dr2_dt * self.dt

        return r2
    ### end GrainCalc ###
    #####################





