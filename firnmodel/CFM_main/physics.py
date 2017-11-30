import numpy as np
from constants import *
from scipy import interpolate
import sys

# The standard parameters that get passed are:
# (iii, steps, gridLen, bdotSec, bdot_mean, bdot_type, Tz, T10m, rho, sigma, dt, Ts, r2, physGrain):
# if you want to add physics that require more parameters, you need to change the 'PhysParams' dictionary in both the spin and nospin classes.
#Hello Emma!
class FirnPhysics:

    def __init__(self,PhysParams):

        '''

        bdot_mean units are m I.E./year
        bdotSec units are m I.E./time step

        '''
        for k,v in list(PhysParams.items()):
            setattr(self,k,v)
        self.RD = {}

    def HL_dynamic(self):
        '''

        Accumulation units are m W.E. per year

        :param steps: # of steps per year
        :param gridLen:
        :param bdotSec:
        :param Tz:
        :param rho:

        :return drho_dt:
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

        if self.bdot_type == 'instant':
            drho_dt[self.rho < RHO_1]     = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * A_instant**aHL * 1000 / S_PER_YEAR
            drho_dt[self.rho >= RHO_1]    = k2 * np.exp(-Q2 / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I_MGM - self.rho[self.rho >= RHO_1] / 1000) * A_instant**bHL * 1000 / S_PER_YEAR

        elif self.bdot_type == 'mean':
            drho_dt[self.rho < RHO_1]     = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * (A_mean[self.rho < RHO_1])**aHL * 1000 / S_PER_YEAR
            drho_dt[self.rho >= RHO_1]    = k2 * np.exp(-Q2 / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I_MGM - self.rho[self.rho >= RHO_1] / 1000) * (A_mean[self.rho >= RHO_1])**bHL * 1000 / S_PER_YEAR

        self.RD['drho_dt'] = drho_dt
        return self.RD
        # return drho_dt

    def HL_Sigfus(self):
        '''

        Accumulation units are m W.E. per year (zone 1); uses stress for zone 2

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param Tz:
        :param rho:
        :param sigma:

        :return drho_dt:
        '''
        Q1  = 10160.0
        Q2  = 21400.0
        k1  = 11.0
        k2  = 575.0
        aHL = 1.0

        A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean = self.bdot_mean * RHO_I_MGM

        drho_dt = np.zeros(self.gridLen)
        f550 = interpolate.interp1d(self.rho, self.sigma)
        sigma550 = f550(RHO_1)
        rhoDiff = (RHO_I_MGM - self.rho / 1000)



        k = np.power(k2 * np.exp(-Q2 / (R * self.Tz[(self.rho >= RHO_1) & (self.rho < RHO_I)])), 2) / S_PER_YEAR
        sigmaDiff = (self.sigma[(self.rho >= RHO_1) & (self.rho < RHO_I)] - sigma550)
        # print('sigmaDiff',len(sigmaDiff))
        if self.bdot_type == 'instant':
            drho_dt[self.rho < RHO_1] = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * A_instant**aHL * 1000 / S_PER_YEAR
        elif self.bdot_type == 'mean':
            drho_dt[self.rho < RHO_1] = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * (A_mean[self.rho < RHO_1])**aHL * 1000 / S_PER_YEAR
        # print('rhoDiff',len(rhoDiff[(self.rho >= RHO_1) & (self.rho < RHO_I)]))
        # print('drho',len(drho_dt[(self.rho >= RHO_1) & (self.rho < RHO_I)]))
        # print('k',len(k))
        drho_dt[(self.rho >= RHO_1) & (self.rho < RHO_I)]  = k * (sigmaDiff * rhoDiff[(self.rho >= RHO_1) & (self.rho < RHO_I)]) / (GRAVITY * np.log((RHO_I_MGM - RHO_1 / 1000) / (rhoDiff[(self.rho >= RHO_1) & (self.rho < RHO_I)])))

        drho_dt[(self.rho >= RHO_1) & (self.rho >= RHO_I)] = 0

        # drho_dt[self.rho >= RHO_1]  = k * (sigmaDiff * rhoDiff[self.rho >= RHO_1]) / (GRAVITY * np.log((RHO_I_MGM - RHO_1 / 1000) / (rhoDiff[self.rho >= RHO_1])))
        
        # drho_dt[(self.rho >= RHO_1) and (self.rho >= RHO_I)] = 0

        # self.viscosity = np.ones(self.gridLen)
        self.RD['drho_dt'] = drho_dt
        return self.RD

        # return drho_dt

    def Li_2004(self):
        '''
        Accumulation units are m W.E. per year (?)
        Equation from Arthern, 2010 (eq. 2): not sure where Rob got that? 
        (Arthern implies accumulation is m I.E./year for bdot; not doing that here.)
        Paper would lead me to believe that it is m W.E.
        
        ***Needs to have the vapor flux coded in if we want to use these physics properly.

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param T10m:
        :param rho:

        :return drho_dt:
        '''

        A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean = self.bdot_mean * RHO_I_MGM

        # dr_dt = np.zeros(self.gridLen)
        
        if self.bdot_type == 'instant':
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Li and Zwally 2004 physics")
            dr_dt = (RHO_I - self.rho) * A_instant * (139.21 - 0.542 * self.T10m) * 8.36 * (K_TO_C - self.Tz) ** -2.061
        elif self.bdot_type == 'mean':
            dr_dt = (RHO_I - self.rho) * A_mean * (139.21 - 0.542 * self.T10m) * 8.36 * (K_TO_C - self.Tz) ** -2.061   

        drho_dt = dr_dt / S_PER_YEAR

        # self.viscosity = np.ones(self.gridLen)
        self.RD['drho_dt'] = drho_dt
        return self.RD
        # return drho_dt

    def Li_2011(self):
        '''
        Accumulation units are m W.E. per year (?)
        Potentially beta should be calculated with the long-term mean accumulation rate

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param bdot_mean:
        :param bdot_type:
        :param Tz:
        :param T10m:
        :param rho:

        :return drho_dt:
        '''

        A_instant     = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean = self.bdot_mean * RHO_I_MGM


        # TmC   = self.T10m - K_TO_C
        TmC   = self.T_mean - K_TO_C

        dr_dt = np.zeros(self.gridLen)

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
        # self.viscosity = np.ones(self.gridLen)
        self.RD['drho_dt'] = drho_dt
        return self.RD

        # return drho_dt

    def Arthern_2010S(self):
        '''
        This is the steady-state solution described in the main text of Arthern et al. (2010)
        Accumulation units are kg/m^2/year

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param bdot_type:
        :param bdot_mean:
        :param Tz:
        :param T10m:
        :param rho:

        :return drho_dt:
        '''

        ar1 = 0.07
        ar2 = 0.03
        Ec  = 60.0e3
        Eg  = 42.4e3

        A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM * 1000
        A_mean_1 = self.bdot_mean[self.rho < RHO_1] * RHO_I_MGM * 1000
        A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * RHO_I_MGM * 1000
        dr_dt = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Arthern 2010 physics")
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T10m))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T10m))
        elif self.bdot_type == 'mean':
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T10m))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T10m))

        drho_dt = dr_dt / S_PER_YEAR
        # self.viscosity = np.ones(self.gridLen)
        self.RD['drho_dt'] = drho_dt
        return self.RD
        # return drho_dt

    def Arthern_2010T(self):
        '''
        This is the transient solution described in the appendix of Arthern et al. (2010)

        Uses stress rather than accumulation rate.

        :param gridLen: 
        :param Tz:
        :param rho:
        :param sigma:
        :param r2:
        :param physGrain:

        :return drho_dt:
        :return viscosity_z1:
        :return viscosity_z2:
        '''

        kc1 = 9.2e-9
        kc2 = 3.7e-9
        Ec  = 60.0e3
        

        if not self.physGrain:
           print("Grain growth should be on for Arthern Transient")
           return


        drho_dt = np.zeros(self.gridLen)
        self.viscosity = np.zeros(self.gridLen)
    
        drho_dt[self.rho < RHO_1]  = kc1 * (RHO_I - self.rho[self.rho < RHO_1]) * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1])) * self.sigma[self.rho < RHO_1] / (self.r2[self.rho < RHO_1])
        # self.viscosity[self.rho < RHO_1] = ((1 / (2*kc1)) * (self.rho[self.rho < RHO_1] / (RHO_I - self.rho[self.rho < RHO_1])) * np.exp (Ec / (R * self.Tz[self.rho < RHO_1])) * (self.r2[self.rho < RHO_1])) / S_PER_YEAR
        
        drho_dt[self.rho >= RHO_1] = kc2 * (RHO_I - self.rho[self.rho >= RHO_1]) * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1])) * self.sigma[self.rho >= RHO_1] / (self.r2[self.rho >= RHO_1])
        # self.viscosity[self.rho >= RHO_1] = ((1 / (2*kc2)) * (self.rho[self.rho >= RHO_1] / (RHO_I - self.rho[self.rho >= RHO_1])) * np.exp (Ec / (R * self.Tz[self.rho >= RHO_1])) * (self.r2[self.rho >= RHO_1])) / S_PER_YEAR
        self.RD['drho_dt'] = drho_dt
        return self.RD
        # return drho_dt

    def Helsen_2008(self):
        '''
        Accumulation units are m W.E. per year (?)
        Equation is from Arthern et al. 2010 (2)
        (Arthern implies units are m I.E.; not doing that here) 


        :param steps:
        :param bdotSec:
        :param bdot_mean:
        :param Tz:
        :param Ts:
        :param rho:
        :param bdot_type:

        :return drho_dt:
        '''

        A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean = self.bdot_mean * RHO_I_MGM


        trigger = False

        if self.bdot_type == 'instant':
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Helsen 2008 physics")            
            dr_dt = (RHO_I - self.rho) * A_instant * (76.138 - 0.28965 * self.T_mean) * 8.36 * (K_TO_C - self.Tz) ** -2.061
        elif self.bdot_type == 'mean':
            dr_dt = (RHO_I - self.rho) * A_mean * (76.138 - 0.28965 * self.T_mean) * 8.36 * (K_TO_C - self.Tz) ** -2.061
            if (dr_dt<0).any():
                print('negative dr_dt at ', self.iii)
                print('depth:', self.z[dr_dt<0])
                trigger = True
        
        if (self.iii>1015 or trigger == True):
            print(self.iii)
            print('A_mean', max(A_mean), min(A_mean))
            print('Tz', max(self.Tz), min(self.Tz))
            print('Tz1st10', self.Tz[0:10])
            yy = np.where(self.Tz==max(self.Tz))[0][0]

            print('maxTz', yy)
            print('mTzdep', self.z[yy])
            print('T_mean', self.T_mean)
            print('rhos', max(self.rho),min(self.rho))
            print('maxdr', max(dr_dt))
            print('mindr', min(dr_dt))
            print('dr_dt', dr_dt[0:20])
            print('!!!!!', self.iii)

        drho_dt = dr_dt / S_PER_YEAR
        # self.viscosity = np.ones(self.gridLen)
        self.RD['drho_dt'] = drho_dt
        return self.RD

        # return drho_dt

    def Simonsen_2013(self):
        '''
        Accumulation units are kg/m^2/year

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param bdot_mean:
        :param bdot_type:
        :param Tz:
        :param T10m:
        :param rho:

        :return drho_dt:
        '''
        ar1 = 0.07
        ar2 = 0.03
        Ec  = 60.0e3
        Eg  = 42.4e3
        F0  = 0.68 #firnmice value?
        F1  = 1.03 #firnmice value?
        # F0=0.8 # Simonsen's recommended (email correspondence)
        # F1=1.25 # Simonsen's recommended (email correspondence)


        A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM * 1000
        A_mean_1 = self.bdot_mean[self.rho < RHO_1] * RHO_I_MGM * 1000
        A_mean_2 = self.bdot_mean[self.rho >= RHO_1]* RHO_I_MGM * 1000
        dr_dt = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Simonsen physics")
            gamma = 61.7 / (A_instant ** (0.5)) * np.exp(-3800. / (R * self.T10m))
            dr_dt[self.rho < RHO_1]  = F0 * (RHO_I - self.rho[self.rho < RHO_1]) * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T10m))
            dr_dt[self.rho >= RHO_1] = F1 * gamma * (RHO_I - self.rho[self.rho >= RHO_1]) * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T10m))
        elif self.bdot_type == 'mean':
            gamma = 61.7 / (A_mean_2 ** (0.5)) * np.exp(-3800.0 / (R * self.T10m))
            dr_dt[self.rho < RHO_1]  = F0 * (RHO_I - self.rho[self.rho < RHO_1]) * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T10m))
            dr_dt[self.rho >= RHO_1] = F1 * gamma * (RHO_I - self.rho[self.rho >= RHO_1]) * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T10m))

        drho_dt = dr_dt / S_PER_YEAR
        # self.viscosity = np.ones(self.gridLen)
        self.RD['drho_dt'] = drho_dt
        return self.RD
        # return drho_dt

    def Ligtenberg_2011(self):
        '''

        Units are mm W.E. per year
        b_dot is meant to be accumulation over a reference period (20 years for spin up, 1 year for regular?) (not mean over the lifetime  of a parcel)

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param bdot_mean:
        :param bdot_type:
        :param Tz:
        :param T10m:
        :param rho:

        :return drho_dt :
        '''
        ar1 = 0.07
        ar2 = 0.03
        Ec = 60.0e3
        Eg = 42.4e3

        


        dr_dt = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Ligtenberg 2011 physics")
            A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM * 1000
            M_0 = 1.435 - 0.151 * np.log(A_instant)
            M_1 = 2.366 - 0.293 * np.log(A_instant)
            M_0 = np.max((0.25,M_0))
            M_1 = np.max((0.25,M_1))
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T10m))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1])+ Eg / (R * self.T10m))
        elif self.bdot_type == 'mean':
            A_mean_1 = self.bdot_mean[self.rho < RHO_1] * RHO_I
            A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * RHO_I
            M_0 = 1.435 - 0.151 * np.log(A_mean_1)
            M_1 = 2.366 - 0.293 * np.log(A_mean_2)
            M_0[M_0<0.25]=0.25
            M_1[M_1<0.25]=0.25
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T10m))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T10m))


            # dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.Tz[self.rho < RHO_1]))
            # dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.Tz[self.rho >= RHO_1]))
        drho_dt = dr_dt / S_PER_YEAR

        # self.viscosity = np.ones(self.gridLen)
        self.RD['drho_dt'] = drho_dt
        return self.RD

        # return drho_dt

    def Barnola_1991(self):
        '''

        uses m W.E. (zone 1) and stress (zone 2)

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param Tz:
        :param rho:
        :param sigma:

        :return drho_dt:
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
        drho_dt = np.zeros(self.gridLen)
        D = self.rho / RHO_I
        nBa = n * np.ones(self.gridLen)
        A0 = A0b * np.ones(self.gridLen) / 1.e18 #this is for the n=3 region.

        sigmaEff = self.sigma

        ### Zone 1 ###
        A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean_1 = self.bdot_mean[self.rho < RHO_1] * RHO_I_MGM

        if self.bdot_type == 'instant':
            drho_dt[self.rho < RHO_1] = dr_dt = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * A_instant ** aHL * 1000 / S_PER_YEAR
        elif self.bdot_type == 'mean':
            drho_dt[self.rho < RHO_1] = dr_dt = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * A_mean_1 ** aHL * 1000 / S_PER_YEAR

        ### Zone 2 ###
        fe = 10.0 ** (alphaBarnola * (self.rho[self.rho <= RHO_2] / 1000) ** 3. + betaBarnola * (self.rho[self.rho <= RHO_2] / 1000) ** 2. + deltaBarnola * self.rho[self.rho <= RHO_2] / 1000 + gammaBarnola)
        drho_dt[self.rho <= RHO_2] = self.rho[self.rho <= RHO_2] * A0[self.rho <= RHO_2] * np.exp(-QBarnola / (R * self.Tz[self.rho <= RHO_2])) * fe * (sigmaEff[self.rho <= RHO_2] ** nBa[self.rho <= RHO_2])

        # zone 3
        fs = (3. / 16.) * (1 - self.rho[self.rho > RHO_2] / RHO_I) / (1 - (1 - self.rho[self.rho > RHO_2] / RHO_I) ** (1. / 3.)) ** 3.
        drho_dt[self.rho > RHO_2] = self.rho[self.rho > RHO_2] * A0[self.rho > RHO_2] * np.exp(-QBarnola / (R * self.Tz[self.rho > RHO_2])) * fs * (sigmaEff[self.rho > RHO_2] ** nBa[self.rho > RHO_2])
        
        # self.viscosity = np.ones(self.gridLen)
        self.RD['drho_dt'] = drho_dt
        return self.RD

        # return drho_dt
    
    def Morris_HL_2014(self):
        '''

        Uses stress instead of accumulation.

        Need to choose physics for zone 2. Herron and Langway here.

        4/7/17: The zone 1 physics need to be examined. We currently ignore the m term (equation 7), assuming that it is zero for a steady-state.
        Dividing drho/dt by 100 for zone 1 gives reasonable numbers, indicating that we might need to check units of everything.

        :param steps:
        :param spin:
        :param iii:
        :param gridLen:
        :param Tz:
        :param dt:
        :param rho:
        :param sigma:
        :param age:
        :param bdotSec:
        :param Hx:

        :return drho_dt:
        '''
        if self.iii ==0:
            print('CAUTION: MORRIS PHYSICS ARE STILL UNDER CODE DEVELOPMENT!')
            print('see physics.py for more details')     

        QMorris = 110.e3
        kMorris = 11.0

        drho_dt = np.zeros(self.gridLen)
        # self.viscosity = np.zeros(self.gridLen)

        
        drho_dt[self.rho < RHO_1] = (kMorris / (RHO_W_KGM * GRAVITY)) * ((RHO_I - self.rho[self.rho < RHO_1])) * (1 / self.Hx[self.rho < RHO_1]) * np.exp(-QMorris / (R * self.Tz[self.rho < RHO_1])) * self.sigma[self.rho < RHO_1] / 100
        
        #Use HL Dynamic for zone 2 b/c Morris does not specify zone 2.
        Q2  = 21400.0
        k2  = 575.0
        bHL = 0.5
        A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM
        A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * RHO_I_MGM
        if self.bdot_type == 'instant':
            drho_dt[self.rho >= RHO_1]   = k2 * np.exp(-Q2 / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I_MGM - self.rho[self.rho >= RHO_1] / 1000) * A_instant ** bHL * 1000 / S_PER_YEAR
        elif self.bdot_type == 'mean':
            drho_dt[self.rho >= RHO_1]   = k2 * np.exp(-Q2 / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I_MGM - self.rho[self.rho >= RHO_1] / 1000) * A_mean_2 ** bHL * 1000 / S_PER_YEAR
        
        self.RD['drho_dt'] = drho_dt
        return self.RD
        # return drho_dt

    def KuipersMunneke_2015(self):
        '''

        Units are mm W.E. per year
        b_dot is meant to be accumulation over a reference period (20 years for spin up, 1 year for regular?) (not mean over the lifetime  of a parcel)

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param bdot_mean:
        :param bdot_type:
        :param Tz:
        :param T10m:
        :param rho:

        :return drho_dt :
        '''
        ar1 = 0.07
        ar2 = 0.03
        Ec = 60.0e3
        Eg = 42.4e3

        # yrmn = 20
        # if self.iii<=yrmn*12:
        #     A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM * 1000
        # else:
        #     A_instant = np.mean(self.bdotSec[self.iii-12*yrmn:self.iii+1]) * self.steps * S_PER_YEAR * RHO_I_MGM * 1000

        dr_dt = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            A_instant = self.bdotSec[self.iii] * self.steps * S_PER_YEAR * RHO_I_MGM * 1000
            if self.iii==0:
                print("It is not recommended to use instant accumulation with Ligtenberg 2011 physics")
            M_0 = 1.042 - 0.0916 * np.log(A_instant)
            M_1 = 1.734 - 0.2039 * np.log(A_instant)
            M_0 = np.max((0.25,M_0))
            M_1 = np.max((0.25,M_1))
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T10m))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1])+ Eg / (R * self.T10m))

        elif self.bdot_type == 'mean':
            A_mean_1 = self.bdot_mean[self.rho < RHO_1] * RHO_I
            A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * RHO_I

            # print "Amean", A_mean_1[0:10]

            M_0 = 1.042 - 0.0916 * np.log(A_mean_1)
            M_1 = 1.734 - 0.2039 * np.log(A_mean_2)
            # print 'M_0', M_0
            # M_0 = np.max((0.25,M_0))
            # M_1 = np.max((0.25,M_1))

            M_0[M_0<0.25]=0.25
            M_1[M_1<0.25]=0.25

            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T10m))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T10m))

        elif self.bdot_type == 'stress':
            # A_mean_1 = self.bdot_mean[self.rho < RHO_1] * RHO_I
            # A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * RHO_I
            A_mean_1 = self.mass[self.rho < RHO_1]*10
            A_mean_2 = self.mass[self.rho >= RHO_1]*10
            # print "Amean", A_mean_1[0:10]

            M_0 = 1.042 - 0.0916 * np.log(A_mean_1)
            # print('Amean',A_mean_1[0:10])
            # print('mass',self.mass[0:10])
            # input()
            M_1 = 1.734 - 0.2039 * np.log(A_mean_2)
            # print 'M_0', M_0
            # M_0 = np.max((0.25,M_0))
            # M_1 = np.max((0.25,M_1))

            M_0[M_0<0.25]=0.25
            M_1[M_1<0.25]=0.25

            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T10m))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T10m))

        drho_dt = dr_dt / S_PER_YEAR
        drho_dt[self.rho>=RHO_I] = 0
        # self.viscosity = np.ones(self.gridLen)
        self.RD['drho_dt'] = drho_dt
        return self.RD
        # return drho_dt

    def Goujon_2003(self):
        '''

        Uses stress

        :return:
        '''

        global Gamma_Gou, Gamma_old_Gou, Gamma_old2_Gou, ind1_old

        atmosP = 101325.0 # Atmospheric Pressure
        dDdt        = np.zeros(self.gridLen) # Capital D is change in relative density

        top2m           = np.nonzero(self.z <= 1.)
        
        self.rho[top2m] = self.rhos0 # top 2 meters of Goujon model are always set to surface density
        
        sigma_MPa   = self.sigma / (1.0e6)
        sigma_bar   = self.sigma / (1.0e5)
        Qgj         = 60.0e3
        n           = 3.0 

        rhoi2cgs    = .9165 * (1.-1.53e-4 * (self.T10m - 273.15)) # Density of ice, temperature dependent, g/cm^3
        rhoi2       = rhoi2cgs * 1000.0 # density of ice, kg/m^3
        
        D           = self.rho / rhoi2 # Relative density
        Dm23        = 0.9 #transition from zone 2 to 3
        rho23       = Dm23 * rhoi2 #density of 2/3 transition
        
        D0          = 0.00226 * self.T10m + 0.03 #D0 is the relative density of transition from zone 1 to 2. Here is from Arnaud et al. (2000) eq. 8, not Goujon (Goujon uses C, not K)

        if D0 > 0.59: #Model requires zone 1/2 transition to be less than D=0.6
            D0      = 0.59
        Dms         = D0 + 0.009 #D0 + epsilon factor, maximum is 0.599
        Dmsrho      = Dms * rhoi2 # density of zone 1/2 transition
        
        ind1        = np.argmax(D >= Dms) #first index where the density is greater than or equal to the transition
        Dm          = D[ind1] #actual transition relative density. Use this so that the transition falls at a node
        Dmrho      = Dm * rhoi2 #density of first node, zone 2

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

        # gamma_An=(5.3*A[ind1] * (Dms**2*Dms)**(1.0/3.0) * (a[ind1]/np.pi)**(1.0/2.0) * (sigmastar[ind1]/3.0)**n) / ((sigma_bar[ind1]/(Dms**2))*(1-(5.0/3.0*Dms))); 
        gamma_An=(5.3*A[ind1] * (Dms**2*D0)**(1.0/3.0) * (a[ind1]/np.pi)**(1.0/2.0) * (sigmastar[ind1]/3.0)**n) / ((sigma_bar[ind1]/(Dms**2))*(1-(5.0/3.0*Dms))); 

        if self.iii == 0 or ind1 != ind1_old:
            Gamma_Gou       = 0.5 / S_PER_YEAR
            Gamma_old_Gou   = Gamma_Gou
        else:
            Gamma_Gou       = Gamma_old_Gou

        # Gamma_Gou = 0.5 / S_PER_YEAR

        dDdt[0:ind1+1]=Gamma_Gou*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
        # dDdt[0:ind1+1]=gamma_An*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
        dDdt[ind1+1:]=5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n         
        gfrac       = 0.03
        gam_div     = 1 + gfrac #change this if want: making it larger will make the code run faster. Must be >=1.
        
        ########## iterate to increase gamma first if not in steady state    
        if self.iii != 0 and dDdt[ind1] <= dDdt[ind1+1] and Gamma_Gou!=Gamma_old2_Gou: #and dDdt_old[ind1]!=dDdt_old[ind1]:
            cc = 1
            while dDdt[ind1] < dDdt[ind1 + 1]:
                Gamma_Gou       = Gamma_Gou * (gam_div)
                dDdt[0:ind1+1]=Gamma_Gou*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
                dDdt[ind1+1:]=5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n

                cc = cc + 1
                if cc>10000:
                    print('Goujon is not converging. exiting')
                    sys.exit()
                # if cc>1000:
                #     print 'cc', cc
                #     print 'dDdt',dDdt[ind1:ind1+2]
                #     dd = dDdt[D<=Dm]
                #     ee = dDdt[D>Dm]
                #     print 'dd',dd[-1]
                #     print 'ee',ee[0]
                #     raw_input()

        ### then iterate to find the maximum value of gamma that will make a continuous drho/dt
        counter = 1
        while dDdt[ind1] >= dDdt[ind1 + 1]:
            # print 'iterating', counter
            Gamma_Gou       = Gamma_Gou / (1 + gfrac/2.0)
            dDdt[0:ind1+1]=Gamma_Gou*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
            dDdt[ind1+1:]=5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n
            counter = counter +1
            if counter>10000:
                print('Goujon is not converging. exiting')
                sys.exit()
            # if counter >100:
            #     print 'dDdt',dDdt[ind1:ind1+2]
            #     print 'counter',counter

        # #dDdt[D<=Dm]=gamma_An*(sigma_bar[D<=Dm])*(1.0-(5.0/3.0)*D[D<=Dm])/((D[D<=Dm])**2.0)
        if self.iii<10:
            print('dDdt',dDdt[ind1:ind1+2])
        Gamma_old2_Gou  = Gamma_old_Gou
        Gamma_old_Gou  = Gamma_Gou
        ind1_old = ind1
        #####################


                ########## iterate to increase gamma first if not in steady state    
        # if dDdt[ind1] <= dDdt[ind1+1]: #and dDdt_old[ind1]!=dDdt_old[ind1]:
        #     cc = 1
        #     while dDdt[ind1] < dDdt[ind1 + 1]:
        #         Gamma_Gou       = Gamma_Gou * (gam_div)
        #         dDdt[0:ind1+1]=Gamma_Gou*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
        #         dDdt[ind1+1:]=5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n

        #         cc = cc + 1
        #         if cc>10000:
        #             print 'Goujon is not converging. exiting'
        #             sys.exit()
        #         # if cc>1000:
        #         #     print 'cc', cc
        #         #     print 'dDdt',dDdt[ind1:ind1+2]
        #         #     dd = dDdt[D<=Dm]
        #         #     ee = dDdt[D>Dm]
        #         #     print 'dd',dd[-1]
        #         #     print 'ee',ee[0]
        #         #     raw_input()

        # ### then iterate to find the maximum value of gamma that will make a continuous drho/dt
        # counter = 1
        # while dDdt[ind1] >= dDdt[ind1 + 1]:
        #     # print 'iterating', counter
        #     Gamma_Gou       = Gamma_Gou / (1 + gfrac/2.0)
        #     dDdt[0:ind1+1]=Gamma_Gou*(sigma_bar[0:ind1+1])*(1.0-(5.0/3.0)*D[0:ind1+1])/((D[0:ind1+1])**2.0)
        #     dDdt[ind1+1:]=5.3*A[ind1+1:]* (((D[ind1+1:]**2.0)*D0)**(1/3.)) * (a[ind1+1:]/np.pi)**(1.0/2.0) * (sigmastar[ind1+1:]/3.0)**n
        #     counter = counter +1
        #     if counter>10000:
        #         print 'Goujon is not converging. exiting'
        #         sys.exit()
        #     # if counter >100:
        #     #     print 'dDdt',dDdt[ind1:ind1+2]
        #     #     print 'counter',counter

        # # dDdt[D<=Dm]=gamma_An*(sigma_bar[D<=Dm])*(1.0-(5.0/3.0)*D[D<=Dm])/((D[D<=Dm])**2.0)
        # if self.iii<10:
        #     print 'dDdt',dDdt[ind1:ind1+2]
        # # Gamma_old2_Gou  = Gamma_old_Gou
        # # Gamma_old_Gou  = Gamma_Gou
        # # ind1_old = ind1
        #####################

        
        rhoC        = RHO_2 #should be Martinerie density
        frho2       = interpolate.interp1d(self.rho,sigma_bar)
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

        dDdt[D>Dm23] = 2.*A[D>Dm23] * ( (D[D>Dm23]*(1-D[D>Dm23])) / (1-(1-D[D>Dm23])**(1/n))**n ) * (2*sigmaEff[D>Dm23]/n)**3.0
        Ad              = 1.2e3 * np.exp(-Qgj / (R * self.Tz)) * 1.0e-1
        T34             = 0.98
        dDdt[D>T34] = 9/4*Ad[D>T34]*(1-D[D>T34])*sigmaEff[D>T34]

        dDdt_old        = dDdt
        drho_dt         = dDdt*rhoi2
        drho_dt[top2m]  = 0.0
        
        self.RD['drho_dt'] = drho_dt
        return self.RD
        # return drho_dt

    def Crocus(self):
        '''

        Uses stress

        :return:
        '''

        f1 = 1.0 # unitless
        f2 = 4.0 # unitless
        nu_0 = 7.62237e6 # kg s^-1
        a_n = 0.1 # K^-1
        b_n = 0.023 # m^3 kg^-1
        c_n = 250 # kg m^-3

        viscosity = f1 * f2 * nu_0 * self.rho / c_n * np.exp(a_n * (273.15 - self.Tz) + b_n * self.rho)
        dr_dt = self.rho * self.sigma / viscosity

        dr_dt[self.rho>=RHO_I] = 0

        drho_dt = dr_dt #/ S_PER_YEAR
        
        self.RD['drho_dt'] = drho_dt
        return self.RD
        # return drho_dt

    def grainGrowth(self):
        '''

        :param Tz:
        :param Ts:
        :param iii:
        :param dt:

        :return r2:
        '''

        kgr = 1.3e-7 															# grain growth rate from Arthern (2010)
        Eg  = 42.4e3
        

        dr2_dt = kgr * np.exp(-Eg / (R * self.Tz))
        r2 = self.r2 + dr2_dt * self.dt
    #     r2_time = np.concatenate(([t * iii + 1], r2))

        if self.calcGrainSize:												# uses surface temperature to get an initial grain size
            r2 = np.concatenate(([-2.42e-9 * self.Ts[self.iii] + 9.46e-7], r2[:-1]))
        else:
            r2 = np.concatenate(([(0.1e-3) ** 2], r2[:-1]))

        return r2

    def THistory(self):
        self.Hx = self.Hx + np.exp(-110.0e3 / (R * self.Tz)) * self.dt
        Hx_new  = np.exp(-110.0e3 / (R * self.Tz[0])) * self.dt 
        self.Hx = np.concatenate(([Hx_new],self.Hx[:-1]))
        return self.Hx






