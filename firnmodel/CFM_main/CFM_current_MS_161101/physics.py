import numpy as np
from constants import *
from scipy import interpolate

# The standard parameters that get passed are:
# (iii, steps, gridLen, bdotSec, bdot_mean, bdot_type, Tz, T_mean, rho, sigma, dt, Ts, r2, physGrain):
# if you want to add physics that require more parameters, you need to change the 'PhysParams' dictionary in both the spin and nospin classes.

class FirnPhysics:
    def __init__(self,PhysParams):
        for k,v in PhysParams.items():
            setattr(self,k,v)

    def HL_dynamic(self):
        '''

        :param steps:
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

        A = self.bdotSec[self.iii] * self.steps * BDOT_TO_A
        drho_dt = np.zeros(self.gridLen)
        dr_dt = drho_dt

        dr_dt[self.rho < RHO_1]   = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * np.power(A, aHL) * 1000 / S_PER_YEAR
        drho_dt[self.rho < RHO_1] = dr_dt[self.rho < RHO_1]

        dr_dt[self.rho >= RHO_1]   = k2 * np.exp(-Q2 / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I_MGM - self.rho[self.rho >= RHO_1] / 1000) * np.power(A, bHL) * 1000 / S_PER_YEAR
        drho_dt[self.rho >= RHO_1] = dr_dt[self.rho >= RHO_1]

        return drho_dt

    def HL_Sigfus(self):
        '''

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

        A = self.bdotSec[self.iii] * self.steps * BDOT_TO_A
        drho_dt = np.zeros(self.gridLen)
        f550 = interpolate.interp1d(self.rho, self.sigma)
        sigma550 = f550(RHO_1)
        rhoDiff = (RHO_I_MGM - self.rho / 1000)

        k = np.power(k2 * np.exp(-Q2 / (R * self.Tz[self.rho >= RHO_1])), 2) / S_PER_YEAR
        sigmaDiff = (self.sigma[self.rho >= RHO_1] - sigma550)

        drho_dt[self.rho < RHO_1] = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * np.power(A, aHL) * 1000 / S_PER_YEAR

        drho_dt[(self.rho >= RHO_1) & (self.rho < RHO_I)]  = k * (sigmaDiff * rhoDiff[self.rho >= RHO_1]) / (GRAVITY * np.log((RHO_I_MGM - RHO_1 / 1000) / (rhoDiff[(self.rho >= RHO_1) & (self.rho < RHO_I)])))
        drho_dt[(self.rho >= RHO_1) & (self.rho >= RHO_I)] = 0

        # self.viscosity = np.ones(self.gridLen)

        return drho_dt

    def Li_2004(self):
        '''

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param T_mean:
        :param rho:

        :return drho_dt:
        '''

        A = self.bdotSec[self.iii] * self.steps * BDOT_TO_A

        # dr_dt = np.zeros(self.gridLen)
        dr_dt = (RHO_I - self.rho) * A * (139.21 - 0.542 * self.T_mean) * 8.36 * (K_TO_C - self.Tz) ** -2.061

        drho_dt = dr_dt / S_PER_YEAR

        # self.viscosity = np.ones(self.gridLen)

        return drho_dt

    def Li_2011(self):
        '''

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param bdot_mean:
        :param bdot_type:
        :param Tz:
        :param T_mean:
        :param rho:

        :return drho_dt:
        '''

        A     = self.bdotSec[self.iii] * self.steps * BDOT_TO_A
        TmC   = self.T_mean - K_TO_C
        beta1 = -9.788 + 8.996 * A - 0.6165 * TmC
        beta2 = beta1 / (-2.0178 + 8.4043 * A - 0.0932 * TmC)
        dr_dt = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            dr_dt[self.rho <= RHO_1] = (RHO_I - self.rho[self.rho <= RHO_1]) * A * beta1 * 8.36 * (K_TO_C - self.Tz[self.rho <= RHO_1]) ** -2.061
            dr_dt[self.rho > RHO_1]  = (RHO_I - self.rho[self.rho > RHO_1]) * A * beta2 * 8.36 * (K_TO_C - self.Tz[self.rho > RHO_1]) ** -2.061
        elif self.bdot_type == 'mean':
            dr_dt[self.rho <= RHO_1] = (RHO_I - self.rho[self.rho <= RHO_1]) * self.bdot_mean[self.rho <= RHO_1] * self.steps * BDOT_TO_A * beta1 * 8.36 * (K_TO_C - self.Tz[self.rho <= RHO_1]) ** -2.061
            dr_dt[self.rho > RHO_1]  = (RHO_I - self.rho[self.rho > RHO_1]) * self.bdot_mean[self.rho > RHO_1] * self.steps * BDOT_TO_A  * beta2 * 8.36 * (K_TO_C - self.Tz[self.rho > RHO_1]) ** -2.061

        drho_dt = dr_dt / S_PER_YEAR
        # self.viscosity = np.ones(self.gridLen)

        return drho_dt

    def Arthern_2010S(self):
        '''

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param bdot_type:
        :param bdot_mean:
        :param Tz:
        :param T_mean:
        :param rho:

        :return drho_dt:
        '''

        ar1 = 0.07
        ar2 = 0.03
        Ec  = 60.0e3
        Eg  = 42.4e3

        A_instant = self.bdotSec[self.iii] * self.steps * BDOT_TO_A * 1000
        A_mean_1 = self.bdot_mean[self.rho < RHO_1] * self.steps * BDOT_TO_A * 1000
        A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * self.steps * BDOT_TO_A * 1000
        dr_dt = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean))
        elif self.bdot_type == 'mean':
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean))

        drho_dt = dr_dt / S_PER_YEAR
        # self.viscosity = np.ones(self.gridLen)

        return drho_dt

    def Arthern_2010T(self):
        '''

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
        

        if not physGrain:
           print "Grain growth should be on for Arthern Transient"
           return


        drho_dt = np.zeros(self.gridLen)
        self.viscosity = np.zeros(self.gridLen)
    
        drho_dt[self.rho < RHO_1]  = kc1 * (RHO_I - self.rho[self.rho < RHO_1]) * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1])) * self.sigma[self.rho < RHO_1] / (self.r2[self.rho < RHO_1])
        # self.viscosity[self.rho < RHO_1] = ((1 / (2*kc1)) * (self.rho[self.rho < RHO_1] / (RHO_I - self.rho[self.rho < RHO_1])) * np.exp (Ec / (R * self.Tz[self.rho < RHO_1])) * (self.r2[self.rho < RHO_1])) / S_PER_YEAR
        
        drho_dt[self.rho >= RHO_1] = kc2 * (RHO_I - self.rho[self.rho >= RHO_1]) * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1])) * self.sigma[self.rho >= RHO_1] / (self.r2[self.rho >= RHO_1])
        # self.viscosity[self.rho >= RHO_1] = ((1 / (2*kc2)) * (self.rho[self.rho >= RHO_1] / (RHO_I - self.rho[self.rho >= RHO_1])) * np.exp (Ec / (R * self.Tz[self.rho >= RHO_1])) * (self.r2[self.rho >= RHO_1])) / S_PER_YEAR

        return drho_dt

    def Helsen_2008(self):
        '''

        :param steps:
        :param bdotSec:
        :param bdot_mean:
        :param Tz:
        :param Ts:
        :param rho:
        :param bdot_type:

        :return drho_dt:
        '''

        A_instant = self.bdotSec[self.iii] * self.steps * BDOT_TO_A
        A_mean = self.bdot_mean * self.steps * BDOT_TO_A

        if self.bdot_type == 'instant':
            dr_dt = (RHO_I - self.rho) * A_instant * (76.138 - 0.28965 * self.Tz) * 8.36 * (K_TO_C - self.Tz) ** -2.061
        elif self.bdot_type == 'mean':
            dr_dt = (RHO_I - self.rho) * A_mean * (76.138 - 0.28965 * self.Tz) * 8.36 * (K_TO_C - self.Tz) ** -2.061

        drho_dt = dr_dt / S_PER_YEAR
        # self.viscosity = np.ones(self.gridLen)

        return drho_dt

    def Simonsen_2013(self):
        '''

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param bdot_mean:
        :param bdot_type:
        :param Tz:
        :param T_mean:
        :param rho:

        :return drho_dt:
        '''
        ar1 = 0.07
        ar2 = 0.03
        Ec  = 60.0e3
        Eg  = 42.4e3
        F0  = 0.68
        F1  = 1.03

        A_instant = self.bdotSec[self.iii] * self.steps * BDOT_TO_A * 1000
        A_mean_1 = self.bdot_mean[self.rho < RHO_1] * self.steps * BDOT_TO_A * 1000
        A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * self.steps * BDOT_TO_A * 1000
        dr_dt = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            gamma = 61.7 / ((self.bdotSec * self.steps * S_PER_YEAR * RHO_I) ** (0.5)) * np.exp(-3800. / (R * self.T_mean))
            dr_dt[self.rho < RHO_1]  = F0 * (RHO_I - self.rho[self.rho < RHO_1]) * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean))
            dr_dt[self.rho >= RHO_1] = F1 * gamma * (RHO_I - self.rho[self.rho >= RHO_1]) * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean))
        elif self.bdot_type == 'mean':
            gamma = 61.7 / ((self.bdot_mean[self.rho >= RHO_1] * self.steps * S_PER_YEAR * RHO_I) ** (0.5)) * np.exp(-3800.0 / (R * self.T_mean))
            dr_dt[self.rho < RHO_1]  = F0 * (RHO_I - self.rho[self.rho < RHO_1]) * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean))
            dr_dt[self.rho >= RHO_1] = F1 * gamma * (RHO_I - self.rho[self.rho >= RHO_1]) * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean))

        drho_dt = dr_dt / S_PER_YEAR
        # self.viscosity = np.ones(self.gridLen)

        return drho_dt

    def Ligtenberg_2011(self):
        '''

        :param steps:
        :param gridLen:
        :param bdotSec:
        :param bdot_mean:
        :param bdot_type:
        :param Tz:
        :param T_mean:
        :param rho:

        :return drho_dt :
        '''
        ar1 = 0.07
        ar2 = 0.03
        Ec = 60.0e3
        Eg = 42.4e3

        A_instant = self.bdotSec * self.steps * BDOT_TO_A * 1000
        A_mean_1 = self.bdot_mean[self.rho < RHO_1] * self.steps * BDOT_TO_A * 1000
        A_mean_2 = self.bdot_mean[self.rho >= RHO_1] * self.steps * BDOT_TO_A * 1000
        dr_dt = np.zeros(self.gridLen)

        if self.bdot_type == 'instant':
            M_0 = 1.435 - 0.151 * np.log(self.bdotSec * S_PER_YEAR * RHO_I)
            M_1 = 2.366 - 0.293 * np.log(self.bdotSec * S_PER_YEAR * RHO_I)
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1])+ Eg / (R * self.T_mean))
        elif self.bdot_type == 'mean':
            M_0 = 1.435 - 0.151 * np.log(self.bdot_mean[self.rho < RHO_1] * S_PER_YEAR * RHO_I)
            M_1 = 2.366 - 0.293 * np.log(self.bdot_mean[self.rho >= RHO_1] * S_PER_YEAR * RHO_I)
            dr_dt[self.rho < RHO_1]  = (RHO_I - self.rho[self.rho < RHO_1]) * M_0 * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho < RHO_1]) + Eg / (R * self.T_mean))
            dr_dt[self.rho >= RHO_1] = (RHO_I - self.rho[self.rho >= RHO_1]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * self.Tz[self.rho >= RHO_1]) + Eg / (R * self.T_mean))

        drho_dt = dr_dt / S_PER_YEAR
        # self.viscosity = np.ones(self.gridLen)

        return drho_dt

    def Barnola_1991(self):
        '''

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

        self.rho[self.rho > RHO_I] = RHO_I #The Barnola model will go a fraction over the ice density (order 10^-3), so this stops that.
        drho_dt = np.zeros(self.gridLen)
        D = self.rho / RHO_I
        nBa = n * np.ones(self.gridLen)
        A0 = A0b * np.ones(self.gridLen) / 1.e18 #this is for the n=3 region.

        sigmaEff = self.sigma

        # zone 1
        A = self.bdotSec * self.steps * BDOT_TO_A
        drho_dt[self.rho < RHO_1] = dr_dt = k1 * np.exp(-Q1 / (R * self.Tz[self.rho < RHO_1])) * (RHO_I_MGM - self.rho[self.rho < RHO_1] / 1000) * np.power(A[self.iii], aHL) * 1000 / S_PER_YEAR

        # zone 2
        fe = 10.0 ** (alphaBarnola * (self.rho[self.rho <= RHO_2] / 1000) ** 3. + betaBarnola * (self.rho[self.rho <= RHO_2] / 1000) ** 2. + deltaBarnola * self.rho[self.rho <= RHO_2] / 1000 + gammaBarnola)
        drho_dt[self.rho <= RHO_2] = self.rho[self.rho <= RHO_2] * A0[self.rho <= RHO_2] * np.exp(-QBarnola / (R * self.Tz[self.rho <= RHO_2])) * fe * (sigmaEff[self.rho <= RHO_2] ** nBa[self.rho <= RHO_2])

        # zone 3
        fs = (3. / 16.) * (1 - self.rho[self.rho > RHO_2] / RHO_I) / (1 - (1 - self.rho[self.rho > RHO_2] / RHO_I) ** (1. / 3.)) ** 3.
        drho_dt[self.rho > RHO_2] = self.rho[self.rho > RHO_2] * A0[self.rho > RHO_2] * np.exp(-QBarnola / (R * self.Tz[self.rho > RHO_2])) * fs * (sigmaEff[self.rho > RHO_2] ** nBa[self.rho > RHO_2])
        
        # self.viscosity = np.ones(self.gridLen)

        return drho_dt


        
    def Morris_HL_2013(self):
        '''

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
        
        QMorris = 110.e3
        kMorris = 11.0
        rhoW    = 1000.0

        drho_dt = np.zeros(self.gridLen)
        # self.viscosity = np.zeros(self.gridLen)
        
        drho_dt[self.rho < RHO_1] = (kMorris / (rhoW * GRAVITY)) * ((RHO_I - self.rho[self.rho < RHO_1]) / self.rho[self.rho < RHO_1]) * (1 / self.Hx[self.rho < RHO_1]) * np.exp(-QMorris / (R * self.Tz[self.rho < RHO_1])) * self.sigma[self.rho < RHO_1]
        # self.viscosity[self.rho < RHO_1] = ((rhoW * GRAVITY) / (2*kMorris)) * (self.rho[self.rho < RHO_1] / (RHO_I - self.rho[self.rho < RHO_1])) * self.Hx [self.rho < RHO_1] * np.exp(-QMorris / (R * self.Tz[self.rho < RHO_1]))

        #Use HL Dynamic for zone 2 b/c Morris does not specify zone 2.
        Q2  = 21400.0
        k2  = 575.0
        bHL = 0.5
        A = self.bdotSec[self.iii] * self.steps * BDOT_TO_A

        drho_dt[self.rho >= RHO_1]   = k2 * np.exp(-Q2 / (R * self.Tz[self.rho >= RHO_1])) * (RHO_I_MGM - self.rho[self.rho >= RHO_1] / 1000) * np.power(A, bHL) * 1000 / S_PER_YEAR
        # self.viscosity[self.rho >= RHO_1] = -(np.power(GRAVITY, (0.5)) / 2*k2) * np.exp(Q2 / (R*self.Tz[self.rho >= RHO_1]))*((self.rho[self.rho >= RHO_1]/(RHO_I_MGM - self.rho[self.rho >= RHO_1]))) * np.power (self.sigma[self.rho >= RHO_1], 0.5) * np.power(self.age[self.rho >= RHO_1], (0.5)) * 1000 / S_PER_YEAR
       
        return drho_dt

    def Spencer_2001(self):
        '''
        :return:
        '''
        pass

    def Goujon_2003(self):
        '''
        :return:
        '''
        pass

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
        self.Hx = self.Hx + np.exp(-110.0e3/(R*self.Tz))*self.dt
        Hx_new = np.exp(-110.0e3/(R*self.Tz[0]))*self.dt 
        self.Hx = np.concatenate(([Hx_new],self.Hx[:-1]))
        return self.Hx






