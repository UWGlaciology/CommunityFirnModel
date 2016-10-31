import numpy as np
from constants import *
from scipy import interpolate

def HL_dynamic(steps, gridLen, bdotSec, Tz, rho):
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

    A = bdotSec * steps * BDOT_TO_A
    drho_dt = np.zeros(gridLen)
    dr_dt = drho_dt

    dr_dt[rho < RHO_1]   = k1 * np.exp(-Q1 / (R * Tz[rho < RHO_1])) * (RHO_I_MGM - rho[rho < RHO_1] / 1000) * np.power(A, aHL) * 1000 / S_PER_YEAR
    drho_dt[rho < RHO_1] = dr_dt[rho < RHO_1]

    dr_dt[rho >= RHO_1]   = k2 * np.exp(-Q2 / (R * Tz[rho >= RHO_1])) * (RHO_I_MGM - rho[rho >= RHO_1] / 1000) * np.power(A, bHL) * 1000 / S_PER_YEAR
    drho_dt[rho >= RHO_1] = dr_dt[rho >= RHO_1]

    return drho_dt

def HL_Sigfus(steps, gridLen, bdotSec, Tz, rho, sigma):
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

    A = bdotSec * steps * BDOT_TO_A
    drho_dt = np.zeros(gridLen)
    f550 = interpolate.interp1d(rho, sigma)
    sigma550 = f550(RHO_1)
    rhoDiff = (RHO_I_MGM - rho / 1000)

    k = np.power(k2 * np.exp(-Q2 / (R * Tz[rho >= RHO_1])), 2) / S_PER_YEAR
    sigmaDiff = (sigma[rho >= RHO_1] - sigma550)

    drho_dt[rho < RHO_1] = k1 * np.exp(-Q1 / (R * Tz[rho < RHO_1])) * (RHO_I_MGM - rho[rho < RHO_1] / 1000) * np.power(A, aHL) * 1000 / S_PER_YEAR

    drho_dt[(rho >= RHO_1) & (rho < RHO_I)]  = k * (sigmaDiff * rhoDiff[rho >= RHO_1]) / (GRAVITY * np.log((RHO_I_MGM - RHO_1 / 1000) / (rhoDiff[(rho >= RHO_1) & (rho < RHO_I)])))
    drho_dt[(rho >= RHO_1) & (rho >= RHO_I)] = 0

    return drho_dt

def Li_2004(steps, gridLen, bdotSec, T_mean, rho):
    '''

    :param steps:
    :param gridLen:
    :param bdotSec:
    :param T_mean:
    :param rho:

    :return drho_dt:
    '''

    A = bdotSec * steps * BDOT_TO_A
    drho_dt = np.zeros(gridLen)
    dr_dt = (RHO_I - rho) * A * (139.21 - 0.542 * T_mean) * 8.36 * (K_TO_C - Tz) ** -2.061

    drho_dt = dr_dt / S_PER_YEAR

    return drho_dt

def Li_2011(steps, gridLen, bdotSec, bdot_mean, bdot_type, Tz, T_mean, rho):
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

    A     = bdotSec * steps * BDOT_TO_A
    TmC   = T_mean - K_TO_C
    beta1 = -9.788 + 8.996 * A - 0.6165 * TmC
    beta2 = beta1 / (-2.0178 + 8.4043 * A - 0.0932 * TmC)
    dr_dt = np.zeros(gridLen)

    if bdot_type == 'instant':
        dr_dt[rho <= RHO_1] = (RHO_I - rho[rho <= RHO_1]) * A * beta1 * 8.36 * (K_TO_C - Tz[rho <= RHO_1]) ** -2.061
        dr_dt[rho > RHO_1]  = (RHO_I - rho[rho > RHO_1]) * A * beta2 * 8.36 * (K_TO_C - Tz[rho > RHO_1]) ** -2.061
    elif bdot_type == 'mean':
        dr_dt[rho <= RHO_1] = (RHO_I - rho[rho <= RHO_1]) * bdot_mean[rho <= RHO_1] * steps * BDOT_TO_A * beta1 * 8.36 * (K_TO_C - Tz[rho <= RHO_1]) ** -2.061
        dr_dt[rho > RHO_1]  = (RHO_I - rho[rho > RHO_1]) * bdot_mean[rho > RHO_1] * steps * BDOT_TO_A  * beta2 * 8.36 * (K_TO_C - Tz[rho > RHO_1]) ** -2.061

    drho_dt = dr_dt / S_PER_YEAR

    return drho_dt

def Arthern_2010S(steps, gridLen, bdotSec, bdot_mean, bdot_type, Tz, T_mean, rho):
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

    A_instant = bdotSec * steps * BDOT_TO_A * 1000
    A_mean_1 = bdot_mean[rho < RHO_1] * steps * BDOT_TO_A * 1000
    A_mean_2 = bdot_mean[rho >= RHO_1] * steps * BDOT_TO_A * 1000
    dr_dt = np.zeros(gridLen)

    if bdot_type == 'instant':
        dr_dt[rho < RHO_1]  = (RHO_I - rho[rho < RHO_1]) * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * Tz[rho < RHO_1]) + Eg / (R * T_mean))
        dr_dt[rho >= RHO_1] = (RHO_I - rho[rho >= RHO_1]) * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * Tz[rho >= RHO_1]) + Eg / (R * T_mean))
    elif bdot_type == 'mean':
        dr_dt[rho < RHO_1]  = (RHO_I - rho[rho < RHO_1]) * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * Tz[rho < RHO_1]) + Eg / (R * T_mean))
        dr_dt[rho >= RHO_1] = (RHO_I - rho[rho >= RHO_1]) * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * Tz[rho >= RHO_1]) + Eg / (R * T_mean))

    drho_dt = dr_dt / S_PER_YEAR

    return drho_dt

def Arthern_2010T(gridLen, Tz, rho, sigma, r2, physGrain):
    '''

    :param gridLen:
    :param Tz:
    :param rho:
    :param sigma:
    :param r2:
    :param physGrain:

    :return drho_dt:
    '''

    kc1 = 9.2e-9
    kc2 = 3.7e-9
    Ec  = 60.0e3

    if not physGrain:
        print "Grain growth should be on for Arthern Transient"
        return

    drho_dt = np.zeros(gridLen)

    drho_dt[rho < RHO_1]  = kc1 * (RHO_I - rho[rho < RHO_1]) * np.exp(-Ec / (R * Tz[rho < RHO_1])) * sigma[rho < RHO_1] / (r2[rho < RHO_1])
    drho_dt[rho >= RHO_1] = kc2 * (RHO_I - rho[rho >= RHO_1]) * np.exp(-Ec / (R * Tz[rho >= RHO_1])) * sigma[rho >= RHO_1] / (r2[rho >= RHO_1])

    return drho_dt

def Helsen_2008(steps, bdotSec, bdot_mean, Tz, Ts, rho, bdot_type):
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

    A_instant = bdotSec * steps * BDOT_TO_A
    A_mean = bdot_mean * steps * BDOT_TO_A

    if bdot_type == 'instant':
        dr_dt = (RHO_I - rho) * A_instant * (76.138 - 0.28965 * Ts) * 8.36 * (K_TO_C - Tz) ** -2.061
    elif bdot_type == 'mean':
        dr_dt = (RHO_I - rho) * A_mean * (76.138 - 0.28965 * Ts) * 8.36 * (K_TO_C - Tz) ** -2.061

    drho_dt = dr_dt / S_PER_YEAR

    return drho_dt

def Simonsen_2013(steps, gridLen, bdotSec, bdot_mean, bdot_type, Tz, T_mean, rho):
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

    A_instant = bdotSec * steps * BDOT_TO_A * 1000
    A_mean_1 = bdot_mean[rho < RHO_1] * steps * BDOT_TO_A * 1000
    A_mean_2 = bdot_mean[rho >= RHO_1] * steps * BDOT_TO_A * 1000
    dr_dt = np.zeros(gridLen)

    if bdot_type == 'instant':
        gamma = 61.7 / ((bdotSec * steps * S_PER_YEAR * RHO_I) ** (0.5)) * np.exp(-3800. / (R * T_mean))
        dr_dt[rho < RHO_1]  = F0 * (RHO_I - rho[rho < RHO_1]) * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * Tz[rho < RHO_1]) + Eg / (R * T_mean))
        dr_dt[rho >= RHO_1] = F1 * gamma * (RHO_I - rho[rho >= RHO_1]) * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * Tz[rho >= RHO_1]) + Eg / (R * T_mean))
    elif bdot_type == 'mean':
        gamma = 61.7 / ((bdot_mean[rho >= RHO_1] * steps * S_PER_YEAR * RHO_I) ** (0.5)) * np.exp(-3800.0 / (R * T_mean))
        dr_dt[rho < RHO_1]  = F0 * (RHO_I - rho[rho < RHO_1]) * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * Tz[rho < RHO_1]) + Eg / (R * T_mean))
        dr_dt[rho >= RHO_1] = F1 * gamma * (RHO_I - rho[rho >= RHO_1]) * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * Tz[rho >= RHO_1]) + Eg / (R * T_mean))

    drho_dt = dr_dt / S_PER_YEAR

    return drho_dt

def Ligtenberg_2011(steps, gridLen, bdotSec, bdot_mean, bdot_type, Tz, T_mean, rho):
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

    A_instant = bdotSec * steps * BDOT_TO_A * 1000
    A_mean_1 = bdot_mean[rho < RHO_1] * steps * BDOT_TO_A * 1000
    A_mean_2 = bdot_mean[rho >= RHO_1] * steps * BDOT_TO_A * 1000
    dr_dt = np.zeros(gridLen)

    if bdot_type == 'instant':
        M_0 = 1.435 - 0.151 * np.log(bdotSec * S_PER_YEAR * RHO_I)
        M_1 = 2.366 - 0.293 * np.log(bdotSec * S_PER_YEAR * RHO_I)
        dr_dt[rho < RHO_1]  = (RHO_I - rho[rho < RHO_1]) * M_0 * ar1 * A_instant * GRAVITY * np.exp(-Ec / (R * Tz[rho < RHO_1]) + Eg / (R * T_mean))
        dr_dt[rho >= RHO_1] = (RHO_I - rho[rho >= RHO_1]) * M_1 * ar2 * A_instant * GRAVITY * np.exp(-Ec / (R * Tz[rho >= RHO_1])+ Eg / (R * T_mean))
    elif bdot_type == 'mean':
        M_0 = 1.435 - 0.151 * np.log(bdot_mean[rho < RHO_1] * S_PER_YEAR * RHO_I)
        M_1 = 2.366 - 0.293 * np.log(bdot_mean[rho >= RHO_1] * S_PER_YEAR * RHO_I)
        dr_dt[rho < RHO_1]  = (RHO_I - rho[rho < RHO_1]) * M_0 * ar1 * A_mean_1 * GRAVITY * np.exp(-Ec / (R * Tz[rho < RHO_1]) + Eg / (R * T_mean))
        dr_dt[rho >= RHO_1] = (RHO_I - rho[rho >= RHO_1]) * M_1 * ar2 * A_mean_2 * GRAVITY * np.exp(-Ec / (R * Tz[rho >= RHO_1]) + Eg / (R * T_mean))

    drho_dt = dr_dt / S_PER_YEAR

    return drho_dt

def Barnola_1991(steps, gridLen, bdotSec, Tz, rho, sigma):
    '''

    :param steps:
    :param gridLen:
    :param bdotSec:
    :param Tz:
    :param rho:
    :param sigma:

    :return drho_dt:
    '''
    Q1 				= 10160.0
    k1 				= 11.0
    aHL 			= 1.0
    alphaBarnola 	= -37.455
    betaBarnola 	= 99.743
    deltaBarnola 	= -95.027
    gammaBarnola 	= 30.673
    A0b 			= 2.54e4
    n 				= 3.0
    QBarnola		= 60.0e3
    closeOff        = 800.0

    rho[rho > RHO_I] = RHO_I #The Barnola model will go a fraction over the ice density (order 10^-3), so this stops that.
    drho_dt = np.zeros(gridLen)
    D = rho / RHO_I
    nBa = n * np.ones(gridLen)
    A0 = A0b * np.ones(gridLen) / 1.e18 #this is for the n=3 region.

    sigmaEff = sigma

    # zone 1
    A = bdotSec * steps * BDOT_TO_A
    drho_dt[rho < RHO_1] = dr_dt = k1 * np.exp(-Q1 / (R * Tz[rho < RHO_1])) * (RHO_I_MGM - rho[rho < RHO_1] / 1000) * np.power(A[ii], aHL) * 1000 / S_PER_YEAR

    # zone 2
    fe = 10.0 ** (alphaBarnola * (rho[rho <= RHO_2] / 1000) ** 3. + betaBarnola * (rho[rho <= RHO_2] / 1000) ** 2. + deltaBarnola * rho[rho <= RHO_2] / 1000 + gammaBarnola)
    drho_dt[rho <= RHO_2] = rho[rho <= RHO_2] * A0[rho <= RHO_2] * np.exp(-QBarnola / (R * Tz[rho <= RHO_2])) * fe * (sigmaEff[rho <= RHO_2] ** nBa[rho <= RHO_2])

    # zone 3
    fs = (3. / 16.) * (1 - rho[rho > RHO_2] / RHO_I) / (1 - (1 - rho[rho > RHO_2] / RHO_I) ** (1. / 3.)) ** 3.
    drho_dt[rho > RHO_2] = rho[rho > RHO_2] * A0[rho > RHO_2] * np.exp(-QBarnola / (R * Tz[rho > RHO_2])) * fs * (sigmaEff[rho > RHO_2] ** nBa[rho > RHO_2])

    return drho_dt

def Morris_HL_2013(steps, spin, iter, gridLen, Tz, dt, rho):
    '''

    :param steps:
    :param spin:
    :param iter:
    :param gridLen:
    :param Tz:
    :param dt:
    :param rho:

    :return drho_dt:
    '''
    QMorris = 110e3
    kMorris = 11.0
    rhoW    = 1000.0

    if spin and iter == 0:
        Hx = np.zeros(gridLen)
        Hx = Hx + np.exp(-QMorris / (R * Tz)) * dt

    drho_dt = np.zeros(gridLen)
    drho_dt[rho <= RHO_1] = (kMorris / (rhoW * GRAVITY)) * ((RHO_I - rho[rho <= RHO_1]) / rho[rho <= RHO_1]) * 1 / Hx[rho <= RHO_1] * np.exp(-QMorris / (R * Tz[rho <= RHO_1])) * sigma[rho <= RHO_1]

    #HL Dynamic
    Q2  = 21400.0
    k2  = 575.0
    bHL = 0.5
    A   = bdotSec * steps * BDOT_TO_A

    drho_dt[rho >= RHO_1] = k2 * np.exp(-Q2 / (R * Tz[rho >= RHO_1])) * (RHO_I_MGM - rho[rho >= RHO_1] / 1000) * np.power(A[iter], bHL) * 1000 / S_PER_YEAR

    return drho_dt

def Spencer_2001():
    '''
    :return:
    '''
    pass

def Goujon_2003():
    '''
    :return:
    '''
    pass

def grain_growth(Tz, Ts, iter, dt):
    '''

    :param Tz:
    :param Ts:
    :param iter:
    :param dt:

    :return r2:
    '''

    kgr = 1.3e-7 															# grain growth rate from Arthern (2010)
    Eg  = 42.4e3

    dr2_dt = kgr * np.exp(-Eg / (R * Tz))
    r2 = r2 + dr2_dt * dt
    r2_time = np.concatenate(([t * iter + 1], r2))

    if calcGrainSize:												# uses surface temperature to get an initial grain size
        r2 = np.concatenate(([-2.42e-9 * Ts +9.46e-7], r2[:-1]))
    else:
        r2 = np.concatenate(([(0.1e-3) ** 2], r2[:-1]))

    return r2