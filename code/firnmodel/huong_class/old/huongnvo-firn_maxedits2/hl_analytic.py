from constants import *
import numpy as np

def hl_analytic(rhos0, h, THL, AHL):
    '''
    Model steady-state firn density and age profiles and bubble close-off, uses m w.e. a^-1

    :param rhos0:
    :param h:
    :param THL:
    :param AHL:

    :return age: age vector of firn column with steady-state dynamics
    :return rho: density vector of firn column with steady state dynamics
    '''

    hSize = np.size(h)
    rhos = rhos0 / 1000.0

    A = AHL * RHO_I_MGM
    k0 = 11.0 * np.exp(-10160 / (R * THL))
    k1 = 575.0 * np.exp(-21400 / (R * THL))

    # depth of critical density, eqn 8 from Herron and Langway
    h0_55 = 1 / (RHO_I_MGM * k0) * (np.log(RHO_1_MGM / (RHO_I_MGM - RHO_1_MGM)) - np.log(rhos / (RHO_I_MGM - rhos)))
    Z0 = np.exp(RHO_I_MGM * k0 * h + np.log(rhos / (RHO_I_MGM - rhos)))

    # The boundary from zone 1 to zone 2 = t0_55
    t0_55 = 1 / (k0 * A) * np.log((RHO_I_MGM - rhos) / (RHO_I_MGM - RHO_1_MGM))
    rho_h0 = (RHO_I_MGM * Z0) / (1 + Z0)
    if np.max(rho_h0) >= RHO_I_MGM:
        t0 = np.zeros(hSize)
        for jj in xrange(hSize):
            if rho_h0[jj] <= RHO_I_MGM - 0.001:
                t0[jj] = (1 / (k0 * A) * np.log((RHO_I_MGM - rhos) / (RHO_I_MGM - rho_h0[jj])))
                jj_max = jj
            else:
                t0[jj] = (t0[jj_max])
    else:
        t0 = 1 / (k0 * A) * np.log((RHO_I_MGM - rhos) / (RHO_I_MGM - rho_h0))

    Z1 = np.exp(RHO_I_MGM * k1 * (h - h0_55) / np.sqrt(A) + np.log(RHO_1_MGM / (RHO_I_MGM - RHO_1_MGM)))
    Z = np.concatenate((Z0[h < h0_55], Z1[h > h0_55]))
    rho_h = (RHO_I_MGM * Z) / (1 + Z)
    tp = np.ones(hSize)
    for j in xrange(hSize):
        if rho_h[j] < RHO_I_MGM - 0.01:
            tp[j] = 1 / (k1 * np.sqrt(A)) * np.log((RHO_I_MGM - RHO_1_MGM) / (RHO_I_MGM - rho_h[j])) + t0_55
            jMax = j
        else:
            tp[j] = tp[jMax]

    # Zone 1 and Zone 2 repsectively
    age = np.concatenate((t0[h < h0_55], tp[h > h0_55])) * S_PER_YEAR
    rho = rho_h * 1000

    return age, rho
