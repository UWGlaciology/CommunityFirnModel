#!/usr/bin/env python

'''
SEB_turbulent.py
'''

import numpy as np 
# import scipy.interpolate as interpolate
from scipy import optimize
from constants import *
# import os
import sys
import math, cmath

FORCING_KEYS = {
    'wind_speed'        : 'WS',
    'specific_humidity' : 'QV2m',
    'surface_pressure'  : 'PS',
    'roughness_m'       : 'z0m',
    'roughness_h'       : 'z0h',
    'roughness_q'       : 'z0q',
}

def q_sat_ice(Ts, PS):
    '''
    Saturation specific humidity over ice, via a Clausius-Clapeyron-
    type fit commonly used in glacier SEB models (a Buck-type
    enhancement over ice). Ts [K], PS [Pa]. Returns q_sat [kg/kg].

    VERIFY: enhancement-factor coefficients (22.452, 0.61) before
    relying on this for science - multiple slightly different
    versions exist in the literature.
    '''
    Tc = Ts - 273.15
    e_sat = 611.15 * np.exp(22.452 * Tc / (Ts - 0.61))  # [Pa]
    q_sat = 0.622 * e_sat / (PS - 0.378 * e_sat)
    return q_sat
### end q_sat_ice ###
######################

def latent_heat_turbulent(Ts):
    '''
    Selects sublimation latent heat for a dry, sub-freezing surface,
    or vaporization latent heat if the surface is at the melting
    point (implying liquid water is present). LV_LIQUID and
    LS_SUBLIMATION should be added to constants.py.
    '''
    if Ts >= T_MELT:
        return LV_LIQUID
    else:
        return LS_SUBLIMATION
### end latent_heat_turbulent ###
##################################

def andreas_roughness(z0m, ustar, nu=1.461e-5):
    '''
    Scalar roughness lengths (z0h, z0q) parameterized from momentum
    roughness z0m and friction velocity, following the piecewise
    regression structure of Andreas (1987) in terms of the roughness
    Reynolds number Re* = ustar*z0m/nu.

    VERIFY: the regression coefficients below are my recollection of
    Andreas (1987) Table 1 and should be checked against the original
    paper (or a trusted secondary source, e.g. Van As et al. 2005)
    before production use - I would not trust these numbers blindly.
    '''
    Re_star = ustar * z0m / nu

    if Re_star <= 0.135:
        b0h, b1h, b2h = 1.250, 0.0, 0.0
        b0q, b1q, b2q = 1.610, 0.0, 0.0
    elif Re_star <= 2.5:
        b0h, b1h, b2h = 0.149, -0.550, 0.0
        b0q, b1q, b2q = 0.351, -0.628, 0.0
    else:
        b0h, b1h, b2h = 0.317, -0.565, -0.183
        b0q, b1q, b2q = 0.396, -0.512, -0.180

    lnRe = np.log(Re_star)
    ln_z0h_z0m = b0h + b1h*lnRe + b2h*lnRe**2
    ln_z0q_z0m = b0q + b1q*lnRe + b2q*lnRe**2

    z0h = z0m * np.exp(-ln_z0h_z0m)
    z0q = z0m * np.exp(-ln_z0q_z0m)

    return z0h, z0q
### end andreas_roughness ###
##############################

class RoughnessLengths:
    '''
    Resolves z0m, z0h, z0q from (in priority order):
      1. forcing-data-supplied values, if present in climateTS
      2. Andreas (1987) parameterization, if SEB_roughness_mode is
         'parameterized'
      3. fixed config/default values otherwise
    '''
    def __init__(self, config, z0m_forcing=None, z0h_forcing=None, z0q_forcing=None):
        self.mode = config.get('SEB_roughness_mode', 'fixed')
        self.z0m_default = config.get('SEB_z0m', 1.0e-3)
        self.z0h_default = config.get('SEB_z0h', self.z0m_default)
        self.z0q_default = config.get('SEB_z0q', self.z0m_default)
        self.z0m_forcing = z0m_forcing
        self.z0h_forcing = z0h_forcing
        self.z0q_forcing = z0q_forcing

    def get_z0m(self, iii):
        if self.z0m_forcing is not None:
            return self.z0m_forcing[iii]
        return self.z0m_default

    def get_z0h_z0q(self, iii, ustar, z0m):
        if (self.z0h_forcing is not None) and (self.z0q_forcing is not None):
            return self.z0h_forcing[iii], self.z0q_forcing[iii]
        if self.mode == 'parameterized':
            return andreas_roughness(z0m, ustar)
        return self.z0h_default, self.z0q_default
### end RoughnessLengths.get_z0h_z0q ###
#########################################

class TurbulentFluxModel:
    '''
    Base class for turbulent flux (QH, QL) calculation strategies.
    Subclasses implement compute(). depends_on_Ts tells the solver
    whether the fast analytic quartic path is valid (False) or
    whether a general iterative solver is required (True).
    '''
    depends_on_Ts = True

    def compute(self, Ts, iii):
        raise NotImplementedError
### end TurbulentFluxModel ###
###############################

class PrescribedTurbulentFlux(TurbulentFluxModel):
    '''
    Current CFM default: QH, QL taken directly from the forcing
    product, independent of the CFM's own Ts.
    '''
    depends_on_Ts = False

    def __init__(self, QH_series, QL_series):
        self.QH_series = QH_series
        self.QL_series = QL_series

    def compute(self, Ts, iii):
        return self.QH_series[iii], self.QL_series[iii]
    ### end PrescribedTurbulentFlux.compute ###
    ############################################


class BulkRichardsonFlux(TurbulentFluxModel):
    '''
    Bulk-aerodynamic turbulent fluxes, bulk Richardson number
    stability correction (Louis 1979-style, as adapted for glacier/
    ice-sheet SEB models e.g. Munro 1989; Braithwaite 1995).

    Arrays passed in are already sliced to match the SEB object's
    start_ind - slicing/FORCING_KEYS lookup happens once, centrally,
    in SurfaceEnergyBudget.__init__.
    '''
    depends_on_Ts = True

    def __init__(self, config, WS, T_a, q_a, PS, roughness_model):
        self.WS  = WS
        self.T_a = T_a
        self.q_a = q_a
        self.PS  = PS
        self.z_ref_wind = config.get('SEB_z_ref_wind', 2.0)
        self.z_ref_T    = config.get('SEB_z_ref_T', 2.0)
        self.z_ref_q    = config.get('SEB_z_ref_q', 2.0)
        self.rough = roughness_model

    def compute(self, Ts, iii):
        k = VON_KARMAN
        g = GRAVITY

        WS_i = np.maximum(self.WS[iii], 0.1)  # avoid singularity at calm winds
        T_a  = self.T_a[iii]
        q_a  = self.q_a[iii]
        PS_i = self.PS[iii]

        rho = PS_i / (R_DRY * T_a)
        q_s = q_sat_ice(Ts, PS_i)

        z0m = self.rough.get_z0m(iii)
        ustar_neutral = k * WS_i / np.log(self.z_ref_wind / z0m)
        z0h, z0q = self.rough.get_z0h_z0q(iii, ustar_neutral, z0m)

        Chn = k**2 / (np.log(self.z_ref_wind/z0m) * np.log(self.z_ref_T/z0h))
        Cqn = k**2 / (np.log(self.z_ref_wind/z0m) * np.log(self.z_ref_q/z0q))

        Tv_a = T_a * (1 + 0.61*q_a)
        Rib = (g * self.z_ref_wind * (T_a - Ts)) / (Tv_a * WS_i**2)

        # VERIFY: Louis (1979) stability-function coefficients (b=5
        # here) and unstable-branch form against source literature
        # before production use.
        b = 5.0
        if Rib >= 0:
            Fh = 1.0 / (1.0 + b*Rib)**2
        else:
            Fh = 1.0 - b*Rib / (1.0 + b*np.abs(Rib))

        Ch = Chn * Fh
        Cq = Cqn * Fh

        QH = rho * CP_AIR * Ch * WS_i * (T_a - Ts)
        Lheat = latent_heat_turbulent(Ts)
        QL = rho * Lheat * Cq * WS_i * (q_a - q_s)

        return QH, QL
    ### end BulkRichardsonFlux.compute ###
    #######################################


class MoninObukhovFlux(TurbulentFluxModel):
    '''
    Bulk-aerodynamic turbulent fluxes via Monin-Obukhov similarity,
    iterating on the Obukhov length L. Unstable stability functions
    follow Paulson (1970); stable functions follow Beljaars &
    Holtslag (1991).

    Arrays passed in are already sliced to match the SEB object's
    start_ind, same convention as BulkRichardsonFlux.
    '''
    depends_on_Ts = True

    def __init__(self, config, WS, T_a, q_a, PS, roughness_model, max_iter=20, tol=1e-3):
        self.WS  = WS
        self.T_a = T_a
        self.q_a = q_a
        self.PS  = PS
        self.z_ref_wind = config.get('SEB_z_ref_wind', 2.0)
        self.z_ref_T    = config.get('SEB_z_ref_T', 2.0)
        self.z_ref_q    = config.get('SEB_z_ref_q', 2.0)
        self.rough = roughness_model
        self.max_iter = max_iter
        self.tol = tol

    @staticmethod
    def _psi_m(zeta):
        '''Integrated stability correction, momentum. VERIFY coefficients.'''
        if zeta >= 0:
            a, b, c, d = 1.0, 0.667, 5.0, 0.35
            return -(a*zeta + b*(zeta - c/d)*np.exp(-d*zeta) + b*c/d)
        else:
            x = (1 - 16*zeta)**0.25
            return (2*np.log((1+x)/2) + np.log((1+x**2)/2)
                    - 2*np.arctan(x) + np.pi/2)
    ### end MoninObukhovFlux._psi_m ###
    ####################################

    @staticmethod
    def _psi_h(zeta):
        '''Integrated stability correction, scalars. VERIFY coefficients.'''
        if zeta >= 0:
            a, b, c, d = 1.0, 0.667, 5.0, 0.35
            return -((1 + 2*a*zeta/3)**1.5 + b*(zeta - c/d)*np.exp(-d*zeta) + b*c/d - 1)
        else:
            x = (1 - 16*zeta)**0.5
            return 2*np.log((1+x)/2)
    ### end MoninObukhovFlux._psi_h ###
    ####################################

    def compute(self, Ts, iii):
        k = VON_KARMAN
        g = GRAVITY

        WS_i = np.maximum(self.WS[iii], 0.1)
        T_a  = self.T_a[iii]
        q_a  = self.q_a[iii]
        PS_i = self.PS[iii]

        rho = PS_i / (R_DRY * T_a)
        q_s = q_sat_ice(Ts, PS_i)
        Tv_a = T_a * (1 + 0.61*q_a)

        z0m = self.rough.get_z0m(iii)
        ustar = k * WS_i / np.log(self.z_ref_wind / z0m)
        z0h, z0q = self.rough.get_z0h_z0q(iii, ustar, z0m)

        L = 1.0e6  # start near-neutral
        thetastar = 0.0
        qstar = 0.0

        for _ in range(self.max_iter):
            zeta_wind = self.z_ref_wind / L
            zeta_T    = self.z_ref_T / L
            zeta_q    = self.z_ref_q / L
            zeta_0m   = z0m / L
            zeta_0h   = z0h / L
            zeta_0q   = z0q / L

            ustar_new = k*WS_i / (np.log(self.z_ref_wind/z0m)
                                   - self._psi_m(zeta_wind) + self._psi_m(zeta_0m))
            thetastar = k*(T_a - Ts) / (np.log(self.z_ref_T/z0h)
                                         - self._psi_h(zeta_T) + self._psi_h(zeta_0h))
            qstar = k*(q_a - q_s) / (np.log(self.z_ref_q/z0q)
                                      - self._psi_h(zeta_q) + self._psi_h(zeta_0q))

            L_new = (ustar_new**2 * Tv_a / (k*g*thetastar)) if thetastar != 0 else 1.0e6

            if np.abs(L_new - L) < self.tol:
                L = L_new
                ustar = ustar_new
                break

            L = L_new
            ustar = ustar_new
            z0h, z0q = self.rough.get_z0h_z0q(iii, ustar, z0m)

        QH = rho * CP_AIR * ustar * thetastar
        Lheat = latent_heat_turbulent(Ts)
        QL = rho * Lheat * ustar * qstar

        return QH, QL
    ### end MoninObukhovFlux.compute ###
    #####################################