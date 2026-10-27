#!/usr/bin/env python
'''
Functions to handle strain effects on firn densification.

Falk Oraschewski, 17.05.2022
'''

import numpy as np
import os

import scipy.interpolate as interpolate

from reader import read_input
from constants import *

#############
def check_strain_settings(self):
    '''
    Check, if strain options are set in the config file. If not, the corresponding modules are deactivated. If 
    deprecated config keys are use, they are updated to the current version.
    
    Falk Oraschewski, 17.05.2022
    
    :param self.c:
    
    :return self.c: 
    '''
    # Check for deprecated config keys.
    if 'strain' in self.c:
        if 'horizontal_divergence' not in self.c:
            print("DEPRECATION WARNING: 'strain' is a deprecated config key. Use 'horizontal_divergence' instead.")
            self.c['horizontal_divergence'] = self.c['strain']
        else:
            print("DEPRECATION WARNING: both 'strain' (deprecated) and 'horizontal_divergence' are present in the config. "
                  "Using 'horizontal_divergence' value; 'strain' is ignored.")

    if 'InputFileNamedudx' in self.c:
        if 'InputFileNameStrain' not in self.c:
            print("DEPRECATION WARNING: 'InputFileNamedudx' is a deprecated config key. Use 'InputFileNameStrain' instead.")
            self.c['InputFileNameStrain'] = self.c['InputFileNamedudx']
        else:
            print("DEPRECATION WARNING: both 'InputFileNamedudx' (deprecated) and 'InputFileNameStrain' are present in the config. "
                  "Using 'InputFileNameStrain' value; 'InputFileNamedudx' is ignored.")

    # Check if strain settings are set in the config. Otherwise deactivate.
    if 'horizontal_divergence' not in self.c:
        self.c['horizontal_divergence'] = False
    if 'strain_softening' not in self.c:
        self.c['strain_softening'] = False
    if 'residual_strain' not in self.c:
        self.c['residual_strain'] = 2e-4
    if 'tuning_bias_correction' not in self.c:
        self.c['tuning_bias_correction'] = False
        
    return self.c

### end check_strain_settings ###
#####################

#############
def _strain_from_mapping(src):
    '''
    Pull horizontal strain rates out of a dictionary in the same format as climateTS: a 'time' array (decimal year) plus
    EITHER 'eps' (1-D: divergence; or 2xN / 3xN: components in rows) OR named arrays: 'eps_divergence' (divergence), or
    'eps_1' and 'eps_2' (principal), or 'eps_xx', 'eps_yy' and 'eps_xy'. Rates in yr^-1, extension positive.

    :param src: dict-like (e.g., climateTS or the strain_rates argument)
    :return (eps, time), or None if src holds no strain data
    '''
    if 'eps' in src:
        eps = np.asarray(src['eps'], dtype=float)
        if eps.ndim == 2 and eps.shape[0] not in (2, 3) and eps.shape[1] in (2, 3):
            eps = eps.T  # accept (N, k) as well as (k, N)
    elif all(k in src for k in ('eps_xx', 'eps_yy', 'eps_xy')):
        eps = np.vstack([np.asarray(src[k], dtype=float) for k in ('eps_xx', 'eps_yy', 'eps_xy')])
    elif all(k in src for k in ('eps_1', 'eps_2')):
        eps = np.vstack([np.asarray(src[k], dtype=float) for k in ('eps_1', 'eps_2')])
    elif 'eps_divergence' in src:
        eps = np.asarray(src['eps_divergence'], dtype=float)
    else:
        return None
    if 'time' not in src:
        raise KeyError("Strain data given as a dictionary needs a 'time' array (decimal years), as in climateTS.")
    time = np.asarray(src['time'], dtype=float)
    if eps.shape[-1] != len(time):
        raise ValueError(f"Strain data has {eps.shape[-1]} values but 'time' has {len(time)}.")
    return eps, time
### end _strain_from_mapping ###
################################

def get_strain_input(self, climateTS=None, strain_rates=None):
    '''
    Find the horizontal strain-rate forcing. Sources, in order of priority:

      1. the strain_rates argument given to FirnDensityNoSpin: a number (constant), a list/array of 1-3 numbers
         (constant), or a dictionary like climateTS (see _strain_from_mapping);
      2. strain keys in climateTS itself (same format; same time axis as the other climate forcing);
      3. the .json key 'strain_rate_const': a number (divergence) or a list of 2 (principal rates) or 3 (eps_xx, eps_yy, eps_xy);
      4. the csv file named by 'InputFileNameStrain' (original behavior).

    Number of values means: 1 = horizontal divergence rate, 2 = principal horizontal strain rates, 3 = eps_xx, eps_yy,
    eps_xy. All in yr^-1, extension positive (so divergence = -(eps_xx + eps_yy)).

    :return eps: constant case: 1-D array of length 1, 2 or 3. Time-series case: 1-D (divergence) or 2-D (k x N)
    :return year: time axis of eps (decimal years), or None for a constant
    :return source: description of where the data came from, for the log
    '''
    if strain_rates is not None:
        if isinstance(strain_rates, dict):
            got = _strain_from_mapping(strain_rates)
            if got is None:
                raise KeyError("strain_rates dictionary needs 'eps', 'eps_divergence', ('eps_1','eps_2') or ('eps_xx','eps_yy','eps_xy').")
            return got[0], got[1], 'strain_rates argument (dictionary)'
        return np.atleast_1d(np.asarray(strain_rates, dtype=float)), None, 'strain_rates argument (constant)'

    if climateTS is not None:
        got = _strain_from_mapping(climateTS)
        if got is not None:
            return got[0], got[1], 'climateTS'

    if self.c.get('strain_rate_const') is not None:
        return np.atleast_1d(np.asarray(self.c['strain_rate_const'], dtype=float)), None, "json key 'strain_rate_const'"

    if 'InputFileNameStrain' not in self.c:
        raise KeyError("horizontal_divergence/strain_softening is on but no strain forcing was found. Provide one of: the "
                       "strain_rates argument, strain keys in climateTS, 'strain_rate_const', or 'InputFileNameStrain' in the .json.")
    eps, year, _, _ = read_input(os.path.join(self.c['InputFileFolder'], self.c['InputFileNameStrain']))
    return eps, year, f"csv file {self.c['InputFileNameStrain']}"
### end get_strain_input ###
############################

def load_strain(self, spin=False, strain_rates=None, climateTS=None):
    '''
    Load strain rate inputs for spin-up and main run (source: see get_strain_input). Depending on the number of values,
    the input is taken as (1-valued) the divergence rate, (2-valued) the principal horizontal strain rates or (3-valued)
    the horizontal strain rate components eps_xx, eps_yy and eps_xy.

    Sets eps_eff_hor_2 and eps_divergence (used by strain_softening and horizontal_divergence) and the components
    eps_xx, eps_yy, eps_xy (given to densification models that return 'eps_zz'; see physics.get_drho_dt), each with
    one value per time step (spin: constant, from the mean of the input or, if spinup_climate_type is 'initial',
    its first value).

    Falk Oraschewski, 17.05.2022; generalized input sources 26/10.

    :param self:
    :param spin: True when called from the spin-up
    :param strain_rates: optional constant or dictionary of strain rates, see get_strain_input
    :param climateTS: optional climateTS dictionary, may contain strain keys, see get_strain_input

    :return self.eps_eff_hor_2, self.eps_divergence, self.c:
    '''

    int_type = self.c['int_type']
    eps, year, source = get_strain_input(self, climateTS=climateTS, strain_rates=strain_rates)
    print(f'Horizontal strain rates from {source}.')

    eps = np.asarray(eps, dtype=float)
    if year is None:    # constant: one value per component
        comps = [np.array([v]) for v in np.atleast_1d(eps)]
    elif eps.ndim <= 1:  # time series of the divergence
        comps = [eps]
    else:               # time series of 2 or 3 components, in rows
        comps = [eps[i, :] for i in range(eps.shape[0])]
    k = len(comps)
    n = self.stp if spin else len(self.modeltime)  # number of time steps

    series = []
    for comp in comps:
        if year is None:  # constant in time
            series.append(float(comp[0]) * np.ones(n))
        elif spin:
            series.append((comp[0] if self.c['spinup_climate_type'] == 'initial' else np.mean(comp)) * np.ones(n))
        else:
            series.append(interpolate.interp1d(year, comp, int_type, fill_value='extrapolate')(self.modeltime))

    if k == 1:
        if not spin:
            print('Single valued strain rate input is taken as horizontal divergence rates. Strain softening is deactivated.')
        self.eps_divergence = series[0]
        self.c['strain_softening'] = False  # Deactivate strain softening
        self.eps_eff_hor_2 = np.zeros_like(self.eps_divergence)  # Create empty effective horizontal strain rate.
        # Components are not known from the divergence alone; assume isotropic horizontal strain (eps_xx = eps_yy, eps_xy = 0).
        self.eps_xx, self.eps_yy, self.eps_xy = - self.eps_divergence / 2, - self.eps_divergence / 2, np.zeros_like(self.eps_divergence)
    elif k == 2:
        eps_1, eps_2 = series
        self.eps_eff_hor_2 = (eps_1 ** 2 + eps_2 ** 2) / 2
        self.eps_divergence = - (eps_1 + eps_2)
        self.eps_xx, self.eps_yy, self.eps_xy = eps_1, eps_2, np.zeros_like(eps_1)  # principal axes taken as x, y
    elif k == 3:
        eps_xx, eps_yy, eps_xy = series
        self.eps_eff_hor_2 = (eps_xx ** 2 + eps_yy ** 2 + 2 * eps_xy ** 2) / 2
        self.eps_divergence = - (eps_xx + eps_yy)
        self.eps_xx, self.eps_yy, self.eps_xy = eps_xx, eps_yy, eps_xy
    else:
        raise ValueError(f'Strain input must have 1, 2 or 3 components; got {k}.')
    return self.eps_eff_hor_2, self.eps_divergence, self.c
### end load_strain ###
#######################

#############
def strain_softening(self, drho_dt, iii):
    '''
    Process strain softening in firn stage 2 as described in (Oraschewski & Grinsted, 2022) by scaling the densification
    rate. A tuning bias correction can be applied, to correct for the impliciet contribution of strain softening 
    captured by the HL model.
    
    Falk Oraschewski, 17.05.2022
    
    :param self: 
    :param drho_dt: 
    :param iii: 
    
    :returns drho_dt, self.viscosity: 
    '''
    
    eps_zz_classic = - drho_dt / self.rho * S_PER_YEAR  # Vertical strain rate from original densification rate.
    eps_zz_classic = eps_zz_classic - np.abs(self.c['residual_strain'])  # Regularisation, eq. 21 in (Oraschewski & Grinsted, 2022)
    eps_eff_classic_2 = 0.5 * (eps_zz_classic ** 2)  # Effective strain rate of original model
    rH2 = self.eps_eff_hor_2[iii] / eps_eff_classic_2  # Square of eq. 19 in (Oraschewski & Grinsted, 2022)
    
    # Solution for n=3
    # K = (13.5 * rH2 + 1.5 * np.sqrt(81 * rH2 ** 2 + 12 * rH2) + 1) ** (1 / 3)
    # scale_factor_softening = (K + 1 / K + 1) / 3
    
    # Solution for n=4
    scale_factor_softening = np.ones_like(rH2)
    rH2_nz = rH2[np.nonzero(rH2)]  # Avoid division by 0.
    
    # Equation 20 in (Oraschewski & Grinsted, 2022):
    a1 = np.cbrt(9 * rH2_nz ** 4 + np.sqrt(81 * rH2_nz ** 8 + 768 * rH2_nz ** 9))
    a2 = np.sqrt(1 + 8 * rH2_nz + np.cbrt(32/9) * a1 - np.cbrt(8192/3) * rH2_nz ** 3 / a1)
    a3 = np.sqrt(1/2 + 4 * rH2_nz - a1 / np.cbrt(18) + np.cbrt(128/3) * rH2_nz ** 3 / a1
                 + (1 + 12 * rH2_nz + 24 * rH2_nz ** 2) / (2 * a2))
    scale_factor_softening[np.nonzero(rH2)] = np.sqrt(1/4 + 1/4 * a2 + 1/2 * a3)
    
    if self.c['tuning_bias_correction']:
    # Repeat computation with the strain rate that corresponds to the tuning bias of the HL model.
        eps_eff_cor = 4.5e-4  # Tuning bias for the HL model. Values for other models need to be determined.
        eps_eff_cor_2 = eps_eff_cor ** 2
    
        rH2 = eps_eff_cor_2 / eps_eff_classic_2
    
        # Solution for n=4
        scale_factor_softening_cor = np.ones_like(rH2)
        rH2_nz = rH2[np.nonzero(rH2)]
    
        a1 = np.cbrt(9 * rH2_nz ** 4 + np.sqrt(81 * rH2_nz ** 8 + 768 * rH2_nz ** 9))
        a2 = np.sqrt(
            1 + 8 * rH2_nz + np.cbrt(32 / 9) * a1 - np.cbrt(8192 / 3) * rH2_nz ** 3 / a1)
        a3 = np.sqrt(1 / 2 + 4 * rH2_nz - a1 / np.cbrt(18) + np.cbrt(128 / 3) * rH2_nz ** 3 / a1
                     + (1 + 12 * rH2_nz + 24 * rH2_nz ** 2) / (2 * a2))
        scale_factor_softening_cor[np.nonzero(rH2)] = np.sqrt(1 / 4 + 1 / 4 * a2 + 1 / 2 * a3)
    
        # Corrected scale factor as in Eq. 23 in (Oraschewski & Grinsted, 2022)
        scale_factor_softening = scale_factor_softening / scale_factor_softening_cor
    
    # Testing effect on a assumed linearly increasing contribution of power-law creep in stage 1:
    # z1mask = (self.rho < RHO_1)
    # stage1_scale = np.linspace(0, 1, sum(z1mask))
    # drho_dt[z1mask] = drho_dt[z1mask] * (1 + (viscosity_scale[z1mask]-1) * stage1_scale)
    
    z2mask = (self.rho >= RHO_1)  # Only apply strain softening in second firn stage, where power-law creep is dominant.
    drho_dt[z2mask]        = drho_dt[z2mask] * scale_factor_softening[z2mask]
    self.viscosity[z2mask] = self.viscosity[z2mask] / scale_factor_softening[z2mask]
    
    return drho_dt, self.viscosity
#############

#############
def horizontal_strain_rates(self, iii):
    '''
    Horizontal strain rates at time step iii, for densification models whose constitutive law depends on them
    (those returning 'eps_zz'; see physics.get_drho_dt). Indexed as in horizontal_divergence. Zero when
    horizontal_divergence is off, consistent with the mass rescaling being off.

    :param self:
    :param iii: time step index
    :return eps_xx, eps_yy, eps_xy: horizontal strain rates [s^-1], extension positive
    '''
    if self.c['horizontal_divergence']:
        return (self.eps_xx[iii] / S_PER_YEAR, self.eps_yy[iii] / S_PER_YEAR, self.eps_xy[iii] / S_PER_YEAR)
    return 0.0, 0.0, 0.0
### end horizontal_strain_rates ###
###################################

#############    
def horizontal_divergence(self,iii):
    '''
    Rescale mass per layer according to the horizontal divergence rate (Horlings et al., 2020)
    
    :param self: 
    :param iii: 
    :return self.mass: 
    '''
    
    self.mass = (1 + self.eps_divergence[iii] / S_PER_YEAR * self.dt[iii]) * self.mass
    
    return self.mass
