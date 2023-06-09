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
        self.c['horizontal_divergence'] = self.c['strain']
    if 'InputFileNamedudx' in self.c:
        self.c['InputFileNameStrain'] = self.c['InputFileNamedudx']

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
#############

#############
def load_strain(self, spin=False):
    '''
    Load strain rate inputs for spin-up and main run. Depending on the number of data rows in the strain input file, the
    input is taken as (1-valued) the divergence rate, (2-valued) the principal horizontal strain rates or (3-valued) the
    horizontal strain rate components eps_xx, eps_yy and eps_xy. 
    
    Falk Oraschewski, 17.05.2022
    
    :param self: 
    :param spin: 
    
    :return self.eps_eff_hor_2, self.eps_divergence, self.c: 
    '''
    
    int_type = self.c['int_type']
    
    input_eps, input_year_eps, input_eps_full, input_year_eps_full = read_input(os.path.join(self.c['InputFileFolder'], self.c['InputFileNameStrain']))
    if np.ndim(input_eps) == 1:
        if spin:
            if self.c['spinup_climate_type'] == 'initial':
                self.eps_divergence = input_eps[0] * np.ones(self.gridLen)
            else:
                self.eps_divergence = np.mean(input_eps) * np.ones(self.gridLen) 
        else:
            print('Single valued strain rate input is taken as horizontal divergence rates. Strain softening is deactivated.')
            dusf = interpolate.interp1d(input_year_eps, input_eps, int_type, fill_value='extrapolate')
            self.eps_divergence = dusf(self.modeltime)
        self.c['strain_softening'] = False  # Deactivate strain softening
        self.eps_eff_hor_2 = np.zeros_like(self.eps_divergence)  # Create empty effective horizontal strain rate.
        
    elif np.min(np.shape(input_eps)) == 2:
        input_eps_1, input_eps_2 = input_eps[0, :], input_eps[1, :]
        if spin:
            if self.c['spinup_climate_type'] == 'initial':
                eps_1 = input_eps_1[0] * np.ones(self.gridLen)
                eps_2 = input_eps_2[0] * np.ones(self.gridLen)
            else:
                print('2-valued strain rate input is taken as principal horizontal strain rates.')
                eps_1 = np.mean(input_eps_1) * np.ones(self.gridLen)
                eps_2 = np.mean(input_eps_2) * np.ones(self.gridLen)
        else:
            d1sf               = interpolate.interp1d(input_year_eps, input_eps_1, int_type, fill_value='extrapolate')
            d2sf               = interpolate.interp1d(input_year_eps, input_eps_2, int_type, fill_value='extrapolate')
            eps_1         = d1sf(self.modeltime)
            eps_2         = d2sf(self.modeltime)
        self.eps_eff_hor_2 = (eps_1 ** 2 + eps_2 ** 2) / 2
        self.eps_divergence = - (eps_1 + eps_2)
        
    elif np.min(np.shape(input_eps)) == 3:
        input_eps_xx, input_eps_yy, input_eps_xy = input_eps[0, :], input_eps[1, :], input_eps[2, :]
        if spin:
            if self.c['spinup_climate_type'] == 'initial':
                eps_xx = input_eps_xx[0] * np.ones(self.gridLen)
                eps_yy = input_eps_yy[0] * np.ones(self.gridLen)
                eps_xy = input_eps_xy[0] * np.ones(self.gridLen)
            else:
                print('3-valued strain rate input is taken as the horizontal components eps_xx, eps_yy and eps_xy in that order.')
                eps_xx = np.mean(input_eps_xx) * np.ones(self.gridLen)
                eps_yy = np.mean(input_eps_yy) * np.ones(self.gridLen)
                eps_xy = np.mean(input_eps_xy) * np.ones(self.gridLen)
        else:
            dxxsf = interpolate.interp1d(input_year_eps, input_eps_xx, int_type, fill_value='extrapolate')
            dyysf = interpolate.interp1d(input_year_eps, input_eps_yy, int_type, fill_value='extrapolate')
            dxysf = interpolate.interp1d(input_year_eps, input_eps_xy, int_type, fill_value='extrapolate')
            eps_xx = dxxsf(self.modeltime)
            eps_yy = dyysf(self.modeltime)
            eps_xy = dxysf(self.modeltime)
        self.eps_eff_hor_2 = (eps_xx ** 2 + eps_yy ** 2 + 2 * eps_xy ** 2) / 2
        self.eps_divergence = - (eps_xx + eps_yy)
    return self.eps_eff_hor_2, self.eps_divergence, self.c
#############

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
def horizontal_divergence(self,iii):
    '''
    Rescale mass per layer according to the horizontal divergence rate (Horlings et al., 2020)
    
    :param self: 
    :param iii: 
    :return self.mass: 
    '''
    
    self.mass = (1 + self.eps_divergence[iii] / S_PER_YEAR * self.dt[iii]) * self.mass
    
    return self.mass
