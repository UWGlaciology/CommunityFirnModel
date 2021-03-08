#!/usr/bin/env python
'''
Code for isotope diffusion.
'''
import numpy as np 
# from solver import solver
from solver import transient_solve_TR
# from Gasses import gasses 
# from Sites import sites 
from reader import read_input, read_init
import json
import scipy.interpolate as interpolate
from constants import *
import os
import sys

class isotopeDiffusion:

    '''
    Isotope diffusion class.
    '''

    def __init__(self,spin,config,isotope,stp,z,modeltime=None):

        '''
        Initialize Isotope diffusion class.
        '''
        self.c = config
        self.isotope = isotope

        try:
            fn = os.path.splitext(self.c['InputFileNameIso'])
            isofile = fn[0] + '_{}'.format(self.isotope) + fn[1]
            print(isofile)
            if isotope=='NoDiffusion':
                isofile = fn[0] + '_dD' + fn[1]
            input_iso, input_year_iso = read_input(os.path.join(self.c['InputFileFolder'],isofile))

            if spin:
                if self.c['spinup_climate_type']=='initial':
                        del_s0  = input_iso[0]
                elif self.c['spinup_climate_type']=='mean':
                        del_s0  = np.mean(input_iso)

                self.del_s  = del_s0 * np.ones(stp)
                self.del_z  = del_s0 * np.ones_like(z)
                self.iso_sig2_s = 0 * np.zeros(stp)
                self.iso_sig2_z = 0 * np.ones_like(z) 

            elif not spin:
                del_z_init          = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'IsoSpin_{}'.format(self.isotope))
                iso_sig2_init        = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'iso_sig2_{}'.format(self.isotope))
                self.del_z = del_z_init[1:]
                self.iso_sig2_z = iso_sig2_init[1:]
                Isf  = interpolate.interp1d(input_year_iso,input_iso,self.c['int_type'],fill_value='extrapolate') # interpolation function
                self.del_s = Isf(modeltime) # isotopes interpolated to modeltime
                self.iso_sig2_s = np.zeros(stp)

        except:
            print('No external file for surface isotope values found ({}), but you specified in the config file that isotope diffusion is on. The model will generate its own synthetic isotope data for you.'.format(self.isotope))
            print('Double check that file name is correct. New module will add d18O or dD to filename for input.')

            if spin:
                del_s0  = -50.0
                print('Currently this is -50 per mil, regardless of isotope you choose.')
                self.del_s  = del_s0 * np.ones(stp)
                self.del_z  = del_s0 * np.ones_like(z)
                self.iso_sig2_s = 0 * np.zeros(stp)
                self.iso_sig2_z = 0 * np.ones_like(z)

            elif not spin:
                del_z_init         = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'IsoSpin_{}'.format(self.isotope))
                iso_sig2_init        = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'iso_sig2_{}'.format(self.isotope))
                self.del_z = del_z_init[1:]
                self.iso_sig2_z = iso_sig2_init[1:]
                ar1             = 0.9   # red noise memory coefficient
                std_rednoise    = 2    # red noise standard deviation
                self.del_s      = std_rednoise*np.random.randn(stp) # white noise
                for x in range(1,stp):
                    self.del_s[x] = self.del_s[x-1]*ar1 + np.random.randn()  # create red noise from white
                self.del_s      = self.del_s - 50
                self.iso_sig2_s = np.zeros(stp)

        if 'site_pressure' not in self.c:
            print('site_pressure is not in .json; defaulting to 1013.25')
            self.c['site_pressure'] = 1013.25


    def isoDiff(self,IsoParams,iii):
        '''
        Isotope diffusion function

        :param iter:
        :param z:
        :param dz:
        :param rho:
        :param iso:

        :returns self.phi_t:
        '''

        for k,v in list(IsoParams.items()):
            setattr(self,k,v)
        
        nz_P        = len(self.z)   # number of nodes in z
        nz_fv       = nz_P - 2      # number of finite volumes in z
        nt          = 1             # number of time steps

        # z_edges_vec = self.z[1:-2] + self.dz[2:-1] / 2                        # uniform edge spacing of volume edges
        # z_edges_vec = np.concatenate(([self.z[0]], z_edges_vec, [self.z[-1]]))
        # z_P_vec   = self.z

        z_edges_vec1 = self.z[0:-1] + np.diff(self.z) / 2
        z_edges_vec = np.concatenate(([self.z[0]], z_edges_vec1, [self.z[-1]]))
        z_P_vec     = self.z

        ### Node positions
        phi_s       = self.del_z[0] # isotope value at surface
        # phi_s       = self.del_z[1] # isotope value at surface
        # if iii==0:
        #     print('Caution! line 121, isotopeDiffusion.py')
        phi_0       = self.del_z # initial isotope profile

        ### Define diffusivity for each isotopic species
        ### Establish values needed for diffusivity calculation
        m           = 0.018 # kg/mol; molar mass of water
        pz          = 3.454e12 * np.exp(-6133 / self.Tz) # Pa; saturation vapor pressure over ice

        alpha_18_z  = 0.9722 * np.exp(11.839 / self.Tz) # fractionation factor for 18_O
        # alpha_18_z    = np.exp(11.839/Tz-28.224*np.power(10,-3))          # alternate formulation from Eric's python code
        alpha_D_z   = 0.9098 * np.exp(16288 / (self.Tz**2)) # fractionation factor for D
        Po          = 1.0 # reference pressure in atm
        # P         = 1.0
        P           = self.c['site_pressure']/1013.25 # Pressure at WAIS from Christo's thesis.

        ### Set diffusivity in air (units of m^2/s)
        Da          = 2.11e-5 * (self.Tz / 273.15)**1.94 * (Po / P)

        ### Calculate tortuosity
        invtau      = np.zeros(int(len(self.dz)))
        b           = 1.3 # Tortuosity parameter (Johnsen 2000)
        invtau[self.rho < RHO_I / np.sqrt(b)]   = 1.0 - (b * (self.rho[self.rho < RHO_I / np.sqrt(b)] / RHO_I)**2)
        # b           = 0.25 # Tortuosity parameter (Emma's value), effective b = 0.0625 (taken outside of squared term)
        # invtau[self.rho < RHO_I / np.sqrt(b)]   = 1.0 - (b * (self.rho[self.rho < RHO_I / np.sqrt(b)] / RHO_I))**2 #Emma
        invtau[self.rho >= RHO_I / np.sqrt(b)]  = 0.0

        ### Set diffusivity for each isotope
        c_vol = np.ones_like(self.rho) # Just a filler here to make the diffusion work.
        if ((self.isotope == '18') or (self.isotope == 'd18O')):
            Da_18   = Da / 1.0285 # account for fractionation factor for 18_O, fixed Johnsen typo
            D       = m * pz * invtau * Da_18 * (1 / self.rho - 1 / RHO_I) / (R * self.Tz * alpha_18_z)
            D       = D + 1.5e-15 # Emma added - not sure why? prevent negative?
            self.del_z  = transient_solve_TR(z_edges_vec, z_P_vec, nt, self.dt, D, phi_0, nz_P, nz_fv, phi_s, self.rho, c_vol)

        elif ((self.isotope == 'D') or (self.isotope == 'dD')):
            Da_D    = Da / 1.0251 # account for fractionation factor for D, fixed Johnsen typo
            D       = m * pz * invtau * Da_D * (1 / self.rho - 1 / RHO_I) / (R * self.Tz * alpha_D_z)
            D[D<=0.0]   = 1.0e-20
            self.del_z  = transient_solve_TR(z_edges_vec, z_P_vec, nt, self.dt, D, phi_0, nz_P, nz_fv, phi_s, self.rho, c_vol)

        elif ((self.isotope == 'NoDiffusion') or (self.isotope == 'ND')):
            D = np.zeros_like(self.z)           

        dsig2_dt = 2 * (-1*self.drho_dt/self.rho) * self.iso_sig2_z + 2 * D
        self.iso_sig2_z = self.iso_sig2_z + dsig2_dt * self.dt
            
        # Advect profile down
        if self.bdot>0.0:
            self.del_z = np.concatenate(([self.del_s[iii]], self.del_z[:-1]))
            self.iso_sig2_z = np.concatenate(([self.iso_sig2_s[iii]],self.iso_sig2_z[:-1]))
        else:
            pass

        # self.del_z[0] = self.del_z[1]
        # print('Caution!!! You are altering the upper isotope value! Line 173, isotopeDiffusion.py')

        return self.del_z, self.iso_sig2_z 
