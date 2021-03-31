#!/usr/bin/env python

'''

Class for spinning up the model.

'''

from diffusion import heatDiff
from hl_analytic import hl_analytic
from reader import read_input
from writer import write_spin_hdf5
from physics import *
from constants import *
from isotopeDiffusion import isotopeDiffusion
import numpy as np
import scipy.interpolate as interpolate
import csv
import json
import sys
import math
from shutil import rmtree
import os
import shutil
import time
import h5py
from regrid import *
try:
    from merge import mergeall
except Exception:
    print('CFMmerge not found; preferential flow will not work')
try:
    import pandas as pd
except:
    print('You do not have the pandas python package installed. It is')
    print('only required if you are importing an initial condition.')

class FirnDensitySpin:
    '''

    Parameters used in the model, for the initialization as well as the time evolution:

    : gridLen: size of grid used in the model run
                (unit: number of boxes, type: int)
    : dx: vector of width of each box, used for stress calculations
                (unit: m, type: array of ints)
    : dz: vector of thickness of each box
                (unit: m, type: float)
    : z:  vector of edge locations of each box (value is the top of the box)
                (unit: m, type: float)
    : dt: number of seconds per time step
                (unit: seconds, type: float)
    : t: number of years per time step
                (unit: years, type: float)
    : modeltime: linearly spaced time vector from indicated start year to indicated end year
                (unit: years, type: array of floats)
    : years: total number of years in the model run
                (unit: years, type: float)
    : stp: total number of steps in the model run
                (unit: number of steps, type: int)
    : T_mean: interpolated temperature vector based on the model time and the initial user temperature data
                (unit: ???, type: array of floats)
    : Ts: interpolated temperature vector based on the model time & the initial user temperature data
                may have a seasonal signal imposed depending on number of years per time step (< 1)
                (unit: ???, type: array of floats)
    : bdot: bdot is meters of ice equivalent/year. multiply by 0.917 for W.E. or 917.0 for kg/year
                (unit: ???, type: )
    : bdotSec: accumulation rate vector at each time step
                (unit: ???, type: array of floats)
    : rhos0: surface accumulate rate vector
                (unit: ???, type: array of floats)

    :returns D_surf: diffusivity tracker
                (unit: ???, type: array of floats)

    '''

    def __init__(self, configName):
        '''

        Sets up the initial spatial grid, time grid, accumulation rate, age, density, mass, stress, and temperature of the model run
        :param configName: name of json config file containing model configurations

        '''

        ### load in json config file and parses the user inputs to a dictionary
        self.spin=True
        with open(configName, "r") as f:
            jsonString  = f.read()
            self.c      = json.loads(jsonString)

        print('Spin run started')
        print("physics are", self.c['physRho'])
        try:
            print('Merging is:',self.c["merging"])
        except Exception:
            print('"merging" missing from .json fields')
            pass

        ### create directory to store results. Deletes if it exists already.
        # Vincent says we do not want to remove existing (preferential flow?) - 4/24/19
        if os.path.exists(self.c['resultsFolder']):
            dir_exts = [os.path.splitext(fname)[1] for fname in os.listdir(self.c['resultsFolder'])]
            dir_unique = list(set(dir_exts))
            print('dir_unique',dir_unique)
            CFM_exts = ['.json','.hdf5']
            if CFM_exts and all(((elem == ".json") or (elem=='.hdf5')) for elem in dir_unique):
                print('rmtree')
                rmtree(self.c['resultsFolder'])
                os.makedirs(self.c['resultsFolder'])
            else:
                print('WARNING: THE DIRECTORY YOU ARE USING CONTAINS NON-CFM FILES')
                print('CFM will delete all files in the results directory with .hdf5 extension')
                files_in_directory = os.listdir(self.c['resultsFolder'])
                filtered_files = [file for file in files_in_directory if file.endswith(".hdf5")]
                for file in filtered_files:
                    path_to_file = os.path.join(self.c['resultsFolder'], file)
                    os.remove(path_to_file)

        else:
            os.makedirs(self.c['resultsFolder'])

        ############################
        ##### load input files #####
        ############################

        ### temperature ###
        input_temp, input_year_temp = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNameTemp']))
        if input_temp[0] < 0.0:
            input_temp              = input_temp + K_TO_C
        try:
            if self.c['spinup_climate_type']=='initial':
                self.temp0                  = input_temp[0]
            elif self.c['spinup_climate_type']=='mean':
                self.temp0                  = np.mean(input_temp)

        except Exception:
            print("You should add key 'spinup_climate_type' to the config .json file")
            print("spinup is based on mean climate of input")
            self.temp0                  = np.mean(input_temp)

        ### accumulation rate ###
        input_bdot, input_year_bdot = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamebdot']))
        try:
            if self.c['spinup_climate_type']=='initial':
                self.bdot0      = input_bdot[0]
            elif self.c['spinup_climate_type']=='mean':
                self.bdot0      = np.mean(input_bdot)
        except:
            self.bdot0      = np.mean(input_bdot)

        if 'manual_climate' in self.c:
            pass
        else:
            print("Add 'manual_climate' to the json to enable specifying long-term bdot and T")
            self.c['manual_climate']=False

        if self.c['manual_climate']: # If we want to use a manually specified climate for spin up (e.g. known long-term values).
            self.temp0 = self.c['deepT'] #specify deep T as mean temperature for spin up calculations (compaction,grain growth)
            self.bdot0 = self.c['bdot_long']# *1e-3/0.917 #specify long term accumulation as mean accumulation for spin up calculations (compaction,grain growth) + conversion from mmWE/yr to mIE/yr
            print('make sure "bdot_long" has units of mIE/yr!')

        print('Spin-up accumulation rate is', self.bdot0)
        print('Spin-up temperature is', self.temp0)
        ### could include others, e.g. surface density
        ############################

        ############################
        ### set up model grid ######
        ############################
        self.gridLen    = int((self.c['H'] - self.c['HbaseSpin']) / (self.bdot0 / self.c['stpsPerYear'])) # number of grid points
        gridHeight      = np.linspace(self.c['H'], self.c['HbaseSpin'], self.gridLen)
        self.z          = self.c['H'] - gridHeight
        self.dz         = np.diff(self.z)
        self.dz         = np.append(self.dz, self.dz[-1])
        self.dx         = np.ones(self.gridLen)
        print('Grid length is', self.gridLen)
        ############################

        ############################
        ### if the regridding module is being used, do the
        ### initial regridding
        ############################
        try:
            self.doublegrid = self.c['doublegrid']
            if self.c['doublegrid']:
                self.nodestocombine, self.z, self.dz, self.gridLen, self.dx, self.gridtrack = init_regrid(self)
        except:
            self.doublegrid = False
            print('you should add "doublegrid" to the json')

        ############################
        ### get an initial depth/density profile based on H&L analytic solution
        ############################
        # if not self.c['initprofile']: #VV
        THL                 = self.temp0
        AHL                 = self.bdot0

        try: #VV use Reeh corrected T
            if self.c['ReehCorrectedT'] and self.c['MELT']:
                input_snowmelt, input_year_snowmelt = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamemelt'])) #VV
                meanmelt = np.mean(input_snowmelt) # mean melt per year [mIE/yr] (units are specified in Reeh 2008)
                meanacc  = self.bdot0 # mean annual accumulation [mIE/yr]
                self.SIR = min(meanmelt,0.6*meanacc) # Reeh 1991 and Reeh 2008 PMAX value is set at 0.6 melt becomes superimposed ice until it reaches 0.6 of annual acc, then runoff
                THL = self.temp0 + 26.6*self.SIR
                THL = min(THL,273.15)
            elif (self.c['ReehCorrectedT'] and not self.c['MELT']):
                print('"ReehCorrectedT" is True but melt is not turned on. That is wierd. Exiting.')
                sys.exit()

        except:
            print('add "ReehCorrectedT" to .json to enable melt-corrected temperature')
            self.c['ReehCorrectedT'] = False
            pass

        self.age, self.rho     = hl_analytic(self.c['rhos0'], self.z, THL, AHL) # self.age is in age in seconds

        # if self.c['initprofile']: # VV filler values to avoid model blow up if THL and AHL are out of HL calibration range
            # self.age = S_PER_YEAR*100*np.ones_like(self.dz) #VV this does not matter as it is rectified when we initialise profie below
            # self.rho = 500*np.ones_like(self.dz)#VV this does not matter as it is rectified when we initialise profile
        ############################

        ############################
        ### set up time stepping
        if self.c['AutoSpinUpTime']: # automatic, based on time that it will take for a parcel to get to 850 kg m^-3
            try:
                zz          = np.min(self.z[self.rho > 850.0])
                self.years  = int(zz / self.bdot0)
            except ValueError:
                print("auto spin up error; using spin up time from json")
                self.years = self.c['yearSpin'] # number of years to spin up for
        else: # based on time taken to spin up in the config file.
            self.years = self.c['yearSpin'] # number of years to spin up for

        dt1         = S_PER_YEAR / self.c['stpsPerYear']
        self.stp    = int(self.years*S_PER_YEAR/dt1)
        self.t      =  1.0 / self.c['stpsPerYear'] # years per time step
        self.dt     = dt1 * np.ones(self.stp)
        ############################

        ############################
        ### Initial and boundary conditions
        ############################
        ### Surface temperature for each time step
        self.Ts         = self.temp0 * np.ones(self.stp)

        if self.c['SeasonalTcycle']: #impose seasonal temperature cycle of amplitude 'TAmp'
            if self.c['SeasonalThemi'] == 'north':
                self.Ts         = self.Ts - self.c['TAmp'] * (np.cos(2 * np.pi * np.linspace(0, self.years, self.stp))) # This is for Greenland

            elif self.c['SeasonalThemi'] == 'south':
                if self.c['coreless']:
                    self.Ts     = self.Ts + self.c['TAmp'] * (np.cos(2 * np.pi * np.linspace(0, self.years, self.stp)) + 0.3 * np.cos(4 * np.pi * np.linspace(0, self.years, self.stp))) # Coreless winter, from Orsi
                else:
                    self.Ts     = self.Ts + self.c['TAmp'] * (np.cos(2 * np.pi * np.linspace(0, self.years, self.stp))) # This is basic for Antarctica
            else:
                print('You have turned on the SeasonalTcycle, but you do not have')
                print('the hemisphere selected. Exiting. (set to south or north')
                sys.exit()

        ### initial temperature profile
        # init_Tz         = self.temp0 * np.ones(self.gridLen)
        self.T_mean     = self.temp0 * np.ones(self.stp)
        self.Tz         = self.temp0 * np.ones(self.gridLen)
        self.T50        = np.mean(self.Tz[self.z<50])
        self.T10m       = np.mean(self.T_mean)
        if self.c['ReehCorrectedT']:
            self.Tz         = self.Tz + 26.6*self.SIR # temperature correction accounting for latent heat, Reeh 1991 (5) Reeh 2008 (16)
            self.Tz         = np.minimum(self.Tz,273.15)
            self.T_mean     = np.mean(self.Tz[self.z<50])
            self.T10m       = self.T_mean

        try:
            ctest = self.c['conductivity']
        except:
            self.c['conductivity'] = 'Anderson'

        ### Accumulation rate for each time step
        self.bdotSec0   = self.bdot0 / S_PER_YEAR / self.c['stpsPerYear'] # accumulation (m I.E. per second)
        self.bdotSec    = self.bdotSec0 * np.ones(self.stp) # vector of accumulation at each time step
        self.bdot_mean = np.ones_like(self.dz)*self.bdot0
        self.bdot_av = self.bdot0 * np.ones(self.stp) #The long-term mean accumulation rate

        ### Surface isotope values for each time step
        if self.c['isoDiff']:
            self.Isotopes   = {} #dictionary of class instances
            self.iso_out    = {} # outputs for each isotope
            self.Isoz       = {} # depth profile of each isotope, at each time step
            self.Iso_sig2_z   = {} # diffusion length profile
            self.iso_sig2_out = {}

            for isotope in self.c['iso']:
                self.Isotopes[isotope] = isotopeDiffusion(self.spin,self.c,isotope,self.stp,self.z)

        ### Surface Density
        self.rhos0      = self.c['rhos0'] * np.ones(self.stp) # could configure this so that user specifies vector of some noise

        ### initial mass, stress, and mean accumulation rate
        self.mass       = self.rho * self.dz
        self.sigma      = self.mass * self.dx * GRAVITY
        self.sigma      = self.sigma.cumsum(axis = 0)
        self.mass_sum   = self.mass.cumsum(axis = 0)

        ### longitudinal strain rates
        if 'strain' not in self.c:
            self.c['strain'] = False
        if 'strain_softening' not in self.c:
            self.c['strain_softening'] = False
        if 'residual_strain' not in self.c:
            self.c['residual_strain'] = 2e-4
        if 'strain_correction' not in self.c:
            self.c['strain_correction'] = False
        if 'strain_heating' not in self.c:
            self.c['strain_heating'] = False

        if self.c['strain']:
            input_eps, input_year_eps = read_input(os.path.join(self.c['InputFileFolder'], self.c['InputFileNameStrain']))
            input_eps_1, input_eps_2 = input_eps[0, :], input_eps[1, :]
            if self.c['spinup_climate_type'] == 'initial':
                self.eps_1 = input_eps_1[0]
                self.eps_2 = input_eps_2[0]
            elif self.c['spinup_climate_type'] == 'mean':
                self.eps_1 = np.mean(input_eps_1)
                self.eps_2 = np.mean(input_eps_2)
            else:
                self.eps_1 = np.mean(input_eps_1)
                self.eps_2 = np.mean(input_eps_2)
            print('Spin-up strain rates are', self.eps_1, 'and', self.eps_2)
            self.eps_eff_hor_2 = (self.eps_1 ** 2 + self.eps_2 ** 2) / 2
            if self.c['strain_correction']:
                eps_eff_hor = np.sqrt(self.eps_eff_hor_2)
                eps_eff_hor = eps_eff_hor - 4e-4
                if eps_eff_hor < 0:
                    eps_eff_hor = 0
                self.eps_eff_hor_2 = eps_eff_hor**2
            self.eps_zz        = - (self.eps_1 + self.eps_2)
        else:
            if self.c['strain_softening']:
                self.c['strain_softening']=False
            if self.c['strain_heating']:
                self.c['strain_heating']=False

        ### initial grain growth (if specified in config file)
        if self.c['physGrain']:
            if self.c['calcGrainSize']:
                r02 = surfacegrain(self,0) #VV
                self.r2 = r02 * np.ones(self.gridLen)
            else:
                self.r2 = np.linspace(self.c['r2s0'], (6 * self.c['r2s0']), self.gridLen)
        else:
            self.r2 = None

        ### "temperature history" if using Morris physics
        if self.c['physRho']=='Morris2014':
            if 'QMorris' not in self.c:
                print('Add "QMorris" to the .json. CFM is setting to 110e3')
                self.c['QMorris'] = 110.0e3
            # initial temperature history function (units seconds)
            self.Hx     = np.exp(-1*self.c['QMorris']/(R*self.Tz))*(self.age+self.dt[0])
            self.THist  = True
        else:
            self.THist  = False

        self.LWC = np.zeros_like(self.z)
        self.MELT = False

        ### values for Goujon physics
        if self.c['physRho']=='Goujon2003':
            self.Gamma_Gou      = 0
            self.Gamma_old_Gou  = 0
            self.Gamma_old2_Gou = 0
            self.ind1_old       = 0
        #######################

        #######################
        try:
            if self.c['no_densification']:
                print('CAUTION: densification if OFF!')
            else:
                pass
        except:
            print('no_densification not in .json; setting to false')
            self.c['no_densification']=False
        #######################

    ############################
    ##### END INIT #############
    ############################

    def time_evolve(self):
        '''
        Evolve the spatial grid, time grid, accumulation rate, age, density, mass, stress, and temperature through time
        based on the user specified number of timesteps in the model run. Updates the firn density using a user specified
        '''
        self.steps = 1 / self.t # this is time steps per year

        ####################################
        ##### START TIME-STEPPING LOOP #####
        ####################################

        for iii in range(self.stp):
            ### create dictionary of the parameters that get passed to physics
            PhysParams = {
                'iii':          iii,
                'steps':        self.steps,
                'gridLen':      self.gridLen,
                'bdotSec':      self.bdotSec,
                'bdot_mean':    self.bdot_mean,
                'bdot_type':    self.c['bdot_type'],
                'Tz':           self.Tz,
                'T_mean':       self.T_mean,
                'T10m':         self.T10m,
                'T50':          self.T50,
                'rho':          self.rho,
                'mass':         self.mass,
                'sigma':        self.sigma,
                'dt':           self.dt[iii],
                'Ts':           self.Ts,
                'r2':           self.r2,
                'age':          self.age,
                'physGrain':    self.c['physGrain'],
                'calcGrainSize':self.c['calcGrainSize'],
                'r2s0':         self.c['r2s0'],
                'GrGrowPhysics':self.c['GrGrowPhysics'],
                'z':            self.z,
                'rhos0':        self.rhos0[iii],
                'dz':           self.dz,
                'LWC':          self.LWC,
                'MELT':         self.MELT,
                'FirnAir':      False,
                'bdot_av':      self.bdot_av
            }

            if self.c['physRho']=='Morris2014':
                PhysParams['Hx'] = self.Hx
                PhysParams['QMorris'] = self.c['QMorris']

            if self.c['physRho']=='Goujon2003':
                PhysParams['Gamma_Gou']      = self.Gamma_Gou
                PhysParams['Gamma_old_Gou']  = self.Gamma_old_Gou
                PhysParams['Gamma_old2_Gou'] = self.Gamma_old2_Gou
                PhysParams['ind1_old']       = self.ind1_old

            ### choose densification-physics based on user input
            physicsd = {
                'HLdynamic':            FirnPhysics(PhysParams).HL_dynamic,
                'HLSigfus':             FirnPhysics(PhysParams).HL_Sigfus,
                'Barnola1991':          FirnPhysics(PhysParams).Barnola_1991,
                'Li2004':               FirnPhysics(PhysParams).Li_2004,
                'Li2011':               FirnPhysics(PhysParams).Li_2011,
                'Li2015':               FirnPhysics(PhysParams).Li_2015,
                'Ligtenberg2011':       FirnPhysics(PhysParams).Ligtenberg_2011,
                'Arthern2010S':         FirnPhysics(PhysParams).Arthern_2010S,
                'Simonsen2013':         FirnPhysics(PhysParams).Simonsen_2013,
                'Morris2014':           FirnPhysics(PhysParams).Morris_HL_2014,
                'Helsen2008':           FirnPhysics(PhysParams).Helsen_2008,
                'Arthern2010T':         FirnPhysics(PhysParams).Arthern_2010T,
                'Goujon2003':           FirnPhysics(PhysParams).Goujon_2003,
                'KuipersMunneke2015':   FirnPhysics(PhysParams).KuipersMunneke_2015,
                'Crocus':               FirnPhysics(PhysParams).Crocus,
                'GSFC2020':             FirnPhysics(PhysParams).GSFC2020,
                'MaxSP':                FirnPhysics(PhysParams).MaxSP
            }

            RD      = physicsd[self.c['physRho']]()
            drho_dt = RD['drho_dt']
            if self.c['no_densification']:
                drho_dt = np.zeros_like(drho_dt)
            self.viscosity = RD['viscosity']

            if self.c['strain']: # consider additional change in box height due to longitudinal strain rate
                eps_zz_classic = - drho_dt / self.rho * S_PER_YEAR  # Densification rate in case of no strain.
                eps_zz_classic = eps_zz_classic - np.abs(self.c['residual_strain'])

                if 'strain_softening' in self.c:
                    if self.c['strain_softening']:
                        eps_eff_classic_2 = 0.5 * (eps_zz_classic ** 2)
                        r_hor2 = self.eps_eff_hor_2 / eps_eff_classic_2
                        K = (13.5 * r_hor2 + 1.5 * np.sqrt(81 * r_hor2 ** 2 + 12 * r_hor2) + 1) ** (1 / 3)
                        correction_factor_viscosity = (K + 1 / K + 1) / 3

                        z2mask = (self.rho >= RHO_1)
                        eps_zz_classic[z2mask] = eps_zz_classic[z2mask] * correction_factor_viscosity[z2mask]
                        drho_dt[z2mask]        = drho_dt[z2mask] * correction_factor_viscosity[z2mask]
                        self.viscosity[z2mask] = self.viscosity[z2mask] / correction_factor_viscosity[z2mask]

                if 'strain_heating' in self.c:
                    if self.c['strain_heating']:
                        self.eps_sum = eps_zz_classic
                        self.eps_eff2 = 0.5 * (self.eps_1**2 + self.eps_2**2 + (eps_zz_classic + self.eps_zz)**2)

                self.mass   = (1 + self.eps_zz / S_PER_YEAR * self.dt[iii]) * self.mass

            if self.c['physRho']=='Goujon2003':
                self.Gamma_Gou      = RD['Gamma_Gou']
                self.Gamma_old_Gou  = RD['Gamma_old_Gou']
                self.Gamma_old2_Gou = RD['Gamma_old2_Gou']
                self.ind1_old       = RD['ind1_old']

            ### update density and age of firn
            self.age = np.concatenate(([0], self.age[:-1])) + self.dt[iii]
            self.rho = self.rho + self.dt[iii] * drho_dt

            if self.THist:
                self.Hx = RD['Hx']

            ### update temperature grid and isotope grid if user specifies
            if self.c['heatDiff']:
                self.Tz, self.T10m = heatDiff(self,iii)

            if self.c['isoDiff']:
                IsoParams = {
                    'Tz':           self.Tz,
                    'rho':          self.rho,
                    'dt':           self.dt[iii],
                    'z':            self.z,
                    'rhos0':        self.rhos0[iii],
                    'dz':           self.dz,
                    'drho_dt':      drho_dt,
                    'bdot':         self.bdotSec[iii]
                }

                for isotope in self.c['iso']:
                    self.Isoz[isotope], self.Iso_sig2_z[isotope] = self.Isotopes[isotope].isoDiff(IsoParams,iii)

            self.T50     = np.mean(self.Tz[self.z<50])

            ### update model grid mass, stress, and mean accumulation rate
            dzNew           = self.bdotSec[iii] * RHO_I / self.rhos0[iii] * S_PER_YEAR
            self.dz_old     = np.copy(self.dz)
            self.dz         = self.mass / self.rho * self.dx
            self.dz         = np.concatenate(([dzNew], self.dz[:-1]))
            self.z          = self.dz.cumsum(axis = 0)
            self.z          = np.concatenate(([0], self.z[:-1]))
            self.rho        = np.concatenate(([self.rhos0[iii]], self.rho[:-1]))

            ### VV corrected temperature profile with latent heat release from meltwater,
            ### following Reeh 1991 parameterisation ##
            if self.c['ReehCorrectedT']:
                self.Tz         = np.concatenate(([self.Ts[iii]]+26.6*self.SIR, self.Tz[:-1]))
            else:
                self.Tz         = np.concatenate(([self.Ts[iii]], self.Tz[:-1]))
            ##

            massNew         = self.bdotSec[iii] * S_PER_YEAR * RHO_I
            self.mass       = np.concatenate(([massNew], self.mass[:-1]))
            self.sigma      = self.mass * self.dx * GRAVITY
            self.sigma      = self.sigma.cumsum(axis = 0)
            self.mass_sum   = self.mass.cumsum(axis = 0)
            self.bdot_mean  = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] * self.t / (self.age[1:] * RHO_I))))*self.c['stpsPerYear']*S_PER_YEAR

            ### Update grain growth #VV ###
            #VV calculate this before accumulation (because the new surface layer should not be subject to grain growth yet
            if self.c['physGrain']:
                self.r2 = FirnPhysics(PhysParams).graincalc(iii)
                r2surface = FirnPhysics(PhysParams).surfacegrain() # This considers whether to use a fixed or calculated surface grain size.
                self.r2 = np.concatenate(([r2surface], self.r2[:-1])) #VV form the new grain size array

            if self.doublegrid:
                self.gridtrack = np.concatenate(([1],self.gridtrack[:-1]))
                if self.gridtrack[-1]==2:
                    self.dz, self.z, self.rho, self.Tz, self.mass, self.sigma, self. mass_sum, self.age, self.bdot_mean, self.LWC, self.gridtrack, self.r2 = regrid(self)

            # write results at the end of the time evolution
            if (iii == (self.stp - 1)):
                if self.c['initprofile']:
                    print('Updating density using init file')
                    initfirn = pd.read_csv(self.c['initfirnFile'],delimiter=',')
                    init_depth      = initfirn['depth'].values
                    self.rho = np.interp(self.z,init_depth,initfirn['density'].values)
                    if 'temperature' in list(initfirn):
                        print('and temperature')
                        init_temp = initfirn['temperature'].values
                        if init_temp[0]<0:
                            init_temp = init_temp + 273.15
                        self.Tz = np.interp(self.z,init_depth,init_temp)
                    if 'age' in list(initfirn):
                        print('and age')
                        self.age = np.interp(self.z,init_depth,initfirn['age'].values*S_PER_YEAR)
                    if 'lwc' in list(initfirn):
                        self.LWC = np.interp(self.z,init_depth,initfirn['lwc'].values)
                    # if 'bdot_mean' in list(initfirn):
                    #     self.write_bdot = True
                    #     self.bdot_mean = np.interp(self.z,init_depth,initfirn['bdot_mean'].values)

                # ### Manual gridding...
                # zold = self.z.copy()
                # dzM = 0.001
                # manualZ = np.arange(0,20+dzM,dzM) #XXX
                # self.z = manualZ
                # self.Tz = initfirn['temperature']
                # self.rho = initfirn['density']
                # self.age = np.interp(self.z,zold,self.age)
                # ###

                self.rho_time        = np.concatenate(([self.t * iii + 1], self.rho))
                self.Tz_time         = np.concatenate(([self.t * iii + 1], self.Tz))
                self.age_time        = np.concatenate(([self.t * iii + 1], self.age))
                self.z_time          = np.concatenate(([self.t * iii + 1], self.z))


                if self.c['physGrain']:
                    self.r2_time     = np.concatenate(([self.t * iii + 1], self.r2))
                else:
                    self.r2_time     = None
                if self.THist:
                    self.Hx_time     = np.concatenate(([self.t * iii + 1], self.Hx))
                else:
                    self.Hx_time     = None
                if self.c['isoDiff']:
                    for isotope in self.c['iso']:
                        self.Iso_sig2_z[isotope] = np.interp(self.z,zold,self.Iso_sig2_z[isotope]) ###XXX


                        self.iso_out[isotope]    = np.concatenate(([self.t * iii + 1], self.Isoz[isotope]))
                        self.iso_sig2_out[isotope] = np.concatenate(([self.t * iii + 1], self.Iso_sig2_z[isotope]))
                        if 'iso{}'.format(isotope) in list(initfirn):
                            print('Interpolating isotope {}'.format(isotope))
                            isoIntFun = interpolate.interp1d(init_depth,initfirn['iso{}'.format(isotope)].values,'nearest',fill_value='extrapolate')
                            self.iso_out[isotope] = np.concatenate(([self.t * iii + 1], isoIntFun(self.z)))


                            # self.iso_out[isotope] = np.interp(self.z,init_depth,initfirn['iso{}'.format(isotope)].values)
                else:
                    self.iso_time    = None
                if self.c['MELT']:
                    self.LWC_time     = np.concatenate(([self.t * iii + 1], self.LWC)) #VV
                else: #VV
                    self.LWC_time     = None #VV
                if self.doublegrid:
                    self.grid_time   = np.concatenate(([self.t * iii + 1], self.gridtrack))
                else:
                    self.grid_time   = None
                # if self.write_bdot:
                    # self.bdot_mean_time = np.concatenate(([self.t * iii + 1], self.bdot_mean))
                # else:
                    # self.bdot_mean_time = None

                write_spin_hdf5(self)

            ####################################
            ##### END TIME-STEPPING LOOP #####
            ####################################
