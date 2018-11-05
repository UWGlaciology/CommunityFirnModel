#!/usr/bin/env python
from diffusion import heatDiff
from diffusion import isoDiff
from hl_analytic import hl_analytic
from reader import read_input
from writer import write_spin_hdf5
from physics import *
from constants import *
import numpy as np
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
# from merge import mergeall
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
        self.spin=False
        with open(configName, "r") as f:
            jsonString  = f.read()
            self.c      = json.loads(jsonString)

        print('Spin run started')
        print("physics are", self.c['physRho'])
        try:
            print('Merging is:',self.c["merging"])
        except Exception:
            pass

        ### create directory to store results. Deletes if it exists already.
        if os.path.exists(self.c['resultsFolder']):
            rmtree(self.c['resultsFolder'])
        os.makedirs(self.c['resultsFolder'])

        ############################
        ##### load input files #####
        ############################
        ### temperature
        input_temp, input_year_temp = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNameTemp']))
        if input_temp[0] < 0.0:
            input_temp              = input_temp + K_TO_C
        try:
            if self.c['spinup_climate_type']=='initial':
                self.temp0                  = input_temp[0]
            elif self.c['spinup_climate_type']=='mean':
                self.temp0                  = np.mean(input_temp)
        except:
            print("You should add key 'spinup_climate_type' to the config .json file")
            print("spinup is based on mean climate of input")
            self.temp0                  = np.mean(input_temp)
        
        # self.temp0 = 270.0
        # self.temp0 = mean(input_temp[0:12]) #Make sure that this is what we want!

        ### accumulation rate
        input_bdot, input_year_bdot = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamebdot']))
        
        try:
            if self.c['spinup_climate_type']=='initial':
                self.bdot0      = input_bdot[0]
            elif self.c['spinup_climate_type']=='mean':
                self.bdot0      = np.mean(input_bdot)
        except:
            self.bdot0      = np.mean(input_bdot)

        try:
            if self.c['initprofile']:
                try:
                    if self.c['singleleap'] or self.c['singlenormal']: #VV If we use the initprofile but the input forcing data does not cover an entire year
                        self.temp0 = self.c['deepT'] #VV use deep T as mean temperature for spin up calculations (compaction,grain growth)
                        self.bdot0 = self.c['bdot_long']*1e-3/0.917 #VV use long term accumulation as mean accumulation for spin up calculations (compaction,grain growth) + conversion from mmWE/yr to mIE/yr
                except Exception:
                    print("you are using an initial profile for spinup,")
                    print("but you do not have 'singleleap' or 'singlenormal' in the .json")
        except Exception:
            pass
        
        print('bdot0', self.bdot0)
        ### could include others, e.g. surface density
        ############################

        ############################
        ### set up model grid ######
        ############################
        self.gridLen    = int((self.c['H'] - self.c['HbaseSpin']) / (self.bdot0 / self.c['stpsPerYearSpin'])) # number of grid points
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
        # if self.c['initprofile'] == False: #VV
        THL                 = self.temp0
        AHL                 = self.bdot0
        self.age, self.rho     = hl_analytic(self.c['rhos0'], self.z, THL, AHL) # self.age is in age in seconds
        # elif self.c['initprofile'] == True: # VV we have an initial temperature and density profile
            # Just avoid the model to blow up because it has no age and rho variables, real values are given below
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
        
        self.dt     = S_PER_YEAR / self.c['stpsPerYearSpin']
        self.stp    = int(self.years*S_PER_YEAR/self.dt)
        self.t      =  1.0 / self.c['stpsPerYearSpin'] # years per time step
        ############################

        ############################
        ### Initial and boundary conditions
        ############################
        ### Surface temperature for each time step
        self.Ts         = self.temp0 * np.ones(self.stp)
        if self.c['SeasonalTcycle']: #impose seasonal temperature cycle of amplitude 'TAmp'
            # self.Ts     = self.Ts + self.c['TAmp'] * (np.cos(2 * np.pi * np.linspace(0, self.years, self.stp )) + 0.3 * np.cos(4 * np.pi * np.linspace(0, self.years, self.stp ))) #Orsi, coreless winter
            self.Ts         = self.Ts - self.c['TAmp'] * (np.cos(2 * np.pi * np.linspace(0, self.years, self.stp))) # This is basic for Greenland (for Antarctica the it should be a plus instead of minus)

        ### initial temperature profile
        # init_Tz       = input_temp[0] * np.ones(self.gridLen)
        init_Tz         = np.mean(self.Ts) * np.ones(self.gridLen)

        ### Accumulation rate for each time step
        self.bdotSec0   = self.bdot0 / S_PER_YEAR / self.c['stpsPerYearSpin'] # accumulation (m I.E. per second)
        self.bdotSec    = self.bdotSec0 * np.ones(self.stp) # vector of accumulation at each time step

        ### Surface isotope values for each time step
        if self.c['isoDiff']:
            try:
                input_iso, input_year_iso = read_input(self.c['InputFileNameIso'])
                del_s0  = input_iso[0]
            except:
                print('No external file for surface isotope values found, but you specified in the config file that isotope diffusion is on. The model will generate its own synthetic isotope data for you.')
                del_s0  = -50.0

            self.del_s  = del_s0 * np.ones(self.stp)
            init_del_z  = del_s0 * np.ones(self.gridLen)
            self.del_z  = init_del_z
        else:
            self.del_s  = None
            init_del_z  = None    
        
        ### Surface Density
        self.rhos0      = self.c['rhos0'] * np.ones(self.stp)
        # could configure this so that user specifies vector of surface elevation
        # could add noise too

        ### initial mass, stress, and mean accumulation rate
        self.mass       = self.rho * self.dz
        self.sigma      = self.mass * self.dx * GRAVITY
        self.sigma      = self.sigma.cumsum(axis = 0)
        self.mass_sum   = self.mass.cumsum(axis = 0)
        #VVself.bdot_mean  = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] / (self.age[1:] * RHO_I / self.t)))) * self.c['stpsPerYear'] * S_PER_YEAR

        ### Mean temperatures and accumulation rates #VV ###
        ### MS: commented out for now 10/4/18
        # if self.c['initprofile']:
        #     self.T_av = self.c['Ts_long'] * np.ones(self.stp) # If we use initial profile, Ts_long is best guess for Tav
        #     self.bdot_av = self.c['bdot_long']*1e-3/0.917 * np.ones(self.stp) # If we use initial profile, bdotlong is best guess for bdot_av + conversion from mmWE/yr to mIE/yr
        #     self.bdot_mean = np.ones_like(self.dz)*self.c['bdot_long']*1e-3/0.917 # Also best guess for bdot_mean
        # else:
        self.T_av = self.temp0 * np.ones(self.stp) #VV This is to be used instead of T10m
        self.bdot_av = self.bdot0* np.ones(self.stp) #VV 
        self.bdot_mean  = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] / (self.age[1:] * RHO_I / self.t)))) * self.c['stpsPerYear'] * S_PER_YEAR #VV

        ### longitudinal strain rate
        if self.c['strain']:
            self.du_dx      = np.zeros(self.gridLen)
            self.du_dx[1:]  = self.c['du_dx']/(S_PER_YEAR)
        
        ### initial temperature grid 
        self.Tz         = init_Tz
        self.T_mean     = np.mean(self.Tz[self.z<50])
        self.T10m       = self.T_mean

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
            # initial temperature history function (units seconds)
            self.Hx     = np.exp(-110.0e3/(R*init_Tz))*(self.age+self.dt)
            self.THist  = True
        else:
            self.THist  = False

        self.LWC = np.zeros_like(self.z)
        self.MELT = False
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
                'rho':          self.rho,
                'mass':         self.mass,
                'sigma':        self.sigma,
                'dt':           self.dt,
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
                'T_av':         self.T_av,
                'bdot_av':      self.bdot_av
            }

            if self.THist:
                PhysParams['Hx'] = self.Hx

            ### choose densification-physics based on user input
            physicsd = {
                'HLdynamic':            FirnPhysics(PhysParams).HL_dynamic,
                'HLSigfus':             FirnPhysics(PhysParams).HL_Sigfus,
                'Barnola1991':          FirnPhysics(PhysParams).Barnola_1991,
                'Li2004':               FirnPhysics(PhysParams).Li_2004,
                'Li2011':               FirnPhysics(PhysParams).Li_2011,
                'Ligtenberg2011':       FirnPhysics(PhysParams).Ligtenberg_2011,
                'Arthern2010S':         FirnPhysics(PhysParams).Arthern_2010S,
                'Simonsen2013':         FirnPhysics(PhysParams).Simonsen_2013,
                'Morris2014':           FirnPhysics(PhysParams).Morris_HL_2014,
                'Helsen2008':           FirnPhysics(PhysParams).Helsen_2008,
                'Arthern2010T':         FirnPhysics(PhysParams).Arthern_2010T,
                'Goujon2003':           FirnPhysics(PhysParams).Goujon_2003,
                'KuipersMunneke2015':   FirnPhysics(PhysParams).KuipersMunneke_2015,
                'Crocus':               FirnPhysics(PhysParams).Crocus,
                'Max2018':              FirnPhysics(PhysParams).Max2018
            }

            RD      = physicsd[self.c['physRho']]()
            drho_dt = RD['drho_dt']

            ### update density and age of firn
            self.age = np.concatenate(([0], self.age[:-1])) + self.dt
            self.rho = self.rho + self.dt * drho_dt
            
            if self.THist:
                self.Hx = FirnPhysics(PhysParams).THistory()

            ### update temperature grid and isotope grid if user specifies
            if self.c['heatDiff']:
                self.Tz, self.T10m = heatDiff(self,iii)

            if self.c['isoDiff']:
                self.del_z  = isoDiff(self,iii)

            self.T_mean     = np.mean(self.Tz[self.z<50])

            if self.c['strain']: # consider additional change in box height due to longitudinal strain rate
                self.dz     = ((-self.du_dx)*self.dt + 1)*self.dz 
                self.mass   = self.mass*((-self.du_dx)*self.dt + 1)

            ### update model grid mass, stress, and mean accumulation rate
            dzNew           = self.bdotSec[iii] * RHO_I / self.rhos0[iii] * S_PER_YEAR
            self.dz         = self.mass / self.rho * self.dx
            self.dz_old     = self.dz    
            self.dz         = np.concatenate(([dzNew], self.dz[:-1]))     
            self.z          = self.dz.cumsum(axis = 0)
            self.z          = np.concatenate(([0], self.z[:-1]))
            self.rho        = np.concatenate(([self.rhos0[iii]], self.rho[:-1]))
            self.Tz         = np.concatenate(([self.Ts[iii]], self.Tz[:-1]))
            massNew         = self.bdotSec[iii] * S_PER_YEAR * RHO_I
            self.mass       = np.concatenate(([massNew], self.mass[:-1]))
            self.sigma      = self.mass * self.dx * GRAVITY
            self.sigma      = self.sigma.cumsum(axis = 0)
            self.mass_sum   = self.mass.cumsum(axis = 0)
            self.bdot_mean  = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] * self.t / (self.age[1:] * RHO_I))))*self.c['stpsPerYear']*S_PER_YEAR
                
            # update grain radius
            if self.c['physGrain']:
                # self.r2, self.dr2_dt    = FirnPhysics(PhysParams).grainGrowth() #old
                ## Vincent's new:
                self.r2 = FirnPhysics(PhysParams).graincalc()
                r2surface = FirnPhysics(PhysParams).surfacegrain()
                self.r2 = np.concatenate(([r2surface], self.r2[:-1])) #VV form the new grain size array

            if self.doublegrid:
                self.gridtrack = np.concatenate(([1],self.gridtrack[:-1]))
                if self.gridtrack[-1]==2:
                    self.dz, self.z, self.rho, self.Tz, self.mass, self.sigma, self. mass_sum, self.age, self.bdot_mean, self.LWC, self.gridtrack, self.r2 = regrid(self)

            # write results at the end of the time evolution
            if (iii == (self.stp - 1)):

                if self.c['initprofile']:
                    initfirn = pd.read_csv(self.c['initfirnFile'],delimiter=',')

                    init_depth      = initfirn['depth'].values
                    self.rho = np.interp(self.z,init_depth,initfirn['density'].values)
                    if 'temperature' in list(initfirn):
                        init_temp = initfirn['temperature'].values
                        if init_temp[0]<0:
                            init_temp = init_temp + 273.15
                        self.Tz = np.interp(self.z,init_depth,init_temp)                      
                    if 'age' in list(initfirn):
                        self.age = np.interp(self.z,init_depth,initfirn['age'].values)
                    if 'lwc' in list(initfirn):
                        self.LWC = np.interp(self.z,init_depth,initfirn['lwc'].values)

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
                    self.iso_time    = np.concatenate(([self.t * iii + 1], self.del_z))
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

                write_spin_hdf5(self)

            ####################################
            ##### END TIME-STEPPING LOOP #####
            ####################################
