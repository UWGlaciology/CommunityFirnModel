from diffusion import Diffusion
from reader import read_temp
from reader import read_bdot
from reader import read_init
from writer import write_nospin
from writer import write_nospin_init
from writer import write_nospin_BCO
from writer import write_nospin_LIZ
from writer import write_nospin_DIP
from physics import *
from constants import *
import numpy as np
import csv
import json
import sys
import math
#from scipy import interpolate
#from scipy.sparse import spdiags
#import scipy.sparse.linalg as splin
#from plot import plotData
from shutil import rmtree
#import matplotlib.pyplot as plt
import os
from string import join
import shutil
import time
#import data_interp as IntpData




class FirnDensityNoSpin:
    def __init__(self, configName):
        '''
        Sets up the initial spatial grid, time grid, accumulation rate, age, density, mass, stress, temperature, and diffusivity of the model run
        :param configName: name of json config file containing model configurations
        '''

        # load in json config file and parses the user inputs to a dictionary
        with open(configName, "r") as f:
            jsonString = f.read()
            self.c = json.loads(jsonString)

        # read in initial depth, age, density, temperature from spin-up results
        initDepth, initAge, initDensity, initTemp = read_init(self.c['resultsFolder'])

        # set up the initial age and density of the firn column
        self.age     = initAge[1:]
        self.rho     = initDensity[1:]

        # set up model grid
        self.z  = initDepth[1:]
        self.dz = np.diff(self.z)
        self.dz = np.append(self.dz, self.dz[-1])

        self.gridLen = np.size(self.z)
        self.dx = np.ones(self.gridLen)

        input_temp, input_year_temp = read_temp(self.c['InputFileNameTemp'])
        input_bdot, input_year_bdot = read_bdot(self.c['InputFileNamebdot'])

        yr_start  = max(input_year_temp[0], input_year_bdot[0])
        yr_end    = min(input_year_temp[-1], input_year_bdot[-1])
        
        self.years     = (yr_end - yr_start) * 1.0
        self.stp       = int(self.years * self.c['stpsPerYear'])
        self.modeltime = np.linspace(yr_start, yr_end, self.stp + 1)

        self.dt = self.years * S_PER_YEAR / self.stp
        self.t  = 1.0 / self.c['stpsPerYear']

        self.TPeriod = self.years
        self.Ts = np.interp(self.modeltime, input_year_temp, input_temp)
        self.T_mean = self.Ts
        if self.t < 1.0:
            self.Ts = self.Ts + self.c['TAmp'] * (np.cos(2 * np.pi * np.linspace(0, self.TPeriod, self.stp + 1)) + 0.3 * np.cos(4 * np.pi * np.linspace(0, self.TPeriod, self.stp + 1)))

        self.bdot = np.interp(self.modeltime, input_year_bdot, input_bdot)
        self.bdotSec = self.bdot / S_PER_YEAR / (self.stp / self.years)

        self.rhos0 = self.c['rhos0'] * np.ones(self.stp)
        self.D_surf = self.c['D_surf'] * np.ones(self.stp)

        # set up the initial grid of diffusivity constant
        self.Dcon = self.c['D_surf'] * np.ones(self.gridLen)

        # set up vector of times data will be written
        self.TWrite = self.modeltime[1::c['TWriteInt']]

        # set up initial mass, stress, and mean accumulation rate
        self.mass = self.rho * self.dz
        self.sigma = self.mass * self.dx * GRAVITY
        self.sigma = self.sigma.cumsum(axis = 0)
        self.mass_sum = self.mass.cumsum(axis = 0)
        self.bdot_mean = np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] / (self.age[1:] * RHO_I / self.t)))

        # set up class to handle heat/isotope diffusion using user provided data for initial temperature vector
        self.diffu = Diffusion(self.z, self.stp, self.gridLen, initTemp[1:])

        # set up initial values for density, temperature, age, depth, diffusivity, Clim??, and accumulation to write
        rho_time  = np.append(self.modeltime[0], self.rho)
        Tz_time   = np.append(self.modeltime[0], self.diffu.Tz)
        age_time  = np.append(self.modeltime[0], self.age)
        z_time    = np.append(self.modeltime[0], self.z)
        D_time    = np.append(self.modeltime[0], self.Dcon)
        Clim_time = np.append(self.modeltime[0], [self.bdot[0], self.Ts[0]])  # not sure if bdot or bdotSec
        bdot_time = np.append(self.modeltime[0], self.bdot_mean)

        # set up initial grain growth (if specified in config file)
        if self.c['physGrain']:
            r2Path   = os.path.join(self.c['resultsFolder'], 'r2Spin.csv')
            initr2 = np.genfromtxt(r2Path, delimiter = ',')
            self.r2 = initr2[1:]
            r20 = self.r2
            r2_time = np.append(self.modeltime[0], self.r2)
        else:
            r2_time = None

        # write initial values to the results folder
        write_nospin_init(self.c['resultsFolder'], self.c['physGrain'], rho_time, Tz_time, age_time, z_time, D_time, Clim_time, bdot_time, r2_time)

        # set up initial values for bubble close-off depth & age, lock-in zone depth & age, and depth integrated porosity
        self.bcoAgeMartAll = []
        self.bcoDepMartAll = []
        self.bcoAge815All  = []
        self.bcoDep815All  = []
        self.LIZAgeAll     = []
        self.LIZDepAll     = []
        self.intPhiAll     = []

        self.update_BCO()
        self.update_LIZ()
        self.update_DIP()

    def time_evolve(self):
        '''
        Evolve the spatial grid, time grid, accumulation rate, age, density, mass, stress, temperature, and diffusivity through time
        based on the user specified number of timesteps in the model run. Updates the firn density using a user specified 
        '''
        self.steps = 1 / self.t
        if not self.c['physGrain']:
            r2_time = None

        start_time=time.time()
        for iii in xrange(self.stp):
            mtime = self.modeltime[iii]

            # the parameters that get passed to physics
            PhysParams = {
                'iii':          iii,
                'steps':        self.steps,
                'gridLen':      self.gridLen,
                'bdotSec':      self.bdotSec,
                'bdot_mean':    self.bdot_mean,
                'bdot_type':    self.c['bdot_type'],
                'Tz':           self.diffu.Tz,
                'T_mean':       self.T_mean,
                'rho':          self.rho,
                'sigma':        self.sigma,
                'dt':           self.dt,
                'Ts':           self.Ts,
                'r2':           self.r2,
                'physGrain':    self.c['physGrain'],
                'calcGrainSize':self.c['calcGrainSize']
            }

            # choose densification-physics based on user input
            physicsd = {
                'HLdynamic':       FirnPhysics(PhysParams).HL_dynamic,
                'HLSigfus':        FirnPhysics(PhysParams).HL_Sigfus,
                'Barnola1991':     FirnPhysics(PhysParams).Barnola_1991,
                'Li2004':          FirnPhysics(PhysParams).Li_2004,
                'Li2011':          FirnPhysics(PhysParams).Li_2011,
                'Ligtenberg2011':  FirnPhysics(PhysParams).Ligtenberg_2011,
                'Arthern2010S':    FirnPhysics(PhysParams).Arthern_2010S,
                'Simonsen2013':    FirnPhysics(PhysParams).Simonsen_2013,
                'Morris2013':      FirnPhysics(PhysParams).Morris_HL_2013,
                'Helsen2008':      FirnPhysics(PhysParams).Helsen_2008,
                'Arthern2010T':    FirnPhysics(PhysParams).Arthern_2010T,
                'Spencer2001':     FirnPhysics(PhysParams).Spencer_2001,
                'Goujon2003':      FirnPhysics(PhysParams).Goujon_2003,
            }

            try:
                drho_dt = physicsd[self.c['physRho']]()
            except KeyError:
                default()

            # update density and age of firn
            self.age = np.concatenate(([0], self.age[:-1])) + self.dt
            self.rho = self.rho + self.dt * drho_dt
            self.rho  = np.concatenate(([self.rhos0[iii]], self.rho[:-1]))
            self.Dcon = np.concatenate(([self.D_surf[iii]], self.Dcon[:-1]))

            # update temperature grid and isotope grid if user specifies
            if self.c['heatDiff']:
                self.diffu.heatDiff(self.z, self.dz, self.Ts[iii], self.rho, self.dt)
            if self.c['isoDiff']:
                self.diffu.isoDiff(iii, self.z, self.dz, self.rho, self.c['iso'], self.gridLen, self.dt)

            # update model grid
            dzNew = self.bdotSec[iii] * RHO_I / self.rhos0[iii] * S_PER_YEAR
            self.dz = self.mass / self.rho * self.dx
            self.dz = np.concatenate(([dzNew], self.dz[:-1]))
            self.z = self.dz.cumsum(axis = 0)
            self.z = np.concatenate(([0], self.z[:-1]))

            # update mass, stress, and mean accumulation rate
            massNew = self.bdotSec[iii] * S_PER_YEAR * RHO_I
            self.mass = np.concatenate(([massNew], self.mass[:-1]))
            self.sigma = self.mass * self.dx * GRAVITY
            self.sigma = self.sigma.cumsum(axis = 0)
            self.mass_sum  = self.mass.cumsum(axis = 0)
            self.bdot_mean = np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] * self.t / (self.age[1:] * RHO_I)))

            # update grain radius
            if self.c['physGrain']:
                self.r2 = FirnPhysics(PhysParams).grainGrowth()

            # write results as often as specified in the init method
            if [True for iii in self.TWrite if iii == mtime] == [True]:
                rho_time  = np.append(mtime, self.rho)
                Tz_time   = np.append(mtime, self.diffu.Tz)
                age_time  = np.append(mtime, self.age)
                z_time    = np.append(mtime, self.z)
                Dcon_time = np.append(mtime, self.Dcon)
                Clim_time = np.append(mtime, [self.bdot[int(iii)], self.Ts[int(iii)]])
                bdot_time = np.append(mtime, self.bdot_mean)
                if self.c['physGrain']:
                    r2_time = np.append(mtime, self.r2)

                write_nospin(self.c['resultsFolder'], self.c['physGrain'], rho_time, Tz_time, age_time, z_time, Dcon_time, Clim_time, bdot_time, r2_time)

                self.update_BCO()
                self.update_LIZ()
                self.update_DIP()
        time_done=time.time()
        print 'time for loop=', time_done-start_time, 'seconds'

        # write BCO, LIZ, DIP at the end of the time evolution
        write_nospin_BCO(self.c['resultsFolder'], self.bcoAgeMartAll, self.bcoDepMartAll, self.bcoAge815All, self.bcoDep815All,self.modeltime,self.TWrite)
        write_nospin_LIZ(self.c['resultsFolder'], self.LIZAgeAll, self.LIZDepAll,self.modeltime,self.TWrite)
        write_nospin_DIP(self.c['resultsFolder'], self.intPhiAll,self.modeltime,self.TWrite)

    # def define_parameters(self):
    #     '''
    #     Return the parameters used in the model, for the initialization as well as the time evolution

    #     :returns gridLen: size of grid used in the model run
    #                      (unit: number of boxes, type: int)
    #     :returns dx: vector of width of each box, used for stress calculations
    #                 (unit: ???, type: array of ints)
    #     :returns dt: number of seconds per time step
    #                 (unit: seconds, type: float)
    #     :returns t: number of years per time step
    #                (unit: years, type: float)
    #     :returns modeltime: linearly spaced time vector from indicated start year to indicated end year
    #                        (unit: years, type: array of floats)
    #     :returns years: total number of years in the model run
    #                    (unit: years, type: float)
    #     :returns stp: total number of steps in the model run
    #                  (unit: number of steps, type: int)
    #     :returns T_mean: interpolated temperature vector based on the model time and the initial user temperature data
    #                     (unit: ???, type: array of floats)
    #     :returns Ts: interpolated temperature vector based on the model time & the initial user temperature data
    #                 may have a seasonal signal imposed depending on number of years per time step (< 1)
    #                 (unit: ???, type: array of floats)
    #     :returns bdot: bdot is meters of ice equivalent/year. multiply by 0.917 for W.E. or 917.0 for kg/year
    #                   (unit: ???, type: )
    #     :returns bdotSec: accumulation rate vector at each time step
    #                      (unit: ???, type: array of floats)
    #     :returns rhos0: surface accumulate rate vector
    #                    (unit: ???, type: array of floats)
    #     :returns D_surf: diffusivity tracker
    #                     (unit: ???, type: array of floats)
    #     '''

    #     gridLen = np.size(self.z)
    #     dx = np.ones(gridLen)

    #     input_temp, input_year_temp = read_temp(self.c['InputFileNameTemp'])
    #     input_bdot, input_year_bdot = read_bdot(self.c['InputFileNamebdot'])

    #     yr_start  = max(input_year_temp[0], input_year_bdot[0])
    #     yr_end    = min(input_year_temp[-1], input_year_bdot[-1])
        
    #     years     = (yr_end - yr_start) * 1.0
    #     stp       = int(years * self.c['stpsPerYear'])
    #     modeltime = np.linspace(yr_start, yr_end, stp + 1)

    #     dt = years * S_PER_YEAR / stp
    #     t  = 1.0 / self.c['stpsPerYear']

    #     TPeriod = years
    #     Ts = np.interp(modeltime, input_year_temp, input_temp)
    #     T_mean = Ts
    #     if t < 1.0:
    #         Ts = Ts + self.c['TAmp'] * (np.cos(2 * np.pi * np.linspace(0, TPeriod, stp + 1)) + 0.3 * np.cos(4 * np.pi * np.linspace(0, TPeriod, stp + 1)))

    #     bdot = np.interp(modeltime, input_year_bdot, input_bdot)
    #     bdotSec = bdot / S_PER_YEAR / (stp / years)

    #     rhos0 = self.c['rhos0'] * np.ones(stp)
    #     D_surf = self.c['D_surf'] * np.ones(stp)

    #     return gridLen, dx, dt, t, modeltime, years, stp, Ts, T_mean, bdot, bdotSec, rhos0, D_surf

    def update_BCO(self):
        '''
        Updates the bubble close-off depth and age based on the Martinerie criteria as well as through assuming the critical density is 815 kg/m^3
        '''

        bcoMartRho = 1 / (1 / (917.0) + self.diffu.T10m * 6.95E-7 - 4.3e-5)  # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
        bcoAgeMart = min(self.age[self.rho >= bcoMartRho]) / S_PER_YEAR  # close-off age from Martinerie
        bcoDepMart = min(self.z[self.rho >= (bcoMartRho)])
        self.bcoAgeMartAll.append(bcoAgeMart)  # age at the 815 density horizon
        self.bcoDepMartAll.append(bcoDepMart)  # this is the 815 close off depth

        # bubble close-off age and depth assuming rho_crit = 815kg/m^3
        bcoAge815 = min(self.age[self.rho >= (RHO_2)]) / S_PER_YEAR  # close-off age where rho = 815 kg m^-3
        bcoDep815 = min(self.z[self.rho >= (RHO_2)])
        self.bcoAge815All.append(bcoAge815)  # age at the 815 density horizon
        self.bcoDep815All.append(bcoDep815)  # this is the 815 close off depth

    def update_LIZ(self):
        '''
        Updates the lock-in zone depth and age
        '''

        bcoMartRho = 1 / (1 / (917.0) + self.diffu.T10m * 6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
        LIZMartRho = bcoMartRho - 14.0  # LIZ depth (Blunier and Schwander, 2000)
        self.LIZAgeMart = min(self.age[self.rho > LIZMartRho]) / S_PER_YEAR  # lock-in age
        self.LIZDepMart = min(self.z[self.rho >= (LIZMartRho)])  # lock in depth
        self.LIZAgeAll.append(self.LIZAgeMart)
        self.LIZDepAll.append(self.LIZDepMart)

    def update_DIP(self):
        '''
        Updates the depth-integrated porosity
        '''

        bcoMartRho = 1 / (1 / (917.0) + self.diffu.T10m * 6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
        phi = 1 - self.rho / RHO_I  # total porosity
        phi[phi <= 0] = 1e-16
        phiC = 1 - bcoMartRho / RHO_I;  # porosity at close off
        phiClosed = 0.37 * phi * (phi / phiC) ** -7.6  # Closed porosity, from Goujon. See Buizert thesis (eq. 2.3) as well

        phiOpen = phi - phiClosed  # open porosity
        phiOpen[phiOpen <= 0] = 1.e-10  # don't want negative porosity.

        intPhi = np.sum(phi * self.dz)  # depth-integrated porosity
        self.intPhiAll.append(intPhi)
