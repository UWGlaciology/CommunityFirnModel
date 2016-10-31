from diffusion import Diffusion
from hl_analytic import hl_analytic
from reader import read_temp
from reader import read_bdot
from writer import write_spin
from physics import *
from constants import *
import numpy as np
#  FOR DEBUGGING
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


class FirnDensitySpin:
    def __init__(self, configName):
        '''
        Sets up the initial spatial grid, time grid, accumulation rate, age, density, mass, stress, and temperature of the model run
        :param configName: name of json config file containing model configurations
        '''

        # load in json config file and parses the user inputs to a dictionary
        with open(configName, "r") as f:
            jsonString = f.read()
            self.c = json.loads(jsonString)

        # create directory to store results
        if os.path.exists(self.c['resultsFolder']):
            rmtree(self.c['resultsFolder'])
        os.makedirs(self.c['resultsFolder'])

        # read in initial temperatures and accumulation rates
        input_temp, input_year_temp = read_temp(self.c['InputFileNameTemp'])
        input_bdot, input_year_bdot = read_bdot(self.c['InputFileNamebdot'])

        # load in model parameters
        self.bdot0 = input_bdot[0]
        self.temp0 = input_temp[0]

        self.gridLen = int((self.c['H'] - self.c['HbaseSpin']) / (self.bdot0 / self.c['stpsPerYearSpin']))
        gridHeight = np.linspace(self.c['H'], self.c['HbaseSpin'], self.gridLen)
        self.z  = self.c['H'] - gridHeight
        THL = input_temp[0]
        AHL = input_bdot[0]
        self.age, self.rho = hl_analytic(self.c['rhos0'], self.z, THL, AHL)
        dx, dt, t, stp, Ts, T_mean, bdotSec, rhos0 = self.define_parameters()

        # set up model grid
        # moved up code --> line 54, 55
        #gridheight = np.linspace(c['H'], c['HbaseSpin'], gridLen)
        #self.z  = self.c['H'] - gridHeight
        self.dz = np.diff(self.z)
        self.dz = np.append(self.dz, self.dz[-1])

        # set up the initial age and density of the firn column using herron & langway analytic
        ### moved up the code --> line 56, 57, 58
        #THL = input_temp[0]
        #AHL = input_bdot[0]
        #self.age, self.rho = hl_analytic(self.c['rhos0'], self.z, THL, AHL)

        # set up initial mass, stress, and mean accumulation rate
        self.mass  = self.rho * self.dz
        self.sigma = self.mass * dx * GRAVITY
        self.sigma = self.sigma.cumsum(axis = 0)
        self.mass_sum = self.mass.cumsum(axis = 0)
        self.bdot_mean = np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] / (self.age[1:] * RHO_I / t)))

        # set up initial temperature grid as well as a class to handle heat/isotope diffusion
        init_Tz = input_temp[0] * np.ones(self.gridLen)
        self.diffu = Diffusion(self.z, stp, self.gridLen, init_Tz)

        # set up initial grain growth (if specified in config file)
        if self.c['physGrain']:
            if self.c['calcGrainSize']:
                r02 = -2.42e-9 * (self.c['Ts0']) + 9.46e-7
                self.r2 = r02 * np.ones(self.gridLen)
            else:
                self.r2 = np.linspace(self.c['r2s0'], (6 * self.c['r2s0']), self.gridLen)

    def time_evolve(self):
        '''
        Evolve the spatial grid, time grid, accumulation rate, age, density, mass, stress, and temperature through time
        based on the user specified number of timesteps in the model run. Updates the firn density using a user specified 
        '''
        
        # load in model parameters
        dx, dt, t, stp, Ts, T_mean, bdotSec, rhos0 = self.define_parameters()
        steps = 1 / t #this is time steps per year
        if not self.c['physGrain']:
            r2_time = None

        for iter in xrange(stp):
#             print 'steps=',steps
#             print 'gridlen=',self.gridLen
#             print 'bdotsec=',len(bdotSec)
#             print 'diffu.Tz=',len(self.diffu.Tz)
#             print 'rho=',len(self.rho)
#             print 'z110=',len(self.z)
            
            # getting the right physics for firn density based on user input
            physics = {
                'HLdynamic':       HL_dynamic,
                'HLSigfus':        HL_Sigfus,
                'Barnola1991':     Barnola_1991,
                'Li2004':          Li_2004,
                'Li2011':          Li_2011,
                'Ligtenberg2011':  Ligtenberg_2011,
                'Arthern2010S':    Arthern_2010S,
                'Simonsen2013':    Simonsen_2013,
                'Morris2013':      Morris_HL_2013,
                'Helsen2008':      Helsen_2008,
                'Arthern2010T':    Arthern_2010T,
                'Spencer2001':     Spencer_2001,
                'Goujon2003':      Goujon_2003,
            }

            parameters = {
                'HLdynamic':      [iter, steps, self.gridLen, bdotSec, self.diffu.Tz, self.rho],
                'HLSigfus':       [iter, steps, self.gridLen, bdotSec, self.diffu.Tz, self.rho, self.sigma],
                'Barnola1991':    [iter, steps, self.gridLen, bdotSec, self.diffu.Tz, self.rho, self.sigma],
                'Li2004':         [iter, steps, self.gridLen, bdotSec, T_mean, self.rho],
                'Li2011':         [iter, steps, self.gridLen, bdotSec, self.bdot_mean, self.c['bdot_type'], self.diffu.Tz, T_mean, self.rho],
                'Ligtenberg2011': [iter, steps, self.gridLen, bdotSec, self.bdot_mean, self.c['bdot_type'], self.diffu.Tz, T_mean, self.rho],
                'Arthern2010S':   [iter, steps, self.gridLen, bdotSec, self.bdot_mean, self.c['bdot_type'], self.diffu.Tz, T_mean, self.rho],
                'Simonsen2013':   [iter, steps, self.gridLen, bdotSec, self.bdot_mean, self.c['bdot_type'], self.diffu.Tz, T_mean, self.rho],
                'Morris2013':     [iter, steps, self.gridLen, self.diffu.Tz, dt, self.rho, True],
                'Helsen2008':     [iter, steps, bdotSec, self.c['bdot_type'], self.bdot_mean, self.diffu.Tz, Ts, self.rho],
                'Arthern2010T':   [iter, self.gridLen, self.diffu.Tz, self.rho, self.sigma, self.r2, self.c['physGrain']],
                'Spencer2001':    [],
                'Goujon2003':     [],
            }
            try:
                drho_dt = physics[self.c['physRho']](*parameters[self.c['physRho']])
            except KeyError:
                default()

            # update density and age of firn
            self.age = np.concatenate(([0], self.age[:-1])) + dt
            self.rho = self.rho + dt * drho_dt
            self.rho  = np.concatenate(([rhos0[iter]], self.rho[:-1]))

            # update temperature grid and isotope grid if user specifies
            if self.c['heatDiff']:
                self.diffu.heatDiff(self.z, self.dz, Ts[iter], self.rho, dt)
            if self.c['isoDiff']:
                self.diffu.isoDiff(iter, self.z, self.dz, self.rho, self.c['iso'], self.gridLen, dt)

            # update model grid
            dzNew = bdotSec[iter] * RHO_I / rhos0[iter] * S_PER_YEAR
            self.dz = self.mass / self.rho * dx
            self.dz = np.concatenate(([dzNew], self.dz[:-1]))
            self.z = self.dz.cumsum(axis = 0)
            self.z = np.concatenate(([0], self.z[:-1]))

            # update mass, stress, and mean accumulation rate
            massNew = bdotSec[iter] * S_PER_YEAR * RHO_I
            self.mass = np.concatenate(([massNew], self.mass[:-1]))
            self.sigma = self.mass * dx * GRAVITY
            self.sigma = self.sigma.cumsum(axis = 0)
            self.mass_sum  = self.mass.cumsum(axis = 0)
            self.bdot_mean = np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] * t / (self.age[1:] * RHO_I)))

            # update grain radius
            if self.c['physGrain']:
                self.r2 = grainGrowth(self.diffu.Tz, Ts, iter, dt, self.r2, self.c['calcGrainSize'])

            # write results at the end of the time evolution
            if (iter == (stp - 1)):
                rho_time = np.concatenate(([t * iter + 1], self.rho))
                Tz_time  = np.concatenate(([t * iter + 1], self.diffu.Tz))
                age_time = np.concatenate(([t * iter + 1], self.age))
                z_time   = np.concatenate(([t * iter + 1], self.z))
                if self.c['physGrain']:
                    r2_time = np.concatenate(([t * iter + 1], self.r2))

                write_spin(self.c['resultsFolder'], self.c['physGrain'], rho_time, Tz_time, age_time, z_time, r2_time)

    def define_parameters(self):
        '''
        Return the parameters used in the model, for the initialization as well as the time evolution

        :returns gridLen: size of grid used in the model run
                         (unit: number of boxes, type: int)
        :returns dx: vector of width of each box, used for stress calculations
                    (unit: ???, type: array of ints)
        :returns dt: number of seconds per time step
                    (unit: seconds, type: float)
        :returns t: number of years per time step
                   (unit: years, type: float)
        :returns stp: total number of steps in the model run
                     (unit: number of steps, type: int)
        :returns T_mean: interpolated temperature vector based on the model time and the initial user temperature data
                        (unit: ???, type: array of floats)
        :returns Ts: interpolated temperature vector based on the model time & the initial user temperature data
                    may have a seasonal signal imposed depending on number of years per time step (< 1)
                    (unit: ???, type: array of floats)
        :returns bdotSec: accumulation rate vector at each time step
                         (unit: ???, type: array of floats)
        :returns rhos0: surface accumulate rate vector
                       (unit: ???, type: array of floats)
        '''

        # TOOK THIS AND MOVED TO INITIALIZATION... NEED IT EARLIER
        #gridLen = int((self.c['H'] - self.c['HbaseSpin']) / (self.bdot0 / self.c['stpsPerYearSpin']))
        dx = np.ones(self.gridLen)

        try:
            zz = np.min(self.z[self.rho > 850.0])
        except ValueError:
            pass
        years = int(zz / self.bdot0)
        stp = int(years * self.c['stpsPerYearSpin'])
        dt = years * S_PER_YEAR / stp
        t  = 1.0 / self.c['stpsPerYearSpin']

        if self.c['stpsPerYearSpin'] or (stp / years) >= 1.:
            Ts = self.temp0 * np.ones(stp)
        else:
            TPeriod = self.c['yearSpin']
            Ts = self.temp0 + self.c['TAmp'] * (np.cos(2 * np.pi * np.linspace(0, TPeriod, stp)) + 0.3 * np.cos(4 * np.pi * np.linspace(0, TPeriod, stp)))
        T_mean = Ts

        bdotSec0 = self.bdot0 / S_PER_YEAR / self.c['stpsPerYearSpin']
        bdotSec = bdotSec0 * np.ones(stp)

        rhos0 = self.c['rhos0'] * np.ones(stp)

        return dx, dt, t, stp, Ts, T_mean, bdotSec, rhos0