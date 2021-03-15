#!/usr/bin/env python

'''
firn_density_nospin.py
======================
The script for the transient model run.

Copyright © 2021 C. Max Stevens

Distributed under terms of the MIT license. 
'''

from diffusion import *
from reader import read_input
from reader import read_init
from writer import write_spin_hdf5
from writer import write_nospin_hdf5
from writer import SpinUpdate
from physics import *
from constants import *
from melt import *
from isotopeDiffusion import isotopeDiffusion
from firn_density_spin import FirnDensitySpin
import numpy as np
import csv
import json
import sys
import math
from shutil import rmtree
import os
import shutil
import time
import inspect
import h5py
import scipy.interpolate as interpolate
from firn_air import FirnAir
from regrid import *
try:
    import pandas as pd
except:
    print('You do not have the pandas python package installed.')
    print('It is used to create a running mean temperature.')
from merge import mergeall #VV
from merge import mergesurf #VV
from merge import mergenotsurf #VV
from re_snowpack import resingledomain #VV
from prefflow_snowpack import prefflow #VV
from sublim import sublim #VV
from ModelOutputs import ModelOutputs

class FirnDensityNoSpin:
    '''
    Class for the main, transient model run. 

    Parameters
    ----------
    : gridLen: size of grid used in the model run
                (unit: number of boxes, type: int)
    : dx: vector of width of each box, used for stress calculations
                (unit: m, type: array of ints)
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
                (unit: K, type: array of floats)
    : Ts: interpolated temperature vector based on the model time & the initial user temperature data
                may have a seasonal signal imposed depending on number of years per time step (< 1)
                (unit: K, type: array of floats)
    : bdot: bdot is meters of ice equivalent/year. multiply by 0.917 for W.E. or 917.0 for kg/year
                (unit: m ice eq. per year, type: array of floats)
    : bdotSec: accumulation rate vector at each time step
                (unit: m ice eq. per second, type: array of floats)
    : rhos0: surface accumulate rate vector
                (unit: kg m^-3, type: array of floats)
    : bdot_mean: mean accumulation over the lifetime of each parcel
                (units are m I.E. per year)
                
    :returns D_surf: diffusivity tracker
                (unit: ???, type: array of floats)
    '''

    def __init__(self, configName, climateTS = None, NewSpin = False):
        '''
        Sets up the initial spatial grid, time grid, accumulation rate, age, 
        density, mass, stress, temperature, and diffusivity of the model run
        
        :param configName: name of json config file containing model configurations       
        '''

        ### load in json config file and parses the user inputs to a dictionary
        
        with open(configName, "r") as f:
            jsonString      = f.read()
            self.c          = json.loads(jsonString)

        spinner = os.path.exists(os.path.join(self.c['resultsFolder'], self.c['spinFileName']))
        
        if ((not spinner) or NewSpin):
            if self.c['timesetup']=='exact':
                if climateTS != None:
                    self.c['stpsPerYear'] = 1/np.mean(np.diff(climateTS['time']))
                    print('stepsperyear:', self.c['stpsPerYear'])
                else:
                    input_bdot, input_year_bdot = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamebdot']))
                    self.c['stpsPerYear'] = 1/np.mean(np.diff(input_year_bdot))

            firnS = FirnDensitySpin(self.c, climateTS = climateTS)
            firnS.time_evolve()
        else:
            pass

        print("Main run starting")
        print("physics are", self.c['physRho'])

        ### read in initial depth, age, density, temperature from spin-up results
        initDepth   = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'depthSpin')
        initAge     = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'ageSpin')
        initDensity = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'densitySpin')
        initTemp    = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'tempSpin')

        try: #VV for reading initial lwc from the spin up file
            initLWC = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'LWCSpin')
            print('Initial LWC provided by spin-up')
        except: 
            pass 

        try:
            self.doublegrid = self.c['doublegrid']
            if self.doublegrid:
                initGrid = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'gridSpin')
                self.gridtrack = initGrid[1:]
                self.nodestocombine = self.c['nodestocombine']

        except:
            self.doublegrid = False
            print('you should add "doublegrid" to the json')

        ### set up the initial age and density of the firn column
        self.age        = initAge[1:]
        self.rho        = initDensity[1:]

        ### set up model grid
        self.z          = initDepth[1:]
        self.dz         = np.diff(self.z)
        self.dz         = np.append(self.dz, self.dz[-1])
        self.gridLen    = np.size(self.z)
        self.dx         = np.ones(self.gridLen)

        ### Feature to update the spin file to not have to repeat a long spin up
        if 'spinUpdate' not in self.c:
            self.c['spinUpdate'] = False
        if self.c['spinUpdate']:
            updatedStartDate = initDepth[0] # if the spin file has been updated, this is the date to start the run (find this time in the forcing data)
        else:
            updatedStartDate = None

        ### get temperature and accumulation rate from input csv file
        try:
            if self.c['manualT']:
                if self.c['timesetup']!= 'exact':
                    print('"timesetup" in .json must be "exact" to use manualT. Exiting.')
                    sys.exit()
                if self.c['heatDiff']:
                    print('"heatDiff" should be false for manualT runs. Run will continue.')
        except:
            # print('"manualT" not in .json; run will continue with original model behavior')
            self.c['manualT'] = False 

        if self.c['manualT']:
            bigTmat = np.loadtxt(os.path.join(self.c['InputFileFolder'],self.c['ManualTFilename']),delimiter = ',')
            self.manualT_time = bigTmat[0,1:]
            self.manualT_dep = bigTmat[1:,0]
            self.manualT_temp = bigTmat[1:,2:] # use 2: because the 1 (first T column) is the init profile
            if self.manualT_temp[0,0]<=0:
                self.manualT_temp = self.manualT_temp+273.15
            input_temp = bigTmat[1,1:] # surface
            input_year_temp = self.manualT_time # time again, just to gel with other code
            init_Tz = bigTmat[1:,1] # temp at zeroeth time step (i.e. at init)

        else:
            if climateTS != None:
                start_ind = np.where(climateTS['time']>=updatedStartDate)[0][0]
                input_temp = climateTS['TSKIN'][start_ind:]
                input_year_temp = climateTS['time'][start_ind:]
                
            else:
                input_temp, input_year_temp = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNameTemp']), updatedStartDate)
            if input_temp[0] < 0.0:
                input_temp      = input_temp + K_TO_C
            input_temp[input_temp>T_MELT] = T_MELT
        #####################

        ### bdot ############
        if climateTS != None:
            input_bdot = climateTS['BDOT'][start_ind:]
            input_year_bdot = climateTS['time'][start_ind:]        
        else:
            input_bdot, input_year_bdot = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamebdot']), updatedStartDate)
        #####################

        #####################
        ### MELT ##############
        if 'MELT' not in self.c:
            print('You should add "MELT" to your .json (True/False)')
            self.c['MELT']      = False
            input_snowmelt      = None
            input_year_snowmelt = None
            self.LWC            = np.zeros_like(self.z)
            self.PLWC_mem       = np.zeros_like(self.z) #VV keep track of water content that was in PFdom
            self.raininput      = False #VV no rain input

        if self.c['MELT']:
            if ((climateTS != None) and ('SMELT' in climateTS.keys())):
                input_snowmelt = climateTS['SMELT'][start_ind:]
                input_year_snowmelt = climateTS['time'][start_ind:]
            else:
                input_snowmelt, input_year_snowmelt = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamemelt']), updatedStartDate)
            self.MELT           = True
            try: 
                self.LWC        = initLWC[1:]
            except:
                self.LWC        = np.zeros_like(self.z)
            self.PLWC_mem   = np.zeros_like(self.z) #VV keep track of water content that was in PFdom

            if 'RAIN' not in self.c:
                self.c['RAIN'] = False
            if self.c['RAIN']:
                if ((climateTS != None) and ('RAIN' in climateTS.keys())):
                    input_rain = climateTS['RAIN'][start_ind:]
                    input_year_rain = climateTS['time'][start_ind:]
                else:
                    input_rain, input_year_rain = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamerain']), updatedStartDate)
            
            if 'liquid' not in self.c:
                print('Melt is on, but you did not specify which perolation scheme in the .json')
                print('Defaulting to original CFM bucket scheme')
                self.c['liquid'] = 'bucketVV'
        else:
            self.MELT           = False
            # print("Melt is not turned on.")
            input_snowmelt      = None
            input_year_snowmelt = None
            self.LWC            = np.zeros_like(self.z)
            self.PLWC_mem       = np.zeros_like(self.z)
            self.c['RAIN'] = False
        #####################

        #####################
        ### time ############
        if 'timesetup' not in self.c:
            # print('You should add "timesetup" to the .json')
            ### could add details here
            self.c['timesetup']='interp'

        # year to start and end, from the input file. If inputs have different start/finish, take only the overlapping times
        if (self.c['timesetup']=='interp' or self.c['SeasonalTcycle']):
            if (self.c['SeasonalTcycle'] and self.c['timesetup']=='exact'):
                print('"Exact" time setup does not work with "SeasonalTcycle" Switching to "interp"')
            yr_start        = max(input_year_temp[0], input_year_bdot[0])   # start year
            yr_end          = min(input_year_temp[-1], input_year_bdot[-1]) # end year
            
            # self.years      = np.ceil((yr_end - yr_start) * 1.0)
            self.years      = (yr_end - yr_start) * 1.0
            self.dt         = S_PER_YEAR / self.c['stpsPerYear'] # seconds per time step
            self.stp        = int(self.years * S_PER_YEAR/self.dt)#-1       # total number of time steps, as integer
            self.dt = self.dt * np.ones(self.stp)
            # self.modeltime  = np.linspace(yr_start, yr_end, self.stp + 1)   # vector of time of each model step
            self.modeltime  = np.linspace(yr_start, yr_end, self.stp)
            # self.t          = 1.0 / self.c['stpsPerYear']                   # years per time step
            self.t = (1.0 / self.c['stpsPerYear']) * np.ones_like(self.dt)
            init_time = -9999.0


        elif self.c['timesetup']=='exact':
            # print('"Exact" time setup will not work properly if input forcing does not all have the same time')
            # yr_start = input_year_temp[0] #previously had [1] - error? 20/3/3
            # yr_end = input_year_temp[-1]
            self.dt = np.diff(input_year_temp)*S_PER_YEAR #[units s]
            # self.dt = np.append(self.dt,self.dt[-1])
            self.stp = len(self.dt)
            self.modeltime = input_year_temp[1:] # this offset because use diff above
            # self.modeltime = input_year_temp[0:-1]
            yr_start = self.modeltime[0]
            yr_end = self.modeltime[-1]
            # self.t = np.mean(np.diff(input_year_temp))
            self.t = np.diff(input_year_temp)
            init_time = input_year_temp[0]

        elif self.c['timesetup']=='retmip': #VV retmip experiments require to match perfectly their 3h time step
            # might be able to just use 'exact'?
            self.years      = (yr_end - yr_start) #VV
            self.dt         = 10800. #VV 3 hours
            self.stp        = len(input_temp)
            self.modeltime = np.zeros(self.stp)
            self.modeltime[0] = yr_start
            if (int(yr_start)%4 == 0):
                leap=1 #leap year
            elif (int(yr_start)%4 != 1):
                leap=0 #normal year
            for ii in range(1,len(self.modeltime)):
                if leap == 1:
                    self.modeltime[ii] = self.modeltime[ii-1]+self.dt/31622400
                    if (int(self.modeltime[ii]+self.dt/31622400)-int(self.modeltime[ii-1]))==1:
                        # At transition between two years, set time exactly at new year (this avoids propagation of small errors)
                        self.modeltime[ii] = int(self.modeltime[ii-1])+1
                elif leap == 0:
                    self.modeltime[ii] = self.modeltime[ii-1]+self.dt/31536000
                    if (int(self.modeltime[ii]+self.dt/31536000)-int(self.modeltime[ii-1]))==1:
                        # At transition between two years, set time exactly at new year (this avoids propagation of small errors)
                       self.modeltime[ii] = int(self.modeltime[ii-1])+1
                if int(self.modeltime[ii])%4 == 0:
                    leap = 1
                elif int(self.modeltime[ii])%4 != 0:
                    leap = 0
            init_time = -9999
            ### next two lines - Vincent's code included, but may have been a bug?
            # self.modeltime  = np.linspace(yr_start, yr_end, self.stp)
            # self.t          = 1.0 / self.c['stpsPerYear']                   # years per time step
        #####################
      
        ###############################
        ### surface boundary conditions
        ### temperature, accumulation, melt, isotopes, surface density
        ###############################
        int_type            = self.c['int_type']
        # print('Climate interpolation method is %s' %int_type)

        ### Temperature #####       
        Tsf                 = interpolate.interp1d(input_year_temp,input_temp,int_type,fill_value='extrapolate') # interpolation function
        self.Ts             = Tsf(self.modeltime) # surface temperature interpolated to model time
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

        if 'conductivity' not in self.c:
            self.c['conductivity'] = 'Calonne2019'
            print('conductivity was not set in .json; Defaulting to Calonne2019')
        #####################

        ### Accumulation ####
        bsf                 = interpolate.interp1d(input_year_bdot,input_bdot,int_type,fill_value='extrapolate') # interpolation function
        self.bdot           = bsf(self.modeltime) # m ice equivalent per year
        self.bdot[self.bdot<1e-6] = 0.0
        
        self.bdotSec        = self.bdot / S_PER_YEAR / (S_PER_YEAR/self.dt) # accumulation at each time step (meters i.e. per second). gets multiplied by S_PER_YEAR later. (sort of hacky, I know)
        
        #old way:    
        # self.bdotSec        = self.bdot / S_PER_YEAR / self.c['stpsPerYear'] # accumulation for each time step (meters i.e. per second)

        try: #Rolling mean average surface temperature and accumulation rate (vector)
        # (i.e. the long-term average climate)
            Nyears = 10 #number of years to average for T_mean
            NN = int(np.mean(S_PER_YEAR/self.dt)*Nyears)
            # NN = int(self.c['stpsPerYear']*Nyears)
            self.T_mean = pd.Series(self.Ts).rolling(window=NN+1,win_type='hamming').mean().values
            self.T_mean[np.isnan(self.T_mean)] = self.T_mean[NN]
            self.bdot_av = pd.Series(self.bdot).rolling(window=NN+1,win_type='hamming').mean().values
            self.bdot_av[np.isnan(self.bdot_av)] = self.bdot_av[NN]
        except Exception:
            self.T_mean = np.mean(self.Ts) * np.ones(self.stp)
            self.bdot_av = np.mean(self.bdot) * np.ones(self.stp)
            print('Error calculating T_mean, using mean surface over all time')

        if self.c['manual_climate']: #in the case of very short runs, you want to set the longer-term climate manually
            self.T_mean = self.c['deepT'] * np.ones(self.stp)
            self.bdot_av = self.c['bdot_long'] * np.ones(self.stp)

        try: 
            if self.c['manual_iceout']:
                self.iceout = self.c['iceout']
                print('Ensure that your iceout value has units m ice eq. per year!')
            else:
                self.iceout = np.mean(self.bdot) # this is the rate of ice flow advecting out of the column, units m I.E. per year.

        except Exception:
            print('add field "manual_iceout" to .json file to set iceout value manually')
            self.iceout = np.mean(self.bdot) # this is the rate of ice flow advecting out of the column, units m I.E. per year.

        self.w_firn         = np.mean(self.bdot) * RHO_I / self.rho 

        if (np.any(self.bdotSec<0.0) and self.c['bdot_type']=='instant'):
            print('ERROR: bdot_type set to "instant" in .json and input')
            print('accumulation has at least one negative value.') 
            print('QUITTING MODEL RUN.')
            sys.exit()
        #####################
        
        ### Melt and rain inputs ############
        if self.MELT:
            ssf                 = interpolate.interp1d(input_year_snowmelt,input_snowmelt,int_type,fill_value='extrapolate')
            self.snowmelt       = ssf(self.modeltime)
            self.snowmeltSec    = self.snowmelt / S_PER_YEAR / (S_PER_YEAR/self.dt) # melt for each time step (meters i.e. per second)
            self.c['LWCheat'] = 'enthalpy' # Filler for future testing.

            if self.c['RAIN'] == True: ##VV use rain climatic input
                rsf             = interpolate.interp1d(input_year_rain,input_rain,int_type,fill_value='extrapolate')
                self.rain       = rsf(self.modeltime) # [mIE/yr]
                self.rainSec    = self.rain / S_PER_YEAR / (S_PER_YEAR/self.dt) # rain for each time step (mIE/s)
            else:
                self.rainSec    = np.zeros(self.stp) #VV to avoid problem in the conditions to call for liquid water routine
        #####################
 
        ### Surface Density ####
        if self.c['variable_srho']:
            if self.c['srho_type']=='userinput':
                input_srho, input_year_srho = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamerho']), updatedStartDate)
                Rsf             = interpolate.interp1d(input_year_srho,input_srho,int_type,fill_value='extrapolate') # interpolation function
                self.rhos0      = Rsf(self.modeltime) # surface temperature interpolated to model time
            elif self.c['srho_type']=='param':
                self.rhos0      = 481.0 + 4.834 * (self.T_av - T_MELT) # Kuipers Munneke, 2015
            elif self.c['srho_type']=='noise':
                rho_stdv        = 25 # the standard deviation of the surface density (I made up 25)
                self.rhos0      = np.random.normal(self.c['rhos0'], rho_stdv, self.stp)
                self.rhos0[self.rhos0>450]=450
                self.rhos0[self.rhos0<250]=250

        else:
            self.rhos0      = self.c['rhos0'] * np.ones(self.stp)       # density at surface
        #####################

        #####
        try:
            if self.c['no_densification']:
                print('CAUTION: densification if OFF!')
            else:
                pass
        except:
            self.c['no_densification']=False

        #####################
        if 'rad_pen' not in self.c:
            self.c['rad_pen']=False
            
        if self.c['rad_pen']: 
            print('RADIATION PENETRATION STILL IN DEVELOPMENT STAGE')   
            input_rad, input_year_rad = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNameSWnetrad']), updatedStartDate)
            Radsf             = interpolate.interp1d(input_year_rad,input_rad,int_type,fill_value='extrapolate')
            self.E_sw = Radsf(self.modeltime)
        #####################        

        ### Layer tracker ###
        self.D_surf     = self.c['D_surf'] * np.ones(self.stp)      # layer traking routine (time vector). 
        self.Dcon       = self.c['D_surf'] * np.zeros(self.gridLen)  # layer tracking routine (initial depth vector)
        self.Dcon = np.flipud(np.arange((-1*(len(self.Dcon))),0))
        #####################

        ###############################
        ### set up vector of times data will be written
        try:
            Tind                = np.nonzero(self.modeltime>=self.c['TWriteStart'])[0][0]
        except Exception:
            Tind                = np.nonzero(self.modeltime>=1980.0)[0][0] # old way
            print('You should add new variable "TWriteStart" to the json file!')
            print('Arbitrarily starting writing at 1980.0. See line 243 in firn_density_nospin.py')
        self.TWrite         = self.modeltime[Tind::self.c['TWriteInt']]
        if self.TWrite[-1]!=self.modeltime[-1]:
            print('Adding last time step to Twrite')
            self.TWrite = np.append(self.TWrite,self.modeltime[-1])
        if self.c['TWriteInt']!=1:
            print('Time writing interval is not 1; dH output will not be accurate.')
        TWlen               = len(self.TWrite) #- 1
        self.WTracker       = 1

        ### you can choose to write certain fields at different times.
        ### (Functionality not fully tested) 
        # Tind2                 = np.nonzero(self.modeltime>=1958.0)[0][0]
        # self.TWrite2      = self.modeltime[Tind2::self.c['TWriteInt']]
        # self.TWrite       = np.append(self.modeltime[10],self.TWrite)
        # self.TWrite2      = self.modeltime[0::self.c['TWriteInt']]
        # self.TWrite_out   = self.TWrite
        # TWlen2            = len(self.TWrite2) #- 1
        # self.WTracker2        = 1
        ###############################

        ### set up initial mass, stress, and mean accumulation rate
        self.mass           = self.rho * self.dz
        self.sigma          = (self.mass + (self.LWC * RHO_W_KGM)) * self.dx * GRAVITY
        self.sigma          = self.sigma.cumsum(axis = 0)
        self.mass_sum       = self.mass.cumsum(axis = 0)
        ### mean accumulation over the lifetime of the parcel:
        spinF = self.c['resultsFolder']+'/'+self.c['spinFileName']
        if 'bdot_meanSpin' in h5py.File(spinF,'r').keys():
            self.bdot_mean = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'bdot_meanSpin')[1:]
        else:
            self.bdot_mean      = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] / (self.age[1:] * RHO_I / np.mean(self.t)))))* np.mean(S_PER_YEAR/self.dt) *S_PER_YEAR
        ### It is the mass of the overlying firn divided by the age of the parcel.
        ### transform mass in meters ice equiv -> divide by age(in sec) [m/s] -> multiply by years per step and by steps per year (cancels) -> multiply by secperyear -> [mIE/yr]
        ### for surf layer -> mass in mIE is only multiplied by steps per year: if 1 stp/yr,mean acc is the mass of surf layer; if 2 stps/yr,mean acc is 2* what has been accumulated over the last step, etc.
        #######################

        ### set up longitudinal strain rate
        if self.c['strain']: # input units are yr^-1
            input_dudx, input_year_dudx = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamedudx']), updatedStartDate)
            dusf                        = interpolate.interp1d(input_year_dudx,input_dudx,int_type,fill_value='extrapolate')           
            self.du_dx      = dusf(self.modeltime)
            self.du_dxSec   = self.du_dx / S_PER_YEAR #/ self.c['stpsPerYear'] # strain rate (s^-1) at each time step
        #######################
        if self.c['manualT']:
            self.Tz = np.interp(self.z, self.manualT_dep, init_Tz)
        else:
            self.Tz             = initTemp[1:]

        self.T50            = np.mean(self.Tz[self.z<50])
        self.T10m           = self.Tz[np.where(self.z>=10.0)[0][0]]

        self.compboxes = len(self.z)
        #######################

        ### Additional melt outputs ###
        ### which of these should go to output?
        if self.MELT:
            self.runoff         = np.array([0.]) #total runoff, single value for the whole firn column
            self.refrozen       = np.zeros_like(self.dz) #vrefreezing in every layer, array of size of our grid
            self.totalrunoff    = np.array([0.]) # Might be useful to have a total final value without having to write every time step
            self.totalrefrozen  = np.zeros_like(self.dz) # Might be useful to have a total final value without having to write every time step
            self.totwatersublim = 0. #VV Total amount of liquid water that get sublimated
            self.lwcerror       = 0. #VV
            self.totallwcerror  = 0. #

        ### initial grain growth (if specified in config file)
        if self.c['physGrain']:
            initr2              = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'r2Spin')
            self.r2             = initr2[1:]
            r20                 = self.r2
            # self.dr2_dt         = np.zeros_like(self.z) # dr2_dt not currently set up    
        else:            
            self.r2             = None
            # self.dr2_dt         = None
        #######################

        ### temperature history for Morris physics
        if self.c['physRho'] == 'Morris2014':
            if 'QMorris' not in self.c:
                self.c['QMorris'] = 110.0e3
            
            self.THist          = True
            initHx              = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'HxSpin') 
            self.Hx             = initHx[1:]
        else:
            self.THist          = False
        #####################

        ### values for Goujon physics
        if self.c['physRho']=='Goujon2003':
            self.Gamma_Gou      = 0 
            self.Gamma_old_Gou  = 0
            self.Gamma_old2_Gou = 0
            self.ind1_old       = 0
        #######################

        ### Isotopes ########
        if self.c['isoDiff']:
            self.spin = False
            print('Isotope Diffusion is initialized')
            if 'site_pressure' not in self.c:
                print('site_pressure not in .json')
                print('Defaulting to 1013.25')

            self.Isotopes     = {} #dictionary of class instances
            self.iso_out      = {} # outputs for each isotope
            self.Isoz         = {} # depth profile of each isotope
            self.Iso_sig2_z   = {} # diffusion length profile
            self.iso_sig2_out = {}

            for isotope in self.c['iso']:
                self.Isotopes[isotope]  = isotopeDiffusion(self.spin,self.c,isotope,self.stp,self.z,self.modeltime)
                # self.iso_out[isotope]       = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
                # self.iso_out[isotope][0,:]  = np.append(init_time, self.Isotopes[isotope].del_z)               
                # self.iso_sig2_out[isotope]  = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
                # self.iso_sig2_out[isotope][0,:]   = np.append(init_time, self.Isotopes[isotope].iso_sig2_z)
        #######################

        #####################
        ##### Firn Air ######
        '''
        each gas of interest gets its own instance of the class, each instance
        is stored in a dictionary
        '''

        if self.c['FirnAir']:
            print('Firn air initialized')
            with open(self.c['AirConfigName'], "r") as f:
                jsonString    = f.read()
                self.cg       = json.loads(jsonString)
            self.FA           = {} # dictionary holding each instance of the firn-air class
            # self.gas_out    = {} # outputs for each gas in the simulation
            self.Gz           = {} # Depth profile of each gas
            self.diffusivity  = np.ones_like(self.rho)
            self.gas_age      = np.zeros_like(self.rho)
            self.w_air        = np.ones_like(self.rho)
            self.w_firn       = np.ones_like(self.rho)

            for gas in self.cg['gaschoice']:
                if (gas=='d15N2' or gas=='d40Ar'):
                    input_year_gas = input_year_temp
                    input_gas = np.ones_like(input_year_temp)
                else:
                    input_gas, input_year_gas = read_input(os.path.join(self.c['InputFileFolder'],'%s.csv' %gas), updatedStartDate)
                Gsf     = interpolate.interp1d(input_year_gas,input_gas,'linear',fill_value='extrapolate')
                Gs      = Gsf(self.modeltime)

                self.FA[gas] = FirnAir(self.cg, Gs, self.z, self.modeltime, self.Tz, self.rho, self.dz, gas, self.bdot)
                self.Gz[gas] = np.ones_like(self.rho)
            
            if self.cg['runtype']=='steady':
                print('Steady-state firn air works only with Herron and Langway physics, instant accumulation mode')
                print('This is automatically changed for you')
                self.bdot           = self.cg['steady_bdot']*np.ones_like(self.bdot)
                self.bdotSec        = self.bdot / S_PER_YEAR / self.c['stpsPerYear'] # accumulation for each time step (meters i.e. per second)
                self.iceout         = np.mean(self.bdot)  # units m I.E. per year.
                self.w_firn         = np.mean(self.bdot) * RHO_I / self.rho 
                self.c['physRho']   = 'HLdynamic'
                self.c['bdot_type'] = 'instant'
        else:
            self.cg = None
        ######################            

        ######################
        if 'merging' not in self.c:
            self.c['merging'] = False
        ######################

        ######################
        ### DIP, DHdt, BCO ###
        bcoAgeMart, bcoDepMart, bcoAge830, bcoDep830, LIZAgeMart, LIZDepMart, bcoAge815, bcoDep815  = self.update_BCO(0)

        if 'DIPhorizon' in self.c:
            self.DIPhorizon = self.c['DIPhorizon']
            if self.DIPhorizon > self.z[-1]:
                print('DIPhorizon is deeper than bottom of domain.')
                reset_HZ = np.floor(self.z[-1]*0.8)
                print('setting to 80% of domain depth ({} m)'.format(reset_HZ))
                self.DIPhorizon = reset_HZ
        else:
            self.DIPhorizon = np.floor(self.z[-1]*0.8)

        self.dHAll      = []
        self.dHAllcorr  = []
        intPhi, self.DIPc, z_co = self.update_DIP()
        ind_z = np.where(self.z>=self.DIPhorizon)[0][0]        
        self.dHAll.append(0)
        self.dHAllcorr.append(0)
        dHOut       = 0 # surface elevation change since last time step
        dHOutC      = 0 # cumulative surface elevation change since start of model run
        compOut     = 0 # compaction of just the firn at each time step; no ice dynamics or accumulation
        dHOutcorr   = 0
        dHOutcorrC  = 0
        DIPhz       = self.DIPhorizon #set the first element in DIPhz to be the horizon depth
        

        self.BCO = np.array([bcoAgeMart, bcoDepMart, bcoAge830, bcoDep830, LIZAgeMart, LIZDepMart, bcoAge815, bcoDep815, z_co])
        self.DIP = np.array([intPhi, dHOut, dHOutC, compOut, dHOutcorr, dHOutcorrC, DIPhz])
        #####################
        self.climate = np.array([self.bdot[0],self.Ts[0]])
        #####################

        ######################################
        ### MODEL OUTPUTS ####################
        ######################################
        if 'grid_outputs' not in self.c:
            self.c['grid_outputs'] = False
        rep_dict = {'density':'rho', 'temperature':'Tz', 'depth':'z', 'dcon':'Dcon', 'temp_Hx':'Hx'} #mapping names from output list to variables in CFM
        self.output_list = [rep_dict.get(n, n) for n in self.c['outputs']]
        if 'grainsize' in self.c['outputs']:
            self.output_list.remove('grainsize')
            self.output_list.append('r2')
            # self.output_list.append('dr2_dt')
        if 'meltoutputs' in self.c['outputs']: # Need to work with Vincent on what should be here.
            self.output_list.remove('meltoutputs')
            self.output_list.append('runoff') # This is cumulative for each time step
            self.output_list.append('refrozen') # Should this be in each grid point, or total?
              
        if ((not self.MELT) and ('LWC' in self.output_list)):
            self.output_list.remove('LWC')
        if ((not self.c['FirnAir']) and ('gasses' in self.output_list)):
            self.output_list.remove('gasses')
        if ((not self.c['isoDiff']) and ('isotopes' in self.output_list)):
            self.output_list.remove('isotopes')
        if ((self.c['grid_outputs']) and ('Dcon' in self.output_list)):
            print('Caution: Dcon with grid_outputs uses nearest interp')
        
        if 'compaction' in self.output_list:
            self.compaction          = np.zeros_like(self.z)
        if 'viscosity' in self.output_list:
            self.viscosity          = np.zeros(self.gridLen)      

        if self.c['FirnAir']:
            grep_dict = {'air_advection_rate': 'w_air', 'firn_advection_rate': 'w_firn',}
            gas_output_list = [grep_dict.get(n, n) for n in self.cg['outputs']]
            self.output_list.extend(gas_output_list)
            try:
                self.output_list.remove('gasses')
            except:
                pass

        MOd = {key:value for key, value in self.__dict__.items() if key in self.output_list}

        if self.c['FirnAir']:    
            for gas in self.cg['gaschoice']:
                self.output_list.append(gas)
                MOd[gas] = self.Gz[gas]

        if self.c['isoDiff']:
            for isotope in self.c['iso']:
                self.output_list.append('isotopes_{}'.format(isotope))
                MOd['isotopes_{}'.format(isotope)] = self.Isotopes[isotope].del_z
                self.output_list.append('iso_sig2_{}'.format(isotope)) 
                MOd['iso_sig2_{}'.format(isotope)] = self.Isotopes[isotope].iso_sig2_z

        self.MOutputs = ModelOutputs(self.c,MOd,TWlen, init_time, len(self.dz))
        ### self.MOutputs is a class instance
        ######################################
        ######################################

    ####################    
    ##### END INIT #####
    ####################

    def time_evolve(self):
        '''
        Evolve the spatial grid, time grid, accumulation rate, age, density, mass, stress, temperature, and diffusivity through time
        based on the user specified number of timesteps in the model run. Updates the firn density using a user specified 
        
        Writes to hdf file at end.
        '''

        self.steps = 1 / np.mean(self.t) # steps per year
        start_time=time.time() # this is a timer to keep track of how long the model run takes.

        if self.c['spinUpdate']:
            indUpdate = np.where(self.modeltime>=self.c['spinUpdateDate'])[0][0]
        
        ####################################
        ##### START TIME-STEPPING LOOP #####
        ####################################

        print('modeltime',self.modeltime[0],self.modeltime[-1])
        # print('stp',self.stp)
        for iii in range(self.stp):
            mtime = self.modeltime[iii]
            self.D_surf[iii] = iii
            if iii==1000:
                ntime = time.time()
                print('estimated model run time (seconds):', self.stp*(ntime-start_time)/1000)

            ### Merging process #VV ###
            if self.c['merging']: # merging may be deprecated (check with VV)
                if ((self.dz[1] < self.c['merge_min']) or (self.dz[0] < 1e-4)): # Start with surface merging
                    self.dz,self.z,self.gridLen,self.dx,self.rho,self.age,self.LWC,self.PLWC_mem,self.mass,self.mass_sum,self.sigma,self.bdot_mean,\
                        self.Dcon,self.T_mean,self.T10m,self.r2 = mergesurf(self,self.c['merge_min'])
                if (np.any(self.dz[2:] < self.c['merge_min'])): # Then merge rest of the firn column                       
                    self.dz,self.z,self.gridLen,self.dx,self.rho,self.age,self.LWC,self.PLWC_mem,self.mass,self.mass_sum,self.sigma,self.bdot_mean,\
                        self.Dcon,self.T_mean,self.T10m,self.r2 = mergenotsurf(self,self.c['merge_min'])

            ### dictionary of the parameters that get passed to physics
            PhysParams = {
                'iii':          iii,
                'steps':        self.steps,
                'gridLen':      self.gridLen,
                'bdotSec':      self.bdotSec,
                'bdot_mean':    self.bdot_mean,
                'bdot_av':      self.bdot_av,
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
                'FirnAir':      self.c['FirnAir']
            }
            
            if self.c['physRho']=='Morris2014':
                PhysParams['Hx'] = self.Hx
                PhysParams['QMorris'] = self.c['QMorris']

            if self.c['FirnAir']:
                PhysParams['AirRunType'] = self.cg['runtype']
                PhysParams['steady_T'] = self.cg['steady_T']

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
                'MaxSP':              FirnPhysics(PhysParams).MaxSP
            }

            RD      = physicsd[self.c['physRho']]()
            drho_dt = RD['drho_dt']
            if self.c['no_densification']:
                drho_dt = np.zeros_like(drho_dt)
            if self.c['physRho']=='Goujon2003':
                self.Gamma_Gou      = RD['Gamma_Gou'] 
                self.Gamma_old_Gou  = RD['Gamma_old_Gou']
                self.Gamma_old2_Gou = RD['Gamma_old2_Gou']
                self.ind1_old       = RD['ind1_old']

            ### update density and age of firn
            self.rho_old    = np.copy(self.rho)
            self.rho        = self.rho + self.dt[iii] * drho_dt
            self.dz_old     = np.copy(self.dz) # model volume thicknesses before the compaction
            self.sdz_old    = np.sum(self.dz) # old total column thickness (s for sum)
            self.z_old      = np.copy(self.z)
            self.dz         = self.mass / self.rho * self.dx # new dz after compaction
            
            if self.THist:
                self.Hx = RD['Hx']

            # self.Tz = np.concatenate(([self.Ts[iii]], self.Tz[1:]))

            ######################
            ### MELT #############
            if self.MELT:
                if self.c['liquid'] == 'prefsnowpack':
                    #VV You can choose a date to switch from bucket to prefflow, this should be easy to use from the json file
                    if self.modeltime[iii] >= 1980: # Apply dualperm from a certain date
                        if ((self.snowmeltSec[iii]>0.) or (np.any(self.LWC > 0.)) or (self.rainSec[iii] > 0.)): #i.e. there is water
                            self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn, self.LWC, self.PLWC_mem, self.r2, self.refrozen, self.runoff = prefflow(self,iii)
                        else: #Dry firn column and no input of meltwater                            
                            self.runoff     = np.array([0.]) #VV no runoff
                            self.refrozen   = np.zeros_like(self.dz) #VV no refreezing
                            self.dzn        = self.dz[0:self.compboxes] # Not sure this is necessary
                    elif self.modeltime[iii] < 1980: # Apply bVV until a certain date
                        if (self.snowmeltSec[iii]>0) or (np.any(self.LWC > 0.) or (self.rainSec[iii] > 0.)): #i.e. there is water
                            self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, self.dzn, self.LWC, self.refrozen, self.runoff, self.lwcerror = bucketVV(self,iii)
                        else:
                            #Dry firn column and no input of meltwater
                            self.runoff     = np.array([0.]) #VV no runoff
                            self.refrozen   = np.zeros_like(self.dz) #VV no refreezing
                            self.dzn        = self.dz[0:self.compboxes] # Not sure this is necessary
                    ### Heat ###        
                    self.Tz, self.T10m  = heatDiff(self,iii)
                ### end prefsnowpack ##################

                elif self.c['liquid'] == 'resingledomain':
                    #VV You can choose a date to switch from bucket to prefflow, this should be easy to use from the json file
                    if self.modeltime[iii] >= 1980: # Apply dualperm from a certain date
                        if ((self.snowmeltSec[iii]>0.) or (np.any(self.LWC > 0.)) or (self.rainSec[iii] > 0.)): #i.e. there is water
                            self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn, self.LWC, self.PLWC_mem, self.r2, self.refrozen, self.runoff = resingledomain(self,iii)
                        else:
                            #Dry firn column and no input of meltwater
                            self.runoff     = np.array([0.]) #VV no runoff
                            self.refrozen   = np.zeros_like(self.dz) #VV no refreezing
                            self.dzn        = self.dz[0:self.compboxes] # Not sure this is necessary
                    elif self.modeltime[iii] < 1980: # Apply bVV until a certain date
                        if (self.snowmeltSec[iii]>0) or (np.any(self.LWC > 0.) or (self.rainSec[iii] > 0.)): #i.e. there is water
                            self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, self.dzn, self.LWC, self.refrozen, self.runoff, self.lwcerror = bucketVV(self,iii)
                        else:
                            #Dry firn column and no input of meltwater
                            self.runoff     = np.array([0.]) #VV no runoff
                            self.refrozen   = np.zeros_like(self.dz) #VV no refreezing
                            self.dzn        = self.dz[0:self.compboxes] # Not sure this is necessary
                    ### Heat ###
                    self.Tz, self.T10m  = heatDiff(self,iii)
                ### end prefsnowpack ##################

                elif self.c['liquid'] == 'bucketVV':
                    if (self.snowmeltSec[iii]>0) or (np.any(self.LWC > 0.)) or (self.rainSec[iii] > 0.): #i.e. there is water
                        #VV 09/12/2020: change output variables of bucketVV
                        self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, self.dzn, self.LWC, meltgridtrack, self.refrozen, self.runoff, self.lwcerror = bucketVV(self,iii)
                        #VV 09/12/2020: account for changes in self.gridtrack
                        if self.doublegrid==True: #gridtrack corrected for melting
                            self.gridtrack = np.copy(meltgridtrack)

                    else:
                        #Dry firn column and no input of meltwater
                        self.runoff     = np.array([0.]) #VV no runoff
                        self.refrozen   = np.zeros_like(self.dz) #VV no refreezing
                        self.dzn        = self.dz[0:self.compboxes] # Not sure this is necessary
                    ### Heat ###
                    if self.c['LWCheat']=='enthalpy':
                        self.Tz, self.T10m, self.rho, self.mass, self.LWC = enthalpyDiff(self,iii)     
                    else:
                        self.Tz, self.T10m  = heatDiff(self,iii)
                              
                    coldlayers = np.where(self.Tz < T_MELT)
                    # if np.any(self.LWC[coldlayers[0]]>0.):
                        # print('Problem: water content in a cold layer')
                ### end bucketVV ##################

                elif self.c['liquid'] == 'percolation_bucket': ### Max's bucket scheme:
                    if iii==0:
                        print('percolation_bucket not recommended at present. Use bucketVV.')
                    if (self.MELT and self.snowmeltSec[iii]>0): #i.e. there is melt
                        self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn, self.LWC = percolation_bucket(self,iii)
                    else: # no melt, dz after compaction
                        self.dzn    = self.dz[0:self.compboxes]
                    ### Heat ###
                    if self.c['LWCheat']=='enthalpy':
                        self.Tz, self.T10m, self.rho, self.mass, self.LWC = enthalpyDiff(self,iii)
                    else: # use heat
                        if self.c['LWCheat']=='effectiveT':
                            self.Tz, self.T10m = effectiveT(self,iii)

                        else:
                            self.Tz, self.T10m  = heatDiff(self,iii)
                            if (self.c['LWCheat']=='LWCcorrect') or (self.c['LWCheat']=='lowK'):
                                if iii == 0:
                                    print('LWCcorrect is on')
                                self.Tz, self.LWC = LWC_correct(self)
                ### end percolation_bucket #########

                if self.LWC[-1] > 0.: #VV we don't want to lose water
                    print('LWC in last layer that is going to be removed, amount is:',self.LWC[-1])
                    # self.LWC[-2] += self.LWC[-1] #VV, flow routine will deal with saturation exceeding 1
                    # This should never happen if bottom of modelled firn column is at rho >= 830
                # self.LWC        = np.concatenate(([0], self.LWC[:-1]))

                if self.PLWC_mem[-1] > 0.: #VV
                    self.PLWC_mem[-2] += self.PLWC_mem[-1] #VV
                self.PLWC_mem    = np.concatenate(([0], self.PLWC_mem[:-1])) #VV
            
            else: # no melt, dz after compaction
                self.dzn    = self.dz[0:self.compboxes]
            ### end MELT #########
            ######################

            ### Surface energy balance #####
            if self.c['rad_pen']:
                if iii==0:
                    E_rp0 = 0.0
                elif (self.Ts[iii] - self.Ts[iii-1])<=0:
                    E_rp0 = 0.0                  
                else:    
                    i6cm = np.where(self.z>=0.06)[0][0]         
                    E_rp0 = self.E_sw[iii] #- (deltaT_s*self.mass[0]*CP_I) # estimate of the amount of energy that penetrates to lower layers, W/m^2
                    # k_ex =  -0.0338*self.rho+33.54
                    k_ex = -0.012 * self.rho+14.706 #extinction coefficient
                    
                    deltaE_layers1 = 0.64 * (1 - 0.8) * E_rp0 * np.exp(-k_ex*self.z) #* self.dz * self.dt[iii]
                    # deltaE_layers1 = 0.64 * E_rp0 * np.exp(-k_ex*self.z) #* self.dz * self.dt[iii]
                    diffE = -1*np.diff(deltaE_layers1)
                    deltaT = diffE[i6cm:]/(self.mass[i6cm:-1]*CP_I) * self.dt[iii]
                    self.Tz[0:i6cm] = self.Ts[iii]
                    self.Tz[i6cm:-1] = self.Tz[i6cm:-1]+deltaT
            ############################

            ### heat diffusion #########

            Ts_old = self.Tz[0]

            if (self.c['heatDiff'] and not self.MELT): # no melt, so use regular heat diffusion
                self.Tz, self.T10m  = heatDiff(self,iii)

            elif self.c['manualT']:
                tif = interpolate.interp1d(self.manualT_dep, self.manualT_temp[:,iii],kind='cubic')
                self.Tz = tif(self.z)

            elif (not self.c['heatDiff'] and not self.MELT): # user says no heat diffusion, so just set the temperature of the new box on top.
                self.Tz = self.Ts[iii]*np.ones_like(self.Tz)
                if iii==0:
                    print('warning:diffusion off, setting temp to Ts[iii]')

            self.T50     = np.mean(self.Tz[self.z<50]) # Temperature at 50 m

            ### Calculation of average surface temperature and accumulation rate VV
            # # Case 1: 1 year has not passed yet -> take averages of 1st year
            # if iii < self.steps: # VV
            #     T10m = np.sum(self.Ts[0:int(self.steps+1)])/self.steps # VV
            # # Case 2: at least 1 year has passed -> take averages of the last year (including the current step)
            # elif iii >= self.steps: # VV
            #     T10m = np.sum(self.Ts[int(iii-(self.steps-1)):iii+1])/self.steps # VV
            
            ### Firn Air ###############
            if self.c['FirnAir']: # Update firn air
                AirParams = {
                    'Tz':           self.Tz,
                    'rho':          self.rho,
                    'dt':           self.dt[iii],
                    'z':            self.z,
                    'rhos0':        self.rhos0[iii],
                    'dz_old':       self.dz_old,
                    'dz':           self.dz,
                    'rho_old':      self.rho_old,
                    'w_firn':       self.w_firn
                }
                for gas in self.cg['gaschoice']:        
                    self.Gz[gas], self.diffusivity, self.w_air, self.gas_age = self.FA[gas].firn_air_diffusion(AirParams,iii)
            ####################

            ### Isotopes #######
            if self.c['isoDiff']: # Update isotopes
                # self.del_z  = isoDiff(self,iii)
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
                ### new box gets added on within isoDiff function
                ####################

            if self.c['strain']: #update horizontal strain
                strain      = (-1 * self.du_dxSec[iii] * self.dt[iii] + 1) * np.ones_like(self.z)
                self.dz     = strain * self.dz
                self.mass   = strain * self.mass

            self.sdz_new    = np.sum(self.dz) #total column thickness after densification, melt, horizontal strain,  before new snow added

            ### Dcon: user-specific code goes here. 
            # self.Dcon[self.LWC>0] = self.Dcon[self.LWC>0] + 1 # for example, keep track of how many times steps the layer has had water

            ### Update grain growth ###
            if self.c['physGrain']: # update grain radius
                self.r2 = FirnPhysics(PhysParams).graincalc(iii) # calculate before accumulation b/c new surface layer should not be subject to grain growth yet
            
            ### update model grid, mass, stress, and mean accumulation rate
            if self.bdotSec[iii]>0: # there is accumulation at this time step
            # MS 2/10/17: should double check that everything occurs in correct order in time step (e.g. adding new box on, calculating dz, etc.)               
                self.age        = np.concatenate(([0], self.age[:-1])) + self.dt[iii]      
                self.dzNew      = self.bdotSec[iii] * RHO_I / self.rhos0[iii] * S_PER_YEAR
                self.dz         = np.concatenate(([self.dzNew], self.dz[:-1]))
                self.z          = self.dz.cumsum(axis = 0)
                znew = np.copy(self.z) 
                self.z          = np.concatenate(([0], self.z[:-1]))
                self.rho        = np.concatenate(([self.rhos0[iii]], self.rho[:-1]))

                if self.c['physGrain']: # update grain radius
                    r2surface       = FirnPhysics(PhysParams).surfacegrain() #grain size for new surface layer
                    self.r2         = np.concatenate(([r2surface], self.r2[:-1]))               
                if not self.c['manualT']:
                    self.Tz         = np.concatenate(([self.Ts[iii]], self.Tz[:-1]))
                
                self.Dcon       = np.concatenate(([self.D_surf[iii]], self.Dcon[:-1]))
                massNew         = self.bdotSec[iii] * S_PER_YEAR * RHO_I
                self.mass       = np.concatenate(([massNew], self.mass[:-1]))
                self.compaction = np.append(0,(self.dz_old[0:self.compboxes-1]-self.dzn[0:self.compboxes-1]))#/self.dt*S_PER_YEAR)
                # self.compaction = self.dz_old[0:self.compboxes] - self.dzn[0:self.compboxes]
                if self.doublegrid:
                    self.gridtrack = np.concatenate(([1],self.gridtrack[:-1]))
                self.LWC        = np.concatenate(([0], self.LWC[:-1]))

            elif self.bdotSec[iii]<0:

                self.mass_sum   = self.mass.cumsum(axis = 0)
                
                self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, self.dzn, self.LWC, self.PLWC_mem, self.totwatersublim, sublgridtrack     = sublim(self,iii) # keeps track of sublimated water for mass conservation
                
                self.compaction = (self.dz_old[0:self.compboxes]-self.dzn)
                self.dzNew      = 0
                if self.doublegrid == True: # gridtrack corrected for sublimation
                    self.gridtrack = np.copy(sublgridtrack)
                znew = np.copy(self.z)

            else: # no accumulation during this time step
                self.age        = self.age + self.dt[iii]
                self.z          = self.dz.cumsum(axis=0)
                self.z          = np.concatenate(([0],self.z[:-1]))
                self.dzNew      = 0
                znew = np.copy(self.z)                             
                self.compaction = (self.dz_old[0:self.compboxes]-self.dzn)
                if not self.c['manualT']:
                    self.Tz = np.concatenate(([self.Ts[iii]], self.Tz[1:]))

            self.w_firn = (znew - self.z_old) / self.dt[iii] # advection rate of the firn, m/s

            self.sigma      = (self.mass + (self.LWC * RHO_W_KGM)) * self.dx * GRAVITY
            self.sigma      = self.sigma.cumsum(axis = 0)
            self.mass_sum   = self.mass.cumsum(axis = 0)
            
            self.bdot_mean  = (np.concatenate(( [self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] * self.t[iii] / (self.age[1:] * RHO_I) ))) * np.mean(S_PER_YEAR/self.dt) * S_PER_YEAR

            ###NOTE: sigma = bdot_mean*GRAVITY*age/S_PER_YEAR*917.0) (or, sigma = bdot*g*tau, steady state conversion.)

            #############################################################
            ### write results as often as specified in the init method ##
            #############################################################
            if mtime in self.TWrite:                
                ind         = np.where(self.TWrite == mtime)[0][0]
                mtime_plus1 = self.TWrite[ind]
                
                if 'viscosity' in self.output_list:
                    self.viscosity = RD['viscosity']

                self.climate = np.array([self.bdot[iii],self.Ts[iii]])
   
                bcoAgeMart, bcoDepMart, bcoAge830, bcoDep830, LIZAgeMart, LIZDepMart, bcoAge815, bcoDep815  = self.update_BCO(iii)

                intPhi, self.DIPc, z_co  = self.update_DIP()
                dHOut, dHOutC, compOut, dHOutcorr, dHOutcorrC  = self.update_dH(iii)
                try:
                    ind_z = np.where(self.z>=self.DIPhorizon)[0][0]
                    DIPhz = self.DIPc[ind_z]
                except Exception:
                    DIPhz = np.nan

                if mtime==self.TWrite[0]:
                    self.dHAll  = 0 * self.dHAll
                    self.dHAllcorr = 0 * self.dHAllcorr
                    dH          = 0.0
                    dHtot       = 0.0
                    comp_firn   = 0.0
                    dHcorr      = 0.0
                    dHtotcorr   = 0.0

                self.BCO  = np.array([bcoAgeMart, bcoDepMart, bcoAge830, bcoDep830, LIZAgeMart, LIZDepMart, bcoAge815, bcoDep815, z_co])
                self.DIP  = np.array([intPhi, dHOut, dHOutC, compOut, dHOutcorr, dHOutcorrC,DIPhz])

                MOd = {key:value for key, value in self.__dict__.items() if key in self.output_list}

                if self.c['FirnAir']:    
                    for gas in self.cg['gaschoice']:
                        MOd[gas] = self.Gz[gas]

                if self.c['isoDiff']:
                    for isotope in self.c['iso']:
                        MOd['isotopes_{}'.format(isotope)] = self.Isotopes[isotope].del_z
                        MOd['iso_sig2_{}'.format(isotope)] = self.Isotopes[isotope].iso_sig2_z

                self.MOutputs.updateMO(MOd,mtime_plus1,self.WTracker)

                self.WTracker = self.WTracker + 1
            ################################
            ### End write ##################
            ################################

            if (self.c['spinUpdate'] and iii==indUpdate):
                SpinUpdate(self,mtime)

            if self.doublegrid:
                #VV changes 09/12/2020
                #if self.gridtrack[-1]==2:
                    ## print('regridding now at ', iii)
                    #self.dz, self.z, self.rho, self.Tz, self.mass, self.sigma, self. mass_sum, self.age, self.bdot_mean, self.LWC, self.gridtrack, self.r2 = regrid(self)
                    #if iii<100:
                        #tdep = np.where(self.gridtrack==1)[0][-1]
                        #print('transition at:', self.z[tdep])
                if self.gridtrack[-1]!=3: #VV works for whatever the gridtrack value we have
                    self.dz, self.z, self.rho, self.Tz, self.mass, self.sigma, self. mass_sum, self.age, self.bdot_mean, self.LWC, self.gridtrack, self.r2 = regrid22(self) #VV regrid22

        ##################################
        ##### END TIME-STEPPING LOOP #####
        ##################################

        # write_nospin_hdf5(self)
        write_nospin_hdf5(self,self.MOutputs.Mout_dict)


    ###########################
    ##### END time_evolve #####
    ###########################

    def update_BCO(self,iii):
        '''
        Updates the bubble close-off depth and age based on the Martinerie criteria as well as through assuming the critical density is 815 kg/m^3
        '''

        try:
            if (self.c['FirnAir'] and self.cg['runtype']=='steady'):
                bcoMartRho  = 1 / (1 / (917.0) + self.cg['steady_T'] * 6.95E-7 - 4.3e-5)  # Martinerie density at close off
            else:
                bcoMartRho  = 1 / (1 / (917.0) + self.T_mean[iii] * 6.95E-7 - 4.3e-5)  # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)

            bcoAgeMart  = min(self.age[self.rho >= bcoMartRho]) / S_PER_YEAR  # close-off age from Martinerie
            bcoDepMart  = min(self.z[self.rho >= (bcoMartRho)])

            # bubble close-off age and depth assuming rho_crit = 815kg/m^3
            bcoAge830   = min(self.age[self.rho >= 830.0]) / S_PER_YEAR  # close-off age where rho = 815 kg m^-3
            bcoDep830   = min(self.z[self.rho >= 830.0])
            bcoAge815   = min(self.age[self.rho >= (RHO_2)]) / S_PER_YEAR  # close-off age where rho = 815 kg m^-3
            bcoDep815   = min(self.z[self.rho >= (RHO_2)])

            LIZMartRho = bcoMartRho - 14.0  # LIZ depth (Blunier and Schwander, 2000)
            LIZAgeMart = min(self.age[self.rho > LIZMartRho]) / S_PER_YEAR  # lock-in age
            LIZDepMart = min(self.z[self.rho >= (LIZMartRho)])  # lock in depth

        except:
            
            bcoAgeMart  = -9999
            bcoDepMart  = -9999
            bcoAge830   = -9999
            bcoDep830   = -9999

            LIZDepMart  = -9999
            LIZAgeMart  = -9999
            bcoAge815   = -9999
            bcoDep815   = -9999

        return bcoAgeMart, bcoDepMart, bcoAge830, bcoDep830, LIZAgeMart, LIZDepMart, bcoAge815, bcoDep815

    ### end update_BCO ########
    ###########################

    def update_DIP(self):
        '''
        Updates the depth-integrated porosity
        '''

        bcoMartRho  = 1 / (1 / (917.0) + self.T50* 6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
        phi         = 1 - self.rho / RHO_I  # total porosity
        phi[phi <= 0] = 1e-16
        phiC        = 1 - bcoMartRho / RHO_I;  # porosity at close off
        phiClosed   = 0.37 * phi * (phi / phiC) ** -7.6  # Closed porosity, from Goujon. See Buizert thesis (eq. 2.3) as well
        phiOpen     = phi - phiClosed  # open porosity
        
        try:
            ind_co  = np.where(phiOpen<=1e-10)[0][0]
            z_co    = self.z[ind_co]
        except:
            z_co    = -9999

        phiOpen[phiOpen <= 0] = 1.e-10  # don't want negative porosity.

        intPhi      = np.sum(phi * self.dz)  # depth-integrated porosity
        intPhi_c    = np.cumsum(phi * self.dz)

        # self.intPhiAll.append(intPhi)

        return intPhi, intPhi_c, z_co
    ### end update_DIP ########
    ###########################

    def update_dH(self,iii):
        '''
        updates the surface elevation change
        '''

        self.dH     = (self.sdz_new - self.sdz_old) + self.dzNew - (self.iceout*self.t[iii]) # iceout has units m ice/year, t is years per time step. 
        # self.dH2 = self.z[-1] - self.z_old[-1] #- (self.iceout*self.t) # alternative method. Should be the same?    
        self.dHAll.append(self.dH)
        self.dHtot  = np.sum(self.dHAll)
        
        ### If the bottom of the domain is not the ice density, there is 
        ### compaction that is not accounted for between the bottom of the 
        ### domain and the 917 density horizon.

        iceout_corr = self.iceout*RHO_I/self.rho[-1]
        self.dHcorr = (self.sdz_new - self.sdz_old) + self.dzNew - (iceout_corr*self.t[iii]) # iceout has units m ice/year, t is years per time step. 
        
        self.dHAllcorr.append(self.dHcorr)
        self.dHtotcorr = np.sum(self.dHAllcorr)
        self.comp_firn = self.sdz_new - self.sdz_old #total compaction of just the firn during the previous time step

        return self.dH, self.dHtot, self.comp_firn, self.dHcorr, self.dHtotcorr

    ###########################

