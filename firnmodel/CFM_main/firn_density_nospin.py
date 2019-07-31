#!/usr/bin/env python

from diffusion import *
from reader import read_input
from reader import read_init
from writer import write_spin_hdf5
from writer import write_nospin_hdf5
from physics import *
from constants import *
from melt import *
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

class FirnDensityNoSpin:
    '''
    Parameters used in the model, for the initialization as well as the time evolution:

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
    : bdot_mean: mean accumulation over the lifetime of each parcel
                (units are m I.E. per year)
    :returns D_surf: diffusivity tracker
                (unit: ???, type: array of floats)

    '''

    def __init__(self, configName):
        '''
        Sets up the initial spatial grid, time grid, accumulation rate, age, density, mass, stress, temperature, and diffusivity of the model run
        :param configName: name of json config file containing model configurations
        '''
        ### load in json config file and parses the user inputs to a dictionary
        self.spin = True
        with open(configName, "r") as f:
            jsonString      = f.read()
            self.c          = json.loads(jsonString)
        print("Main run starting")
        print("physics are", self.c['physRho'])

        ### read in initial depth, age, density, temperature from spin-up results
        initDepth   = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'depthSpin')
        initAge     = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'ageSpin')
        initDensity = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'densitySpin')
        initTemp    = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'tempSpin')

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

        ### get temperature and accumulation rate from input csv file
        input_temp, input_year_temp = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNameTemp']))
        if input_temp[0] < 0.0:
            input_temp      = input_temp + K_TO_C
        input_temp[input_temp>T_MELT] = T_MELT

        input_bdot, input_year_bdot = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamebdot']))

        try:
            if self.c['MELT']:
                input_snowmelt, input_year_snowmelt = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamemelt']))
                self.MELT           = True
                self.LWC            = np.zeros_like(self.z)
                print("Melt is initialized")
            else:
                self.MELT           = False
                print("Melt is not turned on.")
                input_snowmelt      = None
                input_year_snowmelt = None
                self.LWC            = np.zeros_like(self.z)

        except:
            self.MELT               = False
            print("No melt; json does not include a melt field")
            input_snowmelt          = None
            input_year_snowmelt     = None
            self.LWC                = np.zeros_like(self.z)

        #####################
        ### time ############
        if 'timesetup' in self.c:
            pass
        else:
            print('You should add "timesetup" to the .json')
            self.c['timesetup']='interp'

        # year to start and end, from the input file. If inputs have different start/finish, take only the overlapping times
        if (self.c['timesetup']=='interp' or self.c['SeasonalTcycle']):
            if (self.c['SeasonalTcycle'] and self.c['timesetup']=='exact'):
                print('"Exact" time setup does not work with "SeasonalTcycle"')
            yr_start        = max(input_year_temp[0], input_year_bdot[0])   # start year
            yr_end          = min(input_year_temp[-1], input_year_bdot[-1]) # end year
            
            # self.years      = np.ceil((yr_end - yr_start) * 1.0)
            self.years      = (yr_end - yr_start) * 1.0
            self.dt         = S_PER_YEAR / self.c['stpsPerYear'] # seconds per time step
            self.stp        = int(self.years * S_PER_YEAR/self.dt)#-1       # total number of time steps, as integer
            self.dt = self.dt * np.ones(self.stp)
            # self.modeltime  = np.linspace(yr_start, yr_end, self.stp + 1)   # vector of time of each model step
            self.modeltime  = np.linspace(yr_start, yr_end, self.stp)
            self.t          = 1.0 / self.c['stpsPerYear']                   # years per time step
            # print('Model year start:', yr_start)
            # print('Model year end:', yr_end)
            # print('stp',self.stp)
            # print('dt',self.dt)
            # print('years',self.years)
            # print('modeltime',self.modeltime[0:5])
            # print('modeltime',self.modeltime[-5:])
            # print('t',self.t)

        elif self.c['timesetup']=='exact':
            print('"Exact" time setup will not work properly if input forcing does not all have the same time')
            yr_start = input_year_temp[1]
            yr_end = input_year_temp[-1]
            self.dt = np.diff(input_year_temp)*S_PER_YEAR
            self.stp = len(self.dt)
            self.modeltime = input_year_temp[1:]
            self.t = np.mean(np.diff(input_year_temp))
            # print('Model year start:', yr_start)
            # print('Model year end:', yr_end)
            # print('stp',self.stp)
            # print('dt',self.dt)
            # # print('years',self.years)
            # print('modeltime',self.modeltime[0:5])
            # print('modeltime',self.modeltime[-5:])
            # print('t',self.t)

        # sys.exit()
        #####################
      
        ###############################
        ### surface boundary conditions
        ### temperature, accumulation, melt, isotopes, surface density
        ###############################
        int_type            = self.c['int_type']
        print('Climate interpolation method is %s' %int_type)

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
        try: #Rolling mean average surface temperature (vector)
            Nyears = 10 #number of years to average for T_mean
            NN = int(self.c['stpsPerYear']*Nyears)
            self.T_mean = pd.Series(self.Ts).rolling(window=NN+1,win_type='hamming').mean().values
            self.T_mean[np.isnan(self.T_mean)]=self.T_mean[NN]
        except Exception:
            self.T_mean = np.mean(self.Ts) * np.ones(self.stp)
            print('Error calculating T_mean, using mean surface over all time')
        #####################

        ### Accumulation ####
        bsf                 = interpolate.interp1d(input_year_bdot,input_bdot,int_type,fill_value='extrapolate') # interpolation function
        self.bdot           = bsf(self.modeltime) # m ice equivalent per year
        # self.bdotSec      = self.bdot / S_PER_YEAR / (self.stp / self.years) # accumulation rate in per second
        self.bdotSec        = self.bdot / S_PER_YEAR / self.c['stpsPerYear'] # accumulation for each time step (meters i.e. per second)       
        try: 
            if self.c['manual_iceout']:
                self.iceout = self.c['iceout']
                print('Ensure that your iceout value has units m ice eq. per year!')
            else:
                self.iceout = np.mean(self.bdot) # this is the rate of ice flow advecting out of the column, units m I.E. per year.
        
            # print('Iceout is:' self.iceout)
            # print('Iceout per time step:' self.iceout*self.t)
        except Exception:
            print('add field "manual_iceout" to .json file to set iceout value manually')
            self.iceout = np.mean(self.bdot) # this is the rate of ice flow advecting out of the column, units m I.E. per year.

        # print('iceout, used for dH, is %s' %self.iceout)
        self.w_firn         = np.mean(self.bdot) * RHO_I / self.rho 

        if (np.any(self.bdotSec<0.0) and self.c['bdot_type']=='instant'):
            print('ERROR: bdot_type set to "instant" in .json and input')
            print('accumulation has at least one negative value.') 
            print('QUITTING MODEL RUN.')
            sys.exit()
        #####################
        
        ### Melt ############
        if self.MELT:
            ssf                 = interpolate.interp1d(input_year_snowmelt,input_snowmelt,int_type,fill_value='extrapolate')
            self.snowmelt       = ssf(self.modeltime)
            self.snowmeltSec    = self.snowmelt / S_PER_YEAR / self.c['stpsPerYear'] # melt for each time step (meters i.e. per second)
        #####################

        ### Isotopes ########
        if self.c['isoDiff']:
            init_del_z          = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'IsoSpin')
            try:
                input_iso, input_year_iso = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNameIso']))
                self.del_s      = np.interp(self.modeltime, input_year_iso, input_iso)
                # del_s0 = input_iso[0]
            except:
                print('No external file for surface isotope values found, but you specified in the config file that isotope diffusion is on. The model will generate its own synthetic isotope data for you.')
                # del_s0 = -50.0
                ar1             = 0.9   # red noise memory coefficient
                std_rednoise    = 2    # red noise standard deviation
                self.del_s      = std_rednoise*np.random.randn(self.stp)    # white noise

                for x in range(1,self.stp):
                    self.del_s[x] = self.del_s[x-1]*ar1 + np.random.randn()  # create red noise from white

                self.del_s      = self.del_s - 50
                #impose seasonal isotope cycle
                # self.del_s    = self.del_s + 5 * (np.cos(2 * np.pi * np.linspace(0, self.years, self.stp )) + 0.3 * np.cos(4 * np.pi * np.linspace(0, self.years, self.stp )))
        #####################
 
        ### Surface Density #
        if self.c['variable_srho']:
            if self.c['srho_type']=='userinput':
                input_srho, input_year_srho = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamerho']))
                Rsf             = interpolate.interp1d(input_year_srho,input_srho,'nearest',fill_value='extrapolate') # interpolation function
                self.rhos0      = Rsf(self.modeltime) # surface temperature interpolated to model time
                # self.rhos0      = np.interp(self.modeltime, input_year_srho, input_srho)
            elif self.c['srho_type']=='param':
                self.rhos0      = 481.0 + 4.834 * (self.Ts - T_MELT) # Kuipers Munneke, 2015
            elif self.c['srho_type']=='noise':
                rho_stdv        = 50 # the standard deviation of the surface density (I made up 25)
                self.rhos0      = np.random.normal(self.c['rhos0'], rho_stdv, self.stp)
                self.rhos0[self.rhos0>700]=700
                self.rhos0[self.rhos0<300]=300
                print('Max surface density is:', np.max(self.rhos0))
                print('Min surface density is:', np.min(self.rhos0))
        else:
            self.rhos0      = self.c['rhos0'] * np.ones(self.stp)       # density at surface
        #####################

        ### Layer tracker ###
        self.D_surf     = self.c['D_surf'] * np.ones(self.stp)      # layer traking routine (time vector). 
        self.Dcon       = self.c['D_surf'] * np.zeros(self.gridLen)  # layer tracking routine (initial depth vector)
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
        if self.c['TWriteInt']!=1:
            print('Time writing interval is not 1; dH output will not be accurate.')
        # self.TWrite       = np.append(self.modeltime[10],self.TWrite)
        # self.TWrite       = self.modeltime[-2::self.c['TWriteInt']]
        # self.TWrite_out   = self.TWrite
        TWlen               = len(self.TWrite) #- 1
        self.WTracker       = 1

        ### you can choose to write certain fields at different times. (Functionality not fully tested)
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
        self.sigma          = (self.mass + self.LWC * RHO_W_KGM) * self.dx * GRAVITY
        self.sigma          = self.sigma.cumsum(axis = 0)
        self.mass_sum       = self.mass.cumsum(axis = 0)
        ### mean accumulation over the lifetime of the parcel:
        self.bdot_mean      = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] / (self.age[1:] * RHO_I / self.t))))*self.c['stpsPerYear']*S_PER_YEAR
        ### It is the mass of the overlying firn divided by the age of the parcel.
        #######################

        ### set up longitudinal strain rate
        if self.c['strain']:
            # input units are yr^-1
            input_dudx, input_year_dudx = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamedudx']))
            dusf                        = interpolate.interp1d(input_year_dudx,input_dudx,int_type,fill_value='extrapolate')           
            self.du_dx      = dusf(self.modeltime)
            self.du_dxSec   = self.du_dx / S_PER_YEAR / self.c['stpsPerYear'] # strain rate (s^-1) at each time step

            # except:
                # print('strain time series not found; using du_dx from json instead')
                # self.du_dx        = np.ones_like(self.modeltime) * self.c['du_dx']
                # self.du_dxSec     = self.du_dx / S_PER_YEAR / self.c['stpsPerYear'] # strain rate (s^-1) at each time step

        #######################
        
        self.Tz             = initTemp[1:]
        self.T50         = np.mean(self.Tz[self.z<50])
        self.T10m           = self.Tz[np.where(self.z>=10.0)[0][0]]

        # self.compboxes      = len(self.z[self.z<80])
        self.compboxes = len(self.z)
        #######################

        ### model outputs
        self.output_list    = self.c['outputs']
        print(self.output_list)
        if ((not self.MELT) and ('LWC' in self.output_list)):
            self.output_list.remove('LWC')
            print('removed LWC from output list (melt is not on)')
        if ((not self.c['FirnAir']) and ('gasses' in self.output_list)):
            self.output_list.remove('gasses')
            print('removed gasses from output list (firn air is not on)')
        if ((not self.c['isoDiff']) and ('isotopes' in self.output_list)):
            self.output_list.remove('isotopes')
            print('removed isotopes from output list')          
        if 'density' in self.output_list:
            self.rho_out            = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
            self.rho_out[0,:]       = np.append(self.modeltime[0], self.rho)
        if 'temperature' in self.output_list:
            self.Tz_out             = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
            self.Tz_out[0,:]        = np.append(self.modeltime[0], self.Tz)
        if 'age' in self.output_list:
            self.age_out            = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
            self.age_out[0,:]       = np.append(self.modeltime[0], self.age/S_PER_YEAR)
        if 'depth' in self.output_list:
            self.z_out              = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
            self.z_out[0,:]         = np.append(self.modeltime[0], self.z)
        if 'dcon' in self.output_list:
            self.D_out              = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
            self.D_out[0,:]         = np.append(self.modeltime[0], self.Dcon)
        if 'bdot_mean' in self.output_list:
            self.bdot_out           = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
            self.bdot_out[0,:]      = np.append(self.modeltime[0], self.bdot_mean)
        if 'climate' in self.output_list:
            self.Clim_out           = np.zeros((TWlen+1,3),dtype='float32')
            self.Clim_out[0,:]      = np.append(self.modeltime[0], [self.bdot[0], self.Ts[0]])  # not sure if bdot or bdotSec
        if 'compaction' in self.output_list:
            self.comp_out          = np.zeros((TWlen+1,self.compboxes+1),dtype='float32')
            self.comp_out[0,:]     = np.append(self.modeltime[0], np.zeros(self.compboxes))
        if 'LWC' in self.output_list:
            self.LWC_out            = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
            self.LWC_out[0,:]       = np.append(self.modeltime[0], self.LWC)
        if 'viscosity' in self.output_list:
            self.viscosity          = np.zeros(self.gridLen)
            self.viscosity_out      = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
            self.viscosity_out[0,:] = np.append(self.modeltime[0], self.viscosity)
        try:
            print('rho_out size (MB):', self.rho_out.nbytes/1.0e6) # print the size of the output for reference
        except:
            pass
        #####################


        ### Writing also the outputs required for retmip #VV ###
        self.runoff = np.array([0.]) #VV total runoff which is a single value for the whole firn column
        self.refrozen = np.zeros_like(self.dz) #VV refreezing in every layer, array of size of our grid
        self.totalrunoff = np.array([0.]) # Might be useful to have a total final value without having to write every time step
        self.totalrefrozen = np.zeros_like(self.dz) # Might be useful to have a total final value without having to write every time step
        self.totwatersublim = 0. #VV Total amount of liquid water that get sublimated
        self.lwcerror = 0. #VV
        self.totallwcerror =0. #VV


        ### initial grain growth (if specified in config file)
        if self.c['physGrain']:
            initr2                  = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'r2Spin')
            self.r2                 = initr2[1:]
            r20                     = self.r2
            self.dr2_dt             = np.zeros_like(self.z)
            if 'grainsize' in self.output_list:
                self.r2_out         = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
                self.r2_out[0,:]    = np.append(self.modeltime[0], self.r2)
                self.dr2_dt_out     = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
                self.dr2_dt_out[0,:]= np.append(self.modeltime[0], self.dr2_dt)
            else:
                self.r2_out         = None
                self.dr2_dt_out     = None            
        else:            
            self.r2                 = None
            self.dr2_dt             = None
        #######################

        ### temperature history for Morris physics
        if self.c['physRho'] == 'Morris2014':
            self.THist              = True
            initHx                  = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'HxSpin')
            self.Hx                 = initHx[1:]
            if 'temp_Hx' in self.output_list:
                self.Hx_out         = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
                self.Hx_out[0,:]    = np.append(self.modeltime[0], self.Hx)
            else:
                self.Hx_out         = None
        else:
            self.THist              = False
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
            self.del_z              = init_del_z[1:]
            if 'isotopes' in self.output_list:
                self.iso_out        = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
                self.iso_out[0,:]   = np.append(self.modeltime[0], self.del_z)
            else:
                self.iso_out        = None
        #######################

        ### DIP, DHdt, LIZ, BCO ###
        self.dHAll      = []
        self.dHAllcorr  = []
        bcoAgeMart, bcoDepMart, bcoAge830, bcoDep830, LIZAgeMart, LIZDepMart, bcoAge815, bcoDep815  = self.update_BCO()
        intPhi, intPhi_c, z_co = self.update_DIP()
        
        self.dHAll.append(0)
        self.dHAllcorr.append(0)
        dHOut   = 0 # surface elevation change since last time step
        dHOutC  = 0 # cumulative surface elevation change since start of model run
        compOut = 0 # compaction of just the firn at each time step; no ice dynamics or accumulation
        dHOutcorr = 0
        dHOutcorrC = 0

        if 'DIP' in self.output_list:
            self.DIP_out        = np.zeros((TWlen+1,7),dtype='float32')   
            self.DIP_out[0,:]   = np.append(self.modeltime[0], [intPhi, dHOut, dHOutC, compOut, dHOutcorr, dHOutcorrC])
            self.DIPc_out       = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
            self.DIPc_out[0,:]  = np.append(self.modeltime[0], intPhi_c)
        if 'BCO' in self.output_list:
            self.BCO_out        = np.zeros((TWlen+1,10),dtype='float32')
            self.BCO_out[0,:]   = np.append(self.modeltime[0], [bcoAgeMart, bcoDepMart, bcoAge830, bcoDep830, LIZAgeMart, LIZDepMart, bcoAge815, bcoDep815, z_co])
        #####################

        ##### Firn Air ######
        '''
        each gas of interest gets its own instance of the class, each instance
        is stored in a dictionary
        '''
        if self.c['FirnAir']:
            print('Firn air initialized')
            with open(self.c['AirConfigName'], "r") as f:
                jsonString      = f.read()
                self.cg         = json.loads(jsonString)
            self.FA         = {}
            self.gas_out    = {}
            self.Gz         = {}
            for gas in self.cg['gaschoice']:
                if (gas=='d15N2' or gas=='d40Ar'):
                    input_year_gas = input_year_temp
                    input_gas = np.ones_like(input_year_temp)
                else:
                    input_gas, input_year_gas = read_input(os.path.join(self.c['InputFileFolder'],'%s.csv' %gas))

                Gsf     = interpolate.interp1d(input_year_gas,input_gas,'linear',fill_value='extrapolate')
                Gs      = Gsf(self.modeltime)

                self.FA[gas] = FirnAir(self.cg, Gs, self.z, self.modeltime, self.Tz, self.rho, self.dz, gas, self.bdot)

                if "gasses" in self.cg['outputs']:
                    self.gas_out[gas]           = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
                    self.gas_out[gas][0,:]      = np.append(self.modeltime[0], np.ones_like(self.rho))
            if "diffusivity" in self.cg['outputs']:
                self.diffu_out          = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
                self.diffu_out[0,:]     = np.append(self.modeltime[0], np.ones_like(self.rho))
            if "gas_age" in self.cg['outputs']:
                self.gas_age_out            = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
                self.gas_age_out[0,:]       = np.append(self.modeltime[0], np.zeros_like(self.rho))
            if "advection_rate" in self.cg['outputs']:
                self.w_air_out          = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
                self.w_air_out[0,:]     = np.append(self.modeltime[0], np.ones_like(self.rho))              
                self.w_firn_out         = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
                self.w_firn_out[0,:]    = np.append(self.modeltime[0], np.ones_like(self.rho))
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
        #####################

    ####################    
    ##### END INIT #####
    ####################

    def time_evolve(self):
        '''
        Evolve the spatial grid, time grid, accumulation rate, age, density, mass, stress, temperature, and diffusivity through time
        based on the user specified number of timesteps in the model run. Updates the firn density using a user specified 
        '''
        self.steps = 1 / self.t # steps per year
        start_time=time.time() # this is a timer to keep track of how long the model run takes.
        
        ####################################
        ##### START TIME-STEPPING LOOP #####
        ####################################
        
        for iii in range(self.stp):
            mtime = self.modeltime[iii]
            # print('########################')
            # if (iii<5 or iii>self.stp-5):
            #     print('mtime', mtime)
                # print('melt', self.snowmelt[iii])
            self.D_surf[iii] = iii
            ### dictionary of the parameters that get passed to physics
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
            
            if self.THist: #add Hx to dictionary if physics is Morris
                PhysParams['Hx'] = self.Hx

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
                'Max2018':              FirnPhysics(PhysParams).Max2018
            }

            RD      = physicsd[self.c['physRho']]()
            drho_dt = RD['drho_dt']

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
                # self.Hx     = FirnPhysics(PhysParams).THistory()

            if (self.MELT and self.snowmeltSec[iii]>0): #i.e. there is melt            
                ### Max's bucket scheme:
                print(mtime)
                self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn, self.LWC = percolation_bucket(self,iii)

                ### Vincent's bucket scheme:
                # self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, self.dzn, self.LWC, self.refrozen, self.runoff, self.lwcerror = bucketVV(self,iii)

            else: # no melt, dz after compaction
                self.dzn    = self.dz[0:self.compboxes]

            ### heat diffusion
            if (self.c['heatDiff'] and not self.MELT): # no melt, so use regular heat diffusion
                self.Tz, self.T10m  = heatDiff(self,iii)

            elif (self.c['heatDiff'] and self.MELT): # there is melt, so use enthalpy method
                self.Tz, self.T10m, self.rho, self.mass, self.LWC = enthalpyDiff(self,iii)
                
                # self.Tz, self.T10m  = heatDiff(self,iii) #Alternatively, use the heat diffusion scheme.

            else: # user says no heat diffusion, so just set the temperature of the new box on top.
                self.Tz     = np.concatenate(([self.Ts[iii]], self.Tz[:-1]))
                pass # box gets added below

            self.T50     = np.mean(self.Tz[self.z<50])

            '''Calculation of average surface temperature and accumulation rate #VV '''
            # Case 1: 1 year has not passed yet -> take averages of 1st year
            if iii < self.steps: # VV
                T10m = np.sum(self.Ts[0:int(self.steps+1)])/self.steps # VV
            # Case 2: at least 1 year has passed -> take averages of the last year (including the current step)
            elif iii >= self.steps: # VV
                T10m = np.sum(self.Ts[int(iii-(self.steps-1)):iii+1])/self.steps # VV
            
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
                    self.Gz[gas], self.diffu, w_p, gas_age0 = self.FA[gas].firn_air_diffusion(AirParams,iii)

            if self.c['isoDiff']: # Update isotopes
                self.del_z  = isoDiff(self,iii)
                ### new box gets added on within isoDiff function
                
            if self.c['strain']: #update horizontal strain
                strain      = (-1 * self.du_dxSec[iii] * self.dt[iii] + 1) * np.ones_like(self.z)
                self.dz     = strain * self.dz
                self.mass   = strain * self.mass

            self.sdz_new    = np.sum(self.dz) #total column thickness after densification, melt, horizontal strain,  before new snow added

            ### Dcon: user-specific code goes here. 
            # self.Dcon[self.LWC>0] = self.Dcon[self.LWC>0] + 1 # for example, keep track of how many times steps the layer has had water

            ### Update grain growth ###
            if self.c['physGrain']: # update grain radius
                self.r2 = FirnPhysics(PhysParams).graincalc() # calculate before accumulation b/c new surface layer should not be subject to grain growth yet
            
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

                self.LWC        = np.concatenate(([0], self.LWC[:-1]))
                self.Tz         = np.concatenate(([self.Ts[iii]], self.Tz[:-1]))
                self.Dcon       = np.concatenate(([self.D_surf[iii]], self.Dcon[:-1]))
                massNew         = self.bdotSec[iii] * S_PER_YEAR * RHO_I
                self.mass       = np.concatenate(([massNew], self.mass[:-1]))
                self.compaction = np.append(0,(self.dz_old[0:self.compboxes-1]-self.dzn[0:self.compboxes-1]))#/self.dt*S_PER_YEAR)
                # self.compaction = self.dz_old[0:self.compboxes] - self.dzn[0:self.compboxes]
                if self.doublegrid:
                    self.gridtrack = np.concatenate(([1],self.gridtrack[:-1]))

            # elif self.bdotSec[iii]<0: #work in progress from VV
            #     self.mass_sum      = self.mass.cumsum(axis = 0) #VV
            #     self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, self.dzn, self.LWC, self.PLWC_mem, self.totwatersublim = sublim(self,iii) #VV keeps track of sublimated water for mass conservation            

            else: # no accumulation during this time step
                
                self.age        = self.age + self.dt[iii]
                self.z          = self.dz.cumsum(axis=0)
                self.z          = np.concatenate(([0],self.z[:-1]))
                self.dzNew      = 0
                znew = np.copy(self.z)              
                
                self.compaction = (self.dz_old[0:self.compboxes]-self.dzn)
                self.Tz = np.concatenate(([self.Ts[iii]], self.Tz[1:]))

            self.w_firn = (znew - self.z_old) / self.dt[iii] # advection rate of the firn, m/s

            ### find the compaction rate
            ### this should all be old (11/28/17)
            # zdiffnew      = (self.z[1:]-self.z[1])
            # zdiffold      = (self.z_old[0:-1]-self.z_old[0])

            # zdn           = self.z[1:]
            # zdo           = self.z_old[0:-1]
            # self.strain   = np.cumsum(zdo-zdn)
            # self.tstrain  = np.sum(zdo-zdn)
            # self.compaction=np.append((zdiffold-zdiffnew)/self.dt*S_PER_YEAR,self.tstrain) #this is cumulative compaction rate in m/yr from 0 to the node specified in depth
            # if not self.snowmeltSec[iii]>0:
            # self.compaction=np.append(0,np.cumsum((self.dz_old[0:compboxes]-self.dz[1:compboxes+1])/self.dt*S_PER_YEAR))

            self.sigma      = (self.mass + self.LWC * RHO_W_KGM) * self.dx * GRAVITY
            self.sigma      = self.sigma.cumsum(axis = 0)
            self.mass_sum   = self.mass.cumsum(axis = 0)
            
            self.bdot_mean  = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] * self.t / (self.age[1:] * RHO_I))))*self.c['stpsPerYear']*S_PER_YEAR
            
            # if self.c['physGrain']: # update grain radius
            #     self.r2, self.dr2_dt    = FirnPhysics(PhysParams).grainGrowth()

            ### write results as often as specified in the init method
            if mtime in self.TWrite:                
                ind         = np.where(self.TWrite == mtime)[0][0]
                mtime_plus1 = self.TWrite[ind] 

                if 'density' in self.output_list:
                    self.rho_out[self.WTracker,:]   = np.append(mtime_plus1, self.rho)
                if 'temperature' in self.output_list:
                    self.Tz_out[self.WTracker,:]    = np.append(mtime_plus1, self.Tz)
                if 'age' in self.output_list:
                    self.age_out[self.WTracker,:]   = np.append(mtime_plus1, self.age/S_PER_YEAR)
                if 'depth' in self.output_list:
                    self.z_out[self.WTracker,:]     = np.append(mtime_plus1, self.z)
                if 'dcon' in self.output_list:    
                    self.D_out[self.WTracker,:]     = np.append(mtime_plus1, self.Dcon)
                if 'climate' in self.output_list:   
                    self.Clim_out[self.WTracker,:]  = np.append(mtime_plus1, [self.bdot[int(iii)], self.Ts[int(iii)]])
                if 'bdot_mean' in self.output_list:   
                    self.bdot_out[self.WTracker,:]  = np.append(mtime_plus1, self.bdot_mean)
                if 'compaction' in self.output_list:    
                    self.comp_out[self.WTracker,:] = np.append(mtime_plus1, self.compaction)
                if 'LWC' in self.output_list:
                    self.LWC_out[self.WTracker,:]   = np.append(mtime_plus1, self.LWC)
                if 'grainsize' in self.output_list:
                    self.r2_out[self.WTracker,:]    = np.append(mtime_plus1, self.r2)
                    self.dr2_dt_out[self.WTracker,:]= np.append(mtime_plus1, self.dr2_dt)
                if 'temp_Hx' in self.output_list:
                    self.Hx_out[self.WTracker,:]    = np.append(mtime_plus1, self.Hx)
                if 'isotopes' in self.output_list:
                    self.iso_out[self.WTracker,:]   = np.append(mtime_plus1, self.del_z)
                if 'viscosity' in self.output_list:
                    self.viscosity = RD['viscosity']   
                    self.viscosity_out[self.WTracker,:] = np.append(mtime_plus1, self.viscosity)

                bcoAgeMart, bcoDepMart, bcoAge830, bcoDep830, LIZAgeMart, LIZDepMart, bcoAge815, bcoDep815  = self.update_BCO()
                intPhi, intPhi_c, z_co  = self.update_DIP()
                dH, dHtot, comp_firn, dHcorr, dHtotcorr    = self.update_dH()
                if mtime==self.TWrite[0]:
                    self.dHAll  = 0 * self.dHAll
                    self.dHAllcorr = 0 * self.dHAllcorr
                    dH          = 0.0
                    dHtot       = 0.0
                    comp_firn   = 0.0
                    dHcorr      = 0.0
                    dHtotcorr   = 0.0

                if 'BCO' in self.output_list:
                    self.BCO_out[self.WTracker,:]       = np.append(mtime_plus1, [bcoAgeMart, bcoDepMart, bcoAge830, bcoDep830, LIZAgeMart, LIZDepMart, bcoAge815, bcoDep815, z_co])                    
                if 'DIP' in self.output_list:
                    self.DIP_out[self.WTracker,:]       = np.append(mtime_plus1, [intPhi, dH, dHtot, comp_firn, dHcorr, dHtotcorr])
                    self.DIPc_out[self.WTracker,:]      = np.append(mtime_plus1, intPhi_c)

                if self.c['FirnAir']:
                    if "gasses" in self.cg['outputs']:
                        for gas in self.cg['gaschoice']:                        
                            self.gas_out[gas][self.WTracker,:]  = np.append(mtime_plus1, self.Gz[gas])
                    if "diffusivity" in self.cg['outputs']:
                        self.diffu_out[self.WTracker,:] = np.append(mtime_plus1, self.diffu)
                    if "gas_age" in self.cg['outputs']:
                        # gas_age0 = self.age/S_PER_YEAR - LIZAgeMart 
                        self.gas_age_out[self.WTracker,:] = np.append(mtime_plus1, gas_age0)
                    if "advection_rate" in self.cg['outputs']:
                        self.w_air_out[self.WTracker,:] = np.append(mtime_plus1, w_p)
                        self.w_firn_out[self.WTracker,:] = np.append(mtime_plus1, self.w_firn)

                self.WTracker = self.WTracker + 1

            if self.doublegrid:
                if self.gridtrack[-1]==2:
                    # print('regridding now at ', iii)
                    self.dz, self.z, self.rho, self.Tz, self.mass, self.sigma, self. mass_sum, self.age, self.bdot_mean, self.LWC, self.gridtrack, self.r2 = regrid(self)
                    if iii<100:
                        tdep = np.where(self.gridtrack==1)[0][-1]
                        print('transition at:', self.z[tdep])

        ##################################
        ##### END TIME-STEPPING LOOP #####
        ##################################

        write_nospin_hdf5(self)

    ###########################
    ##### END time_evolve #####
    ###########################

    def update_BCO(self):
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
        bcoMartRho = 1 / (1 / (917.0) + self.T50* 6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
        phi = 1 - self.rho / RHO_I  # total porosity
        phi[phi <= 0] = 1e-16
        phiC = 1 - bcoMartRho / RHO_I;  # porosity at close off
        phiClosed = 0.37 * phi * (phi / phiC) ** -7.6  # Closed porosity, from Goujon. See Buizert thesis (eq. 2.3) as well

        phiOpen = phi - phiClosed  # open porosity
        
        try:
            ind_co = np.where(phiOpen<=1e-10)[0][0]
            z_co = self.z[ind_co]
        except:
            z_co = -9999

        phiOpen[phiOpen <= 0] = 1.e-10  # don't want negative porosity.

        intPhi = np.sum(phi * self.dz)  # depth-integrated porosity
        intPhi_c = np.cumsum(phi * self.dz)
        # self.intPhiAll.append(intPhi)

        return intPhi, intPhi_c, z_co
    ### end update_DIP ########
    ###########################

    def update_dH(self):
        '''
        updates the surface elevation change
        '''

        self.dH = (self.sdz_new - self.sdz_old) + self.dzNew - (self.iceout*self.t) # iceout has units m ice/year, t is years per time step. 
        # self.dH2 = self.z[-1] - self.z_old[-1] #- (self.iceout*self.t) # alternative method. Should be the same?    
        self.dHAll.append(self.dH)
        self.dHtot = np.sum(self.dHAll)
        
        ### If the bottom of the domain is not the ice density, there is 
        ### compaction that is not accounted for between the bottom of the 
        ### domain and the 917 density horizon.

        iceout_corr = self.iceout*RHO_I/self.rho[-1]
        self.dHcorr = (self.sdz_new - self.sdz_old) + self.dzNew - (iceout_corr*self.t) # iceout has units m ice/year, t is years per time step. 
        self.dHAllcorr.append(self.dHcorr)
        self.dHtotcorr = np.sum(self.dHAllcorr)

        self.comp_firn = self.sdz_new - self.sdz_old #total compaction of just the firn during the previous time step

        # self.dHOut.append(self.dH)
        # self.dHOutC.append(self.dHtot)


        return self.dH, self.dHtot, self.comp_firn, self.dHcorr, self.dHtotcorr

    ###########################

