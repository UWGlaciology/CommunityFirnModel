# from diffusion import Diffusion
from diffusion import heatDiff
from diffusion import isoDiff
# from reader import read_temp
# from reader import read_bdot
from reader import read_input
from reader import read_init
# from writer import write_nospin
from writer import write_nospin_init
# from writer import write_nospin_BCO
# from writer import write_nospin_LIZ
# from writer import write_nospin_DIP
from writer import write_nospin_hdf5
from physics import *
from constants import *
import numpy as np
import csv
import json
import sys
import math
from shutil import rmtree
import os
from string import join
import shutil
import time
import inspect
import h5py

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
    :returns D_surf: diffusivity tracker
                (unit: ???, type: array of floats)
    '''

    def __init__(self, configName):
        '''
        Sets up the initial spatial grid, time grid, accumulation rate, age, density, mass, stress, temperature, and diffusivity of the model run
        :param configName: name of json config file containing model configurations
        '''

        # load in json config file and parses the user inputs to a dictionary
        with open(configName, "r") as f:
            jsonString      = f.read()
            self.c          = json.loads(jsonString)
        print "Main run starting"
        print "physics are", self.c['physRho']

        ### read in initial depth, age, density, temperature from spin-up results
        initDepth   = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'depthSpin')
        initAge     = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'ageSpin')
        initDensity = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'densitySpin')
        initTemp    = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'tempSpin')

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
        input_temp, input_year_temp = read_input(self.c['InputFileNameTemp'])
        if input_temp[0] < 0.0:
            input_temp = input_temp + K_TO_C
        input_bdot, input_year_bdot = read_input(self.c['InputFileNamebdot'])
        input_bdot[input_bdot<=0.0] = 0.01
        # if self.c['variable_srho']:
        #     input_srho, input_year_srho = read_input(self.c['InputFileNamesrho'])

        ### set up time stepping
        # year to start and end, from the input file. If inputs have different start/finish, take only the overlapping times
        yr_start        = max(input_year_temp[0], input_year_bdot[0])   # start year
        yr_end          = min(input_year_temp[-1], input_year_bdot[-1]) # end year
        
        self.years      = (yr_end - yr_start) * 1.0 
        self.dt         = S_PER_YEAR / self.c['stpsPerYear']
        self.stp        = int(self.years * S_PER_YEAR/self.dt + 1)       # total number of time steps, as integer
        # self.modeltime  = np.linspace(yr_start, yr_end, self.stp + 1)   # vector of time of each model step
        self.modeltime  = np.linspace(yr_start, yr_end, self.stp)

        # self.dt         = self.years * S_PER_YEAR / self.stp            # size of time steps, seconds
        self.t          = 1.0 / self.c['stpsPerYear']                   # years per time step

        ### Temperature
        self.Ts         = np.interp(self.modeltime, input_year_temp, input_temp) # surface temperature interpolated to model time
        # print len(self.Ts)
        # print self.Ts[0:13]
        if self.c['SeasonalTcycle']: #impose seasonal temperature cycle of amplitude 'TAmp'
            self.Ts         = self.Ts + self.c['TAmp'] * (np.cos(2 * np.pi * np.linspace(0, self.years, self.stp)) + 0.3 * np.cos(4 * np.pi * np.linspace(0, self.years, self.stp)))
            # print self.Ts[0:13]

        ### Accumulation
        self.bdot       = np.interp(self.modeltime, input_year_bdot, input_bdot) # interpolate accumulation rate to model time ???Should this be nearest?
        # self.bdotSec    = self.bdot / S_PER_YEAR / (self.stp / self.years) # accumulation rate in per second
        self.bdotSec   = self.bdot / S_PER_YEAR / self.c['stpsPerYear'] # accumulation (per second)
        self.iceout     = np.mean(self.bdot) # this is the rate of ice flow advecting out of the column
        

        ##### Isotopes ###########
        if self.c['isoDiff']:
            init_del_z    = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'IsoSpin')
            try:
                input_iso, input_year_iso = read_input(self.c['InputFileNameIso'])
                self.del_s  = np.interp(self.modeltime, input_year_iso, input_iso)
                # del_s0 = input_iso[0]
            except:
                print 'No external file for surface isotope values found, but you specified in the config file that isotope diffusion is on. The model will generate its own synthetic isotope data for you.'
                # del_s0 = -50.0
                ar1 = 0.9   # red noise memory coefficient
                std_rednoise = 2    # red noise standard deviation
                self.del_s = std_rednoise*np.random.randn(self.stp)    # white noise
                for x in range(1,self.stp):
                    self.del_s[x] = self.del_s[x-1]*ar1 + np.random.randn()  # create red noise from white
                self.del_s = self.del_s - 50
                #impose seasonal isotope cycle
                # self.del_s = self.del_s + 5 * (np.cos(2 * np.pi * np.linspace(0, self.years, self.stp )) + 0.3 * np.cos(4 * np.pi * np.linspace(0, self.years, self.stp )))
        ###########################
 
        try:
            if self.c['variable_srho']:
                self.rhos0      = np.interp(self.modeltime, input_year_srho, input_srho)
            else:
                self.rhos0      = self.c['rhos0'] * np.ones(self.stp)       # density at surface
        except:
            print "you should alter the json to include variable_srho"
            self.rhos0      = self.c['rhos0'] * np.ones(self.stp)       # density at surface


        self.D_surf     = self.c['D_surf'] * np.ones(self.stp)      # layer traking routine (time vector). 

        self.Dcon       = self.c['D_surf'] * np.ones(self.gridLen)  # layer tracking routine (initial depth vector)

        # set up vector of times data will be written
        # Tind = np.nonzero(selfself.modeltime>=1958.0)[0][0]

        self.TWrite     = self.modeltime[0::self.c['TWriteInt']]
        # self.TWrite_out = self.TWrite
        TWlen           = len(self.TWrite) #- 1
        self.WTracker        = 1

        # set up initial mass, stress, and mean accumulation rate
        self.mass       = self.rho * self.dz
        self.sigma      = self.mass * self.dx * GRAVITY
        self.sigma      = self.sigma.cumsum(axis = 0)
        self.mass_sum   = self.mass.cumsum(axis = 0)

        ### mean accumulation over the lifetime of the parcel
        # self.bdot_mean  = np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] / (self.age[1:] * RHO_I / self.t)))
        self.bdot_mean  = np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] / (self.age[1:] * RHO_I / self.t)))

        ### set up longitudinal strain rate
        if self.c['strain']:
            self.du_dx = np.zeros(self.gridLen)
            self.du_dx[1:] = self.c['du_dx']/(S_PER_YEAR)
        
        # set up class to handle heat/isotope diffusion using user provided data for initial temperature vector
        # self.diffu      = Diffusion(self.z, self.stp, self.gridLen, initTemp[1:], init_del_z[1:]) # [1:] because first element is a time stamp
        # self.T_mean     = self.diffu.T10m # initially the mean temp is the same as the surface temperature
        self.Tz         = initTemp[1:]
        self.T_mean     = self.Tz[0]
        self.T10m       = self.T_mean

        # set up initial values for density, temperature, age, depth, diffusivity, model climate, and accumulation to write
        # self.write_hdf5 = self.c['writehdf5']

        # self.RD = {}

        self.rho_out = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
        self.Tz_out = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
        self.age_out = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
        self.z_out = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
        self.D_out = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
        self.bdot_out = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
        self.Clim_out = np.zeros((TWlen+1,3),dtype='float32')
        self.DIP_out = np.zeros((TWlen+1,4),dtype='float32')
        self.LIZ_out = np.zeros((TWlen+1,3),dtype='float32')
        self.BCO_out = np.zeros((TWlen+1,5),dtype='float32')

        # print 'TWlen', TWlen+1

        rho_time        = np.append(self.modeltime[0], self.rho)
        Tz_time         = np.append(self.modeltime[0], self.Tz)
        age_time        = np.append(self.modeltime[0], self.age/S_PER_YEAR)
        z_time          = np.append(self.modeltime[0], self.z)
        D_time          = np.append(self.modeltime[0], self.Dcon)
        Clim_time       = np.append(self.modeltime[0], [self.bdot[0], self.Ts[0]])  # not sure if bdot or bdotSec
        bdot_time       = np.append(self.modeltime[0], self.bdot_mean)

        self.rho_out[0,:]        = rho_time
        self.Tz_out[0,:]         = Tz_time
        self.age_out[0,:]        = age_time
        self.z_out[0,:]          = z_time
        self.D_out[0,:]          = D_time
        self.Clim_out[0,:]       = Clim_time
        self.bdot_out[0,:]       = bdot_time

        print 'rho_out size (MB):', self.rho_out.nbytes/1.0e6

        # self.rho_out[0,:]        = np.append(self.modeltime[0], self.rho)
        # self.Tz_out[0,:]         = np.append(self.modeltime[0], self.diffu.Tz)
        # self.age_out[0,:]        = np.append(self.modeltime[0], self.age)
        # self.z_out[0,:]          = np.append(self.modeltime[0], self.z)
        # self.D_out[0,:]          = np.append(self.modeltime[0], self.Dcon)
        # self.Clim_out[0,:]       = np.append(self.modeltime[0], [self.bdot[0], self.Ts[0]])  # not sure if bdot or bdotSec
        # self.bdot_out[0,:]       = np.append(self.modeltime[0], self.bdot_mean)

        # set up initial grain growth (if specified in config file)
        if self.c['physGrain']:
            initr2              = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'r2Spin')
            self.r2             = initr2[1:]
            r20                 = self.r2
            self.r2_out         = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
            r2_time             = np.append(self.modeltime[0], self.r2)
            self.r2_out[0,:]    = r2_time
            
        else:
            self.r2_out         = None
            self.r2             = None

        if self.c['physRho'] == 'Morris2013':
            self.THist          = True
            initHx              = read_init(self.c['resultsFolder'], self.c['spinFileName'], 'r2Spin')
            self.Hx             = initHx[1:]
            self.Hx_out         = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
            Hx_time             = np.append(self.modeltime[0], self.Hx)
            self.Hx_out[0,:]    = Hx_time       
        else:
            self.THist          = False
            self.Hx_out         = None

        if self.c['isoDiff']:
            self.del_z          = init_del_z[1:]
            self.iso_out        = np.zeros((TWlen+1,len(self.dz)+1),dtype='float32')
            del_z_time          = np.append(self.modeltime[0], self.del_z)
            self.iso_out[0,:]   = del_z_time
        else:
            self.iso_out        = None

        ### write initial values to the results folder
        # if not self.write_hdf5:
            # write_nospin_init(self.c['resultsFolder'], self.c['physGrain'], self.THist, rho_time, Tz_time, age_time, z_time, D_time, Clim_time, bdot_time, r2_time, Hx_time)

        # write_nospin_hdf5(self.c['resultsFolder'], self.c['physGrain'], self.THist, rho_time, Tz_time, age_time, z_time, D_time, Clim_time, bdot_time, r2_time, Hx_time)

        # set up initial values for bubble close-off depth & age, lock-in zone depth & age, and depth integrated porosity
        # self.bcoAgeMartAll = []
        # self.bcoDepMartAll = []
        # self.bcoAge815All  = []
        # self.bcoDep815All  = []
        # self.LIZAgeAll     = []
        # self.LIZDepAll     = []
        # self.intPhiAll     = []
        self.dHAll         = []
        # self.dHOut         = []
        # self.dHOutC        = []

        bcoAgeMart, bcoDepMart, bcoAge815, bcoDep815 = self.update_BCO()
        LIZAgeMart, LIZDepMart = self.update_LIZ()
        intPhi = self.update_DIP()
        
        self.dHAll.append(0)
        dHOut = 0
        dHOutC = 0

        BCO_time        = np.append(self.modeltime[0], [bcoAgeMart, bcoDepMart, bcoAge815, bcoDep815])
        LIZ_time        = np.append(self.modeltime[0], [LIZAgeMart, LIZDepMart])
        DIP_time        = np.append(self.modeltime[0], [intPhi, dHOut, dHOutC])

        self.BCO_out[0,:]       = BCO_time
        self.LIZ_out[0,:]       = LIZ_time
        self.DIP_out[0,:]       = DIP_time

    ####################    
    ##### END INIT #####
    ####################

    def time_evolve(self):
        '''
        Evolve the spatial grid, time grid, accumulation rate, age, density, mass, stress, temperature, and diffusivity through time
        based on the user specified number of timesteps in the model run. Updates the firn density using a user specified 
        '''
        self.steps = 1 / self.t
        # if not self.c['physGrain']:
        #     r2_time = None
        # if not self.THist:
        #     Hx_time = None

        start_time=time.time() # this is a timer to keep track of how long the model run takes.
        ####################################
        ##### START TIME-STEPPING LOOP #####
        ####################################
        for iii in xrange(self.stp):
            mtime = self.modeltime[iii]
            
            # print self.Ts[iii]
            # ### set up longitudinal strain rate
            # self.du_dx = np.zeros(self.gridLen)
            # if mtime  < (2500.):
            #     self.du_dx[0:] = (0.)/(S_PER_YEAR)
            # else:
            #     self.du_dx[0:] = (10**-3.)/(S_PER_YEAR)
               
            # the parameters that get passed to physics
            PhysParams = {
                'iii':          iii,
                'steps':        self.steps,
                'gridLen':      self.gridLen,
                'bdotSec':      self.bdotSec,
                'bdot_mean':    self.bdot_mean,
                'bdot_type':    self.c['bdot_type'],
                'Tz':           self.Tz,
                'T_mean':       self.T_mean,
                'rho':          self.rho,
                'sigma':        self.sigma,
                'dt':           self.dt,
                'Ts':           self.Ts,
                'r2':           self.r2,
                'age':          self.age,
                'physGrain':    self.c['physGrain'],
                'calcGrainSize':self.c['calcGrainSize'],
                'z':            self.z,
                'rhos0':        self.rhos0[iii]
            }

            if self.THist: #add Hx to dictionary if physics is Morris
                PhysParams['Hx']=self.Hx

            # choose densification-physics based on user input
            physicsd = {
                'HLdynamic':            FirnPhysics(PhysParams).HL_dynamic,
                'HLSigfus':             FirnPhysics(PhysParams).HL_Sigfus,
                'Barnola1991':          FirnPhysics(PhysParams).Barnola_1991,
                'Li2004':               FirnPhysics(PhysParams).Li_2004,
                'Li2011':               FirnPhysics(PhysParams).Li_2011,
                'Ligtenberg2011':       FirnPhysics(PhysParams).Ligtenberg_2011,
                'Arthern2010S':         FirnPhysics(PhysParams).Arthern_2010S,
                'Simonsen2013':         FirnPhysics(PhysParams).Simonsen_2013,
                'Morris2013':           FirnPhysics(PhysParams).Morris_HL_2013,
                'Helsen2008':           FirnPhysics(PhysParams).Helsen_2008,
                'Arthern2010T':         FirnPhysics(PhysParams).Arthern_2010T,
                'Spencer2001':          FirnPhysics(PhysParams).Spencer_2001,
                'Goujon2003':           FirnPhysics(PhysParams).Goujon_2003,
                'KuipersMunneke2015':   FirnPhysics(PhysParams).KuipersMunneke_2015,
            }

            try:
                drho_dt = physicsd[self.c['physRho']]()
            except KeyError:
                print "Error at line ", info.lineno

            # update density and age of firn
            self.age = np.concatenate(([0], self.age[:-1])) + self.dt
            self.rho = self.rho + self.dt * drho_dt
            # print self.rho
            
            self.Dcon = np.concatenate(([self.D_surf[iii]], self.Dcon[:-1]))

            if self.THist:
                self.Hx = FirnPhysics(PhysParams).THistory()

            # update temperature grid and isotope grid if user specifies
            if self.c['heatDiff']:
                # self.diffu.heatDiff(self.z, self.dz, self.Ts[iii], self.rho, self.dt)
                self.Tz, self.T10m = heatDiff(self,iii)
            else:
                self.Tz = np.concatenate(([self.Ts[iii]], self.Tz[:-1]))

            if self.c['isoDiff']:
                # self.diffu.isoDiff(self.z, self.dz, self.del_s[iii], self.rho, self.c['iso'], self.dt)
                self.del_z = isoDiff(self,iii)
                # print 'del_z', self.del_z[0:2]
            melt = False

            if melt:
                print "Meltwater percolation is still under development. Run without melt for now."
                sys.exit()

            else:
            # MS 2/10/17: should double check that everything occurs in correct order in time step (e.g. adding new box on, calculating dz, etc.) 
                ##### update model grid
                self.dz_old = self.dz
                self.sdz_old = np.sum(self.dz) # old total column thickness
                self.z_old = self.z
                self.dzNew = self.bdotSec[iii] * RHO_I / self.rhos0[iii] * S_PER_YEAR
                self.dz = self.mass / self.rho * self.dx

                if self.c['strain']:
                    self.dz = ((-self.du_dx)*self.dt + 1)*self.dz
                
                self.sdz_new = np.sum(self.dz) #total column thickness after densification, before new snow added               
                self.dz = np.concatenate(([self.dzNew], self.dz[:-1]))
                self.z = self.dz.cumsum(axis = 0)
                self.z = np.concatenate(([0], self.z[:-1]))
                self.rho  = np.concatenate(([self.rhos0[iii]], self.rho[:-1]))
                ##### update mass, stress, and mean accumulation rate
                if self.c['strain']:
                	self.mass = self.mass*((-self.du_dx)*self.dt + 1)
                massNew = self.bdotSec[iii] * S_PER_YEAR * RHO_I
                self.mass = np.concatenate(([massNew], self.mass[:-1]))


            #### find the compaction rate
            zdiffnew=(self.z[1:]-self.z[1])
            zdiffold=(self.z_old[0:-1]-self.z_old[0])
            zdn=self.z[1:]
            zdo=self.z_old[0:-1]
            self.strain=np.cumsum(zdo-zdn)
            self.tstrain=np.sum(zdo-zdn)
            self.compaction_rate=(zdiffold-zdiffnew)/self.dt*S_PER_YEAR #this is cumulative compaction rate in m/yr from 0 to the node specified in depth
            ####

            self.sigma = self.mass * self.dx * GRAVITY
            self.sigma = self.sigma.cumsum(axis = 0)
            self.mass_sum  = self.mass.cumsum(axis = 0)
            self.bdot_mean = np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] * self.t / (self.age[1:] * RHO_I)))

            # update grain radius
            if self.c['physGrain']:
                self.r2 = FirnPhysics(PhysParams).grainGrowth()

            # print 'mtime', mtime
            # write results as often as specified in the init method
            if mtime in self.TWrite:
                
                ind = np.where(self.TWrite == mtime)[0][0]
                # print ind
                # print self.TWrite[ind]
                mtime_plus1 = self.TWrite[ind] #+ (self.TWrite[1]-self.TWrite[0])
                # print mtime_plus1
                # print mtime_plus1
            # if [True for ii in self.TWrite if ii == mtime] == [True]:
                # print self.WTracker
                self.rho_out[self.WTracker,:] = np.append(mtime_plus1, self.rho)
                # print '!!!', self.rho_out[0:2,0:10]
                self.Tz_out[self.WTracker,:]   = np.append(mtime_plus1, self.Tz)
                self.age_out[self.WTracker,:]  = np.append(mtime_plus1, self.age/S_PER_YEAR)
                self.z_out[self.WTracker,:]    = np.append(mtime_plus1, self.z)
                self.D_out[self.WTracker,:] = np.append(mtime_plus1, self.Dcon)
                self.Clim_out[self.WTracker,:] = np.append(mtime_plus1, [self.bdot[int(iii)], self.Ts[int(iii)]])
                self.bdot_out[self.WTracker,:] = np.append(mtime_plus1, self.bdot_mean)
                if self.c['physGrain']:
                    self.r2_out[self.WTracker,:] = np.append(mtime_plus1, self.r2)
                if self.THist:
                    self.Hx_out[self.WTracker,:] = np.append(mtime_plus1, self.Hx)
                if self.c['isoDiff']:
                    self.iso_out[self.WTracker,:] = np.append(mtime_plus1, self.del_z)


                # self.update_BCO()
                # self.update_LIZ()
                # self.update_DIP()
                # self.update_dH()

                bcoAgeMart, bcoDepMart, bcoAge815, bcoDep815 = self.update_BCO()
                LIZAgeMart, LIZDepMart = self.update_LIZ()
                intPhi = self.update_DIP()
                dH, dHtot = self.update_dH()

                # BCO_time        = np.append(self.modeltime[0], [bcoAgeMart, bcoDepMart, bcoAge815, bcoDep815])
                # LIZ_time        = np.append(self.modeltime[0], [LIZAgeMart, LIZDepMart])
                # DIP_time        = np.append(self.modeltime[0], [intPhi, dH, dHtot])

                self.BCO_out[self.WTracker,:]       = np.append(mtime_plus1, [bcoAgeMart, bcoDepMart, bcoAge815, bcoDep815])
                self.LIZ_out[self.WTracker,:]       = np.append(mtime_plus1, [LIZAgeMart, LIZDepMart])
                self.DIP_out[self.WTracker,:]       = np.append(mtime_plus1, [intPhi, dH, dHtot])

                self.WTracker = self.WTracker + 1

        ##################################
        ##### END TIME-STEPPING LOOP #####
        ##################################
        # print mtime[-5:]
        # print 'Twrite', self.TWrite[-5:]
        # print 'modeltime', self.modeltime[-5:]


        # write_nospin_hdf5(self.c['resultsFolder'], self.c['physGrain'], self.THist, self.rho_out, self.Tz_out, self.age_out, self.z_out, self.D_out, self.Clim_out, self.bdot_out, self.r2_out, self.Hx_out)
        

        write_nospin_hdf5(self)

        # write BCO, LIZ, DIP at the end of the time evolution
        # write_nospin_BCO(self.c['resultsFolder'], self.bcoAgeMartAll, self.bcoDepMartAll, self.bcoAge815All, self.bcoDep815All,self.modeltime,self.TWrite)
        # write_nospin_LIZ(self.c['resultsFolder'], self.LIZAgeAll, self.LIZDepAll,self.modeltime,self.TWrite)
        # write_nospin_DIP(self.c['resultsFolder'], self.intPhiAll, self.dHOut, self.dHOutC, self.modeltime,self.TWrite)
        # write_nospin_dH(self.c['resultsFolder'], self.dsurf, self.dsurftot, self.dsurfAll, self.modeltime,self.TWrite)
    ##### END time_evolve #####

    def update_BCO(self):
        '''
        Updates the bubble close-off depth and age based on the Martinerie criteria as well as through assuming the critical density is 815 kg/m^3
        '''

        bcoMartRho = 1 / (1 / (917.0) + self.T10m * 6.95E-7 - 4.3e-5)  # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
        bcoAgeMart = min(self.age[self.rho >= bcoMartRho]) / S_PER_YEAR  # close-off age from Martinerie
        bcoDepMart = min(self.z[self.rho >= (bcoMartRho)])
        # self.bcoAgeMartAll.append(bcoAgeMart)  # age at the 815 density horizon
        # self.bcoDepMartAll.append(bcoDepMart)  # this is the 815 close off depth

        # bubble close-off age and depth assuming rho_crit = 815kg/m^3
        bcoAge815 = min(self.age[self.rho >= (RHO_2)]) / S_PER_YEAR  # close-off age where rho = 815 kg m^-3
        bcoDep815 = min(self.z[self.rho >= (RHO_2)])
        # self.bcoAge815All.append(bcoAge815)  # age at the 815 density horizon
        # self.bcoDep815All.append(bcoDep815)  # this is the 815 close off depth

        return bcoAgeMart, bcoDepMart, bcoAge815, bcoDep815

    #### end update_BCO

    def update_LIZ(self):
        '''
        Updates the lock-in zone depth and age
        '''

        bcoMartRho = 1 / (1 / (917.0) + self.T10m * 6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
        LIZMartRho = bcoMartRho - 14.0  # LIZ depth (Blunier and Schwander, 2000)
        self.LIZAgeMart = min(self.age[self.rho > LIZMartRho]) / S_PER_YEAR  # lock-in age
        self.LIZDepMart = min(self.z[self.rho >= (LIZMartRho)])  # lock in depth
        # self.LIZAgeAll.append(self.LIZAgeMart)
        # self.LIZDepAll.append(self.LIZDepMart)

        return self.LIZAgeMart, self.LIZDepMart
    #### end update_LIZ 

    def update_DIP(self):
        '''
        Updates the depth-integrated porosity
        '''

        bcoMartRho = 1 / (1 / (917.0) + self.T10m * 6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
        phi = 1 - self.rho / RHO_I  # total porosity
        phi[phi <= 0] = 1e-16
        phiC = 1 - bcoMartRho / RHO_I;  # porosity at close off
        phiClosed = 0.37 * phi * (phi / phiC) ** -7.6  # Closed porosity, from Goujon. See Buizert thesis (eq. 2.3) as well

        phiOpen = phi - phiClosed  # open porosity
        phiOpen[phiOpen <= 0] = 1.e-10  # don't want negative porosity.

        intPhi = np.sum(phi * self.dz)  # depth-integrated porosity
        # self.intPhiAll.append(intPhi)

        return intPhi
    #### end update_DIP

    def update_dH(self):
        '''
        updates the surface elevation change
        '''

        # self.dH = (self.sdz_new-self.sdz_old)+self.dzNew-(self.bdot_mean[0]*S_PER_YEAR) #
        
        self.dH = (self.sdz_new-self.sdz_old)+self.dzNew-(self.iceout/(self.rho[-1]/RHO_I))*self.t #

        self.dHAll.append(self.dH)

        self.dHtot = np.sum(self.dHAll)

        # self.dHOut.append(self.dH)
        # self.dHOutC.append(self.dHtot)

        return self.dH, self.dHtot

