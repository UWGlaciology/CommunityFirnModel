# from diffusion import Diffusion
from diffusion import heatDiff
from diffusion import isoDiff
from hl_analytic import hl_analytic
# from reader import read_temp
from reader import read_input
# from writer import write_spin
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
# from string import join
import shutil
import time
import h5py

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
		with open(configName, "r") as f:
			jsonString 	= f.read()
			self.c 		= json.loads(jsonString)

		print('Spin run started')
		print("physics are", self.c['physRho'])

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
			input_temp 				= input_temp + K_TO_C
		self.temp0 					= np.mean(input_temp) #Make sure that this is what we want!
		# self.temp0 = 270.0
		# self.temp0 = mean(input_temp[0:12]) #Make sure that this is what we want!

		### accumulation rate
		input_bdot, input_year_bdot = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamebdot']))
		self.bdot0 					= np.mean(input_bdot) #Make sure that this is what we want!
		
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
		############################

		############################
		### get an initial depth/density profile based on H&L analytic solution
		############################
		THL 				= input_temp[0]
		AHL 				= input_bdot[0]
		self.age, self.rho 	= hl_analytic(self.c['rhos0'], self.z, THL, AHL) # self.age is in age in seconds
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
		
		self.dt 	= S_PER_YEAR / self.c['stpsPerYearSpin']
		self.stp 	= int(self.years*S_PER_YEAR/self.dt)
		self.t 		=  1.0 / self.c['stpsPerYearSpin'] # years per time step
		############################

		############################
		### Initial and boundary conditions
		############################
		### Surface temperature for each time step
		self.Ts         = self.temp0 * np.ones(self.stp)
		# self.T_mean     = np.mean(self.Ts) # MS 3/7/17: is this what we want?
		# self.T_mean     = np.mean(self.Tz[self.z<50])
		if self.c['SeasonalTcycle']: #impose seasonal temperature cycle of amplitude 'TAmp', including coreless winter (Orsi)
			self.Ts     = self.Ts + self.c['TAmp'] * (np.cos(2 * np.pi * np.linspace(0, self.years, self.stp )) + 0.3 * np.cos(4 * np.pi * np.linspace(0, self.years, self.stp )))

		### initial temperature profile
		# init_Tz 		= input_temp[0] * np.ones(self.gridLen)
		init_Tz         = np.mean(self.Ts) * np.ones(self.gridLen)

		### Accumulation rate for each time step
		self.bdotSec0   = self.bdot0 / S_PER_YEAR / self.c['stpsPerYearSpin'] # accumulation (m I.E. per second)
		self.bdotSec    = self.bdotSec0 * np.ones(self.stp) # vector of accumulation at each time step

		### Surface isotope values for each time step
		if self.c['isoDiff']:
			try:
				input_iso, input_year_iso = read_input(self.c['InputFileNameIso'])
				del_s0 	= input_iso[0]
			except:
				print('No external file for surface isotope values found, but you specified in the config file that isotope diffusion is on. The model will generate its own synthetic isotope data for you.')
				del_s0 	= -50.0

			self.del_s 	= del_s0 * np.ones(self.stp)
			init_del_z	= del_s0 * np.ones(self.gridLen)
			self.del_z 	= init_del_z
		else:
			self.del_s 	= None
			init_del_z 	= None    
		
		### Surface Density
		self.rhos0      = self.c['rhos0'] * np.ones(self.stp)
		# could configure this so that user specifies vector of surface elevation
		# could add noise too

		### initial mass, stress, and mean accumulation rate
		self.mass       = self.rho * self.dz
		self.sigma      = self.mass * self.dx * GRAVITY
		self.sigma      = self.sigma.cumsum(axis = 0)
		self.mass_sum   = self.mass.cumsum(axis = 0)
		self.bdot_mean  = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] / (self.age[1:] * RHO_I / self.t)))) * self.c['stpsPerYear'] * S_PER_YEAR

		### longitudinal strain rate
		if self.c['strain']:
			self.du_dx 		= np.zeros(self.gridLen)
			self.du_dx[1:] 	= self.c['du_dx']/(S_PER_YEAR)
		
		### initial temperature grid 
		self.Tz         = init_Tz
		self.T_mean     = np.mean(self.Tz[self.z<50])
		self.T10m       = self.T_mean

		### initial grain growth (if specified in config file)
		if self.c['physGrain']:
			if self.c['calcGrainSize']:
				r02     = -2.42e-9 * (self.Ts) + 9.46e-7 # where does this equation come from?
				self.r2 = r02 * np.ones(self.gridLen)
			else:
				self.r2 = np.linspace(self.c['r2s0'], (6 * self.c['r2s0']), self.gridLen)
		else:
			self.r2 = None

		### "temperature history" if using Morris physics
		if self.c['physRho']=='Morris2014':
			# initial temperature history function (units seconds)
			self.Hx 	= np.exp(-110.0e3/(R*init_Tz))*(self.age+self.dt)
			self.THist 	= True
		else:
			self.THist 	= False

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
				'z':            self.z,
				'rhos0':        self.rhos0[iii]
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
				'Crocus':               FirnPhysics(PhysParams).Crocus
			}

			RD 		= physicsd[self.c['physRho']]()
			drho_dt = RD['drho_dt']

			### update density and age of firn
			self.age = np.concatenate(([0], self.age[:-1])) + self.dt
			self.rho = self.rho + self.dt * drho_dt
			
			if self.THist:
				self.Hx = FirnPhysics(PhysParams).THistory()

			### update temperature grid and isotope grid if user specifies
			### should not need to use the diffusion here because temperature is steady state (but could be used)
			# if self.c['heatDiff']:
				# self.Tz, self.T10m = heatDiff(self,iii)

			if self.c['isoDiff']:
				self.del_z 	= isoDiff(self,iii)

			self.T_mean     = np.mean(self.Tz[self.z<50])

			if self.c['strain']: # consider additional change in box height due to longitudinal strain rate
				self.dz 	= ((-self.du_dx)*self.dt + 1)*self.dz 
				self.mass 	= self.mass*((-self.du_dx)*self.dt + 1)

			### update model grid mass, stress, and mean accumulation rate
			dzNew 			= self.bdotSec[iii] * RHO_I / self.rhos0[iii] * S_PER_YEAR
			self.dz 		= self.mass / self.rho * self.dx
			self.dz_old 	= self.dz    
			self.dz 		= np.concatenate(([dzNew], self.dz[:-1]))	  
			self.z 			= self.dz.cumsum(axis = 0)
			self.z 			= np.concatenate(([0], self.z[:-1]))
			self.rho  		= np.concatenate(([self.rhos0[iii]], self.rho[:-1]))
			self.Tz 		= np.concatenate(([self.Ts[iii]], self.Tz[:-1]))
			massNew 		= self.bdotSec[iii] * S_PER_YEAR * RHO_I
			self.mass 		= np.concatenate(([massNew], self.mass[:-1]))
			self.sigma 		= self.mass * self.dx * GRAVITY
			self.sigma 		= self.sigma.cumsum(axis = 0)
			self.mass_sum  	= self.mass.cumsum(axis = 0)
			self.bdot_mean 	= (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] * self.t / (self.age[1:] * RHO_I))))*self.c['stpsPerYear']*S_PER_YEAR

			# update grain radius
			if self.c['physGrain']:
				self.r2 	= FirnPhysics(PhysParams).grainGrowth()

			# write results at the end of the time evolution
			if (iii == (self.stp - 1)):

				rho_time        = np.concatenate(([self.t * iii + 1], self.rho))
				Tz_time         = np.concatenate(([self.t * iii + 1], self.Tz))
				age_time        = np.concatenate(([self.t * iii + 1], self.age))
				z_time          = np.concatenate(([self.t * iii + 1], self.z))

				if self.c['physGrain']:
					r2_time     = np.concatenate(([self.t * iii + 1], self.r2))
				else:
					r2_time     = None
				if self.THist:                
					Hx_time     = np.concatenate(([self.t * iii + 1], self.Hx))
				else:
					Hx_time     = None
				if self.c['isoDiff']:
					iso_time    = np.concatenate(([self.t * iii + 1], self.del_z))
				else:
					iso_time    = None

				write_spin_hdf5(self.c['resultsFolder'], self.c['spinFileName'], self.c['physGrain'], self.THist, self.c['isoDiff'], rho_time, Tz_time, age_time, z_time, r2_time, Hx_time, iso_time)

			####################################
			##### END TIME-STEPPING LOOP #####
			####################################
