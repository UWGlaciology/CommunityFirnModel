import numpy as np 
from solver import solver
from solver import transient_solve_TR
# from Gasses import gasses 
# from Sites import sites 
from reader import read_input
import json
import scipy.interpolate as interpolate
from constants import *

class FirnAir:
	def __init__(self,AirconfigName,input_year_gas,z,modeltime,Tz, rho, dz):
		'''
		Initialize Firn Air class
		'''
		with open(AirconfigName, "r") as f:
			jsonString 		= f.read()
			self.cg         = json.loads(jsonString)

		# self.gaschoice_all=self.cg["gaschoice"]
		# nogas=len(gaschoice_all)

		# self.gaschoice = self.cg["gaschoice"]
		# input_gas, input_year_gas = read_input(os.path.join(self.c['InputFileFolder'],self.c['InputFileNamebdot']))
		input_gas = np.ones_like(input_year_gas)
		Gsf = interpolate.interp1d(input_year_gas,input_gas,'nearest',fill_value='extrapolate')
		self.Gs = Gsf(modeltime)
		self.Gz	= self.Gs[0]*np.ones_like(z)
		
		self.gam_x, self.M, self.deltaM, self.d_0, self.omega = gasses(self.cg['gaschoice'], Tz[0], P_0)
		self.czd = 4.0
		self.p_a = 1.0e5
		self.air_pressure = self.p_a * np.ones_like(z)
		self.air_pressure_0 = np.copy(self.air_pressure)
		self.air_volume = (1 - rho / RHO_I) * dz

	def diffusivity(self):

		'''
		D_0 is CO2 in free air.
		gam_x is the diffusivity relative to CO2
		D_x=D_0*gam_x is the free-air (or any) diffusity for that species
		'''
		## Constants
		d_eddy_sc	= self.d_0 #Eddy diffusivity in the convective zone

		d_ind		= np.min(np.where(self.z>self.z_co)) #indices of nodes past close-off depth
		
		##### Parameterizations for diffusivity #####
		if self.cg['Diffu_param'] == "Severinghaus": # Use Severinghaus relationship from Cuffey and Paterson
			
			# d_0			= d_0*1.7
			diffu_full 		= 0.1 * self.gam_x * self.d_0 * ((P_0/self.p_a) * (self.Tz/273.15)**1.85 * (2.00 * (1-(self.rho/RHO_I))-0.167))    
			# diffu_full 	= diffu_full - diffu_full[d_ind]
			# iind 			= np.nonzero(diffu_full == np.min(diffu_full[diffu_full>0.0]))
			# diffu_full[diffu_full<=0] = diffu_full[iind]
			# diffu_full 	= diffu_full * 10.0
		
		elif self.cg['Diffu_param'] == "Schwander": # Use Schwander 1988, Eq. 2 Diffusivity (does not work very well) use 4e2 for d_0
			k_sch 			= P_0 / self.p_a * (self.Tz / 253.16)**1.85 # Constant given in Schwander
			# diffu_full 	= 3.72 * 0.5 * k_sch * (23.7 * por_tot-2.84) * 31.5 # Schwander' diffusivity relationship (for CO2). 31.5 is unit conversion. Added extra 3.72* 9/12/13
			diffu_full 		= 2.0 * k_sch * (23.7 * self.por_tot - 2.84) / (1000**2) # Schwander' diffusivity relationship (for CO2). 1/1000**2 is unit conversion. Added extra 3.72* 9/12/13
			# ind 			= np.nonzero(z>LIZ)
			# diffu_full[ind] = 0.001
			# diffu_full 	= diffu_full - diffu_full[d_ind]
					
		elif self.cg['Diffu_param'] == "Freitag": # Use Freitag, 2002, Eq 15 Diffusivity use 9e2 for d_0
			# d_0 			= d_0*4.9    
			diffu_full 		= 1.0 * self.gam_x * self.d_0 * self.por_op ** 2.1
			# diffu_full 	= diffu_full - diffu_full[d_ind]
			
		elif self.cg['Diffu_param']=="Witrant":    ### Use Witrant, 2012
			diffu_full = self.gam_x * self.d_0 * (2.5 * self.por_op - 0.31) * (self.Tz / 273.15)**(1.8) * P_0 / self.p_a
		
		diffu_full[diffu_full<0] = 1.e-40
				  		
		###### Add in high diffusivity in convective zone and low diffusivity below LIZ
		
		###Convective zone###
		d_eddy 			= np.zeros(np.size(diffu_full))
		ind 			= np.nonzero(self.z<self.czd)
		d_eddy_surf		= 2.426405E-5 #Kawamura, 2006
		# H_scale 		= czd
		d_eddy_up 		= d_eddy_surf * np.exp(-1*self.z/self.czd)
		
		#### Lock-in zone physics###
		if self.cg['lockin']:
			ind 			= np.flatnonzero(self.z>self.LIZ)
			ind2 			= np.flatnonzero(self.z<self.z_co)
			ind3 			= np.intersect1d(ind,ind2)
			d_eddy[ind3] 	= diffu_full[ind] #set eddy diffusivity in LIZ equal to diffusivity at LIZ
			diffu_full[ind]	= 1e-40 #set molecular diffusivity equal to zero for "main" diffusivity after LIZ - eddy diffusivity term drives diffusion below
			d_eddy 			= d_eddy + d_eddy_up #make eddy diffusivity vector have convective and lock-in zone values
			
		else:
			diffu_full[np.flatnonzero(self.z>self.z_co)] = 1.0e-40
			d_eddy 			= d_eddy_up #eddy diffusivity includes only convective zone
		 
		diffu=diffu_full

		
		# dd={}
		# dd['Fre']=diffu_full_fre
		# dd['Sev']=diffu_full_sev
		# dd['Sch']=diffu_full_sch
		# dd['Wit']=diffu_full_wit
		# dd['data']=diffu_full_data
		# dd['eddy']=d_eddy
		# dd['porop']=por_op
		# dd['porcl']=por_cl
		# dd['portot']=por_tot  
		
		return diffu , d_eddy #, dd

	def porosity(self): #,rho,T
		
		self.bcoRho = 1/( 1/(RHO_I) + self.Tz[100]*6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
		self.LIDRho = self.bcoRho - 14.0 #LIZ depth (Blunier and Schwander, 2000)
		## Porosity, from Goujon et al., 2003, equations 9 and 10
		self.por_tot = 1-self.rho/RHO_I # Total porosity

		# rho_co = 910.0
		self.rho_co = self.bcoRho #use Martinerie close-off criteria
		#rho_co = 0.815 # User chosen close off-density (put in site-specific in sites?)

		self.por_co = 1 - self.rho_co/RHO_I # Porosity at close-off
		alpha = 0.37 # constant determined in Goujon
		self.por_cl = alpha*self.por_tot*(self.por_tot/self.por_co)**(-7.6)
		ind=self.por_cl>self.por_tot
		self.por_cl[ind]=self.por_tot[ind]
		
		#por_cl[por_cl>1]=1
		#por_cl = por_cl*por_tot # Final closed porosity
		self.por_op = self.por_tot - self.por_cl # Open Porosity
		self.por_op[self.por_op<=0] = 1.0e-25
		
		return self.rho_co, self.por_co, self.por_tot, self.por_cl, self.por_op, self.bcoRho, self.LIDRho


	def firn_air_diffusion(self,PhysParams,iii):

		for k,v in list(PhysParams.items()):
			setattr(self,k,v)

		nz_P = len(self.z)
		nz_fv = nz_P - 2
		nt = 1

		z_edges_vec = self.z[1:-2] + self.dz[2:-1] / 2
		z_edges_vec = np.concatenate(([self.z[0]], z_edges_vec, [self.z[-1]]))
		z_P_vec = self.z
		# phi_s = self.Ts[iii]
		phi_s = self.Gz[0]
		phi_0 = self.Gz

		# K_ice = 9.828 * np.exp(-0.0057 * phi_0)
		# K_firn = K_ice * (self.rho / 1000) ** (2 - 0.5 * (self.rho / 1000))
		self.rho_co, self.por_co, self.por_tot, self.por_cl, self.por_op, self.bcoRho, self.LIDRho = self.porosity() #self.rho, self.Tz

		
		self.air_pressure_old = np.copy(self.air_pressure)
		self.air_volume_old = np.copy(self.air_volume)

		self.air_volume = self.por_op * self.dz
		self.air_pressure = self.air_pressure_0 * self.air_volume_old / self.air_volume # assume air pressure is atmos in entire column

		self.pressure_grad = np.gradient(air_pressure,self.dz) 

		self.z_co = min(self.z[self.rho>=(self.bcoRho)]) #close-off depth; bcoRho is close off density
		self.LIZ = min(self.z[self.rho>=(self.LIDRho)]) #lock in depth; LIDRho is lock-in density

		if iii==900:
			print(max(self.z[self.rho<(self.LIDRho)]))
			print('z_co',self.z_co)
			print('lockin',self.LIZ)

		self.diffu, self.d_eddy = self.diffusivity()

		airdict = {
			'd_eddy': 		self.d_eddy,
			'por_op': 		self.por_op,
			'Tz': 			self.Tz,
			'deltaM': 		self.deltaM,
			'omega': 		self.omega,
			'dz': 			self.dz,
			'rho':			self.rho,
			'gravity': 		self.cg['gravity'],
			'thermal': 		self.cg['thermal'],
			'air_pressure': self.air_pressure,
			'pressure_grad':self.pressure_grad,
			'z': 			self.z, 
			'dt': 			self.dt,
			'z_co': 		self.z_co
			}

		self.Gz = transient_solve_TR(z_edges_vec, z_P_vec, nt, self.dt, self.diffu, phi_0, nz_P, nz_fv, phi_s, airdict)
		self.Gz = np.concatenate(([self.Gs[iii]], self.Gz[:-1]))
		
		return self.Gz

def gasses(gaschoice, T, p_a):
	

	#d_0 = 5.e2 # Free air diffusivity, CO2, m**2/yr Schwander, 1988 reports 7.24 mm**2/s =379 m**2/yr
	d_0 = 1.6e-5 # m^2/s :wikipedia value. changed 9/27/13  Schwander, 1988 reports 7.24 mm**2/s = 7.24e-6 m**2/yr
	M_air = 28.97e-3 #kg/mol
	D_ref_CO2 = 5.75E-10*T**1.81*(101325/p_a) #Christo Thesis, appendix A3
	print('gas choice is ', gaschoice)
	
	if gaschoice == ['CO2']:
		gam_x = 1. #free-air diffusivity relative to CO2. Unitless (gamma in Buizert thesis, page 13).
		M = 44.01e-3 # molecular mass, kg/mol
		decay = 0.
		omega = 0.0
		
		#if hemisphere == 'SOUTH':
		#    conc1=loadtxt(os.path.join(DataPath,'CO2_NH_history.txt'),skiprows=2) #load data: atmospheric CO2 history.
		#    
		#    firn_meas=loadtxt(os.path.join(DataPath,'CO2_samples_NEEM.txt'),skiprows=2)
		#
		#elif hemisphere == 'NORTH':
		#    conc1=loadtxt(os.path.join(DataPath,'CO2_SH_history.txt'),skiprows=2) #load data: atmospheric CO2 history.
		#    firn_meas=loadtxt(os.path.join(DataPath,'CO2_samples_WAIS.txt'),skiprows=2)   
		#
		#elif hemisphere == 'SCENARIO':
		#    conc1=loadtxt(os.path.join(DataPath,'RampUp2.txt'),skiprows=2) #load data: atmospheric CO2 history.
		#    firn_meas=loadtxt(os.path.join(DataPath,'CO2samples_WAIS.txt'),skiprows=2)   
		#    #conc1=conc1[0:1996,:] # May or may not need this to get time to work...
		
			
	elif gaschoice == ['CH4']:
		gam_x = 1.367
		M = 16.04e-3
		decay = 0.
		omega = 0.

	elif gaschoice == ['d15N2']:
		
		gam_x = 1.275*0.9912227 # not sure of the origin here... Christo's model?
		#gam_x = 
		M = 1.E-3 + M_air
		decay = 0.
		omega = 0.0147/1000

	elif gaschoice == 'SF6':
		gam_x = 0.554
		M = 146.06e-3
		decay = 0.
		omega = 0.
		
	elif gaschoice == 'C14':
		gam_x = 0.991368
		M = 46.01e-3
		decay = 1./8267.
		omega = 0.
		
	elif gaschoice == 'C13':
		gam_x = 0.9955648
		M = 45.01e-3
		decay = 0.
		omega = 0.
		
	elif gaschoice == 'CFC11':
		gam_x = 0.525
		M = 137.37e-3
		decay = 0.

	elif gaschoice == 'CFC12':
		gam_x = 0.596
		M = 120.91e-3
		decay = 0.
		omega = 0.

	elif gaschoice == 'C13_CFC12':
		gam_x = 0.59552
		M = 121.91e-3
		decay = 0.
		omega = 0.

	elif gaschoice == 'CC14':
		gam_x = 0.470
		M = 153.82e-3
		decay = 0.
		omega = 0.

	elif gaschoice == 'CFC113':
		gam_x = 0.453
		M = 187.38e-3
		decay = 0.
		omega = 0.

	elif gaschoice == 'CFC115':
		gam_x = 0.532
		M = 154.47e-3
		decay = 0.
		omega = 0.

	elif gaschoice == 'R134a':
		gam_x = 0.630
		M = 102.03e-3
		decay = 0.
		omega = 0.

	elif gaschoice == 'CH3CCl3':
		gam_x = 0.485
		M = 133.40e-3
		decay = 0.
		omega = 0.

	elif gaschoice == 'HCFC22':
		gam_x = 0.710
		M = 86.47e-3
		decay = 0.
		omega = 0.

	elif gaschoice == 'C13_CH4':
		gam_x = 1.340806
		M = 17.04e-3
		decay = 0.
		omega = 0.
		
	elif gaschoice == 'd40Ar':
		gam_x = 1.21
		M = 4.e-3 + M_air
		decay = 0.
		omega = 0.0985/1000.

	elif gaschoice == 'FOG':
		gam_x = 1.0
		M = 44e-3
		decay = 1./100.
		omega = 0.
		
		
	
				
	
	### Load gas history. The file must be located in the correct folder, and have the correct naming convention.
	### If you want to compare to measured samples, make sure that measurements is on.    
	
#     if loadgas:
#         gas_string=gaschoice+'_history_'+hemisphere+'.txt'
#         meas_string=gaschoice+'_samples_'+sitechoice+'.txt'
	
#         conc1=loadtxt(os.path.join(DataPath,gas_string),skiprows=2) #load data: atmospheric CO2 history.
	
# #     if measurements=='on':
# #         firn_meas=loadtxt(os.path.join(DataPath,meas_string),skiprows=2)
# #     else:
# #         firn_meas='None'
	
#     else:
#         conc1=-9999
			
	deltaM = (M-M_air) #delta molecular mass from CO2.
	#gam_x = D_gas #* D_ref_CO2
	d_0=D_ref_CO2
				   
	return gam_x, M, deltaM, d_0, omega
	### D_x is the free-air diffusivity relative to CO2. 
		




