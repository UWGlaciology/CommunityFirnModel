#!/usr/bin/env python

'''
code to deal with modeling firn-air diffusion
'''


import numpy as np 
from solver import solver
from solver import transient_solve_TR
# from Gasses import gasses 
# from Sites import sites 
from reader import read_input
import json
import scipy.interpolate as interpolate
from constants import *
import os

class FirnAir:
    def __init__(self,air_config,Gs,z,modeltime,Tz, rho, dz, gaschoice, bdot):
        '''
        Initialize Firn Air class
        '''

        self.cg = air_config
        self.Gs = Gs

        self.Gz             = self.Gs[0]*np.ones_like(z)
        self.Tz             = Tz
        self.z              = z
        self.rho            = rho
        self.bdot           = bdot
        self.gaschoice      = gaschoice

        # self.gam_x, self.M, self.deltaM, self.d_0, self.omega = gasses(self.cg['gaschoice'], Tz[0], P_0, M_AIR)
        self.gam_x, self.M, self.deltaM, self.d_0, self.omega = gasses(gaschoice, Tz[0], P_0, M_AIR)

        self.czd            = self.cg['ConvectiveZoneDepth']
        self.p_a            = 7.0e4
        # self.air_pressure     = self.p_a * np.ones_like(z)
        self.air_pressure   = self.p_a * np.exp(M_AIR*GRAVITY*z/(R*self.Tz))
        self.air_pressure_base = np.copy(self.air_pressure)
        # self.air_pressure_0 = np.copy(self.air_pressure)

        self.rho_co, self.por_co, self.por_tot, self.por_cl, self.por_op, self.bcoRho, self.LIDRho = self.porosity()
        
        self.air_volume     = self.por_op * dz

        self.gas_age = np.zeros_like(self.z)

    def diffusivity(self):
        '''
        D_0 is CO2 in free air.
        gam_x is the diffusivity relative to CO2
        D_x=D_0*gam_x is the free-air (or any) diffusity for that species
        diffu_full has units m^2/s
        '''

        ## Constants
        d_eddy_sc   = self.d_0 #Eddy diffusivity in the convective zone
        d_ind       = np.min(np.where(self.z>self.z_co)) #indices of nodes past close-off depth
        if self.cg['runtype']=='steady':
            Tz_d = self.cg['steady_T'] * np.ones_like(self.Tz)
        else:
            Tz_d = self.Tz
        
        ##### Parameterizations for diffusivity #####
        if self.cg['Diffu_param'] == "Severinghaus": # Use Severinghaus relationship from Cuffey and Paterson
            diffu_full      = self.gam_x * self.d_0 * ((P_0/self.p_a) * (Tz_d/273.15)**1.85 * (2.00 * (1-(self.rho/RHO_I))-0.167))  

        elif self.cg['Diffu_param'] == "Schwander": # Use Schwander 1988, Eq. 2 Diffusivity (does not work very well) use 4e2 for d_0
            k_sch           = P_0 / self.p_a * (Tz_d / 253.16)**1.85 # Constant given in Schwander
            diffu_full      = self.gam_x * k_sch * (23.7 * self.por_tot - 2.84) / (1000**2) # 1/1000**2 is unit conversion.                 
        
        elif self.cg['Diffu_param'] == "Freitag": # Use Freitag, 2002, Eq 15 Diffusivity use 9e2 for d_0
            diffu_full      = 1.0 * self.gam_x * self.d_0 * self.por_op ** 2.1  

        elif self.cg['Diffu_param'] == "Adolph": # Use Adolph and Albert, 2014
            perm            = 10**-7.29 * self.por_op**3.71 
            diffu_full      = self.d_0 * self.gam_x * (perm/(10**-6.973))**(1/2.79)         
        
        elif self.cg['Diffu_param']=="Witrant":    ### Use Witrant, 2012
            diffu_full      = self.d_0 * (2.5 * self.por_op - 0.31) * (Tz_d / 273.15)**(1.8) * P_0 / self.p_a
            if self.bdot_t<0.1:
                # d_thresh = 3.16e-8
                d_thresh = 1.
            else:
                # d_thresh = 3.16e-6
                d_thresh = 100.

            thresh_ind = np.where(diffu_full*S_PER_YEAR<=d_thresh)[0]
            z_thresh = self.z[thresh_ind[0]]

            diffu_t = (d_thresh - (1e-2)) / (1 + np.exp(50*(self.z - (self.z_co + z_thresh)/2)/self.z_co)) + (1e-2)

            diffu_full[thresh_ind] = diffu_t[thresh_ind]/S_PER_YEAR

            diffu_full = self.gam_x * diffu_full
        
        elif self.cg['Diffu_param']=="Battle":
            diffu_full      = (self.gam_x * 1.16 * (23.7 * self.por_op - 2.84)) / (24*3600)
        
        elif self.cg['Diffu_param']=="Christo":
            pp              = "/Users/maxstev/Documents/Grad_School/Research/FIRN/CFM/CommunityFirnModel/gasmodel/DataImport"
            diffu_data      = np.loadtxt(os.path.join(pp,'c_diffu_NEEM.txt'))
            h               = diffu_data[:,0]
            diffu_full_data = self.gam_x*self.d_0*diffu_data[:,1]

            diffu_full      = np.interp(self.z,h,diffu_full_data)
        
        
                        
        ###### Add in high diffusivity in convective zone and low diffusivity below LID     
        ### Convective zone###
        d_eddy          = np.zeros(np.size(diffu_full))
        ind             = np.nonzero(self.z<self.czd)
        d_eddy_surf     = 2.426405E-5 #Kawamura, 2006
        # H_scale       = czd
        d_eddy_up       = d_eddy_surf * np.exp(-1*self.z/self.czd)
        
        #### Lock-in zone physics ###
        if self.cg['lockin']:
            LI_ind = np.where(diffu_full<=1e-9)[0][0]

            ind             = np.flatnonzero(self.z>self.LID)
            ind2            = np.flatnonzero(self.z<self.z_co)
            ind3            = np.intersect1d(ind,ind2)
            ind4 = np.where(self.z<=self.LID)[0]
            # ind5 = np.where(self.z>=self.LID)[0][0]

            # diffu_full[ind4] = diffu_full[ind4]-diffu_full[ind4[-1]] #re-scale diffusivity so it becomes zero at LID

            d_eddy[ind3]    = diffu_full[ind3]      # set eddy diffusivity in LID equal to diffusivity at LID
            diffu_full[ind] = 1e-20                 # set molecular diffusivity equal to zero for "main" diffusivity after LID - eddy diffusivity term drives diffusion below
            d_eddy          = d_eddy + d_eddy_up #make eddy diffusivity vector have convective and lock-in zone values
            
        else:
            # diffu_full[np.flatnonzero(self.z>self.z_co)] = 1.0e-40
            d_eddy          = d_eddy_up #eddy diffusivity includes only convective zone
         
        # ij = np.where(self.z>self.LID)[0][0]
        # diffu_full 

        # diffu_full[self.z>self.z_co]=1e-40
        diffu_full[self.rho>self.bcoRho]=1e-40

        diffu_full[diffu_full<=0] = 1.e-40

        i4=np.where(diffu_full==1.e-40)[0][0]
        diffu_full[i4:i4+3]=1e-40 # a single layer of low diffusivity is numerically unstable, so ensure that low-diffusivity is at least 3 layers thick.

        diffu=diffu_full
        
        return diffu , d_eddy #, dd

    def porosity(self): #,rho,T
        '''
        Calculate the total, open, and closed porosity.
        co is close-off, cl is closed
        '''
        
        indT=np.where(self.z>20)[0][0]
        if self.cg['runtype']=='steady':
            self.bcoRho         = 1/( 1/(RHO_I) + self.cg['steady_T'] * 6.95E-7 - 4.3e-5) # Martinerie mean close-off density
            # self.bcoRho = 825.3 
        else:
            self.bcoRho         = 1/( 1/(RHO_I) + self.Tz[indT] * 6.95E-7 - 4.3e-5) # Martinerie mean close-off density 
            # self.bcoRho = 825.3
        
        ### Porosity, from Goujon et al., 2003, equations 9 and 10
        self.por_tot        = 1-self.rho/RHO_I # Total porosity
        self.por_co         = 1 - self.bcoRho/RHO_I # mean close-off porosity
        alpha               = 0.37 # constant determined in Goujon
        self.por_cl         = np.zeros_like(self.por_tot)

        self.por_cl[self.por_tot>0]     = alpha*self.por_tot[self.por_tot>0]*(self.por_tot[self.por_tot>0]/self.por_co)**(-7.6)

        ind                 = self.por_cl>self.por_tot # where (if anywhere) closed porosity exceeds total porosity
        self.por_cl[ind]    = self.por_tot[ind] # fix por_cl in that case so it cannot exceed
        self.por_op         = self.por_tot - self.por_cl # Open Porosity

        co_ind              = np.where(self.por_op<=1e-10)[0][0]
        self.rho_co         = self.rho[co_ind]
        self.LIDRho         = self.rho_co - 14 #LID depth (Blunier and Schwander, 2000)

        op_i = np.where(self.por_op<=0)[0][0]
        co_i = np.where(self.rho>=self.rho_co)[0][0]
        ct_i = np.where(self.por_cl>=self.por_tot)[0][0]
        li_i = np.where(self.rho>=self.LIDRho)[0][0]

        # print('first full close', self.z[ct_i])
        # print('open',self.z[op_i])
        # print('open',self.rho[op_i])
        # print('bco',self.z[co_i])
        # print('bco',self.rho[co_i])
        # print('porosity',self.por_tot[co_i])
        # print('LID z',self.z[li_i])
        # print('LID rho',self.rho[li_i])
        # print('liz width',self.z[co_i]-self.z[li_i])

        self.por_op[self.por_op<=0] = 1.0e-25
        
        return self.rho_co, self.por_co, self.por_tot, self.por_cl, self.por_op, self.bcoRho, self.LIDRho


    def firn_air_diffusion(self,AirParams,iii):

        '''
        Solve the diffusion equation.
        Calls solver.py
        '''


        for k,v in list(AirParams.items()):
            setattr(self,k,v)

        nz_P        = len(self.z)
        nz_fv       = nz_P - 2
        nt          = 1

        z_edges1 = self.z[0:-1] + np.diff(self.z) / 2
        z_edges = np.concatenate(([self.z[0]], z_edges1, [self.z[-1]]))
        z_P_vec     = self.z
        # if (iii>100 and iii<110):
        #     print('z_P',z_P_vec[0:5])
        #     print('len p', len(z_P_vec))
        #     print('edges',z_edges[0:5])
        #     print('len edges',len(z_edges))
        # input('enter')
        # phi_s     = self.Ts[iii]
        phi_s       = self.Gz[0]
        phi_0       = self.Gz

        # K_ice     = 9.828 * np.exp(-0.0057 * phi_0)
        # K_firn    = K_ice * (self.rho / 1000) ** (2 - 0.5 * (self.rho / 1000))

        self.rho_co, self.por_co, self.por_tot, self.por_cl, self.por_op, self.bcoRho, self.LIDRho = self.porosity()
        
        # self.air_pressure_old     = np.copy(self.air_pressure)
        porosity_old            = (RHO_I-self.rho_old)/RHO_I

        por_co1              = 1 - self.rho_co/RHO_I # Porosity at close-off
        alpha               = 0.37 # constant determined in Goujon
        por_cl1              = np.zeros_like(porosity_old)
        por_cl1[porosity_old>0] = alpha*porosity_old[porosity_old>0]*(porosity_old[porosity_old>0]/por_co1)**(-7.6)
        ind                 = por_cl1>porosity_old
        por_cl1[ind]         = porosity_old[ind]
        por_op_old          = porosity_old - por_cl1 # Open Porosity

        self.air_volume_old     = por_op_old * self.dz_old
        # self.air_volume_old   = np.copy(self.air_volume)

        self.air_volume         = self.por_op * self.dz
        volfrac = self.air_volume_old / self.air_volume
        # volfrac = np.concatenate(([volfrac[0]],volfrac))
        self.air_pressure       = (self.p_a*np.exp(M_AIR*GRAVITY*self.z/(R*self.Tz))) * volfrac - (self.p_a*np.exp(M_AIR*GRAVITY*self.z/(R*self.Tz))) # assume air pressure is atmos in entire column


        self.pressure_grad      = np.gradient(self.air_pressure,self.z) 
        self.z_co               = min(self.z[self.rho>=(self.bcoRho)]) #close-off depth; bcoRho is mean close off density
        # self.z_co               = min(self.z[self.rho>=(self.rho_co)]) #close-off depth; bcoRho is close off density
        # self.z_co               = self.z[np.where(self.rho>=self.rho_co)[0][0]] #close-off depth; bcoRho is close off density
        self.LID                = min(self.z[self.rho>=(self.LIDRho)]) #lock in depth; LIDRho is lock-in density
        # print(self.LID)
        self.bdot_t = self.bdot[iii]

        self.diffu, self.d_eddy = self.diffusivity()
        # ij = np.where(self.diffu<=1e-20)[0][0]
        # print('Diffu zero',self.z[ij])

        airdict = {
            'd_eddy':           self.d_eddy,
            'por_op':           self.por_op,
            'Tz':               self.Tz,
            'deltaM':           self.deltaM,
            'omega':            self.omega,
            'dz':               self.dz,
            'rho':              self.rho,
            'gravity':          self.cg['gravity'],
            'thermal':          self.cg['thermal'],
            'air_pressure':     self.air_pressure,
            'pressure_grad':    self.pressure_grad,
            'z':                self.z, 
            'dt':               self.dt,
            'z_co':             self.z_co,
            'p_a':              self.p_a,
            'por_tot':          self.por_tot,
            'por_cl':           self.por_cl,
            'w_firn':           self.w_firn,
            'advection_type':   self.cg['advection_type'],
            'iii': iii,
            'gaschoice': self.gaschoice
            }

        msk = np.where(self.z>self.z_co)[0][0]
        c_vol = np.ones_like(self.rho) # Just a filler here to make the diffusion work.
        
        self.Gz, w_p = transient_solve_TR(z_edges, z_P_vec, nt, self.dt, self.diffu, phi_0, nz_P, nz_fv, phi_s, self.rho, c_vol, airdict)
        self.Gz = np.concatenate(([self.Gs[iii]], self.Gz[:-1]))

        ind_LID = np.where(self.z>=self.LID)[0]
        self.gas_age[self.gas_age>0] = self.gas_age[self.gas_age>0] + self.dt/S_PER_YEAR
        self.gas_age = np.concatenate(([0], self.gas_age[:-1]))
        self.gas_age[ind_LID[0]] = 15.0
        ii2 = np.where((self.z>self.LID) & (self.gas_age==0))[0]
        self.gas_age[ii2] = 15

        # print('diffu',self.diffu[0:10])

        return self.Gz, self.diffu, w_p, self.gas_age

def gasses(gaschoice, T, p_a, M_air):
    
    '''
    Function to set up specifics for different gasses.
    '''

    #d_0 = 5.e2 # Free air diffusivity, CO2, m**2/yr Schwander, 1988 reports 7.24 mm**2/s =379 m**2/yr
    d_0 = 1.6e-5 # m^2/s :wikipedia value. changed 9/27/13  Schwander, 1988 reports 7.24 mm**2/s = 7.24e-6 m**2/yr

    D_ref_CO2 = 5.75E-10*T**1.81*(101325/p_a) #Christo Thesis, appendix A3
    print('gas choice is ', gaschoice)
    
    if gaschoice == 'CO2':
        gam_x = 1. #free-air diffusivity relative to CO2. Unitless (gamma in Buizert thesis, page 13).
        M = 44.01e-3 # molecular mass, kg/mol
        decay = 0.
        omega = 0.0
        
    elif gaschoice == 'CH4':
        gam_x = 1.367
        M = 16.04e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'd15N2':
        
        gam_x = 1.275*0.9912227 # not sure of the origin here... Christo's model?
        #gam_x = 
        M = 1.E-3 + M_air
        decay = 0.
        omega = 0.0147/1000 # Grachev and Severinghaus, 2003 

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
        omega = 0.00985/1000. * 4 # Grachev and Severinghaus, 2003

    elif gaschoice == 'FOG':
        gam_x = 1.0
        M = 44e-3
        decay = 1./100.
        omega = 0.
            
    deltaM = (M - M_air) #delta molecular mass from CO2.
    #gam_x = D_gas #* D_ref_CO2
    d_0=D_ref_CO2
                   
    return gam_x, M, deltaM, d_0, omega
    ### D_x is the free-air diffusivity relative to CO2. 
        




