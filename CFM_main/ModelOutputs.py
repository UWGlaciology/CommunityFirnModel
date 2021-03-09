#!/usr/bin/env python
'''
Code for isotope diffusion.
'''
import numpy as np 
import json
import scipy.interpolate as interpolate
from constants import *
import os
import sys

class ModelOutputs:
    '''
    Class to handle making the model output files
    '''
    def __init__(self, config, MOd, TWlen, init_time, Glen):
        '''
        Initialize the model output class
        Main variable is Mout_dict, which is a dictionary that contains all
        of the variables that will get written to file.
        '''


        self.c = config
        self.Mout_dict = {}      

        if 'output_bits' not in self.c:
            self.c['output_bits']='float32'
        if 'grid_outputs' not in self.c:
            self.c['grid_outputs'] = False
        self.MOgrid = self.c['grid_outputs']

        if self.MOgrid:
            self.grid_out = np.arange(MOd['z'][0], MOd['z'][-1], self.c['grid_output_res'])

        self.output_list = list(MOd.keys())

        for varname in self.output_list:
            if varname == 'Dcon':
                intkind = 'nearest'
            else:
                intkind = 'linear'
            
            if varname == 'DIP':
                self.Mout_dict[varname] = np.zeros((TWlen+1,8), dtype = self.c['output_bits'])
                self.Mout_dict[varname][0,:]  = np.append(init_time, MOd[varname])
            elif varname == 'BCO':
                self.Mout_dict[varname] = np.zeros((TWlen+1,10), dtype = self.c['output_bits'])
                self.Mout_dict[varname][0,:]  = np.append(init_time, MOd[varname])
            elif varname == 'climate':
                self.Mout_dict[varname] = np.zeros((TWlen+1,3), dtype = self.c['output_bits'])
                self.Mout_dict[varname][0,:]  = np.append(init_time, MOd[varname])
            elif varname == 'runoff':
                self.Mout_dict[varname] = np.zeros((TWlen+1,2), dtype = self.c['output_bits'])
                self.Mout_dict[varname][0,:]  = np.append(init_time, MOd[varname])

            else:
                if self.MOgrid: #gridding outputs
                    if varname == 'z':
                        self.Mout_dict[varname] = np.zeros((len(self.grid_out)+1),dtype=self.c['output_bits'])
                        self.Mout_dict[varname] = np.append(init_time,self.grid_out)
                    else:
                        self.Mout_dict[varname] = np.zeros((TWlen+1,len(self.grid_out)+1),dtype=self.c['output_bits'])

                    if varname == 'LWC':
                        self.Mout_dict[varname][0,:] = np.append(init_time,RGfun(MOd['z'], MOd[varname], self.grid_out))
                    elif varname != 'z':
                        Ifun = interpolate.interp1d(MOd['z'], MOd[varname], kind = intkind, fill_value='extrapolate')         
                        self.Mout_dict[varname][0,:] = np.append(init_time,Ifun(self.grid_out))
                
                else: #not gridding outputs
                    self.Mout_dict[varname]       = np.zeros((TWlen+1,Glen+1),dtype=self.c['output_bits'])
                    self.Mout_dict[varname][0,:]  = np.append(init_time, MOd[varname])



    def updateMO(self, MOd, mtime, Wtracker):
        '''
        Function to update the output matrices in Mout_dict
        '''

        for varname in self.output_list:
            if varname == 'Dcon':
                intkind = 'nearest'
            else:
                intkind = 'linear'

            if self.MOgrid:
                if varname == 'LWC':
                    self.Mout_dict[varname][0,:] = np.append(mtime,RGfun(MOd[z], MOd[varname], self.grid_out))
                elif ((varname == 'BCO') or (varname == 'DIP') or (varname == 'climate') or (varname == 'runoff')):
                    self.Mout_dict[varname][Wtracker,:] = np.append(mtime,MOd[varname])
                elif varname == 'z':
                    continue
                else:
                    # try:
                        # Ifun = interpolate.interp1d(MOd['z'], MOd[varname], kind = intkind, fill_value='extrapolate')
                    Ifun = interpolate.interp1d(MOd['z'], MOd[varname], kind = intkind,bounds_error=False,fill_value=np.nan)           
                    self.Mout_dict[varname][Wtracker,:] = np.append(mtime,Ifun(self.grid_out))
                    # except:
                        # print(varname)
                        # print(mtime)
                        # sys.exit()
            else:
                self.Mout_dict[varname][Wtracker,:] = np.append(mtime,MOd[varname])
                

    def RGfun(self, z, var, grid):
        '''
        Function to regrid the variables that can not be linearly interpolated
        e.g. LWC needs to conserve mass. 
        '''

        varC = np.cumsum(var)
        newVar = np.interp(grid, z, var)
        return np.diff(newVar,append = newVar[-1])



        # old stuff below


        # if 'output_bits' not in self.c:
        #     self.c['output_bits']='float32'

        # if 'grid_outputs' not in self.c:
        #     self.c['grid_outputs'] = False

        # if ((not self.MELT) and ('LWC' in self.output_list)):
        #     self.output_list.remove('LWC')
        # if ((not self.c['FirnAir']) and ('gasses' in self.output_list)):
        #     self.output_list.remove('gasses')
        # if ((not self.c['isoDiff']) and ('isotopes' in self.output_list)):
        #     self.output_list.remove('isotopes')
        # if ((self.c['grid_outputs']) and ('Dcon' in self.output_list)):
        #     self.output_list.remove('Dcon')
        # if ((self.c['grid_outputs']) and ('compaction' in self.output_list)):
        #     self.output_list.remove('compaction')
        
        # if self.c['grid_outputs']:
        #     self.grid_out = np.arange(self.z[0],self.z[-1],self.c['grid_output_res'])

        # if 'density' in self.output_list:
        #     if self.c['grid_outputs']:
        #         self.rho_out            = np.zeros((TWlen+1,len(self.grid_out)+1),dtype=self.c['output_bits'])
        #         self.rho_out[0,:]       = np.append(init_time, np.interp(self.grid_out, self.z, self.rho))
        #     else:
        #         self.rho_out            = np.zeros((TWlen+1,len(self.dz)+1),dtype=self.c['output_bits'])
        #         self.rho_out[0,:]       = np.append(init_time, self.rho)
        
        # if 'temperature' in self.output_list:
        #     if self.c['grid_outputs']:
        #         self.Tz_out             = np.zeros((TWlen+1,len(self.grid_out)+1),dtype=self.c['output_bits'])
        #         self.Tz_out[0,:]        = np.append(init_time, np.interp(self.grid_out, self.z, self.Tz))
        #     else:
        #         self.Tz_out             = np.zeros((TWlen+1,len(self.dz)+1),dtype=self.c['output_bits'])
        #         self.Tz_out[0,:]        = np.append(init_time, self.Tz)
        
        # if 'age' in self.output_list:
        #     if self.c['grid_outputs']:
        #         self.age_out            = np.zeros((TWlen+1,len(self.grid_out)+1),dtype=self.c['output_bits'])
        #         self.age_out[0,:]       = np.append(init_time, np.interp(self.grid_out, self.z, self.age)/S_PER_YEAR)
        #     else:
        #         self.age_out            = np.zeros((TWlen+1,len(self.dz)+1),dtype=self.c['output_bits'])
        #         self.age_out[0,:]       = np.append(init_time, self.age/S_PER_YEAR)
        
        # if 'depth' in self.output_list:
        #     if self.c['grid_outputs']:
        #         self.z_out              = np.zeros((TWlen+1,len(self.grid_out)+1),dtype=self.c['output_bits'])
        #         self.z_out[0,:]         = np.append(init_time, self.grid_out)
        #     else:
        #         self.z_out              = np.zeros((TWlen+1,len(self.dz)+1),dtype=self.c['output_bits'])
        #         self.z_out[0,:]         = np.append(init_time, self.z)
        
        # if 'dcon' in self.output_list:
        #         self.D_out              = np.zeros((TWlen+1,len(self.dz)+1),dtype=self.c['output_bits'])
        #         self.D_out[0,:]         = np.append(init_time, self.Dcon)
        
        # if 'bdot_mean' in self.output_list:
        #     if self.c['grid_outputs']:
        #         self.bdot_out           = np.zeros((TWlen+1,len(self.grid_out)+1),dtype=self.c['output_bits'])
        #         self.bdot_out[0,:]      = np.append(init_time, np.interp(self.grid_out, self.z, self.bdot_mean))
        #     else:
        #         self.bdot_out           = np.zeros((TWlen+1,len(self.dz)+1),dtype=self.c['output_bits'])
        #         self.bdot_out[0,:]      = np.append(init_time, self.bdot_mean)
        
        # if 'climate' in self.output_list:
        #     self.Clim_out               = np.zeros((TWlen+1,3),dtype=self.c['output_bits'])
        #     self.Clim_out[0,:]          = np.append(init_time, [self.bdot[0], self.Ts[0]])  # not sure if bdot or bdotSec
        
        # if 'compaction' in self.output_list:
        #     self.comp_out               = np.zeros((TWlen+1,self.compboxes+1),dtype=self.c['output_bits'])
        #     self.comp_out[0,:]          = np.append(init_time, np.zeros(self.compboxes))
        
        # if 'LWC' in self.output_list:
        #     if self.c['grid_outputs']:
        #         self.LWC_out            = np.zeros((TWlen+1,len(self.grid_out)+1),dtype=self.c['output_bits'])
        #         self.LWC_out[0,:]       = np.append(init_time, np.interp(self.grid_out, self.z, self.LWC))
        #     else:
        #         self.LWC_out            = np.zeros((TWlen+1,len(self.dz)+1),dtype=self.c['output_bits'])
        #         self.LWC_out[0,:]       = np.append(init_time, self.LWC)

        # if 'PLWC_mem' in self.output_list:
        #     if self.c['grid_outputs']:
        #         self.PLWC_mem_out       = np.zeros((TWlen+1,len(self.grid_out)+1),dtype=self.c['output_bits'])
        #         self.PLWC_mem_out[0,:]  = np.append(init_time, np.interp(self.grid_out, self.z, self.PLWC_mem))                
        #     else:
        #         self.PLWC_mem_out       = np.zeros((TWlen+1,len(self.dz)+1),dtype=self.c['output_bits']) #VV
        #         self.PLWC_mem_out[0,:]  = np.append(init_time, self.PLWC_mem) #VV

        # if 'viscosity' in self.output_list:
        #     self.viscosity          = np.zeros(self.gridLen)
        #     if self.c['grid_outputs']:
        #         self.viscosity_out      = np.zeros((TWlen+1,len(self.grid_out)+1),dtype=self.c['output_bits'])
        #         self.viscosity_out[0,:] = np.append(init_time, np.interp(self.grid_out, self.z, self.viscosity))
        #     else:
        #         self.viscosity_out      = np.zeros((TWlen+1,len(self.dz)+1),dtype=self.c['output_bits'])
        #         self.viscosity_out[0,:] = np.append(init_time, self.viscosity)




