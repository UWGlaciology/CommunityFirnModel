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
        self.output_list = list(MOd.keys())

        if self.MOgrid:
            self.grid_res = self.c['grid_output_res']
            if 'grid_output_max_depth' not in self.c:
                print('grid_output_max_depth not in .json. Defaulting to 200 m.')
                self.c['grid_output_max_depth'] = 200.0
            self.max_grid_len = int(np.ceil(self.c['grid_output_max_depth'] / self.grid_res))

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
                self.Mout_dict[varname] = np.zeros((TWlen+1,6), dtype = self.c['output_bits'])
                self.Mout_dict[varname][0,:]  = np.append(init_time, MOd[varname])
            elif varname == 'refreeze':
                self.Mout_dict[varname] = np.zeros((TWlen+1,2), dtype = self.c['output_bits'])
                self.Mout_dict[varname][0,:]  = np.append(init_time, MOd[varname])
            elif varname == 'runoff':
                self.Mout_dict[varname] = np.zeros((TWlen+1,2), dtype = self.c['output_bits'])
                self.Mout_dict[varname][0,:]  = np.append(init_time, MOd[varname])
            elif varname == 'meltvol':
                self.Mout_dict[varname] = np.zeros((TWlen+1,2), dtype = self.c['output_bits'])
                self.Mout_dict[varname][0,:]  = np.append(init_time, MOd[varname])

            else:
                if self.MOgrid:  # gridding outputs - grid rebuilt each timestep, see _build_gridded_row
                    self.Mout_dict[varname] = np.zeros((TWlen+1,self.max_grid_len+1),dtype=self.c['output_bits'])
                    self.Mout_dict[varname][0,:] = np.append(init_time, self._build_gridded_row(MOd,varname,intkind))
                else:  # not gridding outputs
                    self.Mout_dict[varname]       = np.zeros((TWlen+1,Glen+1),dtype=self.c['output_bits'])
                    self.Mout_dict[varname][0,:]  = np.append(init_time, MOd[varname])

    ### end __init__    

    def _build_gridded_row(self, MOd, varname, intkind):
            '''
            Build a single padded row (length self.max_grid_len) for a gridded output,
            using the CURRENT domain depth (MOd['z']) to define the grid at this
            timestep. Grid extent/spacing depends on actual depth at write time,
            not a value fixed at __init__ - fixes restart-consistency issue
            (see notes, 2026-09-09).
            '''
            n_pts = int(np.floor((MOd['z'][-1] - MOd['z'][0]) / self.grid_res)) + 1
            grid_out_i = MOd['z'][0] + np.arange(n_pts) * self.grid_res
            n = len(grid_out_i)
            if n > self.max_grid_len:
                # print(f'Warning: domain depth exceeds grid_output_max_depth; truncating {varname}.')
                grid_out_i = grid_out_i[:self.max_grid_len]
                n = self.max_grid_len

            row = np.full(self.max_grid_len, np.nan, dtype=self.c['output_bits'])

            if varname == 'z':
                row[:n] = grid_out_i
            elif varname == 'LWC':
                row[:n] = self.RGfun(MOd['z'], MOd[varname], grid_out_i)
            else:
                Ifun = interpolate.interp1d(MOd['z'], MOd[varname], kind=intkind, fill_value='extrapolate')
                row[:n] = Ifun(grid_out_i)

            return row
        ### end _build_gridded_row ###
        #####################

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
                    if varname in ('BCO','DIP','climate','runoff','refreeze','meltvol'):
                        self.Mout_dict[varname][Wtracker,:] = np.append(mtime,MOd[varname])
                    else:
                        self.Mout_dict[varname][Wtracker,:] = np.append(mtime, self._build_gridded_row(MOd,varname,intkind))
                else:
                    self.Mout_dict[varname][Wtracker,:] = np.append(mtime,MOd[varname])
        ### end updateMO ###
        #####################
                    

    def RGfun(self, z, var, grid):
        '''
        Function to regrid the variables that can not be linearly interpolated
        e.g. LWC needs to conserve mass. 
        '''

        varC = np.cumsum(var)
        newVar = np.interp(grid, z, varC)
        return np.diff(newVar,append = newVar[-1])



