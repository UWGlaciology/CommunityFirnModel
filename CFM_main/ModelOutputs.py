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

### Outputs that are an amount per layer (regridded so that their total is conserved)
EXTENSIVE_OUTPUTS = ('LWC', 'refrozen', 'compaction', 'PLWC_mem')

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
            if 'grid_output_max_depth' not in self.c: # the model sets it from the initial domain; this is for other callers
                self.c['grid_output_max_depth'] = np.ceil(MOd['z'][-1])
                print(f"grid_output_max_depth not in .json. Defaulting to the depth of the model domain, {self.c['grid_output_max_depth']} m.")
            self._truncation_reported = False
            # grid points from the surface down to grid_output_max_depth, both included
            self.max_grid_len = int(np.floor(self.c['grid_output_max_depth'] / self.grid_res + 1e-9)) + 1

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
                if not self._truncation_reported:
                    print(f"Warning: the model domain ({MOd['z'][-1]:.1f} m) is deeper than grid_output_max_depth "
                          f"({self.c['grid_output_max_depth']} m); the gridded outputs are cut off at that depth.")
                    self._truncation_reported = True
                grid_out_i = grid_out_i[:self.max_grid_len]
                n = self.max_grid_len

            row = np.full(self.max_grid_len, np.nan, dtype=self.c['output_bits'])

            if varname == 'z':
                row[:n] = grid_out_i
            elif varname in EXTENSIVE_OUTPUTS:  # amounts per layer: regrid conserving the total
                row[:n] = self.RGfun(MOd['z'], MOd[varname], grid_out_i, to_bottom=(n_pts <= self.max_grid_len))
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
                    

    def RGfun(self, z, var, grid, to_bottom=True):
        '''
        Regrid an amount per layer (not a value that can be interpolated), e.g. LWC
        [m w.e.] or the compaction [m] of each layer, to the output grid, keeping the
        total: the sum of the grid cells between two depths is the amount in that
        depth range, and the sum of all cells is the sum over the layers.

        The cumulative amount is interpolated linearly between the layer edges, so a
        layer's amount is spread evenly over its own thickness (nodes with depth z are
        the tops of the layers; the last layer is taken to be as thick as the one above).

        :param z: depth of the top of each layer [m]
        :param var: amount in each layer
        :param grid: depth of the top of each output cell [m], evenly spaced
        :param to_bottom: if True, the last cell also gets whatever lies below it (the
            output grid ends inside the last layer); if False the grid was truncated, and
            the last cell ends one cell spacing below its top
        :return: amount in each output cell
        '''
        edges = np.append(z, z[-1] + (z[-1] - z[-2]))   # layer tops and the bottom of the last layer
        cumulative = np.concatenate(([0.], np.cumsum(var)))
        if len(grid) > 1:
            spacing = grid[1] - grid[0]
        else:
            spacing = edges[-1] - grid[0]
        cell_edges = np.append(grid, grid[-1] + spacing)
        newC = np.interp(cell_edges, edges, cumulative)
        if to_bottom:
            newC[-1] = cumulative[-1]
        return np.diff(newC)



