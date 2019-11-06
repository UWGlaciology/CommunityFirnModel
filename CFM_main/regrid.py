#!/usr/bin/env python
from constants import *
import numpy as np

def regrid(self):
    '''
    Called in both firn_density_spin and firn_density_nospin

    There are 3 subgrids in the regrid module. Grid 1 is the high resolution grid near the surface. Grid 2 is the lower resolution grid at greater depths; a user-defined number of nodes (self.c['nodestocombine']; refer as NTC here) are combined occasionally (every NTC time steps) to make one new node within grid 2. Grid 3 is at the bottom and has split up one grid 2 node back into a high-resolution grid (1 node into NTC nodes), which can be removed at each time step to keep the model Lagrangian. 

    the variable gridtrack keeps track of which subgrid each node is in.
    '''

    ind10   = np.where(self.gridtrack==1)[0] # all of the nodes in subgrid 1. 
    ind1    = np.where(self.gridtrack==1)[0][-1*self.c['nodestocombine']:] # the last NTC nodes of subgrid 1; will be combined into 1 node within subgrid 2.
    ind1a   = ind1[0]
    ind1b   = ind1[-1]
    ind0    = ind1[0] - 1 # new last node of grid 1

    ### create the properties of the new subgrid 2 node
    g2dz    = np.array([np.sum(self.dz[ind1])])
    g2mass  = np.sum(self.mass[ind1])
    g2rho   = g2mass/g2dz
    g2Tz0   = np.sum(self.Tz[ind1]*self.mass[ind1])
    g2Tz    = np.array([g2Tz0 / g2mass]) # Use a weighted average for temperature (effectively the enthalpy)
    g2gt    = 2 #gridtrack
    g2age   = np.mean(self.age[ind1])
    # g2bm  = np.mean(self.bdot_mean[ind1])
    g2bm0   = np.sum(self.bdot_mean[ind1]*self.mass[ind1])
    g2bm    = np.array([g2bm0 / g2mass])
    g2lwc   = np.sum(self.LWC[ind1])

    ### split up the last node in grid 2 into NTC nodes. Each node retains the density, age, etc of the old subgrid 2 node. 
    g3dz    = self.dz[-1]/self.nodestocombine * np.ones(self.nodestocombine)
    g3rho   = self.rho[-1] * np.ones(self.nodestocombine)
    g3mass  = g3rho * g3dz
    g3gt    = 3 * np.ones(self.nodestocombine)
    g3Tz    = self.Tz[-1]* np.ones(self.nodestocombine)
    g3age   = self.age[-1]*np.ones(self.nodestocombine)
    g3bm    = self.bdot_mean[-1]*np.ones(self.nodestocombine)
    g3lwc   = self.LWC[-1]/self.nodestocombine * np.ones(self.nodestocombine)

    ### combine the new and old nodes into the full grid. 
    self.dz         = np.concatenate((self.dz[0:ind1a],g2dz,self.dz[ind1b+1:-1],g3dz))
    self.z          = self.dz.cumsum(axis=0)
    self.z          = np.concatenate(([0], self.z[:-1]))
    self.rho        = np.concatenate((self.rho[0:ind1a],g2rho,self.rho[ind1b+1:-1],g3rho))
    self.Tz         = np.concatenate((self.Tz[0:ind1a],g2Tz,self.Tz[ind1b+1:-1],g3Tz))
    self.mass       = np.concatenate((self.mass[0:ind1a],[g2mass],self.mass[ind1b+1:-1],g3mass))
    self.sigma      = self.mass * self.dx * GRAVITY
    self.sigma      = self.sigma.cumsum(axis = 0)
    self.mass_sum   = self.mass.cumsum(axis = 0)
    self.age        = np.concatenate((self.age[0:ind1a],[g2age],self.age[ind1b+1:-1],g3age))
    self.bdot_mean  = np.concatenate((self.bdot_mean[0:ind1a],g2bm,self.bdot_mean[ind1b+1:-1],g3bm))
    self.LWC        = np.concatenate((self.LWC[0:ind1a],[g2lwc],self.LWC[ind1b+1:-1],g3lwc))
    self.gridtrack  = np.concatenate((self.gridtrack[0:ind1a],[g2gt],self.gridtrack[ind1b+1:-1],g3gt))

    if self.c['physGrain']:
        #g2r2     = np.array([np.mean(self.r2)])
        g2r2      = np.mean(self.r2[ind1]) # VV added
        g3r2     = self.r2[-1]* np.ones(self.nodestocombine)
        self.r2 = np.concatenate((self.r2[0:ind1a],[g2r2],self.r2[ind1b+1:-1],g3r2)) 

    return self.dz, self.z, self.rho, self.Tz, self.mass, self.sigma, self.mass_sum, self.age, self.bdot_mean, self.LWC, self.gridtrack, self.r2

def init_regrid(self):
    '''
    Used in firn_density_spin for the initial regridding. 
    '''

    grid1b          = self.c['grid1bottom']
    self.nodestocombine = self.c['nodestocombine']
    ind1            = np.where(self.z<grid1b)[0]
    ind2            = np.where(self.z>=grid1b)[0]
    grid1z          = self.z[ind1]
    grid2z          = self.z[ind2[0]::self.nodestocombine]
    self.z          = np.concatenate((grid1z,grid2z))
    grid3z          = self.z[-1] + np.cumsum(self.dz[-1*self.nodestocombine:])
    self.z          = np.concatenate((self.z,grid3z))
    self.dz         = np.diff(self.z)
    self.dz         = np.append(self.dz, self.dz[-1])
    self.gridLen    = len(self.z)
    self.dx         = np.ones(self.gridLen)
    self.gridtrack  = 2 * np.ones(self.gridLen)
    self.gridtrack[ind1] = 1
    self.gridtrack[-1*self.nodestocombine:] = 3

    print('After regrid, grid length is', self.gridLen)

    return self.nodestocombine, self.z, self.dz, self.gridLen, self.dx, self.gridtrack
