#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Script for sublimation.
'''

from constants import *
import numpy as np


def sublim(self,iii):
    '''
    Sublimation of the surface layers, partially based on melt.py
    We don't do anything energy-wise (no modification of Tz)
    Layers are sublimated in turn, starting with the surface layer
    Liquid water is sublimated before the ice matrix
    '''
    
    
    lwc_initial = sum(self.LWC)
    
    sublim_volume_IE      = abs(self.bdotSec[iii]) * S_PER_YEAR #/ self.c['stpsPerYear'] # [m]
    sublim_volume_WE      = sublim_volume_IE * RHO_I_MGM # [m]
    sublim_mass           = sublim_volume_WE * 1000. # [kg]
    ind1a               = np.where((np.cumsum(self.mass)+np.cumsum(1000*self.LWC)) <= sublim_mass)[0] # indices of boxes that will be sublimated away
    num_boxes_sublim    = len(ind1a)+1 # number of boxes that sublimate away, include the box that is partially sublimated
    try:
        ind1                = np.where((np.cumsum(self.mass)+np.cumsum(1000*self.LWC)) > sublim_mass)[0][0] # index which will become the new surface
    except:
        print('error!')
        print(iii)
        print(self.modeltime[iii])
        print('sublim_mass',sublim_mass)
        print('rho',self.rho[0:20])
        print('tz',self.Tz[0:20])
        print('bdotsec',(self.bdotSec[iii] * S_PER_YEAR))
        sys.exit()
 
    # ps is the partial sublimation (the model volume that has a portion sublimated away)   
    if ind1 > 0: # if we sublimate the entire layers we retrieve the mass+lwc of these layers to the sublimated mass of the ps layer
        ps_sublim = sublim_mass - (np.cumsum(self.mass[ind1-1])+1000*np.cumsum(self.LWC[ind1-1])) 
    elif ind1 == 0: # if only the surface layer gets partially sublimated, all sublimation occurs in surface layer (logically)
        ps_sublim = sublim_mass
    ps_sublimlwc = min(ps_sublim,1000*self.LWC[ind1]) # first, we sublimate as much possible of the liquid water
    ps_sublimmass = ps_sublim - ps_sublimlwc # rest of sublimated mass will be ice mass
    ps_lwc = np.maximum(0,(self.LWC[ind1] - ps_sublimlwc/1000)) # new surface layer gets its lwc reduced but not below 0
    ps_plwc = np.maximum(self.PLWC_mem[ind1] - ps_sublimlwc/1000, 0.) # assumes plwc is sublimated before mlwc
    ps_mass = np.maximum(1.0e-6,(self.mass[ind1] - ps_sublimmass)) # finally ice is sublimated if there is still mass to sublimate (<-> if ps_sublim<1000*self.LWC[ind1]), avoid rounding errors to cause negative mass
    ps_dz = ps_mass / self.rho[ind1] # remaining thickness [m]

    ## Sublimated boxes are accomodated by just adding more (new) boxes at the bottom of the column
    ## Beware of this if you are not modeling to firn-ice transition depth.
    divider         = num_boxes_sublim #VV nb of sublimated boxes, including the partially sublimated
    self.rho        = np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_sublim))) #VV add at bottom of column as many layers as were sublimated away
    self.LWC        = np.concatenate((self.LWC[ind1:-1] , np.zeros(num_boxes_sublim))) # This is better but should be equivalent as last layer should be an ice layer -> with 0 LWC
    self.LWC[0]     = ps_lwc #VV LWC calculated for the partially sublimated layer
    
    self.PLWC_mem   = np.concatenate((self.PLWC_mem[ind1:-1] , np.zeros(num_boxes_sublim)))
    self.PLWC_mem[0]= ps_plwc
    # all the water that was in the PFdom of the sublimated layers is also for input
    
    self.age        = np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_sublim))) + self.dt[iii] # age of each layer increases of dt
    # self.dz                  = np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_sublim))) # this splits the last box into many.
    self.dz         = np.concatenate((self.dz[ind1:-1] , self.dz[-1]*np.ones(num_boxes_sublim))) # this adds new boxes at the bottom.
    self.dz[0]      = ps_dz #VV dz calculated for the partially sublimated layer
    self.Dcon       = np.concatenate((self.Dcon[ind1:-1] , self.Dcon[-1]*np.ones(num_boxes_sublim)))
    self.dzn        = np.concatenate((np.zeros(num_boxes_sublim), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
    self.dzn        = self.dzn[0:self.compboxes]
    self.Tz         = np.concatenate((self.Tz[ind1:-1] , self.Tz[-1]*np.ones(num_boxes_sublim)))
    self.Tz[0]      = np.copy(self.Ts[iii])
    if self.c['physGrain']:
        self.r2         = np.concatenate((self.r2[ind1:-1] , self.r2[-1]*np.ones(num_boxes_sublim)))
    self.bdot_mean  = np.concatenate((self.bdot_mean[ind1:-1] , self.bdot_mean[-1]*np.ones(num_boxes_sublim)))
    self.z          = self.dz.cumsum(axis = 0)
    self.z          = np.concatenate(([0] , self.z[:-1]))
    self.mass       = self.rho * self.dz
    ### VV changes 09/12/2020
    if self.doublegrid:
        sublgridtrack = np.append(self.gridtrack[ind1:-1],self.gridtrack[-1]*np.ones(num_boxes_sublim))
    elif self.doublegrid==False:
        sublgridtrack = np.zeros_like(self.dz)
    ###

    ## Grid should be ready

    self.totwatersublim += (lwc_initial-sum(self.LWC))

    if np.any(self.LWC<0.0):
        print('negative LWC after sublim')
        print('setting to zero and continuing')
    self.LWC[self.LWC<0]=0.0

    return self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, self.dzn, self.LWC, self.PLWC_mem, self.totwatersublim,sublgridtrack
    
    