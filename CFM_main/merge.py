#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Script that contains 3 functions. These are to be used if we want to proceed to merging of thin layers of the firn column.
I suggest we specify in json input if merge is true/false and the thickness threshold ('merge_min'):
"merging": true,
"merge_min": 5e-3
mergesurf(): for layer[0] and layer[1], to be used in time_evolve() of firn_density_nospin
mergenotsurf(): for layers[2:], to be used in time_evolve() of firn_density_nospin
mergeall(): for all layers, to be used at the end of firn_density_spin
CAUTION:
- not used for all variables (e.g. du_dx)
- nothing is done considering gas neither for isotopes
@author: verjans
'''

import numpy as np
from constants import *

def mergesurf(self,thickmin,iii):
    '''
    This function is to call during time_evolve function of firn_density_nospin.
    We merge the surface layer[0] with the layer[1] below as long as layer[1] remains under a certain thickness threshold.
    By applying condition on layer[1] instead of layer[0], we avoid merging all newly accumulated layers in the case we use a RCM forcing on a short time scale.
    Thickness threshold must be specified and consistent with the one of mergenotsurf().
    '''
    
    if ((self.dz[1] < thickmin) or (self.dz[0] < 1e-4)): #test
        self.rho[1] = (self.rho[1]*self.dz[1]+self.rho[0]*self.dz[0]) / (self.dz[1]+self.dz[0])
        self.Tz[1] = (self.Tz[1]*self.mass[1]+self.Tz[0]*self.mass[0]) / (self.mass[1]+self.mass[0])
        self.r2[1] = (self.r2[1]*self.mass[1]+self.r2[0]*self.mass[0]) / (self.mass[1]+self.mass[0])
        self.age[1] = self.age[1] # suggestion of Max 28Jun, important if we use bdot_mean
        ### Additive variables: take sum ###
        self.LWC[1] = (self.LWC[0]+self.LWC[1])
        self.PLWC_mem[1] = (self.PLWC_mem[0]+self.PLWC_mem[1])

        self.dz[1] += self.dz[0] # add thickness to underlying layer
        
        ### Remove the thin surface layer ###
        self.rho = np.delete(self.rho,0)
        self.Tz = np.delete(self.Tz,0)
        self.r2 = np.delete(self.r2,0)
        self.age = np.delete(self.age,0)
        self.dz = np.delete(self.dz,0)
        self.LWC = np.delete(self.LWC,0)
        self.PLWC_mem = np.delete(self.PLWC_mem,0)
        ## For Dcon, here we remove the layer that is merged but maybe we want to remove the layer that receives the merging (and keep most recent dcon)##
        self.Dcon = np.delete(self.Dcon,0)
        
        self.rho = np.append(self.rho, self.rho[-1])
        self.Tz = np.append(self.Tz, self.Tz[-1])
        self.r2 = np.append(self.r2, self.r2[-1])
        self.age = np.append(self.age, self.age[-1])
        self.dz = np.append(self.dz, self.dz[-1])
        self.LWC = np.append(self.LWC, 0.)
        self.PLWC_mem = np.append(self.PLWC_mem, 0.)
        self.Dcon = np.append(self.Dcon, self.Dcon[-1])
        
        ### Adjustment of variables ###
        self.z = self.dz.cumsum(axis=0)
        self.z = np.delete(np.append(0,self.z),-1)
        self.gridLen = np.size(self.z)
        self.dx = np.ones(self.gridLen)
        self.mass = self.rho*self.dz
        self.mass_sum = self.mass.cumsum(axis = 0)
        self.sigma = (self.mass + self.LWC * RHO_W_KGM) * self.dx * GRAVITY
        self.sigma = self.sigma.cumsum(axis = 0)
        self.bdot_mean = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] / (self.age[1:] * RHO_I / self.t[iii]))))*self.c['stpsPerYear']*S_PER_YEAR
        #Not sure recalculation of T_mean and T10m are necessary
        # self.T_mean         = np.mean(self.Tz[self.z<50])
        self.T10m           = self.T_mean
        # No change of self.compboxes as it keeps value at end of spin up during entire time evolve (nb of layers above 80m depth at end of spinup)

    return (self.dz,self.z,self.gridLen,self.dx,self.rho,self.age,self.LWC,self.PLWC_mem,self.mass,self.mass_sum,self.sigma,self.bdot_mean,\
                    self.Dcon,self.T_mean,self.T10m,self.r2)


def mergenotsurf(self,thickmin,iii):
    '''
    This function is to call during time_evolve function of firn_density_nospin.
    We merge all the layers below a thickness threshold except the layers of indices 0 and 1.
    This allows layers that became too thin due to compaction to be merged with the layer below.
    Minimum thickness threshold must be specified as thickmin
    We don't do this for surface layer because that would lead to any newly accumulated layer to be merged if RCM forcing is on a short time scale. The surface layer has its own function mergesurf().
    '''

    rmind = np.array([]) # list of indices that will have to be removed
    thinlayers = np.where(self.dz[0:(len(self.dz)-2)]<thickmin)[0] # Potential candidates for being merged,-2 because last layer cannot be merged with any underlying layer
    for index in thinlayers:
        # Specify in the if condition that we don't do this for layers[0] and [1]
        if (index>1) and (self.dz[index] < thickmin): # if a layer is too thin (now we take into account the possibility that a layer from thinlayers has maybe been merged with another one)
            ### Non-additive variables: take arithmetic mean ###
            self.rho[index+1] = (self.rho[index+1]*self.dz[index+1]+self.rho[index]*self.dz[index]) / (self.dz[index+1]+self.dz[index])
            ''' For Tz and r2: use weighted mean according to mass rather than dz!!'''
            self.Tz[index+1] = (self.Tz[index+1]*self.mass[index+1]+self.Tz[index]*self.mass[index]) / (self.mass[index+1]+self.mass[index])
            self.r2[index+1] = (self.r2[index+1]*self.mass[index+1]+self.r2[index]*self.mass[index]) / (self.mass[index+1]+self.mass[index])
            self.age[index+1] = self.age[index+1] # suggestion of Max 28Jun, important if we use bdot_mean
            ### Additive variables: take sum ###
            self.LWC[index+1] = (self.LWC[index]+self.LWC[index+1])
            self.PLWC_mem[index+1] = (self.PLWC_mem[index]+self.PLWC_mem[index+1])

            self.dz[index+1] += self.dz[index] # add thickness to underlying 
            rmind = np.append(rmind,index) # add index to the list to remove
    
    rmind = (rmind.astype(np.int32)).tolist()   
    ### Remove the thin layers ###
    self.rho = np.delete(self.rho,rmind)
    self.Tz = np.delete(self.Tz,rmind)
    self.r2 = np.delete(self.r2,rmind)
    self.age = np.delete(self.age,rmind)
    self.dz = np.delete(self.dz,rmind)
    self.LWC = np.delete(self.LWC,rmind)
    self.PLWC_mem = np.delete(self.PLWC_mem,rmind)
    ## For Dcon, here we remove the layer that is merged but maybe we want to remove the layer that receives the merging (and keep most recent dcon)##
    self.Dcon = np.delete(self.Dcon,rmind)
    # We removed len(rmind) layers
    
    self.rho = np.concatenate((self.rho, self.rho[-1]*np.ones(len(rmind))))
    self.Tz = np.concatenate((self.Tz, self.Tz[-1]*np.ones(len(rmind))))
    self.r2 = np.concatenate((self.r2, self.r2[-1]*np.ones(len(rmind))))
    self.age = np.concatenate((self.age, self.age[-1]*np.ones(len(rmind))))
    self.dz = np.concatenate((self.dz, self.dz[-1]*np.ones(len(rmind))))
    self.LWC = np.concatenate((self.LWC, np.zeros(len(rmind))))
    self.PLWC_mem = np.concatenate((self.PLWC_mem, np.zeros(len(rmind))))
    self.Dcon = np.concatenate((self.Dcon, self.Dcon[-1]*np.ones(len(rmind))))

    ### Adjustment of variables ###
    self.z = self.dz.cumsum(axis=0)
    self.z = np.delete(np.append(0,self.z),-1)
    self.gridLen = np.size(self.z)
    self.dx = np.ones(self.gridLen)
    self.mass = self.rho*self.dz
    self.mass_sum = self.mass.cumsum(axis = 0)
    self.sigma = (self.mass + self.LWC * RHO_W_KGM) * self.dx * GRAVITY
    self.sigma = self.sigma.cumsum(axis = 0)
    self.bdot_mean = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] / (self.age[1:] * RHO_I / self.t[iii]))))*self.c['stpsPerYear']*S_PER_YEAR
    #Not sure recalculation of T_mean and T10m are necessary
    # self.T_mean         = np.mean(self.Tz[self.z<50])
    self.T10m           = self.T_mean
    
    return (self.dz,self.z,self.gridLen,self.dx,self.rho,self.age,self.LWC,self.PLWC_mem,self.mass,self.mass_sum,self.sigma,self.bdot_mean,\
                    self.Dcon,self.T_mean,self.T10m,self.r2)
    

def mergeall(self,thickmin,iii):
    ''' 
    We spot the layers that are under a certain thickness threshold and we merge these with the underlying layer
    This has to be launched at the end of the spinup, we use other functions in firn_density_nospin.
    Here we change self.compboxes as we modify the number of layers above 80m depth.
    '''
    rmind = np.array([]) # list of indices that will have to be removed
    thinlayers = np.where(self.dz[0:(len(self.dz)-2)]<thickmin)[0] # Potential candidates for being merged,-2 because last layer cannot be merged with any underlying layer
    for index in thinlayers:
        if self.dz[index] < thickmin: # if a layer is too thin (now we take into account the possibility that a layer from thinlayers has maybe been merged with another one)
            ### Non-additive variables: take arithmetic mean ###
            self.rho[index+1] = (self.rho[index+1]*self.dz[index+1]+self.rho[index]*self.dz[index]) / (self.dz[index+1]+self.dz[index])
            self.Tz[index+1] = (self.Tz[index+1]*self.dz[index+1]+self.Tz[index]*self.dz[index]) / (self.dz[index+1]+self.dz[index])
            self.r2[index+1] = (self.r2[index+1]*self.dz[index+1]+self.r2[index]*self.dz[index]) / (self.dz[index+1]+self.dz[index])
            self.age[index+1] = self.age[index+1] # suggestion of Max 28Jun, important if we use bdot_mean
            ### Additive variables: take sum ###
            self.LWC[index+1] = (self.LWC[index]+self.LWC[index+1])

            self.dz[index+1] += self.dz[index] # add thickness to underlying layer
            rmind = np.append(rmind,index) # add index to the list to remove
    
    ### Remove the thin layers ###
    self.rho = np.delete(self.rho,rmind)
    self.Tz = np.delete(self.Tz,rmind)
    self.r2 = np.delete(self.r2,rmind)
    self.age = np.delete(self.age,rmind)
    self.dz = np.delete(self.dz,rmind)
    self.LWC = np.delete(self.LWC,rmind)

    ### Adjustment of variables ###
    self.z = self.dz.cumsum(axis=0)
    self.z = np.delete(np.append(0,self.z),-1)
    self.gridLen = np.size(self.z)
    self.dx = np.ones(self.gridLen)
    self.mass = self.rho*self.dz
    self.mass_sum = self.mass.cumsum(axis = 0)
    self.sigma = (self.mass + self.LWC * RHO_W_KGM) * self.dx * GRAVITY
    self.sigma = self.sigma.cumsum(axis = 0)
    self.bdot_mean = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] / (self.age[1:] * RHO_I / self.t[iii]))))*self.c['stpsPerYear']*S_PER_YEAR
    self.compboxes = len(self.z[self.z<80])
    #Not sure recalculation of T_mean and T10m are necessary
    # self.T_mean         = np.mean(self.Tz[self.z<50])
    self.T10m           = self.T_mean
    
    
    return (self.dz,self.z,self.gridLen,self.dx,self.rho,self.age,self.LWC,self.mass,self.mass_sum,self.sigma,self.bdot_mean,\
                    self.compboxes,self.T_mean,self.T10m,self.r2)
       



