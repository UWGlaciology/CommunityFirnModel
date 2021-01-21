#!/usr/bin/env python
'''
Script to change the grid to have different resolutions at different depths.
'''

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


############### VV changes 09/12/2020 ###############

def regrid22(self):
    '''
    Called in both firn_density_spin and firn_density_nospin
    5 grids:
        grid1  -> high resolution determined by accumulation events
        grid2  -> low resolution by merging the batch of lowest nodestocombine layers of grid 1
        grid22 -> very low resolution by merging the batch of lowest multnodestocombine layers of grid 2
        grid23 -> low resolution by splitting the lowest layer of grid22 in multnodestocombine thinner layers
        New layer of grid23 is formed only when their stock is empty
        grid3  -> high resolution by splitting the lowest layer of grid23 in nodestocombine layers

    gridtrack keeps track of which grid each layer is in
    '''


    n1      = self.c['nodestocombine'] # nodes to combine from grid1 to grid2 and to split from grid23 to grid3
    n2      = self.c['multnodestocombine'] # nodes to combine from grid2 to grid22 and to split from grid22 to grid23
    # if self.c['multnodestocombine'] is set to 0 -> process of grid22 is turned off and regrid works as old regrid
    inds1   = np.where(self.gridtrack==1)[0] # layers in grid1 
    inds2   = np.where(self.gridtrack==2)[0] # layers in grid2
    i1_2    = inds1[-1*n1:] # layers to transition from grid1 to grid2
    ind2a   = i1_2[0] # index of the future upper layer of grid2
    ind2b   = inds2[0] # index of old upper layer of grid2
    # Next grids
    inds22  = np.where(self.gridtrack==22)[0] # all nodes in grid22
    inds23  = np.where(self.gridtrack==23)[0] # all nodes in grid23

    ### properties of the new subgrid 2 node
    g2dz    = np.array([np.sum(self.dz[i1_2])]) # sum thickness
    g2mass  = np.sum(self.mass[i1_2]) # sum mass
    g2rho   = g2mass/g2dz
    g2Tz0   = np.sum(self.Tz[i1_2]*self.mass[i1_2])
    g2Tz    = np.array([g2Tz0 / g2mass]) # Use a weighted average for temperature (effectively the enthalpy)
    g2gt    = 2 #gridtrack
    g2age   = np.mean(self.age[i1_2]) # mean age
    #g2age   = np.sum(self.age[i1_2]*self.mass[i1_2])/g2mass #VV test weighted average for age -> change is imperceptible
    g2bdm   = np.mean(self.bdot_mean[i1_2]) #mean bdot_mean
    g2lwc   = np.sum(self.LWC[i1_2]) # sum for lwc
    g2r2    = np.mean(self.r2[i1_2]) # mean for r2
    
    if (len(inds23)==0 and n2>0): # No more layer in grid23 -> we have to split a layer from grid22
        ## First: merge n2 layers from grid2 ##
        i2_22 = inds2[-1*n2:] # layers to transition from grid2 to grid22
        ind22a   = i2_22[0] # index of the future upper layer of grid22
        ind22b   = inds22[0] # current upper node of grid22
        # Properties of the new grid22 layer
        g22dz    = np.array([np.sum(self.dz[i2_22])]) # sum thickness
        g22mass  = np.sum(self.mass[i2_22]) # sum mass
        g22rho   = g22mass/g22dz
        g22Tz0   = np.sum(self.Tz[i2_22]*self.mass[i2_22])
        g22Tz    = np.array([g22Tz0 / g22mass]) # Use a weighted average for temperature (the enthalpy)
        g22gt    = 22 #gridtrack
        g22age   = np.mean(self.age[i2_22])
        g22bdm   = np.mean(self.bdot_mean[i2_22])
        g22lwc   = np.sum(self.LWC[i2_22])
        g22r2    = np.mean(self.r2[i2_22])
        
        ## Second: split the last grid22 layer in n2 layers for grid23
        ind22c   = inds22[-1] # current lower layer of grid22 -> to be split (is also the last layer of the column)
        g23dz    = self.dz[ind22c]/n2 * np.ones(n2)
        g23rho   = self.rho[ind22c] * np.ones(n2)
        g23mass  = g23rho * g23dz
        g23gt    = 23 * np.ones(n2) #gridtrack values
        g23Tz    = self.Tz[ind22c]* np.ones(n2)
        g23age   = self.age[ind22c]*np.ones(n2)
        g23bdm   = self.bdot_mean[ind22c]*np.ones(n2)
        g23lwc   = self.LWC[ind22c]/n2 * np.ones(n2)
        g23r2    = self.r2[ind22c]*np.ones(n2)
        
        ## Third: split the last layer of the new grid23 in n1 layers for grid3
        g3dz    = g23dz[-1]/n1 * np.ones(n1)
        g3rho   = g23rho[-1] * np.ones(n1)
        g3mass  = g3rho * g3dz
        g3gt    = 3 * np.ones(n1)
        g3Tz    = g23Tz[-1]* np.ones(n1)
        g3age   = g23age[-1]*np.ones(n1)
        g3bdm   = g23bdm[-1]*np.ones(n1)
        g3lwc   = g23lwc[-1]/n1 * np.ones(n1)
        g3r2    = g23r2[-1]*np.ones(n1)
        
        ## Fourth: concatenate everything together
        self.dz         = np.concatenate((self.dz[0:ind2a],g2dz,self.dz[ind2b:ind22a],g22dz,self.dz[ind22b:ind22c],g23dz[0:-1],g3dz))
        self.z          = self.dz.cumsum(axis=0)
        self.z          = np.concatenate(([0], self.z[:-1]))
        self.rho        = np.concatenate((self.rho[0:ind2a],g2rho,self.rho[ind2b:ind22a],g22rho,self.rho[ind22b:ind22c],g23rho[0:-1],g3rho))
        self.Tz         = np.concatenate((self.Tz[0:ind2a],g2Tz,self.Tz[ind2b:ind22a],g22Tz,self.Tz[ind22b:ind22c],g23Tz[0:-1],g3Tz))
        self.mass       = np.concatenate((self.mass[0:ind2a],[g2mass],self.mass[ind2b:ind22a],[g22mass],self.mass[ind22b:ind22c],g23mass[0:-1],g3mass))
        self.mass_sum   = self.mass.cumsum(axis = 0)
        self.age        = np.concatenate((self.age[0:ind2a],[g2age],self.age[ind2b:ind22a],[g22age],self.age[ind22b:ind22c],g23age[0:-1],g3age))
        self.bdot_mean  = np.concatenate((self.bdot_mean[0:ind2a],[g2bdm],self.bdot_mean[ind2b:ind22a],[g22bdm],self.bdot_mean[ind22b:ind22c],g23bdm[0:-1],g3bdm))
        self.LWC        = np.concatenate((self.LWC[0:ind2a],[g2lwc],self.LWC[ind2b:ind22a],[g22lwc],self.LWC[ind22b:ind22c],g23lwc[0:-1],g3lwc))
        self.sigma      = (self.mass+self.LWC*RHO_W_KGM)*self.dx*GRAVITY
        self.sigma      = self.sigma.cumsum(axis = 0)
        self.gridtrack  = np.concatenate((self.gridtrack[0:ind2a],[g2gt],self.gridtrack[ind2b:ind22a],[g22gt],self.gridtrack[ind22b:ind22c],g23gt[0:-1],g3gt))
        self.r2         = np.concatenate((self.r2[0:ind2a],[g2r2],self.r2[ind2b:ind22a],[g22r2],self.r2[ind22b:ind22c],g23r2[0:-1],g3r2))
        
    if (len(inds23)>0 or n2==0): # Still some layers in grid23 -> no need to split a layer from grid22
        ## Split the last layer of grid23 (layer [-1]) in n1 layers for grid3
        g3dz    = self.dz[-1]/n1 * np.ones(n1)
        g3rho   = self.rho[-1] * np.ones(n1)
        g3mass  = g3rho * g3dz
        g3gt    = 3 * np.ones(n1)
        g3Tz    = self.Tz[-1]* np.ones(n1)
        g3age   = self.age[-1]*np.ones(n1)
        g3bdm   = self.bdot_mean[-1]*np.ones(n1)
        g3lwc   = self.LWC[-1]/n1 * np.ones(n1)
        g3r2    = self.r2[-1]*np.ones(n1)
        
        ## Concatenate everything together
        self.dz         = np.concatenate((self.dz[0:ind2a],g2dz,self.dz[ind2b:-1],g3dz))
        self.z          = self.dz.cumsum(axis=0)
        self.z          = np.concatenate(([0], self.z[:-1]))
        #print("self.z[-1]:",self.z[-1])
        self.rho        = np.concatenate((self.rho[0:ind2a],g2rho,self.rho[ind2b:-1],g3rho))
        self.Tz         = np.concatenate((self.Tz[0:ind2a],g2Tz,self.Tz[ind2b:-1],g3Tz))
        self.mass       = np.concatenate((self.mass[0:ind2a],[g2mass],self.mass[ind2b:-1],g3mass))
        self.mass_sum   = self.mass.cumsum(axis = 0)
        self.age        = np.concatenate((self.age[0:ind2a],[g2age],self.age[ind2b:-1],g3age))
        self.bdot_mean  = np.concatenate((self.bdot_mean[0:ind2a],[g2bdm],self.bdot_mean[ind2b:-1],g3bdm))
        self.LWC        = np.concatenate((self.LWC[0:ind2a],[g2lwc],self.LWC[ind2b:-1],g3lwc))
        self.sigma      = (self.mass+self.LWC*RHO_W_KGM)*self.dx*GRAVITY
        self.sigma      = self.sigma.cumsum(axis = 0)
        self.gridtrack  = np.concatenate((self.gridtrack[0:ind2a],[g2gt],self.gridtrack[ind2b:-1],g3gt))
        self.r2 = np.concatenate((self.r2[0:ind2a],[g2r2],self.r2[ind2b:-1],g3r2)) 

    #self.bdot_mean  = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] * self.t / (self.age[1:] * RHO_I))))*self.c['stpsPerYear']*S_PER_YEAR
    #print('sum(self.dz):',sum(self.dz))

    return self.dz, self.z, self.rho, self.Tz, self.mass, self.sigma, self.mass_sum, self.age, self.bdot_mean, self.LWC, self.gridtrack, self.r2



def regrid22_reciprocal(self):

    '''
    Reciprocal of regrid22: must be called if we accumulate too many grid3 nodes because of
    heavy melting of surface nodes (problematic in ablation area)
    -> merge k batches of n1 grid3 nodes into k grid23 nodes (k is maximum nb of batches of n1 grid3 nodes available)
    if nb of layers in grid2 is <k:
        -> calculate the nb of supplementary grid2 layers required
        -> calculate xx: the number of grid22 layers that must be split to provide the supplementary grid2 layers
        -> merge xx batches of n2 grid23 layers into xx grid22 layer
        -> divide xx grid22 layers into xx*n2 grid2 layers
    -> divide k grid2 layer into n1 grid1 layer
    
    5 grids:
        grid1  -> high resolution determined by accumulation events
        grid2  -> low resolution
        grid22 -> very low resolution
        grid23 -> low resolution
        grid3  -> high resolution

    gridtrack keeps track of which grid each layer is in
    '''

    n1      = self.c['nodestocombine'] # nodes to combine from grid3 to grid23 and to split from grid2 to grid1
    n2      = self.c['multnodestocombine'] # nodes to combine from grid23 to grid22 and to split from grid22 to grid2
    inds1   = np.where(self.gridtrack==1)[0] # all nodes in grid1
    inds2   = np.where(self.gridtrack==2)[0] # all nodes in grid2
    inds22  = np.where(self.gridtrack==22)[0] # all nodes in grid22
    inds23   = np.where(self.gridtrack==23)[0] # layers in grid23
    inds3   = np.where(self.gridtrack==3)[0] # layers in grid3 
    
    # Create list of batches of grid3 nodes that will be merged in grid23 nodes #
    i3_23    = [np.arange(i3,i3+n1) for i3 in range(inds3[0],inds3[-n1],n1)] # layers to transition from grid3 to grid23

    # Initialize the 10 lists of the new grid23 nodes #
    g23dz,g23mass,g23rho,g23Tz,g23gt,g23age,g23bdm,g23lwc,g23r2,g23inds = ([] for ii in range(10))
    # Fill in the lists with the batch properties #
    for bb,batch in enumerate(i3_23):
        # Properties of the new grid23 nodes #
        g23dz.append(np.sum(self.dz[batch])) # sum thickness
        g23mass.append(np.sum(self.mass[batch])) # sum mass
        g23rho.append(g23mass[bb]/g23dz[bb])
        bb_Tz0   = np.sum(self.Tz[batch]*self.mass[batch])
        g23Tz.append(bb_Tz0/g23mass[bb]) # Use a weighted average for temperature (effectively the enthalpy)
        g23gt.append(23) #gridtrack
        g23age.append(np.mean(self.age[batch])) # mean age
        g23bdm.append(np.mean(self.bdot_mean[batch])) #mean bdotmean
        g23lwc.append(np.sum(self.LWC[batch])) # sum for lwc
        g23r2.append(np.mean(self.r2[batch])) # mean for r2
        g23inds.append(batch) #old indices of the nodes merged into grid23
    
    nfl = np.size(i3_23) #nb of fine nodes lost
    nmg = int(np.size(i3_23)/n1) #nb of medium nodes gained
    
    ##### Not enough nodes in grid2 -> we have to split nodes from grid22 (avoid emptying grid2) #####
    if len(inds2)-nmg<=0:
        ncl = np.ceil((nmg-len(inds2)+1)/n2) #nb of nodes to split from grid22 to grid2 (nb of coarse layers lost)
        ncl = ncl.astype(int) #convert to int for indexing
        
        ### First: merge ncl times batches of n2 layers from grid23 to grid22 ###
        ## Enough nodes in initial grid23 to form the ncl coarse nodes ##
        if ncl*n2<=len(inds23):
            # Create list of batches of grid23 nodes that will be merged in grid22 nodes #
            i23_22 = [np.arange(i23,i23+n2) for i23 in range(inds23[0],inds23[int(ncl*n2-1)],n2)] # layers to transition from grid23 to grid22
            hi23 = np.size(i23_22) #highest node in the self.dz[inds23] array that will still be part of grid23
            # Initialize the 9 lists of the new grid22 nodes #
            g22dz,g22mass,g22rho,g22Tz,g22gt,g22age,g22bdm,g22lwc,g22r2 = ([] for ii in range(9))
            # Fill in the lists with the batch properties #
            for bb,batch in enumerate(i23_22):
                # Properties of the new grid22 nodes #
                g22dz.append(np.sum(self.dz[batch])) # sum thickness
                g22mass.append(np.sum(self.mass[batch])) # sum mass
                g22rho.append(g22mass[bb]/g22dz[bb])
                bb_Tz0   = np.sum(self.Tz[batch]*self.mass[batch])
                g22Tz.append(bb_Tz0/g22mass[bb]) # Use a weighted average for temperature (effectively the enthalpy)
                g22gt.append(22) #gridtrack
                g22age.append(np.mean(self.age[batch])) # mean age
                g22bdm.append(np.mean(self.bdot_mean[batch])) #mean bdotmean
                g22lwc.append(np.sum(self.LWC[batch])) # sum for lwc
                g22r2.append(np.mean(self.r2[batch])) # mean for r2
        ## Not enough nodes in initial grid23 to form the ncl coarse nodes -> also merge g23 nodes ##
        elif ncl*n2>len(inds23):
            hi23  = len(inds23) #no nodes of the self.dz[inds23] will still be part of grid23
            ng23l = ncl*n2-len(inds23) #nb of nodes of g23 that will contribute to the merging
            # Create list of batches of the grid23 nodes to be merged in grid22 nodes #
            i23_22 = [np.arange(i23,i23+n2) for i23 in range(inds23[0],inds23[-n2],n2)] # layers to transition from grid23 to grid22
            rem23  = np.arange(i23_22[-1][-1]+1,inds23[-1]+1) #remaining nodes that did not create an entire grid22 node
            i0g23  = n2-len(rem23) #index until which we have to take g23 nodes to compensate for rem23 not having enough nodes
            for sublist in g23inds[0:i0g23]:
                rem23 = np.append(rem23,sublist) #progressively append the g23inds sublists to fill in rem23
            i23_22.append(rem23) #append the batch overlapping grid23 and g23
            for lsi in range(i0g23,ng23l,n2):
                g23i_toap = [ii for sublist in g23inds[lsi:lsi+n2] for ii in sublist] #the g23 indices that form a single batch for the new grid22 layers
                i23_22.append(np.array(g23i_toap)) #append to the list of all indices contributing to the grid22 new nodes' formation
            # Remove the g23 nodes that are going to grid22 from the g23 lists #
            g23dz   = g23dz[ng23l:]
            g23mass = g23mass[ng23l:]
            g23rho  = g23rho[ng23l:]
            g23Tz   = g23Tz[ng23l:]
            g23gt   = g23gt[ng23l:]
            g23age  = g23age[ng23l:]
            g23bdm  = g23bdm[ng23l:]
            g23lwc  = g23lwc[ng23l:]
            g23r2   = g23r2[ng23l:]
            g23inds = g23inds[ng23l:]            
            # Initialize the 9 lists of the new grid22 nodes #
            g22dz,g22mass,g22rho,g22Tz,g22gt,g22age,g22bdm,g22lwc,g22r2 = ([] for ii in range(9))
            # Fill in the lists with the batch properties #
            for bb,batch in enumerate(i23_22):
                # Properties of the new grid22 nodes #
                g22dz.append(np.sum(self.dz[batch])) # sum thickness
                g22mass.append(np.sum(self.mass[batch])) # sum mass
                g22rho.append(g22mass[bb]/g22dz[bb])
                bb_Tz0   = np.sum(self.Tz[batch]*self.mass[batch])
                g22Tz.append(bb_Tz0/g22mass[bb]) # Use a weighted average for temperature (effectively the enthalpy)
                g22gt.append(22) #gridtrack
                g22age.append(np.mean(self.age[batch])) # mean age
                g22bdm.append(np.mean(self.bdot_mean[batch])) #mean bdotmean
                g22lwc.append(np.sum(self.LWC[batch])) # sum for lwc
                g22r2.append(np.mean(self.r2[batch])) # mean for r2

        ### Second: split the ncl highest grid22 nodes in ncl*n2 grid2 nodes ###
        i22_2 = inds22[0:ncl] #nodes to transition from grid22 to grid 2
        # Initialize the 9 lists of the new grid2 nodes #
        g2dz,g2mass,g2rho,g2Tz,g2gt,g2age,g2bdm,g2lwc,g2r2 = ([] for ii in range(9))
        # Fill in the lists with the nodes' properties #
        for i22 in i22_2:
            # Properties of the new grid2 nodes #
            g2dz   = np.append(g2dz,self.dz[i22]/n2*np.ones(n2))
            g2rho  = np.append(g2rho,self.rho[i22]*np.ones(n2))
            g2Tz   = np.append(g2Tz,self.Tz[i22]*np.ones(n2))
            g2gt   = np.append(g2gt,2*np.ones(n2))
            g2age  = np.append(g2age,np.linspace(self.age[i22],self.age[i22+1],n2)) #assume linearly increasing age until layer below
            g2bdm  = np.append(g2age,self.bdot_mean[i22]*np.ones(n2)) 
            g2lwc  = np.append(g2lwc,self.LWC[i22]/n2*np.ones(n2))
            g2r2   = np.append(g2r2,self.r2[i22]*np.ones(n2))
        g2mass = g2dz*g2rho
            
        ### Now there are enough layers in grid2 combined with g2 to form nfl new nodes in grid1 ###
        i2_1 = inds2[0:] #all nodes of grid2 will be split into grid1 nodes
        #ng2l = len(inds2)-nmg #nb of nodes from g2 that will also be split in grid1 nodes
        ng2l = nmg-len(inds2) #nb of nodes from g2 that will also be split in grid1 nodes
        if ng2l==0: #Case where there are just enough nodes from grid2 for the splitting (g2 created for non-empty grid2)
            ig2_1 = np.array([]) #no nodes of g2 will be split
        else:
            ig2_1 = np.arange(0,ng2l) #nodes of g2 that will also be split in grid1 nodes (indices are on the g2 grid!)
        # Initialize the 9 lists of the new grid1 nodes #
        g1dz,g1mass,g1rho,g1Tz,g1gt,g1age,g1bdm,g1lwc,g1r2 = ([] for ii in range(9))
        # Fill in the lists with the nodes' properties #
        for i2 in i2_1: #Proceed first to the splitting of the grid2 nodes
            # Properties of the new grid1 nodes #
            g1dz   = np.append(g1dz,self.dz[i2]/n1*np.ones(n1))
            g1rho  = np.append(g1rho,self.rho[i2]*np.ones(n1))
            g1Tz   = np.append(g1Tz,self.Tz[i2]*np.ones(n1))
            g1gt   = np.append(g1gt,1*np.ones(n1))
            g1age  = np.append(g1age,np.linspace(self.age[i2],self.age[i2+1],n1)) #assume linearly increasing age until layer below
            g1bdm  = np.append(g1bdm,self.bdot_mean[i2]*np.ones(n1)) 
            g1lwc  = np.append(g1lwc,self.LWC[i2]/n1*np.ones(n1))
            g1r2   = np.append(g1r2,self.r2[i2]*np.ones(n1))
        for ig2 in ig2_1: #Then proceed to the splitting of g2 nodes (if necessary, otherwise ig2_1 is empty)
            # Properties of the new grid1 nodes #
            g1dz   = np.append(g1dz,g2dz[ig2]/n1*np.ones(n1))
            g1rho  = np.append(g1rho,g2rho[ig2]*np.ones(n1))
            g1Tz   = np.append(g1Tz,g2Tz[ig2]*np.ones(n1))
            g1gt   = np.append(g1gt,1*np.ones(n1))
            g1age  = np.append(g1age,np.linspace(g2age[ig2],g2age[ig2+1],n1)) #assume linearly increasing age until layer below
            g1bdm  = np.append(g1bdm,g2bdm[ig2]*np.ones(n1))
            g1lwc  = np.append(g1lwc,g2lwc[ig2]/n1*np.ones(n1))
            g1r2   = np.append(g1r2,g2r2[ig2]*np.ones(n1))
        g1mass = g1dz*g1rho
        # Remove the g2 nodes that are going to grid1 from the g1 lists #
        g2dz   = g2dz[ng2l:]
        g2mass = g2mass[ng2l:]
        g2rho  = g2rho[ng2l:]
        g2Tz   = g2Tz[ng2l:]
        g2gt   = g2gt[ng2l:]
        g2age  = g2age[ng2l:]
        g2bdm  = g2bdm[ng2l:]
        g2lwc  = g2lwc[ng2l:]
        g2r2   = g2r2[ng2l:]       
                    
    ##### Enough nodes in grid2 -> simply split nodes from grid2 to grid1 #####
    elif len(inds2)-nmg>0:
        ncl  = 0 #no node from grid22 has been split to g2
        hi23 = 0 #all nodes of the self.dz[inds23] array will still belong to grid23 (since none has been merged into grid22)
        # No new grid22 nodes #
        g22dz,g22mass,g22rho,g22Tz,g22gt,g22age,g22bdm,g22lwc,g22r2 = ([] for ii in range(9))
        # No new grid2 nodes #
        g2dz,g2mass,g2rho,g2Tz,g2gt,g2age,g2bdm,g2lwc,g2r2 = ([] for ii in range(9))
        ### Find nodes in grid2 that will form the nfl new nodes in grid1 ###
        i2_1 = inds2[0:nmg] #highest nodes of grid2 will be split into grid1 nodes
        # Initialize the 9 lists of the new grid1 nodes #
        g1dz,g1mass,g1rho,g1Tz,g1gt,g1age,g1bdm,g1lwc,g1r2 = ([] for ii in range(9))
        # Fill in the lists with the nodes' properties #
        for i2 in i2_1: #Proceed first to the splitting of the grid2 nodes
            # Properties of the new grid1 nodes #
            g1dz   = np.append(g1dz,self.dz[i2]/n1*np.ones(n1))
            g1rho  = np.append(g1rho,self.rho[i2]*np.ones(n1))
            g1Tz   = np.append(g1Tz,self.Tz[i2]*np.ones(n1))
            g1gt   = np.append(g1gt,1*np.ones(n1))
            g1age  = np.append(g1age,np.linspace(self.age[i2],self.age[i2+1],n1)) #assume linearly increasing age until layer below
            g1bdm  = np.append(g1bdm,self.bdot_mean[i2]*np.ones(n1))
            g1lwc  = np.append(g1lwc,self.LWC[i2]/n1*np.ones(n1))
            g1r2   = np.append(g1r2,self.r2[i2]*np.ones(n1))
        g1mass = g1dz*g1rho
    
    ##### Concatenate everything together #####
    self.dz         = np.concatenate((self.dz[inds1],g1dz,self.dz[inds2][nmg:],g2dz,self.dz[inds22][ncl:],g22dz,self.dz[inds23][hi23:],g23dz,self.dz[inds3][nfl:]))
    self.z          = self.dz.cumsum(axis=0)
    self.z          = np.concatenate(([0], self.z[:-1]))
    self.rho        = np.concatenate((self.rho[inds1],g1rho,self.rho[inds2][nmg:],g2rho,self.rho[inds22][ncl:],g22rho,self.rho[inds23][hi23:],g23rho,self.rho[inds3][nfl:]))
    self.Tz         = np.concatenate((self.Tz[inds1],g1Tz,self.Tz[inds2][nmg:],g2Tz,self.Tz[inds22][ncl:],g22Tz,self.Tz[inds23][hi23:],g23Tz,self.Tz[inds3][nfl:]))
    self.mass       = np.concatenate((self.mass[inds1],g1mass,self.mass[inds2][nmg:],g2mass,self.mass[inds22][ncl:],g22mass,self.mass[inds23][hi23:],g23mass,self.mass[inds3][nfl:]))
    self.mass_sum   = self.mass.cumsum(axis = 0)
    self.age        = np.concatenate((self.age[inds1],g1age,self.age[inds2][nmg:],g2age,self.age[inds22][ncl:],g22age,self.age[inds23][hi23:],g23age,self.age[inds3][nfl:]))
    self.bdot_mean  = np.concatenate((self.bdot_mean[inds1],g1bdm,self.bdot_mean[inds2][nmg:],g2bdm,self.bdot_mean[inds22][ncl:],g22bdm,self.bdot_mean[inds23][hi23:],g23bdm,self.bdot_mean[inds3][nfl:]))
    self.LWC        = np.concatenate((self.LWC[inds1],g1lwc,self.LWC[inds2][nmg:],g2lwc,self.LWC[inds22][ncl:],g22lwc,self.LWC[inds23][hi23:],g23lwc,self.LWC[inds3][nfl:]))
    self.gridtrack  = np.concatenate((self.gridtrack[inds1],g1gt,self.gridtrack[inds2][nmg:],g2gt,self.gridtrack[inds22][ncl:],g22gt,self.gridtrack[inds23][hi23:],g23gt,self.gridtrack[inds3][nfl:]))
    self.r2         = np.concatenate((self.r2[inds1],g1r2,self.r2[inds2][nmg:],g2r2,self.r2[inds22][ncl:],g22r2,self.r2[inds23][hi23:],g23r2,self.r2[inds3][nfl:]))

    self.sigma      = (self.mass+self.LWC*RHO_W_KGM)*self.dx*GRAVITY
    self.sigma      = self.sigma.cumsum(axis = 0)

    #self.bdot_mean  = (np.concatenate(([self.mass_sum[0] / (RHO_I * S_PER_YEAR)], self.mass_sum[1:] * self.t / (self.age[1:] * RHO_I))))*self.c['stpsPerYear']*S_PER_YEAR

    return self.dz, self.z, self.rho, self.Tz, self.mass, self.sigma, self.mass_sum, self.age, self.bdot_mean, self.LWC, self.gridtrack, self.r2



def init_regrid22(self):
    '''
    Splits the column in 5 grids: grid1(high res)-grid2(low res)-grid22(v. low res)-grid23(low res)-grid3(high res)
    Used in firn_density_spin for the initial regridding. 
    '''

    grid1b          = self.c['grid1bottom'] # bottom of grid1
    grid2b          = self.c['grid2bottom'] # bottom of grid2
    self.nodestocombine = self.c['nodestocombine'] # nb layers to combine from grid1 to grid2
    self.nodestocombine2 = self.c['multnodestocombine']*self.c['nodestocombine'] # nb layers to combine from grid1 to grid22
    ind1            = np.where(self.z<grid1b)[0] # layers of grid1
    if self.nodestocombine2 > 0: # grid22 process is turned on
        ind2            = np.intersect1d(np.where(self.z>=grid1b)[0],np.where(self.z<grid2b)[0]) # layers of grid2
        ind22           = np.where(self.z>=grid2b)[0] # layers of grid22
        ind23           = ind22[-self.nodestocombine2:] # layers of grid23
        ind22           = ind22[0:-self.nodestocombine2] # remove layers of grid23 from grid22
        grid1z          = self.z[ind1] # z values grid1
        grid2z          = self.z[0:ind22[0]][ind2[0]::self.nodestocombine] # z values grid2 after merging of batches of nodestocombine
        grid22z         = self.z[0:ind23[0]][ind22[0]::self.nodestocombine2] #VV z values of grid22 after merging of batches of nodestocombine2
        grid23z         = self.z[ind23[0]::self.nodestocombine] # layers in grid23 have thickness determined by nodestocombine
        self.z          = np.concatenate((grid1z,grid2z,grid22z,grid23z))
    elif self.nodestocombine2 == 0: # in effect: grid22 process is turned off
        ind2            = np.where(self.z>=grid1b)[0] # layers of grid22
        grid1z          = self.z[ind1] # z values grid1
        grid2z          = self.z[ind2[0]::self.nodestocombine] # z values grid2 after merging of batches of nodestocombine
        self.z          = np.concatenate((grid1z,grid2z))
    grid3z          = self.z[-1] + np.cumsum(self.dz[-1*self.nodestocombine:]) # create grid3
    self.z          = np.concatenate((self.z,grid3z))
    self.dz         = np.diff(self.z)
    self.dz         = np.append(self.dz, self.dz[-1])
    self.gridLen    = len(self.z)
    self.dx         = np.ones(self.gridLen)
    # Define self.gridtrack
    if self.nodestocombine2 > 0: # grid22 process is turned on
        self.gridtrack = np.concatenate((np.ones_like(grid1z),2*np.ones_like(grid2z),22*np.ones_like(grid22z),23*np.ones_like(grid23z),3*np.ones_like(grid3z)))
    elif self.nodestocombine2 == 0: # in effect: grid22 process is turned off
        self.gridtrack = np.concatenate((np.ones_like(grid1z),2*np.ones_like(grid2z),3*np.ones_like(grid3z)))
    
    #print('After regrid, grid length is', self.gridLen)

    return self.z, self.dz, self.gridLen, self.dx, self.gridtrack


