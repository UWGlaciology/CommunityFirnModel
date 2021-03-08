#!/usr/bin/env python
from constants import *
import numpy as np
from diffusion import heatDiff

'''
Functions to handle meltwater percolation.
Bucket schemes only in this file.
'''

def percolation_bucket(self, iii):

    '''
    This is the bucket scheme that allows liquid water to persist in the firn.
    It includes consideration of irreducible liquid water content (LWC) an maximum
    LWC. Water that encounters a slab of a certain density (impermeable_rho) will
    not percolate through.

    LWC is in volume (m^3), and since we are working in one dimension we assume that it is m^3/m^2.
    
    @author: maxstev
    '''

    # maxpore_f                 = 2.0   # factor by which the maximum filled porespace can exceed the irreducible saturation.

    try:
        impermeable_rho = self.c['impermeable_rho']
    except:
        impermeable_rho         = 830.  # impermeable lens density.

    if np.any(self.LWC<0):
        print('ERROR: negative LWC')
        print('(model will continue to run)')

    melt_volume_IE          = self.snowmeltSec[iii] * S_PER_YEAR    # meters
    melt_volume_WE          = melt_volume_IE * RHO_I_MGM            # meters
    melt_mass               = melt_volume_WE * 1000.                # kg
    heat_to_freeze          = melt_mass * LF_I                      # amount of heat needed to refreeze the melt (J)
    ind1a                   = np.where(self.mass_sum <= melt_mass)[0]   # indicies of boxes that will be melted away
    num_boxes_melted        = len(ind1a)+1                              # number of boxes that melt away, include the box that is partially melted
    ind1                    = np.where(self.mass_sum > melt_mass)[0][0] # index which will become the new surface

    ### pm is the partial melt (the model volume that has a portion melted away)
    pm_mass                 = self.mass_sum[ind1] - melt_mass       # the remaining mass of the PM box
    pm_dz                   = pm_mass / self.rho[ind1]              # remaining thickness
    pm_porespace            = (1 - self.rho[ind1]/RHO_I) * pm_dz    # porespace in the PM box
    pm_rho                  = self.rho[ind1]                        # density of the PM box
    pm_lwc                  = self.LWC[ind1]/self.dz[ind1] * pm_dz  # LWC of the PM box

    melt_boxes_LWC_vol      = np.sum(self.LWC[0:ind1+1]) - pm_lwc #include the water mass from the volumes/nodes that melt (currently does not include from the partial melt box)
    melt_boxes_LWC_mass     = melt_boxes_LWC_vol * RHO_W_KGM
    melt_mass_a             = melt_mass + melt_boxes_LWC_mass
    melt_vol_a              = melt_mass_a / RHO_W_KGM

    ###################################
    ### Regrid after melt
    ### Melted boxes are accomodated by just adding more (new) boxes at the bottom of the column
    ### Beware of this if you are not modeling to firn-ice transition depth.
    divider                 = num_boxes_melted
    self.rho                = np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted)))
    self.LWC                = np.concatenate((self.LWC[ind1:-1] , self.LWC[-1]*np.ones(num_boxes_melted)))
    self.LWC[0]             = pm_lwc

    self.age                = np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_melted)))
    # self.dz               = np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted))) # this splits the last box into many.
    self.dz                 = np.concatenate((self.dz[ind1:-1] , self.dz[-1]*np.ones(num_boxes_melted))) # this adds new boxes at the bottom.
    self.dz[0]              = pm_dz
    self.Dcon               = np.concatenate((self.Dcon[ind1:-1] , self.Dcon[-1]*np.ones(num_boxes_melted)))
    self.dzn                = np.concatenate((np.zeros(num_boxes_melted-1), self.dz[0:])) # this will fail in the case that there is no PM box, i.e. the melt mass is exactly equal to the mass of one or several boxes.
    ### old version:
    # self.dzn                = np.concatenate((np.zeros(num_boxes_melted-1), self.dz[0:])) #this is not quite right because is assumes compaction for the pm box is zero.
    ###
    self.dzn                = self.dzn[0:self.compboxes]
    self.Tz                 = np.concatenate((self.Tz[ind1:-1] , self.Tz[-1]*np.ones(num_boxes_melted)))
    self.bdot_mean          = np.concatenate((self.bdot_mean[ind1:-1] , self.bdot_mean[-1]*np.ones(num_boxes_melted)))
    self.z                  = self.dz.cumsum(axis = 0)
    self.z                  = np.concatenate(([0] , self.z[:-1]))
    self.mass               = self.rho * self.dz
    ###################################

    ##########################################
    ### now working all with the new grid ####
    ##########################################
    porosity                = 1 - self.rho / RHO_I      # porosity (unitless)
    porespace_vol           = porosity * self.dz        # pore space volume (meters) of each box - volume of air + water
    porespace_air           = porespace_vol - self.LWC  # pore space that is filled with air (meters)

    cold_content            = CP_I * self.mass * (T_MELT - self.Tz) # cold content of each box, i.e. how much heat to bring it to 273K (J)
    cold_content_sum        = cold_content.cumsum(axis=0)
    refreeze_mass_pot       = cold_content / LF_I   # how much meltwater mass could be refrozen in each node due to cold content (kg)
    refreeze_mass_pot_sum   = refreeze_mass_pot.cumsum(axis=0) # total meltwater mass could be refrozen due to cold content (kg)

    ### calculate what the values will be after refreeze happens (pot stands for potential)
    rho_pot                 = (self.mass + refreeze_mass_pot) / self.dz # what the density of each node would be if it refroze the maximum possible
    porosity_pot            = 1 - rho_pot / RHO_I # potential porosity of each node if it refroze the maximum possible
    porespace_vol_pot       = porosity_pot * self.dz # volume of potential porespace
    # porespace_air_pot       = porespace_vol_pot - self.LWC # I don't think this gets used.

    ### Irreducible calculations
    Wmi                     = 0.057 * (RHO_I - rho_pot) / rho_pot + 0.017 # water per snow-plus-water mass irreducible liquid water content, Langen eqn 3 unitless)
    Swi                     = Wmi / (1 - Wmi) * (rho_pot * RHO_I) / (1000 * (RHO_I - rho_pot))  #irreducible water saturation, volume of water per porespace volume (unitless), Colbeck 1972
    # Swi  = 0.02 * np.ones_like(self.dz)
    irreducible_mass_pot    = Swi * porespace_vol_pot * RHO_W_KGM # [kg] mass of irreducible water for each node (potential - does not consider how much is already present in each node)
    irreducible_vol_pot     = irreducible_mass_pot / RHO_W_KGM # [m] volume of irreducible water for each node

    overIRRinds = np.where(self.LWC>irreducible_vol_pot)[0] # Indicies of nodes where the current LWC volume is more than irreducible
    # if overIRRinds.size > 0: # there are nodes where there is extra water

    # maxpore                 = 0.1 # upper limit on what percentage of the porosity can be filled with water.
    maxpore = 0.3

    # maxLWC1                 = porespace_vol * maxpore   # maximum volume of water that can be stored in each node (meters)
    maxLWC1                 = porespace_air * maxpore # [m]
    maxLWC2                 = ((917.0 * self.dz) - self.mass) / RHO_W_KGM # double check that the LWC does not get too large. 
    maxLWC                  = np.minimum(maxLWC1 , maxLWC2) # [m] maximum volume of water that can exist in each node

    maxLWC[self.rho>impermeable_rho] = 0
    maxLWC_mass             = maxLWC * RHO_W_KGM        # [kg] mass of the maximum volume of water
    maxLWC1_pot             = porespace_vol_pot * maxpore   # [m] maximum volume of water that can be stored in each node (meters)
    maxLWC2_pot             = ((917.0 * self.dz) - (self.mass + refreeze_mass_pot)) / RHO_W_KGM # [m], volume; double check that the LWC does not get too large (i.e. that it could not push the density beyond ice density). 
    maxLWC_pot              = np.minimum(maxLWC1_pot , maxLWC2_pot) # [m], volume
    # maxLWC_pot[rho_pot>impermeable_rho] = 0
    maxLWC_mass_pot         = maxLWC_pot * RHO_W_KGM        # [kg] mass of the maximum volume of water allowed in each node


    liquid_storage_vol_pot  = irreducible_vol_pot - self.LWC # [m] how much additional water can be stored in each node
    liquid_storage_mass_pot = liquid_storage_vol_pot * RHO_W_KGM # [kg] how much additional water can be stored in each node

    extra_liquid_mass       = np.sum(self.LWC[self.LWC > irreducible_vol_pot] * RHO_W_KGM - irreducible_mass_pot[self.LWC > irreducible_vol_pot])
    storage_mass_pot        = liquid_storage_mass_pot + refreeze_mass_pot #how much can be refrozen plus how much will stick around due to capillary
    storage_mass_pot_sum    = storage_mass_pot.cumsum(axis=0)
    total_liquid_mass       = melt_mass_a + extra_liquid_mass
    
    try:
        ind_p   = np.where(storage_mass_pot_sum >= total_liquid_mass)[0][0] # the layer that water will percolate to
    except: # all of the liquid is runoff.
        ind_p   = 0
    ###################################

    ### if there is an impermeable layer, block water from getting through
    if np.any(self.rho[0:ind_p+1] >= impermeable_rho):

        ind_p                   = np.where(self.rho >= impermeable_rho)[0][0] #- 1 # the index of the node that has density greater than the impermeable density
        id1                     = np.where(self.LWC >  irreducible_vol_pot)[0] # indices where the liquid water content is greater than the irreducible
        id2                     = id1[id1<ind_p] 

        extra_liquid_mass       = np.sum(self.LWC[id2] * RHO_W_KGM) - np.sum(irreducible_mass_pot[id2])
        storage_mass_pot        = liquid_storage_mass_pot[0:ind_p] + refreeze_mass_pot[0:ind_p] #how much can be refrozen plus how much will stick around due to capillary
        storage_mass_pot_sum    = storage_mass_pot.cumsum(axis=0)
        total_liquid_mass       = (melt_mass_a + extra_liquid_mass) #*0.5

        ### first, refreeze where possible
        self.mass[0:ind_p]      = self.mass[0:ind_p] + refreeze_mass_pot[0:ind_p]
        self.rho[0:ind_p]       = self.mass[0:ind_p] / self.dz[0:ind_p]
        self.Tz[0:ind_p]        = T_MELT

        # nomoreliquid = np.where(self.rho>830)[0]
        # maxLWC_mass_pot[nomoreliquid] = 0
        # maxLWC_pot[nomoreliquid] = 0


        mass_frozen             = np.sum(refreeze_mass_pot[0:ind_p])
        if mass_frozen >= total_liquid_mass:
            total_liquid_mass   = 0
        else:
            total_liquid_mass   = total_liquid_mass - mass_frozen

        ### then, fill up the nodes above the ice slab
        maxLWC_mass_pot_f       = np.flipud(maxLWC_mass_pot[0:ind_p])
        maxLWC_mass_pot_f_sum   = maxLWC_mass_pot_f.cumsum(axis=0)

        if total_liquid_mass >= np.sum(maxLWC_mass_pot_f): # all porespace gets filled and there is runoff      
            self.LWC[0:ind_p]       = maxLWC_pot[0:ind_p] # each node gets the maximum allowed
            # stored_water_vol      = np.sum(self.LWC[0:ind_p]) # can calculate how much runoff there is, need to consider how much LWC there was previously
        
        else: # fill up however much porespace is needed to accomodate the meltwater

            ind_f                   = np.where(maxLWC_mass_pot_f_sum > total_liquid_mass)[0][0] #index on the flipped grid
            ind_g                   = ind_p - 1 - ind_f #index on the real grid.
            self.LWC[ind_g+1:ind_p] = maxLWC_mass_pot[ind_g + 1:ind_p] / RHO_W_KGM # fill the indices up with the maximum allowed water
            lv_mass                 = total_liquid_mass - np.sum(maxLWC_mass_pot[ind_g + 1:ind_p])  # leftover volume
            self.LWC[ind_g]         = lv_mass / RHO_W_KGM                       # put that into the ind_g node
    ###################################
   
    ### there is not an impermeable layer, water goes to layer ind_p
    elif ind_p>0: 

        ### first, up to ind_p (not inclusive)
        self.mass[0:ind_p]      = self.mass[0:ind_p] + refreeze_mass_pot[0:ind_p]
        self.rho[0:ind_p]       = self.mass[0:ind_p] / self.dz[0:ind_p]
        lwc_old                 = np.copy(self.LWC)
        self.LWC[0:ind_p]       = irreducible_mass_pot[0:ind_p] / RHO_W_KGM
        self.Tz[0:ind_p]        = T_MELT
        lw_mass_retained        = np.sum(refreeze_mass_pot[0:ind_p]) + np.sum(irreducible_mass_pot[0:ind_p]) - np.sum(lwc_old[0:ind_p] * RHO_W_KGM)
        lw_mass_remaining       = total_liquid_mass - lw_mass_retained # mass left that will go into the ind_p node

        ### now deal with the very last node where there may be just freezing or both freezing and some amount of retention
        if lw_mass_remaining <= refreeze_mass_pot[ind_p]: # all remaining water freezes
            latent_heat_released    = lw_mass_remaining * LF_I
            self.Tz[ind_p]          = self.Tz[ind_p] + latent_heat_released / (CP_I * self.mass[ind_p])
            self.mass[ind_p]        = self.mass[ind_p] + lw_mass_remaining
            self.rho[ind_p]         = self.mass[ind_p] / self.dz[ind_p]
            self.LWC[ind_p]         = 0
            
        else:   # some refreeze, some sticks around 
            self.mass[ind_p]        = self.mass[ind_p] + refreeze_mass_pot[ind_p]
            self.rho[ind_p]         = self.mass[ind_p] / self.dz[ind_p]
            self.LWC[ind_p]         = (lw_mass_remaining - refreeze_mass_pot[ind_p]) / RHO_W_KGM
            self.Tz[ind_p]          = T_MELT
    ###################################

    self.LWC[self.LWC<0] = 0

    if np.any(self.rho>917.0):
        print('high rho in melt.py', np.max(self.rho))
        self.rho[self.rho>917.0]=917.0
    

    return self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn, self.LWC
###################################################
### End percolation_bucket ########################
###################################################


def bucketVV(self, iii):
    '''
    Bucket scheme for meltwater percolation, retention, refreezing, and runoff

    @author: verjans
    Used in Verjans et al. (2019)
    '''

    rhoimp      = 830. #impermeable density, threshold to generate runoff, can be changed
    irr         = 0.02 * np.ones_like(self.dz) # Irreducible water content, this is a proportion of the available pore space
    #irr = 0.06 * np.ones_like(self.dz) # Crocus value (Reijmer 2012)
    CLparam     = 0 # Set this to 1 to use Coleou and Lesaffre 1998 parameterisation
    if CLparam == 1:
        irr     = np.zeros_like(self.dz) #calculated below (twice: once before freezing and once after freezing)
    percbottom = 0 #if set to 1: allows percolation until the depth where all nodes have rho>=rhoimp (i.e. allows percolation through ice lenses)
    ##### First: melting of the surface layers, taken from melt.py #####
    melt_volume_IE      = self.snowmeltSec[iii] * S_PER_YEAR # This still has to be checked by Max (division by self.c['stpsPerYear']?) [m]
    melt_volume_WE      = melt_volume_IE * RHO_I_MGM # [m]
    melt_mass           = melt_volume_WE * 1000. # [kg]   
    initial_lwc         = 1 * self.LWC
    
#    heat_to_freeze             = melt_mass * LF_I                         # amount of heat needed to refreeze the melt (J)
    ind1a               = np.where(self.mass_sum <= melt_mass)[0] # indices of boxes that will be melted away
    num_boxes_melted    = len(ind1a)+1 # number of boxes that melt away, include the box that is partially melted
    ind1                = np.where(self.mass_sum > melt_mass)[0][0] # index which will become the new surface
    # pm is the partial melt (the model volume that has a portion melted away)
    pm_mass             = self.mass_sum[ind1] - melt_mass # the remaining mass of the PM box [kg]
    pm_dz               = pm_mass / self.rho[ind1] # remaining thickness [m]
#    pm_porespace             = (1 - self.rho[ind1]/RHO_I) * pm_dz # porespace in the PM box
    pm_rho              = self.rho[ind1] # density of the PM box [kg/m3]
    pm_lwc              = self.LWC[ind1]/self.dz[ind1] * pm_dz # LWC of the PM box [m]
    
    melt_boxes_LWC_vol  = np.sum(self.LWC[0:ind1+1]) - pm_lwc #include the LWC from the boxes that melt (currently does not include from the partial melt box) [m]
    melt_boxes_LWC_mass = melt_boxes_LWC_vol * RHO_W_KGM #include the mass of LWC from the boxes that melt (currently does not include from the partial melt box) [kg]
    melt_mass_a         = melt_mass + melt_boxes_LWC_mass #total liq water from melted boxes(due to melting + LWC at previous time step) [kg]
    melt_vol_a          = melt_mass_a / RHO_W_KGM #total liq water from melted boxes(due to melting + LWC at previous time step) [m]
    # I think, melt_vol_a is the total liquid water input for this time step of the CFM
#    pm_plwc = self.PLWC_mem[ind1]/self.dz[ind1] * pm_dz

    ## Melted boxes are accomodated by just adding more (new) boxes at the bottom of the column
    ## Beware of this if you are not modeling to firn-ice transition depth.
    divider         = num_boxes_melted #VV nb of melted boxes, including the partially melted
    self.rho        = np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted))) #VV add at bottom of column as many layers as were melted away
    self.LWC        = np.concatenate((self.LWC[ind1:-1] , self.LWC[-1]*np.ones(num_boxes_melted)))
    # self.LWC        = np.concatenate((self.LWC[ind1:-1] , np.zeros(num_boxes_melted))) # This is better but should be equivalent as last layer should be an ice layer -> with 0 LWC
    self.LWC[0]     = pm_lwc #VV LWC calculated for the partially melted layer
        
    self.age        = np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_melted)))
    # self.dz                  = np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted))) # this splits the last box into many.
    self.dz         = np.concatenate((self.dz[ind1:-1] , self.dz[-1]*np.ones(num_boxes_melted))) # this adds new boxes at the bottom.
    self.dz[0]      = pm_dz # VV dz calculated for the partially melted layer
    
    if np.any(self.dz<1e-6): # VV change 09/12/2020
        self.dz = np.maximum(self.dz,1e-6) #avoids dz[0] to be extremely small if melt amount makes pm_dz very close to 0

    self.Dcon       = np.concatenate((self.Dcon[ind1:-1] , self.Dcon[-1]*np.ones(num_boxes_melted)))
    self.dzn        = np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
    self.dzn        = self.dzn[0:self.compboxes]
    self.Tz         = np.concatenate((self.Tz[ind1:-1] , self.Tz[-1]*np.ones(num_boxes_melted)))
    self.r2         = np.concatenate((self.r2[ind1:-1] , self.r2[-1]*np.ones(num_boxes_melted)))
    self.bdot_mean  = np.concatenate((self.bdot_mean[ind1:-1] , self.bdot_mean[-1]*np.ones(num_boxes_melted)))
    self.z          = self.dz.cumsum(axis = 0)
    self.z          = np.concatenate(([0] , self.z[:-1]))
    self.mass       = self.rho * self.dz
    
    if self.doublegrid: # if we have doublegrid -> need to adjust gridtrack
        meltgridtrack  = np.concatenate((self.gridtrack[ind1:-1],self.gridtrack[-1]*np.ones(num_boxes_melted)))
    elif self.doublegrid==False:
        meltgridtrack = np.zeros_like(self.dz) # just return a zero array
    ## Grid should be ready

    porosity           = 1 - self.rho/RHO_I # Definition of porosity [/]
    porespace_vol      = porosity * self.dz # Pore space of each layer [m]
    porosity_refr      = porosity*RHO_I/RHO_W_KGM # space available for liq water volume once refrozen, Wever 2014 (9) [/]
    porespace_refr_vol = porosity_refr*self.dz # Available pore space of each layer [m]
    
    cold_content            = CP_I * self.mass * (T_MELT - self.Tz) # cold content of each box, i.e. how much heat to bring it to 273K [J]
    refreeze_mass_pot       = cold_content / LF_I # how much mass of the meltwater could be refrozen due to cold content [kg]
    refreeze_vol_pot        = refreeze_mass_pot/1000. # how much meters of the meltwater could be refrozen due to cold content [m]
    
    self.refrozen = np.zeros_like(self.dz)
    runoff = 0

    ### First refreeze if there is any lwc already in the column ###
    if np.any(self.LWC>0)==True:
        lwc_before_freeze = 1*self.LWC
        lwcpres = np.where(self.LWC > 0)[0] # layers with existing lwc
        coldlay = np.where(self.Tz<273.15)[0]
        frl = np.intersect1d(lwcpres,coldlay) # layers with refreezing
        # if frl.size !=0:
            # print('frl',frl)
        latheat = np.zeros_like(self.dz)
        freeze = np.zeros_like(self.dz)
        freeze[frl] = np.minimum(refreeze_vol_pot[frl],porespace_refr_vol[frl])
        freeze[frl] = np.minimum(self.LWC[frl],freeze[frl])
        self.LWC[frl] -= freeze[frl]
        self.refrozen[frl] += freeze[frl]
        self.mass[frl] += freeze[frl]*1000
        self.rho[frl] = self.mass[frl]/self.dz[frl]
        latheat[frl] = freeze[frl]*1000*LF_I
        cold_content[frl] -= latheat[frl]
        self.Tz[frl] = T_MELT-cold_content[frl]/(CP_I*self.mass[frl])    
        porosity[frl]           = 1 - self.rho[frl]/RHO_I # Definition of porosity [/]
        porespace_vol[frl]      = porosity[frl] * self.dz[frl] # Pore space of each layer [m]
        porosity_refr[frl]      = porosity[frl]*RHO_I/RHO_W_KGM # space available for liq water volume once refrozen, Wever 2014 (9) [/]
        porespace_refr_vol[frl] = porosity_refr[frl]*self.dz[frl] # Available pore space of each layer [m]    
        cold_content[frl]       = CP_I * self.mass[frl] * (T_MELT - self.Tz[frl]) # cold content of each box, i.e. how much heat to bring it to 273K [J]
        refreeze_mass_pot[frl]  = cold_content[frl] / LF_I # how much mass of the meltwater could be refrozen due to cold content [kg]
        refreeze_vol_pot[frl]   = refreeze_mass_pot[frl]/1000. # how much meters of the meltwater could be refrozen due to cold content [m]
    
      
    if CLparam == 1:
        liqmassprop = 0.057*(porosity/(1-porosity)) + 0.017
        irr = liqmassprop/(1-liqmassprop) * (self.rho*RHO_I) / (RHO_W_KGM*(RHO_I-self.rho))
        
    ### Then let the non refrozen water of the column percolate ###
    ### This is only used at the very first time step of the model run (if there is a prescribed LWC from an initial profile) ###
#    print('run time before second for loop =' , time.time()-tic2 , 'seconds')
    if iii == 0: # if water in initial profile, let it percolate according to bucket scheme
        print('Percolation loop')
        irr_limit = np.maximum(irr*porespace_refr_vol,0.) # max amount of lwc that a layer can hold in its pore space
        perclayers = np.where(self.LWC>irr_limit) # layers with lwc in excess of irr
        lwc_above_irr = self.LWC[perclayers[0]]-irr_limit[perclayers[0]]
        if len(lwc_above_irr) > 0: #if there are some layers exceeding irreducible water content
        # Densification (in physics) slightly decreases pore space and thus irr_limit at every time step -> we would have to proceed to entire routine for all layers at irr_limit again and again at every time step
            if ((max(lwc_above_irr) < 1e-6) and (sum(lwc_above_irr)<1e-3)): # if the lwc above irr_limit is very small (due to slight densification): cheat -> add it to runoff
                runoff += sum(lwc_above_irr)
                self.LWC[perclayers[0]] = irr_limit[perclayers[0]]
            elif ((max(lwc_above_irr) >= 1e-6) or (sum(lwc_above_irr)>=1e-3)):
                #print('max(lwc_above_irr), sum(lwc_above_irr) are:',max(lwc_above_irr),sum(lwc_above_irr))
                for index in perclayers[0]:
                    #print('Now we deal with index:',index)
                    ii = 1*index # start from the layer that has lwc in excess of irr
                    toperc = self.LWC[index]-irr_limit[index] # amount of water that has to percolate
                    self.LWC[index] -= toperc # this toperc amount will be lost by the layer
                    p_frozen = np.zeros_like(self.dz) # volume of refrozen water (from toperc) in each layer
                    p_retained = np.zeros_like(self.dz) # volume of retained water (from toperc) in each layer
                    if self.rho[ii] >= rhoimp: #if the layer is impermeable: runoff of all toperc
                        runoff += 1*toperc
                        toperc = 0.
                        
                    while toperc > 0:
                        if self.Tz[ii] < 273.15: # there can be some refreezing
                            p_frozen[ii] = min(refreeze_vol_pot[ii],porespace_refr_vol[ii]) # max freezing possible [mWE] 
                            p_frozen[ii] = min(toperc,p_frozen[ii]) # don't exceed water available
                            self.refrozen[ii] += p_frozen[ii]
                            self.mass[ii] += p_frozen[ii]*1000
                            self.rho[ii] = self.mass[ii]/self.dz[ii]
                            latheat = p_frozen[ii]*1000*LF_I # latent heat released due to the refreezing [J]
                            cold_content[ii] -= latheat # remaining cold content
                            self.Tz[ii] = T_MELT - cold_content[ii]/(CP_I*self.mass[ii]) # temperature is changed accordingly
                            porosity[ii]           = 1 - self.rho[ii]/RHO_I # Definition of porosity [/]
                            porespace_vol[ii]      = porosity[ii] * self.dz[ii] # Pore space of each layer [m]
                            porosity_refr[ii]      = porosity[ii]*RHO_I/RHO_W_KGM # space available for liq water volume once refrozen, Wever 2014 (9) [/]
                            porespace_refr_vol[ii] = porosity_refr[ii]*self.dz[ii] # Available pore space of each layer [m]
                            cold_content[ii]       = CP_I * self.mass[ii] * (T_MELT - self.Tz[ii]) # cold content of each box, i.e. how much heat to bring it to 273K [J]
                            refreeze_mass_pot[ii]  = cold_content[ii] / LF_I # how much mass of the meltwater could be refrozen due to cold content [kg]
                            refreeze_vol_pot[ii]   = refreeze_mass_pot[ii]/1000. # how much meters of the meltwater could be refrozen due to cold content [m]
                            toperc -= p_frozen[ii] # total water toperc is decreased
                            if (self.rho[ii] >= rhoimp): # if we are in an impermeable layer
                                runoff += 1*toperc # runoff the rest of the water
                                toperc = 0.
                        if ((self.rho[ii] < rhoimp) and (irr[ii]*porespace_refr_vol[ii]-self.LWC[ii]) > 0) : # water can be retained in the layer (this happens after possible refreezing)
                            p_retained[ii] = min(irr[ii]*porespace_refr_vol[ii]-self.LWC[ii],toperc)
                            p_retained[ii] = max(p_retained[ii],0)
                            toperc -= p_retained[ii]
                            self.LWC[ii] += p_retained[ii]
                        if self.rho[ii+1] >= rhoimp: # if next layer is impermeable: runoff
                            runoff += 1*toperc
                            toperc = 0.
                        if ii+1 == len(self.dz): # if next layer is end of the column: runoff
                            runoff += 1*toperc
                            toperc = 0.
                        ii += 1 # go to next layer (if no water is left in toperc, while loop will stop)
        if CLparam == 1:
            liqmassprop = 0.057*(porosity/(1-porosity)) + 0.017
            irr = liqmassprop/(1-liqmassprop) * (self.rho*RHO_I) / (RHO_W_KGM*(RHO_I-self.rho))
    
    ### Make water in excess of irreducible water content runoff ###
    irr_limit   = np.maximum(irr*porespace_refr_vol,0.) # max amount of lwc that a layer can hold in its pore space
    rofflayers  = np.where(self.LWC>irr_limit)[0] # layers with lwc in excess of irr
    
    if len(rofflayers)>0: # if there are layers where lwc exceeds irreducible water content
        lwc2 = self.LWC.copy()
        ## Vincent's old code, can be vectorized
        for ll in rofflayers:
            runoff       += self.LWC[ll] - irr_limit[ll] # excess is added to runoff
            self.LWC[ll] = 1*irr_limit[ll] # LWC is reduced 
                  
#         runoff = sum(self.LWC[rofflayers] - irr_limit[rofflayers]) # VV 09/12/2020 fix because runoff should be a scalar       
#         self.LWC[rofflayers] = 1*irr_limit[rofflayers]
        #print('Water above irrlimit runs off, max(self.LWC/(porosity_refr*self.dz)):',max(self.LWC/(porosity_refr*self.dz)))
        
    ### Runoff water present in layers with density exceeding impermeability threshold ###
    if np.any(self.LWC[self.rho>=rhoimp]!=0.):
        runoff += sum(self.LWC[self.rho>=rhoimp])
        self.LWC[self.rho>=rhoimp] = 0.

    ### Percolation + Refreezing of the surface input ###
    try:
       raintoadd = self.rainSec[iii] * S_PER_YEAR * RHO_I_MGM # [m]
       #print('raintoadd is:',raintoadd)
    except:
        raintoadd = 0.

    tofreeze    = 1*melt_vol_a + raintoadd # m water equivalent, the input melt water that has to be refrozen
    frozen      = np.zeros_like(self.dz) # array for refrozen volume in every layer
    retained    = np.zeros_like(self.dz) # array for retained volume in every layer
    
    if tofreeze>0:
        
        if (refreeze_vol_pot[0]>tofreeze) and (porespace_refr_vol[0]>tofreeze): # if all meltwater refreezes in layer[0] -> no need to do entire computation
            frozen[0] += tofreeze
            ipl       = 0
            #print('Full freezing in surface layer')
        elif (refreeze_vol_pot[0]<tofreeze) or (porespace_refr_vol[0]<tofreeze): # otherwise entire computation
            frcap   = np.zeros_like(self.dz) # array for refreezing capacity of the layers
            retcap  = np.zeros_like(self.dz) # array for retention capacity of the layers
            
            if np.any(self.rho>=rhoimp): # check if there is an impermeable layer in the domain
                ipl = np.where(self.rho>=rhoimp)[0][0] # first impermeable layer -> bottom limit for downward percolation
            elif np.all(self.rho<rhoimp): # if no impermeable layer in the domain
                ipl = len(self.dz) # percolation will occur over the entire column potentially
            if percbottom==1:
                if np.any(self.rho<rhoimp): # check if there is an impermeable layer in the domain
                    ipl = np.where(self.rho<rhoimp)[0][-1]+1 # first impermeable layer -> bottom limit for downward percolation
                elif np.all(self.rho>=rhoimp): # if no impermeable layer in the domain
                    ipl = 0 # percolation cannot happen
            frcap[0:ipl] = np.minimum(refreeze_vol_pot[0:ipl],porespace_refr_vol[0:ipl]) # max freezing possible [mWE]
            frcap[self.rho>=rhoimp]  = 0 #avoids refreezing in impermeable layers if we allow water to bypass ice lenses (if percbottom==1)
            rhopot = (self.mass+frcap*1000)/self.dz # potential density if all the refreezing capacity was to be consumed
            porositypot = 1-rhopot/RHO_I # potential porosity
            porosity_refrpot = porositypot*RHO_I/RHO_W_KGM # potential space available for liq water volume once refrozen
            porespace_refr_volpot = porosity_refrpot*self.dz # potential available pore space of each layer [m]
            if CLparam == 1:
                liqmassproppot = 0.057*(porositypot/(1-porositypot)) + 0.017
                irr = liqmassproppot/(1-liqmassproppot) * (rhopot*RHO_I) / (RHO_W_KGM*(RHO_I-rhopot))
            retcap[0:ipl] = irr[0:ipl]*porespace_refr_volpot[0:ipl]-self.LWC[0:ipl] # retention capacity of every layer after refreezing [mWE]
            retcap[self.rho>=rhoimp] = 0 #avoids retention in impermeable layers if we allow water to bypass ice lenses (if percbottom==1)
            if sum(frcap+retcap)>=tofreeze: # enough retention possible to accomodate all liquid water
                ist = np.where(np.cumsum(frcap+retcap)>=tofreeze)[0][0] # layer where all water accomodated -> percolation stops there
                frozen[0:ist]   = np.copy(frcap[0:ist]) # frozen water in layers above ist
                retained[0:ist] = np.copy(retcap[0:ist]) # retained water in layers above ist
                tofreeze        -= sum(frozen+retained) # water left to be frozen/retained in layer ist
                frozen[ist]     = min(tofreeze,frcap[ist]) # freezing in layer ist
                tofreeze        -= frozen[ist] # water left to be retained in layer ist
                retained[ist]   = min(tofreeze,retcap[ist]) # retention in layer ist
                tofreeze        -= retained[ist] # no water should be left in tofreeze
            elif sum(frcap+retcap)<tofreeze: # not possible to accomodate all liquid water -> runoff will occur
                frozen[0:ipl]   = np.copy(frcap[0:ipl]) # frozen water in layers above ipl
                retained[0:ipl] = np.copy(retcap[0:ipl]) # retained water in layers above ipl
                tofreeze -= sum(frozen+retained) # water left for runoff
                runoff += 1*tofreeze # runoff the rest of the water
                tofreeze = 0.
        # print(f'Frozen:  {max(frozen[self.rho>=rhoimp])}')
        # print(f'Retained:  {max(retained[self.rho>=rhoimp])}')
        self.refrozen[0:ipl+1] += frozen[0:ipl+1] # add the freezing to self.refrozen
        self.mass[0:ipl+1] += frozen[0:ipl+1]*1000 # adjust the mass of layers
        self.rho[0:ipl+1] = self.mass[0:ipl+1]/self.dz[0:ipl+1] # adjust density of layers
        latheat = frozen*1000*LF_I # latent heat released due to the refreezing [J]
        cold_content[0:ipl+1] -= latheat[0:ipl+1] # remaining cold content
        self.Tz[0:ipl+1] = T_MELT - cold_content[0:ipl+1]/(CP_I*self.mass[0:ipl+1]) # temperature is changed accordingly
        self.LWC[0:ipl+1] += retained[0:ipl+1] # retained water added to LWC    
        self.LWC[abs(self.LWC)<1e-12] = 0. # adjust for numerical errors
    
    self.runoff = runoff
    self.lwcerror += sum(self.LWC)+sum(self.refrozen)+self.runoff - (melt_volume_WE+sum(initial_lwc)+raintoadd)
    if np.any(self.LWC<0): #VV fixing negative LWC values of the order of 1e-22
        print(self.modeltime[iii])
        print('Inds and Vals LWC<0:',np.where(self.LWC<0)[0],self.LWC[np.where(self.LWC<0)[0]])
        self.LWC = np.maximum(self.LWC,0)
        
    ## Sanity checks
    if abs(sum(self.LWC)+sum(self.refrozen)+self.runoff - (melt_volume_WE+sum(initial_lwc)+raintoadd)) > 1e-12: #check for water balance
        print('Liquid water loss/gain, amount:',sum(self.LWC)+sum(self.refrozen)+self.runoff - (melt_volume_WE+sum(initial_lwc)+raintoadd))
    
    if max(self.Tz>273.16):
        print('Max Tz:',max(self.Tz))
    
    #### Check cold layers are dry ###
    # coldlayers = np.where(self.Tz < T_MELT)
    # if np.any(self.LWC[coldlayers[0]]>0.):
    #     print('Problem: water content in a cold layer')
        
    # return self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, self.dzn, self.LWC, meltgridtrack, self.refrozen, self.runoff, self.lwcerror 

    return self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, self.dzn, self.LWC, meltgridtrack, self.refrozen, self.runoff, self.lwcerror 

def LWC_correct(self):
    '''
    If there is LWC in a layer after temperature diffusion and the temperature
    is less than zero, one option is to just balance the energy to increase the
    temperature and lower the LWC. It isn't the best way to solve the problem 
    but it is one way. 

    This should be vectorized but that is not a priority.
    '''

    ind_wetcold = np.where((self.Tz<T_MELT) & (self.LWC>0))[0]
    if ind_wetcold.size!=0:
        cold_content = CP_I * self.mass * (T_MELT - self.Tz)
        heattofreeze = self.LWC*1000*LF_I
        for kk in ind_wetcold:
            if cold_content[kk] < heattofreeze[kk]:
                # not enough cold content
                # temp needs to be melt
                # some water refreeze to bring T to T_melt

                self.Tz[kk] = T_MELT
                self.LWC[kk] = self.LWC[kk] - (cold_content[kk]/1000/LF_I)
                # self.LWC[kk] = self.LWC[kk] - (cold_content[kk]/1000/LF_I)
            else: #enough cold content, all LWC refreezes
                # Temperature is raised from refreezing
                self.LWC[kk] = 0
                self.Tz[kk] = self.Tz[kk] + heattofreeze[kk]/CP_I/self.mass[kk]
                # self.Tz[kk] = self.Tz[kk] + (heattofreeze[kk]/1000/LF_I)
        if np.any(self.LWC<0):
            print("negative LWC from correction")
            self.LWC[self.LWC<0] = 0
        if np.any(self.Tz > T_MELT):
            print("temps above T_MELT from correction")
            self.Tz[self.Tz>T_MELT] = T_MELT

    return self.Tz, self.LWC

def effectiveT(self,iii):
    '''
    trying to see what happens if we raise the temperature to an 'effective temperature' that is the temperature plus the latent heat from the liquid. Puts the firn above T_melt.
    Potential issue with this method: there might be effective diffusion of mass because we are now diffusing with volumes warmer than T_melt --> adjacent volumes that were dry might end up as having liquid.
    '''
    Q = LF_I * self.LWC * 1000
    deltaT = Q / (self.mass*CP_I)
    Tz_eff = self.Tz + deltaT
    self.Tz = Tz_eff
    T_eff_new, self.T10m = heatDiff(self,iii)
    excessT = np.maximum(0.0,(T_eff_new - T_MELT))
    LWC_new = (excessT * self.mass * CP_I)/ (LF_I * 1000)
    self.LWC = LWC_new
    T_eff_new[self.LWC>0] = T_MELT
    self.Tz = T_eff_new

    return self.Tz, self.T10m