#!/usr/bin/env python
from constants import *
import numpy as np

def percolation_noliquid(self, iii):

    '''
    This is a simple (and very wrong, probably) percolation scheme that does not
    allow any liquid water to exist. Any water that can be refrozed by the 
    cold content of the firn will be. There is no limit of how much of the
    porosity can be filled with that water (which will refreeze). Any water that
    is not refrozed is considered runoff.
    '''
    
    porosity            = 1 - self.rho/RHO_I # porosity
    dz_old              = self.dz
    porespace_vol       = porosity * self.dz #porosity in meters of each box
    melt_volume_IE      = self.snowmeltSec[iii] * S_PER_YEAR #meters
    melt_volume_WE      = melt_volume_IE * RHO_I_MGM #meters
    melt_mass           = melt_volume_WE * 1000. #kg
    heat_to_freeze      = melt_mass * LF_I #amount of heat needed to refreeze the melt (J)

    if (self.mass_sum==melt_mass).any():
        exactmass = True
    else:
        exactmass = False

    ind1a               = np.where(self.mass_sum<=melt_mass)[0] # indicies of boxes that will be melted away
    num_boxes_melted    = len(ind1a)+1 #number of boxes that melt away, include the box that is partially melted
    ind1                = np.where(self.mass_sum>melt_mass)[0][0]   # index which will become the new surface

    ### pm is the partial melt (the box/volume that has a portion melted away)
    pm_mass             = self.mass_sum[ind1] - melt_mass   # the remaining mass of the PM box
    pm_dz               = pm_mass / self.rho[ind1] #remaining thickness
    pm_porespace        = (1 - self.rho[ind1]/RHO_I) * pm_dz #porespace in the PM box
    pm_rho              = self.rho[ind1] #density of the PM box

    cold_content_0      = CP_I * self.mass * (T_MELT - self.Tz) #cold content of each box, i.e. how much heat to bring it to 273K
    cold_content_0_sum  = cold_content_0.cumsum(axis=0)
    cold_content        = cold_content_0[ind1:] #just the boxes that don't melt away
    cold_content[0]     = CP_I * pm_mass * (T_MELT - self.Tz[ind1]) # the partial melt box has its cold content reassigned.
    cold_content_sum    = cold_content.cumsum(axis=0)

    ind2_rel            = np.where(cold_content_sum>heat_to_freeze)[0][0] #freeze horizon index (where the wetting front freezes), index relative to ind1
    ind2                = ind2_rel + ind1 #absolute index on real grid (we have not removed the melted boxes yet)

    if (self.rho[ind1:ind2+1]>830.0).any(): #if there is an ice lens somewhere between the new surface and where freezing should occur
        ind2_rel                = np.where(self.rho[ind1:]>=830.0)[0][0]
        ind2                    = ind2_rel + ind1 #the recalculated freezing front (water can not go past this level)
        cold_content_lens       = cold_content[0:ind2_rel].sum() #cold content that is is available in space between the surface and the lens   
        refreeze_mass           = cold_content_lens / LF_I # this is how much should be able to refreeze in the available pore space above the ice lens.
        melt_volume_WE          = refreeze_mass / RHO_W_KGM
        melt_volume_IE          = melt_volume_WE / RHO_I_MGM #the volume of melt that is not runoff
        runoff_volume_duetolens = (melt_mass - refreeze_mass) * RHO_I_KGM

    else:
        runoff_volume_duetolens = 0.0

    pore_indices        = np.arange(ind1,ind2+1) # indicies of the boxes that are available to fill with water
    pore_indices_flip   = np.flipud(pore_indices)
    porespace_vol[ind1] = pm_porespace
    porespace_0_sum     = porespace_vol.cumsum(axis=0)
    porespace           = porespace_vol[ind1+1:ind2+1] #space available for the water
    porespace_sum       = porespace.cumsum(axis=0)
    porespace_sum_flip  = (np.flipud(porespace)).cumsum(axis=0)
    available_space     = porespace_0_sum[ind2]-porespace_0_sum[ind1]

    if available_space < melt_volume_IE: # melt volume has already been recalculated based on how much can freeze with the cold content

        runoff_volume_duetolimitedporespace = (melt_volume_IE - porespace_sum) * RHO_I_MGM

        self.rho[ind1:ind2+1]   = 830.0 #fill all of the boxes with water.
        self.Tz[ind1:ind2+1]    = T_MELT

        # split up last box into several
        divider         = num_boxes_melted
        self.rho        = np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted)))
        self.age        = np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_melted)))
        self.dz         = np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted)))
        self.dz[0]      = pm_dz
        self.dzn        = np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
        self.dzn        = self.dzn[0:self.compboxes]
        self.Tz         = np.concatenate((self.Tz[ind1:-1],self.Tz[-1]*np.ones(num_boxes_melted)))
        self.bdot_mean  = np.concatenate((self.bdot_mean[ind1:-1],self.bdot_mean[-1]*np.ones(num_boxes_melted)))
        self.z          = self.dz.cumsum(axis = 0)
        self.z          = np.concatenate(([0], self.z[:-1]))
        self.mass       = self.rho*self.dz
    
    elif available_space == 0.0: #the top layer is an ice lens, so the melt runs off

        # split up last box into several
        divider         = num_boxes_melted # ind3 should be removed and replaced with 2 new boxes.
        self.rho        = np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted)))
        self.age        = np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_melted)))
        self.dz         = np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted)))
        self.dz[0]      = pm_dz
        self.dzn        = np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
        self.dzn        = self.dzn[0:self.compboxes]
        self.Tz         = np.concatenate((self.Tz[ind1:-1],self.Tz[-1]*np.ones(num_boxes_melted)))
        self.bdot_mean  = np.concatenate((self.bdot_mean[ind1:-1],self.bdot_mean[-1]*np.ones(num_boxes_melted)))
        self.z          = self.dz.cumsum(axis = 0)
        self.z          = np.concatenate(([0], self.z[:-1]))
        self.mass       = self.rho*self.dz

    else:
        runoff_volume_duetolimitedporespace = 0
        ind3a               = np.where(porespace_sum_flip>melt_volume_IE)[0][0]
        ind3                = ind2 - ind3a #the index of the node that is partially filled with water
        partial_volume      = melt_volume_IE - np.sum(porespace_vol[ind3+1:ind2+1]) # pore space filled in the box that is partially filled
        leftover_porespace  = porespace_vol[ind3]-partial_volume #open pore space in the the partially-filled box

        new_node_1_rho      = self.rho[ind3] #split up the partial box into 2 parts
        new_node_2_rho      = 870.0
        new_node_1_dz       = leftover_porespace / (1 - self.rho[ind3]/RHO_I)
        new_node_2_dz       = self.dz[ind3] - new_node_1_dz

        self.rho[ind3+1:ind2+1] = 870.0
        self.Tz[ind1:ind2+1]    = T_MELT

        # split up last box into several
        divider         = num_boxes_melted # ind3 should be removed and replaced with 2 new boxes.
        self.rho        = np.concatenate((self.rho[ind1:ind3] , [new_node_1_rho,new_node_2_rho] , self.rho[ind3+1:-1] , self.rho[-1]*np.ones(num_boxes_melted-1)))
        self.age        = np.concatenate((self.age[ind1:ind3] , [self.age[ind3],self.age[ind3]] , self.age[ind3+1:-1] , self.age[-1]*np.ones(num_boxes_melted-1)))
        dzhold          = self.dz[ind1+1:ind3]
        dzhold2         = self.dz[ind3+1:-1]
        self.dz         = np.concatenate((self.dz[ind1:ind3] , [new_node_1_dz,new_node_2_dz] , self.dz[ind3+1:-1] ,self.dz[-1]/divider*np.ones(num_boxes_melted-1)))
        self.dzn        = np.concatenate((np.zeros(num_boxes_melted), np.append(dzhold, new_node_1_dz+new_node_2_dz), dzhold2))
        self.dzn        = self.dzn[0:self.compboxes]
        self.dz[0]      = pm_dz
        self.Tz         = np.concatenate((self.Tz[ind1:ind3] , [self.Tz[ind3],self.Tz[ind3]] , self.Tz[ind3+1:-1] , self.Tz[-1]*np.ones(num_boxes_melted-1)))
        self.bdot_mean  = np.concatenate((self.bdot_mean[ind1:ind3] , [self.bdot_mean[ind3],self.bdot_mean[ind3]] , self.bdot_mean[ind3+1:-1] , self.bdot_mean[-1]*np.ones(num_boxes_melted-1)))
        self.z          = self.dz.cumsum(axis = 0)
        self.z          = np.concatenate(([0], self.z[:-1]))
        self.mass       = self.rho*self.dz

    return self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn

def percolation_bucket(self, iii):

    '''
    This is the bucket scheme that allows liquid water to persist in the firn.
    It includes consideration of irreducible liquid water content (LWC) an maximum
    LWC. Water that encounters a slab of a certain density (impermeable_rho) will
    not percolate through.

    LWC is in volume (m^3), and since we are working in one dimension we assume
    that 
    '''

    # maxpore_f                 = 2.0   # factor by which the maximum filled porespace can exceed the irreducible saturation.
    impermeable_rho         = 800.  # impermeable lens density.

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
    self.dzn                = np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
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
    refreeze_mass_pot       = cold_content / LF_I   # how much mass of the meltwater could be refrozen due to cold content (kg)
    refreeze_mass_pot_sum   = refreeze_mass_pot.cumsum(axis=0) 

    ### calculate what the values will be after refreeze happens (pot stands for potential)
    rho_pot                 = (self.mass + refreeze_mass_pot) / self.dz # what the mass of the boxes would be if the refreezemass refroze
    porosity_pot            = 1 - rho_pot / RHO_I
    porespace_vol_pot       = porosity_pot * self.dz
    porespace_air_pot       = porespace_vol_pot - self.LWC

    Wmi                     = 0.057 * (RHO_I - rho_pot) / rho_pot + 0.017 # water per snow-plus- water mass irreducible liquid water content, Langen eqn 3 unitless)
    Swi                     = Wmi / (1 - Wmi) * (rho_pot * RHO_I) / (1000 * (RHO_I - rho_pot))  #irreducible water saturation, volume of water per porespace volume (unitless), Colbeck 1972

    maxpore                 = 0.9 # upper limit on what percentage of the porosity can be filled with water.

    maxLWC1                 = porespace_vol * maxpore   # maximum volume of water that can be stored in each node (meters)
    maxLWC2                 = ((917.0 * self.dz) - self.mass) / RHO_W_KGM # double check that the LWC does not get too large. 
    maxLWC                  = np.minimum(maxLWC1 , maxLWC2)
    maxLWC[self.rho>impermeable_rho] = 0
    maxLWC_mass             = maxLWC * RHO_W_KGM        # mass of the maximum volume of water
    maxLWC1_pot             = porespace_vol_pot * maxpore   # maximum volume of water that can be stored in each node (meters)
    maxLWC2_pot             = ((917.0 * self.dz) - (self.mass + refreeze_mass_pot)) / RHO_W_KGM # double check that the LWC does not get too large. 
    maxLWC_pot              = np.minimum(maxLWC1_pot , maxLWC2_pot)
    # maxLWC_pot[rho_pot>impermeable_rho] = 0
    maxLWC_mass_pot         = maxLWC_pot * RHO_W_KGM        # mass of the maximum volume of water

    irreducible_mass_pot    = Swi * porespace_vol_pot * RHO_W_KGM # mass of irreducible water for each volume (potential - does not separate how much is already there)
    irreducible_vol_pot     = irreducible_mass_pot / RHO_W_KGM
    liquid_storage_vol_pot  = irreducible_vol_pot - self.LWC
    liquid_storage_mass_pot = liquid_storage_vol_pot * RHO_W_KGM

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
        total_liquid_mass       = melt_mass_a + extra_liquid_mass

        ### first, refreeze where possible
        self.mass[0:ind_p]      = self.mass[0:ind_p] + refreeze_mass_pot[0:ind_p]
        self.rho[0:ind_p]       = self.mass[0:ind_p] / self.dz[0:ind_p]
        self.Tz[0:ind_p]        = T_MELT

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
    # print('lwc:',np.sum(self.LWC))
    # print('lwc:',self.LWC)
    

    return self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn, self.LWC

def bucketVV(self, iii):
    
#    tic2=time.time()
    irr = 0.02 # Irreducible water content, this is a proportion of the available pore space
    # We use 830 as threshold to generate runoff, can be changed
    
    ##### First: melting of the surface layers, taken from melt.py #####
    melt_volume_IE      = self.snowmeltSec[iii] * S_PER_YEAR # This still has to be checked by Max (division by self.c['stpsPerYear']?) [m]
    melt_volume_WE      = melt_volume_IE * RHO_I_MGM # [m]
    melt_mass           = melt_volume_WE * 1000. # [kg]
    
    initial_lwc = 1*self.LWC
    
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
    self.dz[0]      = pm_dz #VV dz calculated for the partially melted layer
    self.Dcon       = np.concatenate((self.Dcon[ind1:-1] , self.Dcon[-1]*np.ones(num_boxes_melted)))
    self.dzn        = np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
    self.dzn        = self.dzn[0:self.compboxes]
    self.Tz         = np.concatenate((self.Tz[ind1:-1] , self.Tz[-1]*np.ones(num_boxes_melted)))
    self.r2         = np.concatenate((self.r2[ind1:-1] , self.r2[-1]*np.ones(num_boxes_melted)))
    self.bdot_mean  = np.concatenate((self.bdot_mean[ind1:-1] , self.bdot_mean[-1]*np.ones(num_boxes_melted)))
    self.z          = self.dz.cumsum(axis = 0)
    self.z          = np.concatenate(([0] , self.z[:-1]))
    self.mass       = self.rho * self.dz
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
    lwc_before_freeze = 1*self.LWC
    lwcpres = np.where(self.LWC > 0)[0] # layers with existing lwc
    coldlay = np.where(self.Tz<273.15)[0]
    freezinglayers = np.intersect1d(lwcpres,coldlay)
    latheat = np.zeros_like(self.dz)
    freeze = np.zeros_like(self.dz)
#    print('run time before first for loop =' , time.time()-tic2 , 'seconds')
    for kk in freezinglayers:
        #if self.Tz[kk] < T_MELT: # we can proceed to the refreezing of the lwc if the layer is below melting point
        freeze[kk] = min(refreeze_vol_pot[kk],porespace_refr_vol[kk]) # m water equivalent
        freeze[kk] = min(self.LWC[kk],freeze[kk])
        self.LWC[kk] -= freeze[kk]
        self.refrozen[kk] += freeze[kk]
        self.mass[kk] += freeze[kk]*1000
        self.rho[kk] = self.mass[kk]/self.dz[kk]
        latheat[kk] = freeze[kk]*1000*LF_I
        cold_content[kk] -= latheat[kk]
        self.Tz[kk] = T_MELT - cold_content[kk]/(CP_I*self.mass[kk])
        if ((self.Tz[kk] < T_MELT) and (self.LWC[kk]>0)):
            print('Water remains but cold content still available')
        
        porosity[kk]           = 1 - self.rho[kk]/RHO_I # Definition of porosity [/]
        porespace_vol[kk]      = porosity[kk] * self.dz[kk] # Pore space of each layer [m]
        porosity_refr[kk]      = porosity[kk]*RHO_I/RHO_W_KGM # space available for liq water volume once refrozen, Wever 2014 (9) [/]
        porespace_refr_vol[kk] = porosity_refr[kk]*self.dz[kk] # Available pore space of each layer [m]    
        cold_content[kk]          = CP_I * self.mass[kk] * (T_MELT - self.Tz[kk]) # cold content of each box, i.e. how much heat to bring it to 273K [J]
        refreeze_mass_pot[kk]       = cold_content[kk] / LF_I # how much mass of the meltwater could be refrozen due to cold content [kg]
        refreeze_vol_pot[kk]        = refreeze_mass_pot[kk]/1000. # how much meters of the meltwater could be refrozen due to cold content [m]
        
    ''' The part below is commented out for 2 reasons: very slow and not sure this is part of existing bucket schemes'''    
    ### Then let the non refrozen water of the column percolate ###
    # This is slow! Test not to use it
#    print('run time before second for loop =' , time.time()-tic2 , 'seconds')
#    irr_limit = np.maximum(irr*porespace_refr_vol,0.) # max amount of lwc that a layer can hold in its pore space
#    perclayers = np.where(self.LWC>irr_limit) # layers with lwc in excess of irr
#    lwc_above_irr = self.LWC[perclayers[0]]-irr_limit[perclayers[0]]
#    if len(lwc_above_irr) > 0: #if there are some layers exceeding irreducible water content
#    # Densification (in physics) slightly decreases pore space and thus irr_limit at every time step -> we would have to proceed to entire routine for all layers at irr_limit again and again at every time step
#        if ((max(lwc_above_irr) < 1e-6) and (sum(lwc_above_irr)<1e-3)): # if the lwc above irr_limit is very small (due to slight densification): cheat -> add it to runoff
#            runoff += sum(lwc_above_irr)
#            self.LWC[perclayers[0]] = irr_limit[perclayers[0]]
#        elif ((max(lwc_above_irr) >= 1e-6) or (sum(lwc_above_irr)>=1e-3)):
#            #print('max(lwc_above_irr), sum(lwc_above_irr) are:',max(lwc_above_irr),sum(lwc_above_irr))
#            for index in perclayers[0]:
#                #print('Now we deal with index:',index)
#                ii = index # start from the layer that has lwc in excess of irr
#                toperc = self.LWC[index]-irr_limit[index] # amount of water that has to percolate
#                self.LWC[index] -= toperc # this toperc amount will be lost by the layer
#                p_frozen = np.zeros_like(self.dz) # volume of refrozen water (from toperc) in each layer
#                p_retained = np.zeros_like(self.dz) # volume of retained water (from toperc) in each layer
#                if self.rho[ii] >= 830: #if the layer is impermeable: runoff of all toperc
#                    runoff += 1*toperc
#                    toperc = 0.
#                    
#                while toperc > 0:
#                    if self.Tz[ii] < 273.15: # there can be some refreezing
#                        p_frozen[ii] = min(refreeze_vol_pot[ii],porespace_refr_vol[ii]) # max freezing possible [mWE] 
#                        p_frozen[ii] = min(toperc,p_frozen[ii]) # don't exceed water available
#                        self.refrozen[ii] += p_frozen[ii]
#                        self.mass[ii] += p_frozen[ii]*1000
#                        self.rho[ii] = self.mass[ii]/self.dz[ii]
#                        latheat = p_frozen[ii]*1000*LF_I # latent heat released due to the refreezing [J]
#                        cold_content[ii] -= latheat # remaining cold content
#                        self.Tz[ii] = T_MELT - cold_content[ii]/(CP_I*self.mass[ii]) # temperature is changed accordingly
#                        porosity[ii]           = 1 - self.rho[ii]/RHO_I # Definition of porosity [/]
#                        porespace_vol[ii]      = porosity[ii] * self.dz[ii] # Pore space of each layer [m]
#                        porosity_refr[ii]      = porosity[ii]*RHO_I/RHO_W_KGM # space available for liq water volume once refrozen, Wever 2014 (9) [/]
#                        porespace_refr_vol[ii] = porosity_refr[ii]*self.dz[ii] # Available pore space of each layer [m]
#                        cold_content[ii]       = CP_I * self.mass[ii] * (T_MELT - self.Tz[ii]) # cold content of each box, i.e. how much heat to bring it to 273K [J]
#                        refreeze_mass_pot[ii]  = cold_content[ii] / LF_I # how much mass of the meltwater could be refrozen due to cold content [kg]
#                        refreeze_vol_pot[ii]   = refreeze_mass_pot[ii]/1000. # how much meters of the meltwater could be refrozen due to cold content [m]
#                        toperc -= p_frozen[ii] # total water toperc is decreased
#                        if (self.rho[ii] >= 830): # if we are in an impermeable layer
#                            runoff += 1*toperc # runoff the rest of the water
#                            toperc = 0.
#                    if ((self.rho[ii] < 830) and (irr*porespace_refr_vol[ii]-self.LWC[ii]) > 0) : # water can be retained in the layer (this happens after possible refreezing)
#                        p_retained[ii] = min(irr*porespace_refr_vol[ii]-self.LWC[ii],toperc)
#                        p_retained[ii] = max(p_retained[ii],0)
#                        toperc -= p_retained[ii]
#                        self.LWC[ii] += p_retained[ii]
#                    if self.rho[ii+1] >= 830: # if next layer is impermeable: runoff
#                        runoff += 1*toperc
#                        toperc = 0.
#                    if ii+1 == len(self.dz): # if next layer is end of the column: runoff
#                        runoff += 1*toperc
#                        toperc = 0.
#                    ii += 1 # go to next layer (if no water is left in toperc, while loop will stop)

    ### Percolation + Refreezing of the surface input ###
    tofreeze = 1*melt_vol_a # m water equivalent, the input melt water that has to be refrozen
    jj = 0 # start the refreezing attempt from the surface layer
    #runoff = 0 # Initialised in percolating water loop
    frozen = np.zeros_like(self.dz) # array for refrozen volume in every layer
    retained = np.zeros_like(self.dz) # array for retained volume in every layer
#    print('run time before while loop =' , time.time()-tic2 , 'seconds')
    while tofreeze > 0:
        frozen[jj] = min(refreeze_vol_pot[jj],porespace_refr_vol[jj]) # max freezing possible [mWE]
        frozen[jj] = min(tofreeze,frozen[jj]) # don't exceed water available
        self.refrozen[jj] += frozen[jj]
        self.mass[jj] += frozen[jj]*1000
        self.rho[jj] = self.mass[jj]/self.dz[jj]
        latheat = frozen[jj]*1000*LF_I # latent heat released due to the refreezing [J]
        cold_content[jj] -= latheat # remaining cold content
        self.Tz[jj] = T_MELT - cold_content[jj]/(CP_I*self.mass[jj]) # temperature is changed accordingly
        porosity[jj]           = 1 - self.rho[jj]/RHO_I # Definition of porosity [/]
        porespace_vol[jj]      = porosity[jj] * self.dz[jj] # Pore space of each layer [m]
        porosity_refr[jj]      = porosity[jj]*RHO_I/RHO_W_KGM # space available for liq water volume once refrozen, Wever 2014 (9) [/]
        porespace_refr_vol[jj] = porosity_refr[jj]*self.dz[jj] # Available pore space of each layer [m]
        tofreeze -= frozen[jj]
        #if ((frozen[jj] == porespace_refr_vol[jj]) or (self.rho[jj] >= 830)):
        if (self.rho[jj] >= 830): # if we are in an impermeable layer
            runoff += 1*tofreeze # runoff the rest of the water
            tofreeze = 0.
        #elif frozen[jj] == refreeze_vol_pot[jj]:
        elif self.rho[jj] < 830 : #test
            retained[jj] = min(irr*porespace_refr_vol[jj]-self.LWC[jj],tofreeze)
            retained[jj] = max(retained[jj],0)
            tofreeze -= retained[jj]
            self.LWC[jj] += retained[jj]
        if self.rho[jj+1] >= 830: # if next layer is impermeable: runoff
            runoff += 1*tofreeze
            tofreeze = 0.
        if jj+1 == len(self.dz): # if next layer is end of the column: runoff
            runoff += 1*tofreeze
            tofreeze = 0.
        jj += 1 # go to next layer (if no water is left in tofreeze, while loop will stop)
     
#    print('run time after while loop =' , time.time()-tic2 , 'seconds')
    
    self.runoff = runoff
    self.lwcerror += sum(self.LWC)+sum(self.refrozen)+self.runoff - (melt_volume_WE+sum(initial_lwc))
    
    if abs(sum(self.LWC)+sum(self.refrozen)+self.runoff - (melt_volume_WE+sum(initial_lwc))) > 1e-15: #check for water balance
        print('Liquid water loss/gain, amount:',sum(self.LWC)+sum(self.refrozen)+self.runoff - (melt_volume_WE+sum(initial_lwc)))
    
#    if abs(liquid_bucket + sum(self.LWC) - initial_liquid) > 1e-15:
#        print('Liquid water loss/gain, amount:',liquid_bucket + sum(self.LWC) - initial_liquid)
    
    ### Check cold layers are dry ###
    coldlayers = np.where(self.Tz < T_MELT)
    if np.any(self.LWC[coldlayers[0]]>0.):
        print('Problem: water content in a cold layer')
        
        
    return self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, self.dzn, self.LWC, self.refrozen, self.runoff, self.lwcerror 

