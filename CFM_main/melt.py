#!/usr/bin/env python
from constants import *
import numpy as np
import time
import sys

from diffusion import heatDiff

from darcy_funcs import hydrconducsat_Calonne
from darcy_funcs import vG_Yama_params
from darcy_funcs import phead_vG
from darcy_funcs import krel_vG
from darcy_funcs import thetae_update
from darcy_funcs import thetaeff_equaliser
from darcy_funcs import dfdg_derivative
from darcy_funcs import runoffZuoOerlemans
from darcy_funcs import runoffDarcy
from darcy_funcs import flux_bisection
from darcy_funcs import flux_newtonraphson

'''
Functions to handle meltwater percolation.
'''
#############
def bucket(self,iii):   
    '''
    Percolation bucket scheme, with edits by Max
    Several parameters can be set by the user (see below ### USER CHOICES ###)
    
    Coded by Vincent Verjans

    '''

    ####################
    ### USER CHOICES ###
    try:
        ColeouLesaffre     = self.c['ColeouLesaffre']  # parameterising irreducible water content following Coléou and Lesaffre (1998) formulation [True/False]
        if ColeouLesaffre == False:
            IrrVal         = self.c['IrrVal']   # [%] irreducible water content: proportion of pore space that holds irreducible water
        RhoImp             = self.c['RhoImp']   # density threshold for nodes to be considered as ice lens [kg m-3]
        DownToIce          = self.c['DownToIce']  # allows water to bypass all ice lenses until ice sheet is reached (depth where RhoImp density is definitely reached)
        if DownToIce == False:
            ThickImp       = self.c['ThickImp']    # thickness threshold for ice lens to be impermeable (all ice layers are impermeable if set to 0m) [m] # Using this is slow
        Ponding            = self.c['Ponding']  # allowing LWC ponding above impermeable ice lenses [True/False]
        DirectRunoff       = self.c['DirectRunoff']    # (applicable if Ponding==True) fraction of excess LWC not considered for ponding but running off directly [between 0 and 1]
        RunoffZuoOerlemans = self.c['RunoffZuoOerlemans']  # (applicable if Ponding==True) computing lateral runoff following Zuo and Oerlemans (1996) Eqs.(21,22) [True/False]
        Slope              = self.c['Slope']     # (used only if RunoffZuoOerlemans==True) slope value used in Zuo and Oerlemans (1996) Eq.(22) [/]

    except:
        print('You should add the new melt variables to your .json See melt.py and example.json')
        ColeouLesaffre     = True  # parameterising irreducible water content following Coléou and Lesaffre (1998) formulation [True/False]
        if ColeouLesaffre == False:
            IrrVal         = 0.02   # [%] irreducible water content: proportion of pore space that holds irreducible water
        RhoImp             = 830.   # density threshold for nodes to be considered as ice lens [kg m-3]
        DownToIce          = False  # allows water to bypass all ice lenses until ice sheet is reached (depth where RhoImp density is definitely reached)
        if DownToIce == False:
            ThickImp       = 0.1    # thickness threshold for ice lens to be impermeable (all ice layers are impermeable if set to 0m) [m] # Using this is slow
        Ponding            = False  # allowing LWC ponding above impermeable ice lenses [True/False]
        DirectRunoff       = 0.0    # (applicable if Ponding==True) fraction of excess LWC not considered for ponding but running off directly [between 0 and 1]
        RunoffZuoOerlemans = False  # (applicable if Ponding==True) computing lateral runoff following Zuo and Oerlemans (1996) Eqs.(21,22) [True/False]
        Slope              = 0.1     # (used only if RunoffZuoOerlemans==True) slope value used in Zuo and Oerlemans (1996) Eq.(22) [/]
    ### END USER CHOICES ###
    ########################

    ### Determine mass of melted firn ###
    T_init = self.Tz.copy()

    melt_volume_IE      = self.snowmeltSec[iii]*S_PER_YEAR # [m ie]
    melt_volume_WE      = melt_volume_IE*RHO_I_MGM         # [m] 
    melt_mass           = melt_volume_WE*RHO_W_KGM         # [kg]
    
    ### Define last variables needed for the routine ###
    nnd        = len(self.z)   # number of nodes
    # rhoi       = 917.00001       # avoids numerical errors due to porosity being strictly 0
    rhoi = RHO_I
    runofftot  = 0.            # initialise runoff [m we]
    LWCblocked = np.zeros(nnd) # initialise LWC blocked by impermeable barriers, susceptible to ponding [m]

    if ColeouLesaffre==True:
        IrrVal = 0. # IrrVal is not used in calculations if ColeouLesaffre==True

    self.mass_sum = np.cumsum(self.mass) #cumulative mass [kg]
    
    ### Melting of surface nodes ###
    ind1     = np.where(self.mass_sum>melt_mass)[0][0] #index which will become the new surface
    n_melted = ind1+1 # number of nodes melted

    ### Partially melted node properties ###
    pm_mass = self.mass_sum[ind1]-melt_mass #remaining mass
    pm_dz   = pm_mass/self.rho[ind1] #remaining thickness
    pm_rho  = self.rho[ind1] #density of the pm node
    pm_lwc  = self.LWC[ind1]/self.dz[ind1]*pm_dz #LWC of the pm node
    pm_Tz   = T_MELT

    ### Liquid water input at the surface ###
    liq_in_mass = max(melt_mass + (np.sum(self.LWC[0:ind1+1]) - pm_lwc) * RHO_W_KGM, 0) #avoid negative lwcinput due to numerical round-off errors
    liq_in_vol  = liq_in_mass/RHO_W_KGM

    try: #add rain input if it was provided
       liq_in_vol = liq_in_vol+self.rainSec[iii]*S_PER_YEAR*RHO_I_MGM #[m]
    except:
        pass

    liqmcinit  = pm_lwc+sum(self.LWC[ind1+1:])+liq_in_vol #mass conservation checks

    ### Regridding ###
    if ind1>0:
        self.rho       = np.concatenate((self.rho[ind1:-1],self.rho[-1]*np.ones(n_melted)))
        # self.Tz        = np.concatenate((self.Tz[ind1:-1],self.Tz[-1]*np.ones(n_melted)))
        self.r2        = np.concatenate((self.r2[ind1:-1],self.r2[-1]*np.ones(n_melted)))
        self.bdot_mean = np.concatenate((self.bdot_mean[ind1:-1],self.bdot_mean[-1]*np.ones(n_melted)))
        self.age       = np.concatenate((self.age[ind1:-1],self.age[-1]*np.ones(n_melted))) 
        self.Dcon      = np.concatenate((self.Dcon[ind1:-1],self.Dcon[-1]*np.ones(n_melted)))
        self.dzn       = np.concatenate((np.zeros(n_melted),self.dz[1:]))
        self.dzn       = self.dzn[0:self.compboxes]
    else:
        self.dzn       = self.dz[0:self.compboxes] #VV avoids bug due to undefined self.dzn

    self.LWC       = np.concatenate(([pm_lwc],self.LWC[ind1+1:-1],self.LWC[-1]*np.ones(n_melted)))
    self.dz        = np.concatenate(([pm_dz],self.dz[ind1+1:-1],self.dz[-1]*np.ones(n_melted)))
    self.Tz        = np.concatenate(([pm_Tz],self.Tz[ind1+1:-1],self.Tz[-1]*np.ones(n_melted))) # PM layer should have temp=T_MELT
    self.z         = self.dz.cumsum(axis=0)
    self.z         = np.concatenate(([0],self.z[:-1]))
    self.mass      = self.rho*self.dz
    if self.doublegrid: # if we have doublegrid: need to adjust gridtrack
        meltgridtrack  = np.concatenate((self.gridtrack[ind1:-1],self.gridtrack[-1]*np.ones(n_melted)))
    elif self.doublegrid==False:
        meltgridtrack = np.zeros(nnd) #just return a zero array
    ### end regridding ###

    ### Calculate excessive LWC (above irreducible holding capacity) ###
    phi         = (rhoi - self.rho) / rhoi      # porosity [/]
    phivol      = phi * self.dz                 # pore space [m]
    phivol_av   = phivol * (RHO_I / RHO_W_KGM)  # saturated water content [m], i.e. tot. pot pore space avlbl for any LWC (Eq.9, Wever(2014); Discussion in Yamaguchi(2010))
    ilim        = np.where(self.rho + phivol_av * RHO_W_KGM / self.dz > RHO_I)[0] # nodes where saturation could lead to density>917

    if len(ilim) > 0:     # limit pore space availability for storage in ilim nodes
        phivol_av[ilim] = np.maximum(self.dz[ilim]*(916.99-self.rho[ilim])/RHO_W_KGM,0.)
    
    LWCirr   = IrrVal * phivol_av      # volume of LWC that can be held as irreducible water [m]
    LWCirr[self.rho>=RhoImp] = 0.      # set irreducible water to zero in nodes exceeding impermeability threshold

    if ColeouLesaffre:
        wmi                     = 0.057 * (rhoi - self.rho) / self.rho + 0.017 # irred. water mass per mass of (water+firn) [/] (Coleou and Lesaffre (1998); Eq.3,Langen (2017))
        wmi[self.rho>=RhoImp]   = 0. # set 0 irreducible water in nodes exceeding impermeability threshold
        swi                     = np.zeros_like(wmi)
        imsk                    = self.rho<rhoi
        swi[imsk]               = wmi[imsk]/(1-wmi[imsk]) * rhoi * self.rho[imsk] / (RHO_W_KGM * (rhoi - self.rho[imsk])) # irreducible LWC per porosity space [/] (Eq.4 in Langen (2017))
        LWCirr                  = phivol_av * swi # maximum LWC that can be held as irreducible water [m]
   
    LWC_excess   = np.maximum(self.LWC-LWCirr,0) # LWC in excess of irreducible water content (vector)
    
    ### Calculation of storage capacity in each node ###
    ### storage capacity defined as (refreezing  + irreducible water retention) capacities 
    ### 'excess' LWC means in excess of the irreducible value
    ### 'additional' LWC means additional LWC beyond what is currently in the node, but still less than or equal to the irreducible LWC
    cold_content    = CP_I * self.mass * (T_MELT - self.Tz)           # cold content [J]
    
    refr_cap0       = (cold_content / LF_I) / RHO_W_KGM               # refreezing capacity due to cold content [m we]
    refr_cap0_supp  = np.maximum(0, refr_cap0 - self.LWC)             # refreezing capacity available for any liquid beyond what is currently present [m we]
    refr_cap        = np.minimum(refr_cap0, phivol_av)                # total (existing LWC plus any more) refreezing capacity [m we]
    # (In theory, refr_cap should always be 0 for layers that have any LWC, except the upper most layer.)

    LWC_to_ice      = RHO_W_KGM / rhoi * self.LWC                     # volume that the existing LWC would take if it froze
    refr_vol_supp   = np.maximum(0, phivol_av - LWC_to_ice)           # total volume available for additional LWC [m]
    refr_cap_supp   = np.minimum(refr_cap0_supp, refr_vol_supp)       # refreezing capacity for additional LWC [m we]
    
    rho_pot         = (self.mass + refr_cap * RHO_W_KGM) / self.dz    # potential density after refreezing the maximum possible [kg m-3]
    # rho_pot[rho_pot>=rhoi] = rhoi #avoid rounding errors
    phi_pot         = np.zeros_like(rho_pot)
    con1            = rho_pot<rhoi
    phi_pot[con1]   = (rhoi - rho_pot[con1]) / rhoi                   # potential porosity after refreezing [/]
    phivol_pot      = phi_pot * self.dz                               # potential pore space after refreezing [m]
    phivol_av_pot   = phivol_pot * (RHO_I / RHO_W_KGM)                # total pot. pore space avbl for storage aftr refreezing [m]; Eq.9, Wever(2014); Discussion: Yamaguchi(2010)
    ilim            = np.where(rho_pot + phivol_av_pot * RHO_W_KGM / self.dz > RHO_I)[0] # nodes potentially exceeding 917 density
    
    if len(ilim)>0: # limit pore space availability for storage in ilim nodes
        phivol_av_pot[ilim] = np.maximum(self.dz[ilim] * (916.99 - rho_pot[ilim]) / RHO_W_KGM,0.)
    
    ######
    LWCirr_pot    = IrrVal * phivol_av_pot # LWC that can be held as irreducible water after refreezing occurs[m]
    ######

    if ColeouLesaffre:
        wmi_pot                      = 0.057 * (rhoi - rho_pot) / rho_pot + 0.017 # irred. water mass per mass of (water+firn) [/] (Coleou and Lesaffre (1998); Eq.3,Langen (2017))
        wmi_pot[rho_pot >= RhoImp]   = 0. # set 0 irreducible water in nodes exceeding impermeability threshold
        swi_pot                      = np.zeros_like(wmi_pot)
        imskp                        = rho_pot < rhoi
        swi_pot[imskp]               = wmi_pot[imskp] / (1 - wmi_pot[imskp]) * rhoi * rho_pot[imskp] / (RHO_W_KGM * (rhoi - rho_pot[imskp])) # irreducible LWC per porosity space [/] (Eq.4 in Langen (2017))
        LWCirr_pot                   = phivol_av_pot * swi_pot # maximum LWC that can be held as irreducible water [m]
    
    LWCirr_pot[rho_pot >= RhoImp] = 0. # set 0 irreducible water in nodes exceeding impermeability threshold
    LWCunf                        = np.maximum(0,self.LWC - refr_cap) # unfrozen LWC that will remain in each node after refreeze [m]
    retcap_supp                   = np.maximum(0,LWCirr_pot-LWCunf) # retention capacity for additional LWC (assuming refreezing occurs first), potential minus how much is there already.
    
    # Define storage capacity #
    stcap = refr_cap_supp + retcap_supp #total storage capacity of each node for additional LWC [m]

    ### Ice lens algorithm: find ice lenses satisfying density and thickness criteria ###
    if DownToIce: # only ice sheet nodes are considered impermeable
        if np.any(self.rho < RhoImp): # there is at least one node below density threshold
            imp = np.arange(np.where(self.rho < RhoImp)[0][-1]+1,nnd,1).astype(int) # the nodes below last rho<RhoImp are impermeable
        else: # all nodes above impermeabilty threshold
            imp = np.arange(0,nnd,1).astype(int) # all nodes are impermeable
    
    elif DownToIce==False:
        if ThickImp > 0:
            lens0 = np.array([ii for ii in range(1,nnd) if (self.rho[ii]>=RhoImp and self.rho[ii-1]<RhoImp)]) # top index of each ice lens
            lens1 = np.array([ii for ii in range(0,nnd-1) if (self.rho[ii]>=RhoImp and self.rho[ii+1]<RhoImp)]) # bottom index of each ice lens
            
            if self.rho[0] >= RhoImp: # if surface node is an ice lens
                lens0 = np.append(0,lens0).astype(int) # add to list
            
            lens1 = np.append(lens1,nnd-1).astype(int) # bottom node is always end of the bottom ice lens
            imp   = np.array([]) # prepare vector of impermeable nodes
            for ii in range(len(lens0)):
                lensdz = (sum(self.dz[lens0[ii]:lens1[ii]+1])) # thickness of the ice lens
                if lensdz>=ThickImp or (ii==len(lens0)-1): # impermeability if thickimp reached (bottom of domain is always impermeable)
                    imp = np.append(imp,np.arange(lens0[ii],lens1[ii]+1,1)).astype(int) #impermeable nodes
        else:
            imp = np.where(self.rho>=RhoImp)[0] # all nodes exceeding RhoImp are considered impermeable

    imp = imp.astype(int)

    if len(imp) == 0:
        imp = [(len(self.rho) - 1)] # MS addition: If domain does not extend to full ice density, this will ensure the melt routine works (hack solution?)
        # might create issues if ponding allowed?

    stcap[imp] = 0. # set 0 storage capacity for impermeable nodes
    stcap_cum  = np.cumsum(stcap) # cumulative storage capacity, refreezing + irreducible

    ### Store surface melt according to stcap of nodes from surface to bottom ###
    LWCblocked  = np.zeros(nnd)

    if liq_in_vol > 0: # there is some liquid water input from the surface
        if stcap_cum[-1] >= liq_in_vol: # there is enough storage capacity for all liquid input
            ii0 = np.where(stcap_cum >= liq_in_vol)[0][0] # bottom most node storing surface melt
        else: # not enough storage capacity for liq input
            ii0 = nnd - 1 # set ii0 to bottom node

        if ii0 >= imp[0]: # impermeable barrier or not enough pore space prevents full distribution of meltinput
            ii0             = max(0,imp[0]-1) # ii0 limited to node above impermeable barrier
            storageinp      = np.concatenate((stcap[0:ii0+1],np.zeros(nnd-ii0-1))) # each node above the barrier gets filled with its stcap
            LWCblocked[ii0] = liq_in_vol - sum(storageinp) # volume of water that is excess due to blockage
        else: # no imperbeamble barrier and there is adequate pore space: meltinput is distributed according to storage capacity
            if ii0 == 0: # all water input stored in surface node
                storageinp  = np.concatenate(([liq_in_vol],np.zeros(nnd-1)))
            else: # water input stored in several nodes
                storageinp  = np.concatenate((stcap[0:ii0],[liq_in_vol-stcap_cum[ii0-1]],np.zeros(nnd-ii0-1)))

    elif liq_in_vol == 0: #no liquid water input
        storageinp = np.zeros(nnd) #no input water storage

    stcap1 = stcap - storageinp #update storage capcity

    ### Set LWC_excess in impermeable nodes as blocked LWC ###
    indsblc             = np.intersect1d(np.where(LWC_excess > 0)[0],imp) # imp nodes with some LWC_excess
    LWCblocked[indsblc] = LWCblocked[indsblc] + LWC_excess[indsblc]       # LWC_excess of insdblc assumed blocked
    self.LWC[indsblc]   = self.LWC[indsblc] - LWC_excess[indsblc]           # update LWC
    LWC_excess[indsblc] = 0.                                          # update LWC_excess
    
    ### Distribute LWC_excess in the nodes supporting storage and/or in LWCblocked ###
    LWC1     = np.copy(self.LWC) # LWC will be modified by LWC_excess transfers
    storage1 = np.zeros(nnd)     # LWC stored in the different nodes

    if np.any(LWC_excess)>0: #if there is some excess LWC
        tostore = 0                     # LWC stock that must be stored
        indsexc = np.where(LWC_excess>0)[0] # indices of nodes with excess LWC
        indb1   = indsexc[-1]           # bottom most node where LWC_excess exists
        jj0     = indsexc[0]            # start from most upper node with excess LWC
        if np.any(stcap1 >0):           # Max moved indb2 definition to here, before start of if 
            indb2 = np.where(stcap1 > 0)[0][-1] #bottom most node where LWC_excess can be stored
        else:
            indb2 = 0

        # if np.any(stcap1>0): # there is some storage capacity in the firn column (This is VV original)
        if ((np.any(stcap1[1:]>0)) and (indb2>jj0)): #there is some storage capacity in the firn column, and it is deeper than jj0

            while ((jj0 <= indb1) or (tostore > 0)):
                if (np.where(stcap1[jj0:]>0)[0]).size > 0:    
                    jj1 = jj0+np.where(stcap1[jj0:]>0)[0][0] # next node that can store some of the LWC_excess
                else: # all nodes with positive stcap1 is shallower than jj0, emulate the no storage capacity routine below
                    # though there might be stcap1 in shallower - could route water that way?
                    indsexc_deep = indsexc[indsexc>=jj0]
                    for jj2 in indsexc_deep: # find underlying impermeable barrier for each node with some LWC_excess
                        jj1             = imp[np.where(imp>=jj2)[0][0]]-1 # jj1 becomes index of node above the impermeable barrier
                        LWCblocked[jj1] += LWC_excess[jj2]          # LWC_excess is blocked above the barrier
                        LWC1[jj2]       = LWCirr[jj2]               # LWC of jj0 is reduced to irreducible water content
                    break # Exit the while loop

                if imp[np.where(imp >= jj0)[0][0]] > jj1:       # jj0 and jj1 nodes not separated by an impermeable barrier
                    tostore         += sum(LWC_excess[jj0:jj1+1])   # all LWC_excess from jj0 to jj1 are subject to storage
                    LWC1[jj0:jj1+1] = np.minimum(LWC1[jj0:jj1+1],LWCirr[jj0:jj1+1]) # LWC_excess is evacuated
                    storage1[jj1]   = min(stcap1[jj1],tostore)  # jj1 node stores as much as possible
                    tostore         -= storage1[jj1]            # tostore is reduced, jj1 is filled
                    jj0             = jj1+1                     # go to next node with possible storage capacity
                    if jj0 >= indb2:                            # no possible storage of LWC_excess 
                        jj1 = imp[np.where(imp>=jj0)[0][0]] - 1 # find the next impermeable barrier
                        LWCblocked[jj1] += tostore              # all LWC to be stored is blocked above the barrier
                        tostore = 0.                            # tostore is set to 0

                else: # impermeable barrier between jj0 and jj1
                    jj1             = imp[np.where(imp>=jj0)[0][0]]-1   # jj1 becomes index of node above the impermeable barrier
                    tostore         += sum(LWC_excess[jj0:jj1+1])       # all LWC_excess from jj0 to jj1 are subject to be blocked above the barrier
                    LWC1[jj0:jj1+1] = np.minimum(LWC1[jj0:jj1+1],LWCirr[jj0:jj1+1]) # LWC_excess is evacuated
                    LWCblocked[jj1] += tostore                        # all LWC to be stored is blocked above the barrier
                    tostore         = 0.                              # tostore is set to 0
                    
                    if jj1 < indb1: # still nodes with LWC_excess to be treated
                        jj0 = indsexc[np.where(indsexc>jj1)[0][0]]  # go to next node with LWC_excess>0
                    else: # all nodes with LWC_excess have been treated
                        jj0 = indb1+1                               # terminate the while loop
        
        else: # no storage capacity in the firn column
            for jj0 in indsexc: # find underlying impermeable barrier for each node with some LWC_excess
                jj1 = imp[np.where(imp>=jj0)[0][0]]-1   # jj1 becomes index of node above the impermeable barrier
                LWCblocked[jj1] += LWC_excess[jj0]      # LWC_excess is blocked above the barrier
                LWC1[jj0] = LWCirr[jj0]                 # LWC of jj0 is reduced to irreducible water content
                
    storagetot = storageinp+storage1    # total storage in each node
    LWC1       = LWC1+storagetot        # redistributed LWC

    ### Refreezing ###
    freeze      = np.minimum(LWC1,refr_cap)     # refreezing in each individual node [m we]
    self.mass   = self.mass + RHO_W_KGM*freeze  # update mass [kg]
    self.LWC    = LWC1 - freeze                 # update LWC
    self.rho    = self.mass/self.dz             # update density [kg m-3]
    latheat     = freeze*RHO_W_KGM*LF_I         # latent heat released due to the refreezing [J]
    cold_content    -= latheat                  # remaining cold content [J]   
    refrozentot = sum(freeze)                   # total refrozen water [m we]
    self.Tz[freeze>0] = T_MELT - cold_content[freeze>0]/(CP_I*self.mass[freeze>0]) # update Tz [K]
    
    ### Store LWC blocked ###
    runofftot  = runofftot + DirectRunoff*np.sum(LWCblocked) #Direct runoff of part of the blocked LWC (user choice)
    LWCblocked = (1 - DirectRunoff)*LWCblocked #corresponding decrease of LWCblocked
    if np.any(LWCblocked > 0):
        if Ponding == True: #ponding is allowed
            LWCold = self.LWC.copy()
            rhofinal        = self.rho
            phiempty        = self.dz * (rhoi - rhofinal) / RHO_W_KGM - self.LWC # updte tot pot pore space avbl for LWC [m] (Eq.9; Wever(2014); Discussion in Yamaguchi(2010))
            phiempty[imp] = 0. # set 0 LWC ponding in impermeable nodes
                
            for kk in np.flip(np.where(LWCblocked > 0)[0]): 
                phiempty_cumf = np.cumsum(np.flip(phiempty[0:kk+1])) #cumulative empty porespace until kk included, flipped
                
                if phiempty_cumf[-1] >= LWCblocked[kk]: # enough porosity to accomodate ponding LWC
                    ifill = np.where(phiempty_cumf > LWCblocked[kk])[0][0] # [kk-ifill] is most upper node that accomodates LWCblocked[kk]                 
                else:
                    ifill = kk # ponding until surface node
                    runofftot = runofftot + LWCblocked[kk] - phiempty_cumf[kk] # remove LWC that cannot be accomodated as runoff
                    LWCblocked[kk] = phiempty_cumf[kk] # MS added: need to also remove that volume from LWCblocked

                if ifill == 0: # The excess water can be contained in the kk node 
                    self.LWC[kk] = self.LWC[kk]+LWCblocked[kk] # update LWC
                    phiempty[kk] = phiempty[kk]-LWCblocked[kk] # update phiempty
                else:
                    LWCfinal                    = self.LWC
                    LWCfinal[kk-ifill+1:kk+1]   = LWCfinal[kk-ifill+1:kk+1] + phiempty[kk-ifill+1:kk+1] # fill nodes from kk-ifill (not included)
                    LWCblocked[kk]              = LWCblocked[kk] - np.sum(phiempty[kk-ifill+1:kk+1])  # remaining LWC in LWCblocked[kk]
                    phiempty[kk-ifill+1:kk+1]   = 0. # update phiempty
                    self.LWC[kk-ifill]          = self.LWC[kk-ifill] + LWCblocked[kk] # node[kk-ifill] accomodates remaining of LWCblocked[kk]
                    phiempty[kk-ifill]          = phiempty[kk-ifill] - LWCblocked[kk] # update phiempty
                LWCblocked[kk]                  = 0. # LWCblocked[kk] has been accomodated
        
        elif Ponding == False: #no ponding
            runofftot   = runofftot+np.sum(LWCblocked) #set all LWCblocked as runoff
            LWCblocked  = 0*LWCblocked #LWCblocked is empty
   
    ### Zuo and Oerlemans (1996) runoff routine ###
    if RunoffZuoOerlemans == True: # Calculations with post-refreezing values       
        phi         = (rhoi - self.rho) / rhoi #porosity [/]
        phivol      = phi * self.dz #pore space [m]
        phivol_av   = phivol * (RHO_I / RHO_W_KGM) #total potential pore space available for refreezing [m] (Eq.9 in Wever (2014) and Discussion in Yamaguchi (2010))
        LWCirr      = IrrVal * phivol_av #maximum LWC that can be held as irreducible water [m]
        
        if ColeouLesaffre:
            wmi                   = 0.057 * (rhoi - self.rho) / self.rho + 0.017 # irred. wtr mass per mass of (water+firn) [/] (Coleou and Lesaffre (1998); Eq.3 in Langen (2017))
            wmi[self.rho>=RhoImp] = 0. # set 0 irreducible water in nodes exceeding impermeability threshold
            swi                   = np.zeros_like(wmi)
            imsk                  = self.rho<rhoi
            swi[imsk]                  = wmi[imsk] / (1 - wmi[imsk]) * rhoi * self.rho[imsk] / (RHO_W_KGM * (rhoi - self.rho[imsk])) # irreducible LWC per porosity space [/] (Eq.4, Langen(2017))
            LWCirr                = phivol_av * swi # maximum LWC that can be held as irreducible water [m]
        
        LWCirr[self.rho >= RhoImp] = 0.             # set 0 irreducible water in nodes exceeding impermeability threshold
        LWC_rfZO = np.maximum(0,self.LWC-LWCirr)    # LWC subject to Zuo and Oerlemans runoff [m]
        
        if np.any(LWC_rfZO > 0):
            indsrfZO        = np.where(LWC_rfZO > 0)[0] # nodes subject to Zuo and Oerlemans runoff
            c1zuo           = 1.5*24*3600               # constant from Zuo and Oerlemans (1996), converted in [s]
            c2zuo           = 25.*24*3600               # constant from Zuo and Oerlemans (1996), converted in [s]
            c3zuo           = 140.                      # constant from Zuo and Oerlemans (1996) [/]
            tstar           = c1zuo + c2zuo * np.exp(-1 * c3zuo * Slope) # Eq.(22) Zuo and Oerlemans 1996 [s]
            rfZO            = np.zeros(nnd)             # initialise runoff Zuo and Oerlemans
            rfZO[indsrfZO]  = self.dt[iii] * LWC_rfZO[indsrfZO] / tstar # from Eq.(21) Zuo and Oerlemans 1996 [m]
            self.LWC        = self.LWC - rfZO           # decrease LWC
            runofftot       = runofftot + np.sum(rfZO)  # add the calculated runoff to the total runoff
            
    ### Mass conservation check ###
    liqmcfinal = sum(self.LWC) + refrozentot + runofftot
    if abs(liqmcfinal - liqmcinit) > 1e-3:
        print(f'Mass conservation error (melt.py) at step {iii}\n    Init: {liqmcinit} m\n    Final: {liqmcfinal} m')
    
    ### Dry cold firn check ###
    coldlayers = np.where(self.Tz < T_MELT)[0]
    if np.all(self.LWC[coldlayers] < 1e-9):
        self.LWC[coldlayers] = 0.
    if np.any(self.LWC[coldlayers] > 0.):
        print('Problem: water content in a cold layer')

    return self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, self.dzn, self.LWC, meltgridtrack, refrozentot, runofftot

#############

##########################
def darcyscheme(self,iii):
    '''
    Threshold thickness for ice lenses to be impermeable
    Approach: modify density for calculation of Darcy variables (upper threshold set to 910)
    Liquid water input is distributed in upper nodes rather than using a flux
    boundary condition at the surface node (causes high runoff)
    Input is accomodated until impermeable node reached: remaining input runs off
    Melting and refreezing occur at each Darcy step
    '''

    ticdarcy = time.time()
    timetot  = self.dt[iii] #total duration to be covered by the Darcy routine
    dtsub    = 60 # [s] duration of Darcy time steps, adjusted iteratively
    ### User choices ###
    dtmin    = 60 # [s] minimal time step for the Darcy routine
    dtmax    = 3600 # [s] maximal time step for the Darcy routine
    dtsub    = 60. #starting time step used in the Darcy scheme (default value: 60 sec)
    RhoImp     = 873. #density from which a firn layer is considered impermeable to incoming flux
    ThickImp   = 0.5 #minimum thickness of ice lens to be considered impermeable [m]
    lat_runoff = True #to compute lateral runoff above impermeable ice layers
    slope_proxy = 0.02 #proxy for the slope (used to compute the lateral runoff fluxes) [dimensionless]
    rholim_dcy  = 910. #max density used in Darcy calculations to allow for some permeability of ice lenses [kg m-3]
    
    ### Surface fluxes ###
    melt_mass_tot = self.snowmeltSec[iii]*S_PER_YEAR*917 #total melt over the Darcy routine [kg]
    meltflux_mass = melt_mass_tot/timetot #melt flux throughout the Darcy routine [kg s-1]
    try:
        rain_vol_tot = self.rainSec[iii]*S_PER_YEAR*0.917 #total rain [m we]
    except:
        rain_vol_tot = 0.
    rainflux_vol = rain_vol_tot/timetot #rain flux throughout the Darcy routine [mwe s-1]

    eps_cvg  = 0.1e-3 #convergence criterion for equilibrium head when solving for qlim of Hirashima et al. (2010) [m]

    ### Firn variables ###
    ncv            = len(self.z) #number of nodes (number of control volumes)
    initial_lwc    = np.copy(self.LWC) #LWC before Darcy scheme
    phi            = (917-self.rho)/917 #update porosity
    phi            = np.maximum(0,phi) #avoid numerical errors of very small negative phi
    cp_i           = 152.5+7.122*self.Tz #specific heat of ice [J kg-1 K-1] Cuffey and Paterson 2010 (9.1)
    rg             = np.sqrt(self.r2) #grain radius [m]
    dltz           = np.append(self.dz[0:-1]/2+self.dz[1:]/2,self.dz[-1]/2) #distance between centres of nodes
    runofftot      = 0. #total runoff over the entire Darcy routine
    refr_tot        = 0. #total refreezing over the entire Darcy routine
    if self.doublegrid: #if we have doublegrid: need to adjust gridtrack
        meltgridtrack  = np.copy(self.gridtrack) #prepare gridtrack adjusted for melting
    elif self.doublegrid==False:
        meltgridtrack = np.zeros_like(self.dz) # just return a zero array

    timer     = 0 #timer of the Darcy routine
    rho_lens0 = np.copy(self.rho) #initialise rho used for lenses at the Darcy time step
    glwflux_d2 = np.zeros(ncv-1) #glw flux computed at Darcy step -2
    glwflux_d1 = np.zeros(ncv-1) #glw flux computed at Darcy step -1
    # if self.rho[-1]<RhoImp:
    #     print('Bottom domain does not reach impermable density in Darcy scheme, exiting')
    #     sys.exit()
    if RhoImp>rholim_dcy:
        print('RhoImp must be below rholim_dcy in Darcy scheme, exiting')
        sys.exit()
    while timer<timetot:
        ### Melt upper nodes ###
        if meltflux_mass>0:
            meltstep_mass  = meltflux_mass*dtsub #melt mass over dtsub time interval [kg]
            self.mass_sum  = np.cumsum(self.mass) #depth-cumulated mass
            ind1           = np.where(self.mass_sum>=meltstep_mass)[0][0] #bottom most node affected by melt
            pm_mass        = self.mass_sum[ind1]-meltstep_mass #mass of partially melted node
            pm_dz          = pm_mass/self.rho[ind1] #thickness of partially melted node
            pm_lwc         = self.LWC[ind1]/self.dz[ind1] * pm_dz #lwc of partially melted node
            lwc_p          = sum(self.LWC[0:ind1+1])-pm_lwc #lwc of the melted part of the firn contributing to percolation
            n_mlt          = ind1+1 #number of nodes melted, including the partially melted node
            self.dz        = np.concatenate(([pm_dz],self.dz[ind1+1:-1],self.dz[-1]*np.ones(n_mlt))) #update dz
            self.dzn       = np.concatenate((np.zeros(n_mlt),self.dz[1:])) #taken from the code of Max
            self.dzn       = self.dzn[0:self.compboxes] #taken from the code of Max
            self.z         = np.cumsum(self.dz) #update z
            self.z         = np.append(0,self.z[0:-1]) #update z
            dltz           = np.append(self.dz[0:-1]/2+self.dz[1:]/2,self.dz[-1]/2) #distance between centres of nodes
            self.rho       = np.append(self.rho[ind1:-1],self.rho[-1]*np.ones(n_mlt)) #update rho
            phi            = (917-self.rho)/917 #update porosity
            phi            = np.maximum(0,phi) #avoid numerical errors of very small negative phi
            self.mass      = self.dz*self.rho #update mass
            self.age       = np.append(self.age[ind1:-1],self.age[-1]*np.ones(n_mlt)) #update age
            self.LWC       = np.concatenate(([pm_lwc],self.LWC[ind1+1:-1],self.LWC[-1]*np.ones(n_mlt))) #update LWC
            self.Tz        = np.append(self.Tz[ind1:-1],self.Tz[-1]*np.ones(n_mlt)) #update temperature
            cp_i           = 152.5+7.122*self.Tz #specific heat of ice [J kg-1 K-1] Cuffey and Paterson 2010 (9.1)
            self.r2        = np.append(self.r2[ind1:-1],self.r2[-1]*np.ones(n_mlt)) #update squared grain radius
            rg             = np.sqrt(self.r2) #grain radius [m]
            if self.doublegrid: #if we have doublegrid: need to adjust gridtrack
                meltgridtrack  = np.concatenate((meltgridtrack[ind1:-1],meltgridtrack[-1]*np.ones(n_mlt)))

            # Liquid water input at surface node #
            liq_input  = meltstep_mass/1000+lwc_p+rainflux_vol*dtsub #[m we]
        
        else: #only rain as liquid water input
            n_mlt    = 0 #no melted node
            self.dzn = self.dz[0:self.compboxes] #taken from the code of Max
            liq_input   = rainflux_vol*dtsub #[m we]
            meltstep_mass = 0.
         
        ### Spot ice lenses and evaluate their thickness ###
        if timer==0 or (n_mlt>1 or self.rho[0]>=RhoImp or len(rho_lens0[rho_lens0>=RhoImp])!=len(self.rho[self.rho>=RhoImp])):
            # Ice lens algorithm only at time step 0 or if the lens distribution has changed #
            lens0 = np.array([ii for ii in range(1,ncv) if (self.rho[ii]>=RhoImp and self.rho[ii-1]<RhoImp)]) #top index of each ice lens
            lens1 = np.array([ii for ii in range(0,ncv-1) if (self.rho[ii]>=RhoImp and self.rho[ii+1]<RhoImp)]) #bottom index of each ice lens
            if self.rho[0]>=RhoImp: #if surface node is an ice lens
                lens0 = np.append(0,lens0).astype(int) #add to list
            lens1 = np.append(lens1,ncv-1).astype(int) #bottom node is always end of the bottom ice lens
            imp    = np.array([]) #prepare vector of impermeable nodes
            imptop = np.array([]) #prepare vector of nodes at the top of the impermeable ice lenses
            for ii in range(len(lens0)):
                lensdz = (sum(self.dz[lens0[ii]:lens1[ii]+1])) #thickness of the ice lens
                if lensdz>ThickImp or (ii==len(lens0)-1): #impermeability if ThickImp reached (bottom of domain is always impermeable)
                    imp = np.append(imp,np.arange(lens0[ii],lens1[ii]+1,1)).astype(int) #impermeable nodes
                    imptop = np.append(imptop,lens0[ii]).astype(int) #surface node of the impermeable ice lens
            
            if imp.size==0:
                imp = np.array([len(self.rho)-1])
                imptop = imp.copy()
            # print('imp',imp)
            imp     = imp.astype(int)
        rho_lens0 = np.copy(self.rho) #rho used for lenses at this Darcy time step

        ### Define all Darcy-scheme variables ###
        rhodcy  = np.minimum(self.rho,rholim_dcy) #put limit on density to allow percolation through ice lenses
        rhomod  = np.where(self.rho!=rhodcy)[0] #all nodes of which density is modified
        phidcy  = (917-rhodcy)/917 #porosity with modified rhodcy
        theta_w   = self.LWC/self.dz #volumetric water content
        theta_s   = phidcy*0.917 #theta at saturation, Yamaguchi 2010 obs of 10% pore space filled with air
        theta_s   = theta_s-1e-6 #numerical adjustment to avoid rho slightly above 917
        theta_i   = np.minimum(0.02,(theta_s-1e-9)) #irreducible water content, Yamaguchi 2010, limited to availabe porosity
        theta_i[self.rho>=RhoImp] = 0 #set 0 irreducible water content for ice lenses    
        theta_s[theta_s<1e-6] = 1e-9 #0 porosity of ice layers can cause numerical problems
        LWC_i     = theta_i*self.dz #corresponding irreducible water amount
    
        LWCav     = np.maximum(0,self.LWC-LWC_i) #LWC available for water flow
        LWCacm    = self.dz*theta_s-self.LWC #extra LWC that can be accomodated in each layer
        LWCacm[imp] = 0. #do not allow for water storage in impermeable ice lenses
        i_inp     = np.where(np.cumsum(LWCacm)>=liq_input)[0][0] #lowest node that receives liq water input
        if i_inp>=imptop[0]: #liq water input depth reaches an impermeable lens
            i_inp = imptop[0]-1 #do not allow for input water at and below the impermeable ice lens
        if i_inp==-1: #impermeable lens at the surface
            runoff0 = 1*liq_input #all liquid water input runs off
        elif i_inp==0: #only surface node receives input
            self.LWC[0] = self.LWC[0]+min(liq_input,LWCacm[0]) #input can be limited by LWCacm
            runoff0 = liq_input-min(liq_input,LWCacm[0]) #input not accomodated by the surface node runs off
        elif i_inp>0:
            liq_input1           = liq_input-np.cumsum(LWCacm)[i_inp-1] #liq water input remaining after filling nodes[0:i_inp]
            self.LWC[0:i_inp] = self.LWC[0:i_inp]+LWCacm[0:i_inp] #update LWC of nodes[0:i_inp]
            self.LWC[i_inp]   = self.LWC[i_inp]+min(liq_input1,LWCacm[i_inp]) #node[i_inp] receives remaining liq_input1
            runoff0           = liq_input1-min(liq_input1,LWCacm[i_inp])
        theta_w   = self.LWC/self.dz #volumetric water content
        LWCav     = np.maximum(0,self.LWC-LWC_i) #LWC available for water flow
        LWCacm    = self.dz*theta_s-self.LWC #extra LWC that can be accomodated in each layer
        LWCacm[imp] = 0. #do not allow for water storage in impermeable ice lenses
        
        theta_e   = (theta_w-theta_i)/(theta_s-theta_i) #effective water saturatin, Hirashima 2010 (5)
        stab_e    = 1e-9 #stabilisation theta_e
        theta_e   = np.maximum(stab_e,theta_e) #avoid non-positive effective saturation
        theta_e   = np.minimum(1-stab_e,theta_e) #avoid effective saturation equal to 1
        avG,nvG,mvG = vG_Yama_params(rg,rhodcy) #van Genuchten parameters
        bigk_s      = hydrconducsat_Calonne(rg,rhodcy) #hydraulic conductivity at saturation [m s-1]
        bigk_r      = krel_vG(mvG,theta_e) #relative hydraulic conductivity
        bigk        = bigk_r*bigk_s #hydraulic conductivity Hirashima 2010 (11) [m s-1]
        bigk_d      = np.append(bigk[1:],0) #hydraulic conductivity staggered down
        hd          = phead_vG(avG,nvG,mvG,theta_e) #pressure head [m]
        dlthd       = np.append(np.diff(hd),0) #difference between node pressure head and underlying node head
        dhdz        = dlthd/dltz #ratio delta(head)/delta(z) (>0 implies that absolute head increases downwards)
        # Ratio dhdz determines the conductivity taken at lower interface #
        bigkedg   = np.zeros(ncv)
        # Note: in Szymkiewicz 2009, potential head is taken negative (what matters is sign of q in (1))
        bigkedg[dhdz+1>=0] = bigk[dhdz+1>=0] #Szymkiewicz 2009, Eq. (3c) downward gradient in total head (1b)
        bigkedg[dhdz+1<0]  = bigk_d[dhdz+1<0] #Szymkiewicz 2009, Eq. (3c) upward gradient in total head (1b)

        ### Determine flow at interfaces of volumes following Hirashima et al. (2010) method ###
        ## Detect all nodes that cannot have outflow ##
        indsdry = np.where(theta_e<=1e-3)[0] #dry nodes
        if np.any(theta_e>1e-3):
            i0dry   = np.where(theta_e>1e-3)[0][-1]+1 #node below which all nodes are dry
        else:
            i0dry = 1 #only surface influx in the surface node
        imp_d    = imp-1 #nodes with an impermeable volume below
        ## Find all nodes above an impermeable boundary that are fully saturated ##
        satnofl = [] #continuum of saturated nodes above an impermeable boundary
        for jj1 in imptop: #check for all impermeable boundaries
            if jj1<=i0dry: #calculation is only necessary in part of the domain where water is percolating
                jj1u = jj1-1 #upperlying node
                if jj1u>=0: #check if surface node is reached
                    srf = False
                elif jj1u<0: #check if surface node is reached
                    srf = True
                # Detect continuum of saturated nodes above jj1 #
                while (srf==False and theta_e[jj1u]>=0.95):
                    satnofl.append(jj1u) #saturated
                    jj1u = jj1u-1 #go to upper node
                    if jj1u<0: #surface node has been processed
                        srf,jj1u = True,0 #break the while loop
        satnofl   = np.array(satnofl) #convert to numpy array
        satnofl_d = satnofl-1 #volumes with a satnofl volume below
        # All i00 indices must have 0 outflow #
        i00      = np.unique(np.concatenate((indsdry,imp,imp_d,satnofl_d,[ncv-1]))).astype(int) #add imp, imp_d, satnofl_d and lower boundary volume to volumes without outflow
        # All i11 indices can have outflow >0 #
        i11      = [ii for ii in range(i0dry) if ((ii in i00)==False)]

        ### Iterative guesses to determine qlim of Hirashima et al. (2010) ###
        glw      = np.zeros(ncv-1) #initialise guess of downward transported lwc (equivalent of qlim in Hirashima 2010)
        glwc     = np.copy(self.LWC) #initialise updated lwc after glw has been transported
        glwcacm  = np.copy(LWCacm) #initialise updated lwcacm after glw has been transported
        gtheta_e = np.copy(theta_e) #initialise updated theta_e after glw has been transported
        ghd      = np.zeros(ncv-1) #initialise updated head pressure after glw has been transported
        
        for j1 in np.flip(i11): #find qlim at each interface in turn starting from bottom-most interface 
            # qlim must equalise hd[j1] and hd[j1+1] assuming influx from j1-1 did not occur yet #
            inds11 = np.array([j1,j1+1]) #indices at each side of the interface of interest
            # First guess value for glw[j1] #
            if timer==0 or abs(glwflux_d1[j1]-abs(glwflux_d2[j1])>0.1*glwflux_d1[j1]):
                # glw flux is not in steady state: first guess equalises saturation[j1] and saturation[j1+1] #
                glw[j1] = thetaeff_equaliser(theta_i[inds11],theta_s[inds11],glwc[inds11],self.dz[inds11])
                glw[j1] = max(glw[j1],0) #avoid negative guess values
                glw[j1] = min(glw[j1],min(LWCav[j1],glwcacm[j1+1])) #guess values cannot exceed available LWC
            else:
                # glw flux is stable: initial guess keeps previous value #
                glw[j1] = dtsub*glwflux_d1[j1] #glw set to flux corresponding to glw calculated at previous Darcy step          
                glw[j1] = min(glw[j1],min(LWCav[j1],glwcacm[j1+1])) #guess values cannot exceed available LWC

            # Update gtheta_e and ghd #
            gtheta_e[inds11] = thetae_update(glw[j1],theta_i[inds11],theta_s[inds11],glwc[inds11],self.dz[inds11])
            ghd[inds11] = phead_vG(avG[inds11],nvG[inds11],mvG[inds11],gtheta_e[inds11]) #pressure head [m]
            # Calculate equilibrium #
            f_eq = ghd[j1]-ghd[j1+1]-dltz[j1] #Hirashima 2010 Eq.(20) evaluated at interface

            if abs(f_eq>1.):
                # Initial guess gives f_eq far from 0: Bisection algorithm #
                glw[j1] = flux_bisection(glw[j1],LWCav,glwcacm,theta_i[inds11],theta_s[inds11],glwc[inds11],self.dz[inds11],avG[inds11],nvG[inds11],mvG[inds11],eps_cvg)
            else:
                # Initial guess gives f_eq close to 0: Newton-Raphson algorithm #
                glw[j1] = flux_newtonraphson(glw[j1],LWCav,glwcacm,theta_i[inds11],theta_s[inds11],glwc[inds11],self.dz[inds11],avG[inds11],nvG[inds11],mvG[inds11],eps_cvg)
                             
            # Update glwc and glwcacm according to updated flow (to compute gtheta_e when dealing with upper nodes) #
            glwc         = self.LWC+np.append(0,glw)-np.append(glw,0) #update glwc 
            glwcacm      = self.dz*theta_s-glwc #update glwcacm
            glwcacm[imp] = 0. #do not allow for water storage in impermeable ice lenses

        glwflux_d2 = np.copy(glwflux_d1) #save glw fluxes computed at Darcy step -2
        glwflux_d1 = glw/dtsub #save glw fluxes computed at Darcy step -1
        ### Compute the water fluxes ###
        qlim  = np.append(glw,0) #qlim satisfying Hirashima 2010 Eq.(20) [m] (0 outflow from last node)
        q0    = bigkedg*(dhdz+1) #Darcy flux with initial step conditions [m s-1] Hirashima 2010 Eq.(1)
        ifl   = np.logical_and(qlim>1e-20,q0>1e-20) #indices with non-zero flow    
        
        qstep = np.zeros(ncv) #initialise qstep: total flux over dtsub time step
        qstep[ifl] = qlim[ifl]*(1-np.exp(dtsub*(-1)*q0[ifl]/qlim[ifl])) #Hirashima 2010 Eq.(23) [m we]
        qstep[-1] = 0 #lower boundary condition: no outflow
        qstep  = np.minimum(qstep,LWCav) #make sure that total flow does not exceed LWCav
        
        lwcin  = np.append(0.,qstep[0:-1]) #inflow for each node (with no upper boundary flux)
        lwcout = np.copy(qstep) #outflow for each node [m we]
        if (lat_runoff==True): #calulate lateral runoff in nodes that have no downwards outflux
            ii_rf = np.concatenate((imp_d,satnofl_d)).astype(int)
            ii_rf = np.intersect1d(ii_rf,np.where(LWCav>0)[0])
            # Choose ZuoOerlemans or Darcy runoff #
            #rfout = runoffZuoOerlemans(dtsub,slope_proxy,LWCav,ii_rf) #compute runoff with Zuo and Oerlemans 1996 parameterisation
            rfout = runoffDarcy(dtsub,slope_proxy,bigk,ii_rf) #compute runoff with Darcy parameterisation
            rfout = np.minimum(rfout,LWCav-lwcout) #make sure that runoff does not cause total flow to exceed LWCav
            lwcout += rfout #add to the LWC loss of the nodes
            runoff1 = sum(rfout) #total non-surface runoff [m we]
        else: #no lateral runoff apart from the saturated surface node
            runoff1 = 0.
        self.LWC  = self.LWC+lwcin-lwcout #compute LWC
        if timer+dtsub==timetot: #at end of Darcy routine: avoid LWC>0 in nodes of Darcy-modified density
            runoff1 = runoff1+sum(self.LWC[rhomod]) #LWC of Darcy-modified density assumed to run off
            self.LWC[rhomod] = 0. #Darcy-modified density nodes have 0 LWC
        runofftot = runofftot+runoff0+runoff1 #add contributions to total runoff

        ### Refreezing ###
        cold_content = cp_i*self.mass*(273.15-self.Tz) #cold content of nodes [J]
        refr_pot_ht   = cold_content/LF_I #refreezing potential of nodes from cold content [kg]
        refr_pot_v    = 0.917*phi*self.dz*1000-1e-6 #refreezing potential of nodes from available volume [kg] (use numerical safety)
        refr_pot      = np.minimum(refr_pot_ht,refr_pot_v) #refreezing potential of nodes [kg]
        refr         = np.minimum(self.LWC*1000,refr_pot) #refreezing in each node [kg]
        refr[refr<0] = 0 #avoid negative refreezing
        self.LWC     = np.maximum(0,self.LWC-refr*1e-3) #liquid mass loss, avoids numerical rounding errors
        self.mass    = self.mass+refr #solid mass gain
        self.rho     = self.mass/self.dz #update density
        phi          = (917-self.rho)/917 #update porosity
        phi          = np.maximum(0,phi) #avoid numerical errors of very small negative phi
        latheat      = refr*LF_I # latent heat released due to the refreezing [J]
        cold_content = cold_content-latheat # remaining cold content [J]
        self.Tz      = 273.15 - cold_content/(cp_i*self.mass) #updated temperatures of nodes
        cp_i         = 152.5+7.122*self.Tz #specific heat of ice [J kg-1 K-1] Cuffey and Paterson 2010 (9.1)
        refr_tot      = refr_tot+sum(refr)/1000 #add contributions to total refreezing [m we]

        timer = timer+dtsub #increase timer
        
        ### Adapt time step ###
        if np.any(ifl==True):
            qratio = max(qstep[ifl]/qlim[ifl]) #max ratio qstep/qlim
        else: #no single node had outflow>0
            qratio = 0.
        if qratio<0.25:
            dtsub = 1.25*dtsub #increase time step
        if qratio>0.75:
            dtsub = 0.75*dtsub #decrease time step
        dtsub = max(dtsub,dtmin) #keep time step bounded
        dtsub = min(dtsub,dtmax) #keep time step bounded
        dtsub = min(dtsub,timetot-timer) #ensures to end exactly at timetot

        
    ## Sanity checks
    if abs(sum(self.LWC)+refr_tot+runofftot-(melt_mass_tot/1000+sum(initial_lwc)+rain_vol_tot)) > 1e-12: #check for water balance
        print('Liquid water loss/gain, amount:',sum(self.LWC)+sum(self.refrozen)+self.runoff-(melt_mass_tot/1000+sum(initial_lwc)+rain_vol_tot))
    
    if max(self.Tz>273.15):
        print('Max Tz:',max(self.Tz))

    # Check cold layers are dry #
    coldlayers = np.where(self.Tz<273.15)[0]
    if np.all(self.LWC[coldlayers]<1e-9):
        self.LWC[coldlayers] = 0.
    if np.any(self.LWC[coldlayers]>0.):
        print('Problem: water content in a cold layer')

    if time.time()-ticdarcy>=10:
        print(f'{iii} CFM time: {self.modeltime[iii]}')
        print(f'Darcy scheme run time: {np.around(time.time()-ticdarcy,2)}')    

    return self.rho,self.age,self.dz,self.Tz,self.r2,self.z,self.mass,self.dzn,self.LWC,meltgridtrack,refr_tot,runofftot


def LWC_correct(self):
    '''
    *** TEST FUNCTION ***
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
    *** TEST FUNCTION ***
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













