
'''
This is a code applying water flow in firn with single domain approach:
- Richards Equation in Matrix Flow domain, no Preferential Flow domain

Essentially a simplification of the dual domain flow scheme
'''
import matplotlib.pyplot as plt
import numpy as np
import math
import time

from constants import *
from fcts_snowpackflow import TDMAsolver
from fcts_snowpackflow import splitCFM
from fcts_snowpackflow import combineCFM
from fcts_snowpackflow import restrictdom
from fcts_snowpackflow import lengthendom
from fcts_snowpackflow import Mrefreezing
from fcts_snowpackflow import Msatexcess
from fcts_snowpackflow import runoff
# from fcts_snowpackflow import initialrepartition
from fcts_snowpackflow import Micedryer
# from fcts_snowpackflow import distribute_tostore_single

from merge import mergesurf

def resingledomain(self,iii):
    
    ##### First: melting of the surface layers, taken from melt.py #####
    melt_volume_IE      = self.snowmeltSec[iii] * S_PER_YEAR # This still has to be checked by Max (division by self.c['stpsPerYear']?) [m]
    melt_volume_WE      = melt_volume_IE * RHO_I_MGM # [m]
    melt_mass           = melt_volume_WE * 1000. # [kg]
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
    pm_plwc = self.PLWC_mem[ind1]/self.dz[ind1] * pm_dz 

    ## Melted boxes are accomodated by just adding more (new) boxes at the bottom of the column
    ## Beware of this if you are not modeling to firn-ice transition depth.
    divider         = num_boxes_melted #VV nb of melted boxes, including the partially melted
    self.rho        = np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted))) #VV add at bottom of column as many layers as were melted away
    self.LWC        = np.concatenate((self.LWC[ind1:-1] , self.LWC[-1]*np.ones(num_boxes_melted)))
    # self.LWC        = np.concatenate((self.LWC[ind1:-1] , np.zeros(num_boxes_melted))) # This is better but should be equivalent as last layer should be an ice layer -> with 0 LWC
    self.LWC[0]     = pm_lwc #VV LWC calculated for the partially melted layer
    
    self.PLWC_mem   = np.concatenate((self.PLWC_mem[ind1:-1] , np.zeros(num_boxes_melted)))
    self.PLWC_mem[0]= pm_plwc
    # all the water that was in the PFdom of the melted layers is also for input
    
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
    
    ### Avoid extremely thin surface layer ###
    if self.c['merging'] == True:
        #if self.dz[0] < 0.1*self.c['merge_min']:
        if self.dz[0] < 1e-4:
            self.dz,self.z,self.gridLen,self.dx,self.rho,self.age,self.LWC,self.PLWC_mem,self.mass,self.mass_sum,self.sigma,self.bdot_mean,\
                        self.Dcon,self.T_mean,self.T10m,self.r2 = mergesurf(self,self.c['merge_min'])
    
    ##### Proceed to split of the CFM layers to have layers of max thickness #####
    '''We should split the lwc variable of the CFM in MLWC and PLWC to allow repartition of water between both domains
    Or implement a self.PLWCmemory variable that remembers what was in PFdom, even though we transfer everything in MFdom (maybe this
    is easier for other subroutines), next time we call for PrefFlow, we can simply remove PLWC_mem from LWC, give LWC to MFdom and give
    PLWC_mem to PFdom, see Plwc_mem in commented lines of this code'''
    vert_res = 0.2 # maximal layer thickness allowed, with upstream weighted mean for MK at interfaces, we can take larger value
    rhoC,dzC,TzC,massC,lwcC,Plwc_memC,r2C = restrictdom(self)
    split_list,rho,dz,Tz,mass,lwc,Plwc_mem,r2 = splitCFM(rhoC,dzC,TzC,massC,lwcC,Plwc_memC,r2C,vert_res)    #print('len(lwc) after restrictdom is:',len(lwc))
    
    ### Parameters we might want to change according to different runs ###
    tstep = 900. # Frequency of exchange betwwen domains and of refreezing [s]
    slope = 3.4907e-3 # slope that appears in  lateral runoff formula of Zuo and Oerlemans 1996, possibly in RE
    cst_melt = 1 # constant meltwater input: set this to 1
    sinus_melt = 0 # meltwater input according to sinusoidal daily cycle: set this to 1 !! Make sure dtCFM is 24h !!
    bigF = 1e-12*np.ones_like(dz) # extremely low value since there is no preferential flow domain, keeping it allows to use same functions as dual domain scheme
    impermeability = 1. # 1 to activate impermeability of thick ice layers, 0 otherwise
    impthick = 0.0 # minimum thickness of ice layers to be impermeable [m]
    rhoimp = 810.
    rofflim = 0.95 # effective saturation of MFdom above which we launch extra runoff to have MeffSat == rofflim

    grain = r2 ** (1/2) # grain size [m]
    
    Mlwc = lwc-Plwc_mem
    Plwc = Plwc_mem
    
    depth = np.cumsum(dz)-dz/2 # depth of the center of every layer [m]
    depth_imin1 = np.delete(np.append(0,depth),-1) # depth staggered 1 level below [m]
    zstep = depth-depth_imin1 # distance between centers of adjacent layers (between center[i] and center[i-1]), will be used for the PFdom [m]
    zstep_iplus1 = np.append(np.delete(zstep,0),dz[-1]/2) # zstep vector staggered 1 level lower
    
    ### Time parameters ###
    dtCFM = self.dt #duratin of timesteps in CFM, constant value, [s]
    Mdeltat = 300. #duration of timesteps for this script, choose an arbitrary starting value (<= dtCFM) [s]
    Mdeltat_new = 300. #further, we adjust deltat_new and then assign this value to deltat, start with deltat_new == deltat [s]
    deltat_max = 1*tstep #maximum time step allowed
    deltat_min = 1e-20 #minimum time step allowed (1e-10 in D'Amboise 2017, not specified in Wever 2014) [s]
    Mtime_counter = 0 # time for REMF, will be reset at 0 every tstep seconds
    Pdeltat = 300. #duration of timesteps for flow in PFdom, this is adapted according to max value of Ptheta [s]
    Pdeltat_new = 300. 
    Ptime_counter = 0 # time for PF, will be reset at 0 every tstep seconds
    timer = 0 # this tracks the time passing --> remains between 0 and dtCFM [s]

    ### Calculate pore space available in every layer --> water content at saturation 
    porosity           = 1 - rho/RHO_I # Definition of porosity [/]
    porespace_vol      = porosity * dz # Pore space of each layer [m]
    porosity_refr      = porosity*RHO_I/RHO_W_KGM # space available for liq water volume once refrozen, Wever 2014 (9) [/]
    #porosity_refr      = np.maximum(porosity_refr,1e-4) # allow space for minimum water content required in both domains for numerical stability, 1e-4 is equivalent to 916.9 density
    porosity_refr      = np.maximum(porosity_refr,17e-3) # allow space for minimum water content required in both domains for numerical stability, 17e-3 is equivalent to 900.0 density
    porespace_refr_vol = porosity_refr*dz # Available pore space of each layer [m]
    
    theta_sat = porosity_refr # value of volumetric water content in saturated conditions [/]
    Mtheta_sat = (1-bigF)*theta_sat
    Ptheta_sat = bigF*theta_sat
    totrunoff = 0.

    if iii == 0 and np.any(lwc>0): # avoid problems if initial profile with lwc
        lwc,totrunoff = initialrepartition(dz,rho,lwc,theta_sat,rhoimp,totrunoff)
        Mlwc = (1-bigF)*lwc # ditribute water between MFdom and PFdom

    ### Convergence of the solution between successive iterations ###
    crtn_head = 1e-3 #Criterion on head pressure, see Wever 2014 Appendix
    crtn_theta = 1e-5 #Criterion on water content, see Wever 2014 Appendix
    crtn_MB = 1e-8 # Mass balance criterion, as in SNOWPACK ReSolver line 576
    
    theta_min_fr = 1e-4 # Threshold value to allow freezing, Wever 2014 Appendix
    lwc_min_fr = theta_min_fr*dz # corresponding lwc threshold
    
    ### Use rain climatic inputy: still requires confirmation of Max!!! ###
    try:
       raintoadd = self.rainSec[iii] * S_PER_YEAR * RHO_I_MGM # [m]
    except:
        raintoadd = 0.
    melt_vol_a += raintoadd # we add the rain to the melt volume, everything is in m we (total value for this time step)

    totrefrozen_lwc = 0 # will calculate total lwc refrozen [mWE]
    refrozenlay = np.zeros_like(dz) # will calculate total lwc refrozen in every layer [mWE]
    
    surfrunoff = 0. # all lateral runoff of this CFMstep    
    if rho[0] >= rhoimp:
        surfrunoff += melt_vol_a
        totinflux  = 0.
    elif rho[0] <= rhoimp:
        totinflux = melt_vol_a # total surface input for this CFM timestep (should be double checked) [m]
        
    mean_influx  = totinflux/dtCFM # mean value of the liquid water flux at surface [m/s]
    liqout = 0. # water flowing out at the bottom of the domain [m/s]
    
     ### To calculate head pressure: van Genuchten parameters and correction factor
    alpha_vG = 4.4e6*(rho/(2*grain))**(-0.98) # Hirashima 2014 (5) ; typical value ~35.0 
    n_vG     = 1 + 2.7e-3*(rho/(2*grain))**(0.61) # Hirashima 2014 (6) ; typical value ~4.0
    m_vG     = 1 - 1/n_vG # Wever 2014 (8) ; typical value ~0.75
    h_e   = 5.8e-3 # air entry pressure for pore size of 5mm at 273K [m], Wever 2014 
    Sc    = (1 + (alpha_vG*h_e)**n_vG)**(-m_vG) #Saturation at cut-off point [/], see Ippisch et al., 2006 eq(11)

    h_we = 0.0437/(2*grain*1e3) + 0.01074 # water entry suction for snow, Hirashima 2014 (15) [m]
    MSat_we = (1+(alpha_vG*abs(h_we))**n_vG)**-m_vG / Sc # Effective saturation corresponding to the water entry suction
    mu   = 0.001792 # Dynamic viscosity of water [kg m-1 s-1] , can be added to constants list
    Ksat = RHO_W_KGM*GRAVITY/mu * 3.0*(grain)**2*np.exp(-0.013*rho) # Hydraulic conductivity at saturation (>0) [m s-1], Formula of Calonne et al. 2012, see Wever 2015 (7) and D'Amboise 2017 (10)

    theta_min = crtn_theta/10*np.ones_like(dz) # minimum initial value of theta used for numerical stability
    ### MFdom variables ###
    Mtheta = Mlwc/dz # Volumetric liquid water content [/]
    Mprewetting = dz*(theta_min-Mtheta) # amount of water required for prewetting of MFdom [m]
    Mprewetting = np.maximum(0,Mprewetting) # exclude negative values (where Mtheta already higher than theta_min)
    Mtheta = np.maximum(Mtheta,theta_min) # modify theta
    Mlwc = Mtheta*dz # modify lwc
    ## Define residual water content as in Wever 2014 ##
    Mthetar = np.minimum((np.ones_like(Mtheta)*0.02),0.9*Mtheta) # initial residual water content [/], Wever 2014 (10)
    ## Calculate effective saturation ##
    MeffSat = (Mtheta-Mthetar)/(Mtheta_sat-Mthetar) # effective saturation of the MFdom layers []
    
    ### PFdom variables ###
    Ptheta = Plwc/dz # Volumetric liquid water content [/]
    ## Calculate effective saturation ##
    PeffSat = Ptheta/Ptheta_sat # effective saturation of the MFdom layers []

    totflux = 0
    massconsalarm = 0
    provrunoff = 0.

    ### Move water from layers where effSat is above 1 ###
    # Should be executed before head calculations, otherwise RuntimeWarning
    if np.any(Mtheta > Mtheta_sat): # Evacuate water in layers where we exceed maximal saturation towards other layers
        Mtheta,Mthetar,MeffSat,Mlwc,totrunoff = Msatexcess(dz,rho,Mtheta,Mtheta_sat,crtn_theta,rhoimp,totrunoff)
    
    if np.any(Mtheta[rho>rhoimp]>theta_min[rho>rhoimp]):
        Mtheta,Mthetar,MeffSat,Mlwc,totrunoff = Micedryer(dz,rho,Mtheta,Mtheta_sat,crtn_theta,rhoimp,totrunoff) # Evacuate water from layers of density > rhoimp
      
    ## Spot the aquifer #
    ice1 = np.zeros_like(dz) # bit array for ice layers
    ice1[np.where(rho>=rhoimp)[0]] = 1
    if np.any(ice1==0):
        noflow = np.where(ice1==0)[0][-1] + 1 # this layer and all layers below are saturated or are ice layers
    elif np.all(ice1>0): # If all the domain is ice or saturated
        noflow = 0 # all layers are saturated or are ice layers
     
    indaqM = 1*noflow
    tostore = 0
    if np.any(MeffSat>0.95):
        icebottom = np.zeros_like(dz)
        icebottom[np.where(rho>=rhoimp)[0]] = 1 # ice layers
        Msatbottom = np.zeros_like(dz)
        Msatbottom[np.where(MeffSat>=0.95)[0]] = 1 # saturated layers
        Msaticebottom = Msatbottom+icebottom
        if np.any(Msaticebottom==0):
            indaqM = np.where(Msaticebottom==0)[0][-1]+1 #from this layer, all layers below are saturated or ice
            if indaqM > 0:
                if Mtheta[indaqM-1] > theta_min_fr and Mtheta[indaqM-2] <= theta_min_fr: # case where rest of the storage could not fill entirely a layer -> consider this layer as part of the aquifer
                    indaqM -= 1
        elif np.all(Msaticebottom>0):
            indaqM = 0
        tostore += sum((Mtheta[indaqM:]-theta_min[indaqM:])*dz[indaqM:]) # store all that water, which will be redistributed at end of flow routine
        ### Update all variables
        Mtheta[indaqM:] = 1*theta_min[indaqM:]
        Mlwc[indaqM:] = Mtheta[indaqM:]*dz[indaqM:]
        Mthetar[indaqM:] = np.minimum((np.ones_like(Mtheta[indaqM:])*0.02),0.9*Mtheta[indaqM:])
        MeffSat[indaqM:] = (Mtheta[indaqM:]-Mthetar[indaqM:])/(Mtheta_sat[indaqM:]-Mthetar[indaqM:])
    if tostore > 0:
        ### Consider the top of the aquifer as the lowest layer where we can have incoming flow
        noflow = 1*indaqM # top of the aquifer is the highest layer saturated in both PFdom and MFdom
    
    ### Head pressure ###
    Mhead  = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
    
    ## Avoid letting surface influx come in when the entire firn column is saturated
    Mthetar_old = Mthetar
    #Pbottom = 0
    bottom = len(dz)-1
    if impermeability == 1: #if we use impermeability of thick ice layers
        ice1 = np.zeros_like(dz) # bit array for ice layers
        sat1 = np.zeros_like(dz) # bit array for saturated layers
        ice1[np.where(rho>=rhoimp)[0]] = 1
        sat1[np.where(MeffSat>=0.95)[0]] = 1
        icesat = ice1+sat1 # all layers that have value 0 are neither ice nor saturated
        if np.any(icesat==0):
            bottom = np.where(icesat == 0)[0][-1] # New bottom of the domain: lowest layer that will receive water flow
            if bottom == 0: # if only the surface layer can accomodate some water
                ## Add all possible input in surface layer, rest will be treated as runoff
                Mlwcsat0 = dz[0]*Mtheta_sat[0]*0.95 # maximum amount of water that can be held in surface layer
                inflow = min(totinflux,(Mlwcsat0-Mlwc[0])) # input of water in the surface layer is limited by lwcsat0 or by total input during the time step
                totinflux -= inflow # only the inflow part of the total influx can be accomodated by the firn column (inflow <= totinflux)
                surfrunoff += totinflux # if column is entirely saturated totinflux > 0 // if column is not entirely saturated totinflux == 0 (and runoff is not increased)
                totinflux = 0. # influx has been partitioned in surface layer and runoff (if surface layer now saturated)
                mean_influx = 0.
                Mlwc[0] += inflow # inflow has been added to surface lwc
                Mtheta[0] = Mlwc[0]/dz[0] # update Mtheta
                Mthetar[0] = min((0.02),0.9*Mtheta[0])
                if Mtheta[0]<Mthetar[0]+1e-6:
                    if Mtheta[0]>crtn_theta/10:
                        Mthetar[0] = Mtheta[0] - crtn_theta/10
                    if Mtheta[0]<=crtn_theta/10:
                        Mthetar[0] = 0
                MeffSat[0] = (Mtheta[0]-Mthetar[0])/(Mtheta_sat[0]-Mthetar[0]) # update MeffSat
        elif np.all(icesat>0): # If all the domain is ice or saturated
            bottom = 0 # we will not proceed to RE
            if surfrunoff == 0.: # if the input has not been attributed to surface runoff (because surface layer is not at rhoimp)
                surfrunoff += mean_influx*dtCFM # [m]
                totinflux  = 0.
                mean_influx = 0.
                print('All melt in runoff because all domain is ice or saturated, totinflux is:',totinflux)         
    
    totlwc0 = sum(Mlwc) + sum(Plwc)
    
    #Assign old values
    Mtheta_old = 1*Mtheta
    Mthetar_old = 1*Mthetar
    MeffSat_old = 1*MeffSat
    Mlwc_old = 1*Mlwc
    Mhead_old = 1*Mhead
    Mhead_old_imin1 = np.append(0,np.delete(Mhead,-1))
    Mdtheta_dh = np.zeros_like(dz)
    Mtheta_previous = np.zeros_like(dz)
    Mhead_previous = np.zeros_like(dz)
    delta_Mtheta = np.zeros_like(dz)
    delta_Mhead2 = np.zeros_like(dz)
    
    ##### Start the Global time loop #####
    while timer < dtCFM:
        
        if timer + tstep >= dtCFM: # make sure we end precisely at the time step of the CFM
            tstep = dtCFM - timer
        
        bottom = max(0,min(bottom,noflow-1)) # don't proceed to RE in the aquifer
        sat_layers = np.array([]) # we don't want inflow in fully saturated layers
        while Mtime_counter < tstep:
            Mflag = 0 # Too many iterations or other cases where finer time step is required
            provrunoff = 0.
            
            if np.any(Mtheta<1e-10): # For some unknown reason, this happens on some rare occasions (only in the surface layer when layer below is an ice layer I think)
                instb = np.where(Mtheta<1e-10)[0]
                for ii in instb:
                    Mtheta[ii] = theta_min[ii] # set theta back to numerical stability threshold and update all MFdom variables
                    Mthetar[ii] = 0.
                    MeffSat[ii] = (Mtheta[ii]-Mthetar[ii])/(Mtheta_sat[ii]-Mthetar[ii])
                    Mlwc[ii] = Mtheta[ii]*dz[ii]
                    Mhead[ii] = -1*1/alpha_vG[ii] * ((Sc[ii] * MeffSat[ii])**(-1/m_vG[ii])-1)**(1/n_vG[ii]) # [m] Wever 2014 (3)
                    print('We had to fix an instability in MFdom, layers were:',instb)
            
            ### Liquid water input of the time step ###
            ## Constant melt ##
            if cst_melt == 1:
                liqin = 1*mean_influx #This works because we assume cst influx
                Mliq_in = 1*liqin # All influx in MFdom, Wever 2016 [m/s]
                if Mliq_in*Mdeltat_new >= (rofflim*(Mtheta_sat[0]-Mthetar[0])+Mthetar[0]-Mtheta[0])*dz[0] and bottom>0:
                    if Mdeltat_new >= 60. and rho[0]<rhoimp: #We do this as long as time step is not too small and surface not impermeable
                        Mdeltat_new = Mdeltat_new/2
                        continue
                    elif Mdeltat_new < 60. or rho[0]>=rhoimp: #If time step is very small or surface impermeable
                        Mliq_in = 0. #Set influx to 0 and add this water to the provisional runoff
                        provrunoff = mean_influx*Mdeltat_new # provisional runoff (if backstep it is reset to 0 and not accounted (yet) in totrunoff)
                        #print('About to use provrunoff, Mdeltat and rho[0] are:',Mdeltat,rho[0])
                        
                Mliq_out = liqout # Not sure this works if we use liq_out different than 0
            ## Sinus function for daily melt ##
            if sinus_melt == 1:
                if timer < 86400:
                    flux_sinus = (1/Mdeltat)*totinflux/2 * (np.cos(np.pi/dtCFM * (timer+Mtime_counter))-np.cos(np.pi/dtCFM * (timer+Mtime_counter+Mdeltat))) # integrating the sinus function, [m/s]
                    Mliq_in = flux_sinus # All influx in MFdom, Wever 2016 [m/s]
                if timer > 86400:
                    Mliq_in = 0.
                Mliq_out = liqout # Not sure this works if we use liq_out different than 0
            
            Mdeltat = 1*Mdeltat_new # Use the dynamically adjusted time step
            
            ### Assign old values ###
            Mtheta_old = 1*Mtheta
            Mthetar_old = 1*Mthetar
            MeffSat_old = 1*MeffSat
            Mlwc_old = 1*Mlwc
            totMlwc_old = sum(Mlwc_old)
            Mhead_old = 1*Mhead
            Mhead_old_imin1 = np.append(0,np.delete(Mhead,-1))
            
            Mthetar[0:(bottom+1)] = np.minimum((np.ones_like(Mtheta[0:(bottom+1)])*0.02),0.9*Mtheta[0:(bottom+1)])
            ### Update of effSat and head (because Mthetar might have changed) ##
            MeffSat = (Mtheta-Mthetar)/(Mtheta_sat-Mthetar) # []
            Mhead  = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
    
            ### This is to reduce calculation time: we don't solve RE for the lowest parts of the domain if these are dry or only filled due to emptying of Pend_lwc ###
            if Mtheta[bottom-1] <= theta_min_fr + crtn_theta: # if layer above lowest saturated layer is dry then water at the bottom does not come from matrix flow
                # We are in the case layers are getting filled by Pend_lwc -> we don't need to model these
                wetlayers = np.where(Mtheta[0:bottom] >= (np.maximum(0.9*0.02,Mthetar[0:bottom])))[0] # all the layers having a significant theta value
                if np.size(wetlayers) == 0.: # in the case all our layers are dry
                    ic = 0 # last wet layer can be considered the surface layer as there will be some meltwater input
                elif np.size(wetlayers) > 0.: # in the case we have some layers already wet
                    ic = wetlayers[-1] # ic is the index of the deepest wet layer
                lowestdepth = depth[ic] + 0.5 # we take 0.5 meter below the lowest wet layer as the domain to perform calculations
                if lowestdepth < depth[-1]: # we want to define the maximum depth of the calculation domain
                    lc = np.where(depth >= lowestdepth)[0][0] # lc is the index of the last layer of the calculation domain
                elif lowestdepth >= depth[-1]: # the maximum depth of the calculation domain must not be greater than the maximum depth of the real domain
                    lc = len(dz)-1
                bottom = min(lc,bottom) # replace bottom by lc except in the case where we have ice layers at bottom that are then taken as limit of domain even if the limit just calculated goes further below
            bottom = max(0,min(bottom,noflow-1)) # test
    
            ### Similar to aquifer calculation, because MeffSat might have increased over tstep up to saturation due to inflow ###
            icedom = np.zeros(bottom) # bit array for ice layers
            satdom = np.zeros(bottom) # bit array for saturated layers in MFdom
            icedom[np.where(rho[0:bottom]>=rhoimp)[0]] = 1 # ice layers
            satdom[np.where(MeffSat[0:bottom]>=0.95)[0]] = 1 # saturated layers in MFdom
            icesatdom = icedom+satdom # lowest layer with 0 value will be last layer taken into account in RE
            if np.all(icesatdom>0): # if all layers are saturated or ice
                bottom = 0 # don't proceed to RE
            
            if bottom == 0 and mean_influx>0: # Case where bottom was not set to 0 before entering global time loop, but since then, surface layers became saturated
                    Mliq_in = 0.
                    provrunoff = mean_influx*Mdeltat # firn cannot accomodate input -> in provrunoff
             
            if Mliq_in>0 and rho[1]>=rhoimp: #If there is input but nothing can flow out of layer[0] because layer[1] is impermeable
                Mlwc[0] += Mliq_in*Mdeltat # We immediately give all the input of the time step to the surface layer
                Mtheta[0] = Mlwc[0]/dz[0]
                Mthetar[0] = np.minimum(0.02,0.9*Mtheta[0])
                MeffSat[0] = (Mtheta[0]-Mthetar[0])/(Mtheta_sat[0]-Mthetar[0])
                Mhead[0]  = -1*1/alpha_vG[0] * ((Sc[0] * MeffSat[0])**(-1/m_vG[0])-1)**(1/n_vG[0]) # [m] Wever 2014 (3)
                Mliq_in = 0. # There is no input to give anymore
                if ic == 0: # If all the layers are dry, we don't need to proceed to RE (water only in layer[0] and layer[1] is impermeable)
                    bottom = 0
                Mtheta_old[0] = 1*Mtheta[0]
                Mthetar_old[0] = 1*Mthetar[0]
                MeffSat_old[0] = 1*MeffSat[0]
                Mlwc_old[0] = 1*Mlwc[0]
                Mhead_old[0] = 1*Mhead[0]
        
            Miteration = 0 # Let's count the number of iterations required
                    
            delta_Mhead  = np.ones_like(dz) # Difference in head pressure between two iterative steps, calculated by TDMA
            delta_Mhead2 = np.ones_like(dz) # Difference in head pressure between two iterative steps, directly calculated (in case head is artificially adjusted)
            delta_Mtheta = np.ones_like(dz) # Difference in volumetric liquid water content between two iterative steps
            MB_Merror = np.ones_like(dz) # Mass balance error between discretised flux (with K values of last iteration) and changes calculated by the Picard scheme
            
            ### Spot the ice layers ###
            saturated = np.where(MeffSat >= 0.95)[0]
            sat_layers = np.unique(np.append(sat_layers,saturated))
            ice_layers = np.where(rho>=rhoimp) # density threshold should be used here
                
            if impthick == 0.:
                thick_ice = ice_layers[0]
                
            elif impthick > 0.:            
                if len(ice_layers[0]) > 0: # If there are ice layers
                    thick_ice = np.array([]) # this will contain all the indices of layers that are part of thick ice layers (thickness threshold is to be fixed)
                    start_ice = np.array([]) # this will contain all the indices where an ice layer starts
                    end_ice = np.array([]) # this will contain all the indices where an ice layer ends
                    cc = ice_layers[0][0] # first ice layer
                    while cc <= ice_layers[0][-1]: # loop through all ice layers
                        if cc in ice_layers[0]:
                            start = cc # index of the start of the ice layer
                            while (cc+1 in ice_layers[0]) == True: # goes through all the layers belonging to an individual ice layers
                                cc += 1
                            end = cc # index of the end of the ice layer
                            start_ice = np.append(start_ice,start) # add the start of the layer to the start list
                            end_ice = np.append(end_ice,end) # add the end of the layer to the end list
                        cc += 1
                    for dd in range(len(start_ice)): # For every ice layer
                        thickness = sum(dz[int(start_ice[dd]):int(end_ice[dd])+1]) # calculate thickness of every ice layer
                        if thickness >= impthick: # if the thickness is above the threshold thickness for impermeability
                            for ee in range(int(start_ice[dd]),int(end_ice[dd])+1):
                                thick_ice = np.append(thick_ice,ee) # write all the indices of the thick ice layer in this array
    
            while (np.max(np.append(0,abs(delta_Mtheta[0:(bottom+1)][MeffSat[0:(bottom+1)]<0.99]))) > crtn_theta or Miteration<2 or np.max(abs(MB_Merror[0:(bottom+1)]))>crtn_MB or np.max(np.append(0,abs(delta_Mhead2[0:(bottom+1)][MeffSat[0:(bottom+1)]>=0.99]))) > crtn_head) and Mflag==0 and bottom>0:
            # we append a 0 value for the criterion check because this avoids error using max in case no single layer 
            ##### Iterative loop for RE in MFdom #####
            # 3 conditions: difference below convergence criterion, min 2 iterations, special case flag is not on
                ### "center of layers = nodes" approach, as Wever 2014 and D'Amboise 2017 --> staggering ###
                dzdom     = dz[0:(bottom+1)]
                depthdom = np.cumsum(dzdom)-dzdom/2 # depth of the center of every layer [m]
                depthdom_imin1 = np.delete(np.append(0,depthdom),-1) # depth staggered 1 level below [m]
                zstepdom = depthdom-depthdom_imin1 # distance between centers of adjacent layers (between center[i] and center[i-1]), will be used for the PFdom [m]
                zstepdom_iplus1 = np.append(np.delete(zstepdom,0),dzdom[-1]/2) # zstep vector staggered 1 level lower         
                
                Mheaddown      = np.delete(Mhead,0) # Mhead values except last layer -> 1 element shorter than Mhead
                Mheadup        = np.delete(Mhead,-1) # Mhead values except first layer -> 1 element shorter than Mhead
                Mhead_imin1 = np.append(0,Mheadup) # Mhead vector staggered 1 level lower
                Mhead_iplus1   = np.append(Mheaddown,0) # Mhead vector staggered 1 level higher
                
                ### Hydraulic conductivity of every layer [Wever 2014 (11)] ###
                MKdom = Ksat[0:(bottom+1)] * MeffSat[0:(bottom+1)]**0.5 * (1-(1-MeffSat[0:(bottom+1)]**(1/m_vG[0:(bottom+1)]))**m_vG[0:(bottom+1)])**2 # (>0) [m s-1]
                #MKdom = 1*MK[0:(bottom+1)] # hydraulic saturation only on our calculation domain
                
                if impermeability == 1 and len(ice_layers[0]) > 0 : #if we use impermeability of thick ice layers
                    for iifl in thick_ice: # for all the layers that are part of thick ice layers
                        ii = int(iifl) # convert index to integer
                        if ii <= bottom:
                            MKdom[ii] = 0. # set very low value of hydraulic conductivity
                            if ii > 0:
                                MKdom[ii-1] = 0. # required for low value of MKtop[ii]
                                
                if len(sat_layers) > 0 : # make fully saturated layers impermeable
                    for iifl in sat_layers: # for all fully saturated layers
                        ii = int(iifl) # convert index to integer
                        if ii <= bottom:
                            MKdom[ii] = 0. # set very low value of hydraulic conductivity
                            if ii > 0:
                                MKdom[ii-1] = 0. # required for low value of MKtop[ii]  
                
                ### Hydraulic conductivity at the interfaces ###                                
                MKup         = np.delete(MKdom,-1) # K values except last layer -> 1 element shorter than K
                MKdown       = np.delete(MKdom,0) # K values except first layer -> 1 element shorter than K
                ## Upstream weighted mean for MK at interfaces: Forsyth 1995 (15), Szymkiewicz 2009 (3c)
                MK_inter     = MKup # Take value of layer above
                otherMK = np.where((Mhead_old[1:(bottom+1)]-Mhead_old_imin1[1:(bottom+1)])/zstepdom[1:] -1 >0) # Szymkiewicz 2009 (3c)
                MK_inter[otherMK[0]] = MKdown[otherMK[0]] # in these case, take value of layer below
                
                MKtop        = np.append(0,MK_inter) # Ktop will determine influx in every layer Ktop[0] is at surface layer (set to 0)
                MKbot        = np.append(MK_inter,0) # Kbot will determine outflux in every layer Kbot[-1] is at bottom layer (set to 0)
                
                ### Calculate analytically derivative of theta with respect to head (Ci in Zarba 1988 (3.10) and Celia 1990 (17)) ###
                Mdtheta_dh[0:(bottom+1)] = -1* (Mtheta_sat[0:(bottom+1)]-Mthetar[0:(bottom+1)])/Sc[0:(bottom+1)] * \
                    (-m_vG[0:(bottom+1)]*n_vG[0:(bottom+1)]*alpha_vG[0:(bottom+1)]**n_vG[0:(bottom+1)]) * (abs(Mhead[0:(bottom+1)]))**(n_vG[0:(bottom+1)]-1) / (1+(alpha_vG[0:(bottom+1)]*(abs(Mhead[0:(bottom+1)])))**n_vG[0:(bottom+1)])**(m_vG[0:(bottom+1)]+1)
                ##-1* as we use absolute value of head which is negative
                
                ### Flux calculations (part 1) ###
                Mqbot = MKbot*((Mhead[0:(bottom+1)]-Mhead_iplus1[0:(bottom+1)])/zstepdom_iplus1+1) # flux at bottom of every layer [m/s]
                Mqtop = MKtop*((Mhead_imin1[0:(bottom+1)]-Mhead[0:(bottom+1)])/zstepdom+1) # flux at top of every layer (simply staggered with respect to Mqbot) [m/s]
                Mqtop[0] = Mliq_in # flux boundary condition [m/s]
                Mqbot[-1] = Mliq_out # flux boundary condition [m/s]
                
                ### Parameters for the Tridiagonal Matrix Algorithm, see Zarba (1988) ###
                #Caution: for Zarba: i=0 at bottom, here: i=0 at surface --> be consistent in signs and with staggered grids
                #         alpha, beta, gamma not necessarily same (spatially speaking) as in Zarba 1988 but signs compensate
                #         ex: alpha(CFM) = beta(ZARBA) (spatially) but head_stagdown(CFM) = hi+1(ZARBA)
        
                # Intermediary calculations #
                Malpha = MKtop/(dzdom*zstepdom)
                Mbeta  = MKbot/(dzdom*zstepdom_iplus1)
                Mgamma = (MKtop-MKbot)/dzdom
                Mphi   = (Mdtheta_dh[0:(bottom+1)])/Mdeltat
            
                # Diagonals of the matrix: a is sub, b is principal, c is sup #
                a_diag = -1*Malpha
                a_diag = np.delete(a_diag,0)
                b_diag = Malpha+Mbeta+Mphi 
                c_diag = -1*Mbeta
                c_diag = np.delete(c_diag,-1)
                
                # Calculation of the residual (right-hand side of the equation) #
                Rd = (Mtheta_old[0:(bottom+1)]-Mtheta[0:(bottom+1)])/Mdeltat - Malpha*(Mhead[0:(bottom+1)]-Mhead_imin1[0:(bottom+1)]) - Mbeta*(Mhead[0:(bottom+1)]-Mhead_iplus1[0:(bottom+1)]) + Mgamma 
                
                # Surface boundary condition: constant flux [Zarba (1988) Table 3.2]
                Rd[0] = (Mtheta_old[0]-Mtheta[0])/Mdeltat - Mbeta[0]*(Mhead[0]-Mhead_iplus1[0]) - MKbot[0]/dz[0] + Mliq_in/dz[0]
                # Bottom boundary condition: constant flux
                Rd[bottom] = (Mtheta_old[bottom]-Mtheta[bottom])/Mdeltat - Malpha[bottom]*(Mhead[bottom]-Mhead_imin1[bottom]) + MKtop[bottom]/dzdom[bottom] - Mliq_out/dzdom[bottom]
                
                ### Solve tridiagonal matrix ###
                #delta_Mhead = NPtrid(a_diag,b_diag,c_diag,Rd)
                delta_Mhead = TDMAsolver(a_diag,b_diag,c_diag,Rd)
                
                ### Save values of head and theta at the last iterative step ###
                Mtheta_previous[0:(bottom+1)] = Mtheta[0:(bottom+1)]
                Mhead_previous[0:(bottom+1)] = Mhead[0:(bottom+1)]
                
                ### Assign new value to pressure head ###
                Mhead[0:(bottom+1)] = Mhead[0:(bottom+1)] + delta_Mhead #[m]
                Mtheta[0:(bottom+1)] = Mthetar[0:(bottom+1)] + (Mtheta_sat[0:(bottom+1)]-Mthetar[0:(bottom+1)])*(1/Sc[0:(bottom+1)])*(1+(alpha_vG[0:(bottom+1)]*(abs(Mhead[0:(bottom+1)])))**n_vG[0:(bottom+1)])**(-m_vG[0:(bottom+1)]) # []
                
                
                MeffSat[0:(bottom+1)] = (Mtheta[0:(bottom+1)]-Mthetar[0:(bottom+1)])/(Mtheta_sat[0:(bottom+1)]-Mthetar[0:(bottom+1)]) # []
                Mlwc[0:(bottom+1)]   = Mtheta[0:(bottom+1)]*dz[0:(bottom+1)] # [m]            
                
                delta_Mtheta[0:(bottom+1)] = Mtheta[0:(bottom+1)]-Mtheta_previous[0:(bottom+1)]
                delta_Mhead2[0:(bottom+1)] = Mhead[0:(bottom+1)]-Mhead_previous[0:(bottom+1)] # Useful if head artificially adjusted (eg to minimum value), normally delta_Mhead2==delta_Mhead
                
                ### Flux calculations (part 2) ###
                change_Mlwc = Mlwc-Mlwc_old # change in lwc over the time step [m]
                MB_Merror = Mqbot*Mdeltat+change_Mlwc[0:(bottom+1)]-Mqtop*Mdeltat # change in lwc not explained by flux discretisation [m]
                totMB_Merror = sum(MB_Merror) # sum of the errors in whole column [m]
        
                Miteration += 1 # count iterations
                
                ### Back step from Wever 2014 (don't use back step for MB criterion, let model iterate further!)
                if Miteration == 16 or np.any(Mhead[0:(bottom+1)]>0) or np.any(abs(delta_Mhead2[0:(bottom+1)])>1000) or np.max(MeffSat[0:(bottom+1)])>1.0:#test
                    
                    # Assign old values and we will start again with a shorter time step #
                    Mtheta[0:(bottom+1)] = 1*Mtheta_old[0:(bottom+1)]
                    Mthetar[0:(bottom+1)] = 1*Mthetar_old[0:(bottom+1)]
                    MeffSat[0:(bottom+1)] = 1*MeffSat_old[0:(bottom+1)]
                    Mlwc[0:(bottom+1)] = 1*Mlwc_old[0:(bottom+1)]
                    Mhead[0:(bottom+1)] = 1*Mhead_old[0:(bottom+1)]
                    Mdeltat_new = Mdeltat/2
                    Mflag = 1
                                    
            ### According to number of iterations, adjust the next time step
            if Miteration <= 5 and Mflag == 0: # Increase next time step but don't exceed max time step value
                Mdeltat_new = 1.25*Mdeltat
                Mdeltat_new = min(deltat_max,Mdeltat_new)
            if Miteration > 5 and Miteration <= 10 and Mflag == 0: # Keep same time step
                Mdeltat_new = 1*Mdeltat
            if Miteration > 10 and Miteration <= 15 and Mflag == 0: # Deccrease next time step but don't exceed min time step value
                Mdeltat_new = Mdeltat/2    
                Mdeltat_new = max(deltat_min,Mdeltat_new)   
            if Mflag == 1: #Now we are in the time loop --> don't go until the end of time loop and restart for the process with finer time step
                Mdeltat_new = Mdeltat/2    
                Mdeltat_new = max(deltat_min,Mdeltat_new)
                if Mdeltat_new <= 10*deltat_min:
                    print('Too small Mdeltat, make the simulation crash, self.modeltime[iii] is:',self.modeltime[iii])
                    Mtheta_sat[0:10] = -1e9
                continue # go back to top of Matrix time loop
            
            if provrunoff > 0: # Now that Mdeltat definitely passed, we can add provrunoff to totrunoff
                totrunoff += provrunoff
                provrunoff = 0.
            
    
            totflux += mean_influx*Mdeltat # meltwater influx so far  
            Mtime_counter += Mdeltat # increase the Mtime_counter

            if bottom == noflow-1: # If bottom is at the limit where we allow water to flow
                if Mtheta[bottom] > theta_min_fr: # If water reached that limit -> store it in the aquifer water
                    # store the water to make effSat equal the effSat of overlying layer (no sharp wetting front)
                    Mwaterout = dz[bottom]*(Mtheta[bottom]-(MeffSat[bottom-1]*(Mtheta_sat[bottom]-Mthetar[bottom])+Mthetar[bottom]))
                    Mwaterout = max(Mwaterout,0)
                    Mwaterout = min(Mwaterout,dz[bottom]*(Mtheta[bottom]-theta_min[bottom]))
                    if Mwaterout>0:
                        tostore += Mwaterout
                        Mtheta[bottom] -= Mwaterout/dz[bottom]
                        ### Update all variables
                        Mlwc[bottom] = Mtheta[bottom]*dz[bottom]
                        Mthetar[bottom] = np.minimum((np.ones_like(Mtheta[bottom])*0.02),0.9*Mtheta[bottom])
                        MeffSat[bottom] = (Mtheta[bottom]-Mthetar[bottom])/(Mtheta_sat[bottom]-Mthetar[bottom])
                        Mhead[bottom]  = -1*1/alpha_vG[bottom] * ((Sc[bottom] * MeffSat[bottom])**(-1/m_vG[bottom])-1)**(1/n_vG[bottom])


            ### We want to end precisely at the end of tstep ###
            if Mtime_counter + Mdeltat_new >= tstep:
                Mdeltat_new = tstep - Mtime_counter

        ##### Back in the global time loop #####
        timer += tstep # tstep seconds have passed
            
        ### Reset time counters and initial time steps (otherwise, might ==0 because we make sure we finish at tstep precisely) ###
        Mtime_counter = 0.
        Ptime_counter = 0.
        Mdeltat = 300.
        Mdeltat_new = 1*Mdeltat
        Pdeltat = 300.
        Pdeltat_new = 1*Pdeltat
         
        ##### Additional functions that are called at end of each tstep period #####

        ##### Proceed to refreezing (only in MFdom) #####
        # refreezing process according to cold content of every layer
        if timer == dtCFM: # If we are at the end of the flow routine, we allow for full refreezing
            lwc_min_fr = theta_min*dz # Thus we set minimal water amount that must remain after refreezing to theta_min
   
        rho,Tz,Mhead,Mtheta,Mthetar,MeffSat,Mlwc,Mtheta_sat,Ptheta,PeffSat,Plwc,Ptheta_sat,Ksat,theta_sat,alpha_vG,n_vG,m_vG,Sc,totrefrozen_lwc,refrozenlay,totrunoff = \
            Mrefreezing(dz,zstep,rho,grain,Tz,Mthetar_old,Mlwc,lwc_min_fr,Ptheta,PeffSat,Plwc,h_e,bigF,mu,crtn_theta,rhoimp,totrefrozen_lwc,refrozenlay,totrunoff)
    
        
        ## Runoff from Zuo and Oerlemans formula ##
        Mtheta,Mthetar,MeffSat,Mlwc,totrunoff = \
            runoff(dz,rho,Mhead,Mtheta,Mthetar,Mthetar_old,MeffSat,Mlwc,Mtheta_sat,theta_min_fr,crtn_theta,slope,rhoimp,noflow,tstep,totrunoff)
        Mhead  = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
    
    if tostore > 0: #If there is some water to store in the bottom layers (aquifer)
        Mlwc,Plwc,totrunoff = distribute_tostore_single(dz,rho,tostore,Mlwc,Plwc,rhoimp,bigF,totrunoff) # redistribute the water stored
        ### Update all variables (required for the refreezing) ###
        Mtheta = Mlwc/dz
        Mthetar = np.minimum((np.ones_like(Mtheta)*0.02),0.9*Mtheta)
        MeffSat = (Mtheta-Mthetar)/(Mtheta_sat-Mthetar)
        Mhead = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
        ### Proceed to refreezing ###
        ##Note that this has to be done after the runoff to avoid emptying aquifer
        rho,Tz,Mhead,Mtheta,Mthetar,MeffSat,Mlwc,Mtheta_sat,Ptheta,PeffSat,Plwc,Ptheta_sat,Ksat,theta_sat,alpha_vG,n_vG,m_vG,Sc,totrefrozen_lwc,refrozenlay,totrunoff = \
            Mrefreezing(dz,zstep,rho,grain,Tz,Mthetar_old,Mlwc,lwc_min_fr,Ptheta,PeffSat,Plwc,h_e,bigF,mu,crtn_theta,rhoimp,totrefrozen_lwc,refrozenlay,totrunoff)
        
    ### Check for mass conservation ###
    lwcerror_abs = sum(Mlwc)+sum(Plwc)+totrefrozen_lwc-totlwc0-totflux+totrunoff #total error of lwc with respect to initial state and fluxes
    lwcerror_rel = (lwcerror_abs)/(sum(Mlwc)+sum(Plwc)) #relative error of lwc with respect to initial state and fluxes
    if ((lwcerror_rel > 1e-4) and (lwcerror_abs > 1e-5)):
        massconsalarm = 1
    if ((lwcerror_rel > 1e-4) and (lwcerror_abs > 1e-5)):
        print('Spot-> Year, relative error and absolute error of lwc are:',self.modeltime[iii],lwcerror_rel,lwcerror_abs)
    
    ### Add the surface runoff to the total runoff (at end of routine so it does not interfere with mass balance calculations ###)
    totrunoff += surfrunoff
    
    ### Remove the prewetting ###
    Mlwc -= np.minimum(Mlwc,Mprewetting) # all Mprewetting are at most theta_min but make sure to avoid negative Mlwc
    
    r2 = grain**2 # reconversion of grain radius to squared grain radius
    Plwc_mem = 1*Plwc

    lwctot = Mlwc+Plwc

    if ((max(lwctot/dz)<crtn_theta) and (max(lwctot)>0)):
        layers_to_dry = np.where(lwctot>0)
        #print('Extradrying of [mWE]:',sum(lwctot[layers_to_dry[0]]))
        for index in layers_to_dry:
            lwctot[index] = 0.
            Plwc_mem[index] = 0.
    
    totlwccheck = sum(lwctot)
    
    ##### Set the variables back on the entire column #####
    rhoC,dzC,TzC,massC,lwcC,Plwc_memC,r2C,refrozenC = combineCFM(split_list,rho,dz,Tz,mass,lwctot,Plwc_mem,r2,refrozenlay)
    self.rho,self.dz,self.Tz,self.mass,self.LWC,self.PLWC_mem,self.r2 = lengthendom(self,rhoC,dzC,TzC,massC,lwcC,Plwc_memC,r2C)
    self.refrozen = np.append(refrozenC,np.zeros(len(self.dz)-len(refrozenC))) # in deep layers not involved in flow routine, no refreezing occured

    return self.rho,self.age,self.dz,self.Tz,self.z,self.mass,self.dzn,self.LWC,self.PLWC_mem,self.r2,self.refrozen,totrunoff
    
    
    
    
    
    


    
