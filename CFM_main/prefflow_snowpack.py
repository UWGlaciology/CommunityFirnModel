
'''
This is a code imitating the preferential flow scheme developed in Wever et al. (2016)
Water flow in firn with dual domain approach:
- Richards Equation in Matrix Flow domain and Preferential Flow domain

Few differences with Wever 2016 as:
-Use of a constant bigF value (part of the pore space allocated to each domain)
-Runoff function from Zuo and Oerlemans 1996
-Possible to apply PFfreezing if cold wave penetrates from the surface
-Possible to build up aquifer at end of domain
-Use of upstream weighted mean to determine hydraulic conductivity at interfaces avoids oscillations for large mesh size
-We use a changing bottom boundary if saturated layers accumulate: don't solve RE for the saturated layers at end of the firn column    
- don't solve RE in dry part of the domain

Density of last layer of the domain should always be >= 830
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
from fcts_snowpackflow import entrysuction
from fcts_snowpackflow import layerequaliser_eq
from fcts_snowpackflow import PFleave
from fcts_snowpackflow import PFleaveheat
from fcts_snowpackflow import Mrefreezing
from fcts_snowpackflow import Prefreezing
from fcts_snowpackflow import Psatexcess
from fcts_snowpackflow import Msatexcess
from fcts_snowpackflow import runoff
from fcts_snowpackflow import Micedryer
from fcts_snowpackflow import Picedryer
from fcts_snowpackflow import distribute_tostore

from merge import mergesurf

def prefflow(self,iii):
        
    ##### First: melting of the surface layers, taken from melt.py #####
    melt_volume_IE      = self.snowmeltSec[iii] * S_PER_YEAR # This still has to be checked by Max (division by self.c['stpsPerYear']?) [m]
    melt_volume_WE      = melt_volume_IE * RHO_I_MGM # [m]
    melt_mass           = melt_volume_WE * 1000. # [kg]
    
    liquid_befmelt = melt_volume_WE+sum(self.LWC) #VV tracking disappearing water
    
    ind1a               = np.where(self.mass_sum <= melt_mass)[0] # indices of boxes that will be melted away
    num_boxes_melted    = len(ind1a)+1 # number of boxes that melt away, include the box that is partially melted
    ind1                = np.where(self.mass_sum > melt_mass)[0][0] # index which will become the new surface
    # pm is the partial melt (the model volume that has a portion melted away)
    pm_mass             = self.mass_sum[ind1] - melt_mass # the remaining mass of the PM box [kg]
    pm_dz               = pm_mass / self.rho[ind1] # remaining thickness [m]
    pm_rho              = self.rho[ind1] # density of the PM box [kg/m3]
    pm_lwc              = self.LWC[ind1]/self.dz[ind1] * pm_dz # LWC of the PM box [m]
    melt_boxes_LWC_vol  = np.sum(self.LWC[0:ind1+1]) - pm_lwc #include the LWC from the boxes that melt (currently does not include from the partial melt box) [m]
    melt_boxes_LWC_mass = melt_boxes_LWC_vol * RHO_W_KGM #include the mass of LWC from the boxes that melt (currently does not include from the partial melt box) [kg]
    melt_mass_a         = melt_mass + melt_boxes_LWC_mass #total liq water from melted boxes(due to melting + LWC at previous time step) [kg]
    melt_vol_a          = melt_mass_a / RHO_W_KGM #total liq water from melted boxes(due to melting + LWC at previous time step) [m]
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
    split_list,rho,dz,Tz,mass,lwc,Plwc_mem,r2 = splitCFM(rhoC,dzC,TzC,massC,lwcC,Plwc_memC,r2C,vert_res)
    #print('We are at year:',self.modeltime[iii])
    
    ### Parameters we might want to change according to different runs ###
    tstep       = 900. # Frequency of exchange betwwen domains and of refreezing [s]
    slope       = 3.4907e-3 # slope that appears in  lateral runoff formula of Zuo and Oerlemans 1996, possibly in RE
    cst_melt    = 1 # constant meltwater input: set this to 1
    sinus_melt  = 0 # meltwater input according to sinusoidal daily cycle: set this to 1 !! Make sure dtCFM is 24h !!
    bigF        = 0.2*np.ones_like(dz) # close to observed values of Williams 2010 and references cited in there
    impermeability = 1. # 1 to activate impermeability of thick ice layers, 0 otherwise
    impthick    = 0.0 # minimum thickness of ice layers to be impermeable [m]
    rhoimp      = 810.
    rofflim     = 0.95 # effective saturation of MFdom above which we launch extra runoff to have MeffSat == rofflim
    PSatlim     = 0.1
    kth         = 0.021 + 2.5 * (rho/1000.)**2 # Thermal conductivity [W m-1 K-1] Make sure to use same as in CFM for consistency
    bigN        = 0.2 # tuning parameter of SNOWPACK where it is set to 0, nb of flow paths per m2 Wever 2016 (7)
    rhoPdr      = 1*rhoimp # density from which we dry out PFdom layers if we have a drying front at end of run to avoid keeping tiny amounts of water forever

    grain       = r2 ** (1/2) # grain size [m]
    
    Mlwc        = lwc-Plwc_mem
    Plwc        = Plwc_mem
    
    depth           = np.cumsum(dz)-dz/2 # depth of the center of every layer [m]
    depth_imin1     = np.delete(np.append(0,depth),-1) # depth staggered 1 level below [m]
    zstep           = depth-depth_imin1 # distance between centers of adjacent layers (between center[i] and center[i-1]), will be used for the PFdom [m]
    zstep_iplus1    = np.append(np.delete(zstep,0),dz[-1]/2) # zstep vector staggered 1 level lower
    
    ### Time parameters ###
    dtCFM           = self.dt #duratin of timesteps in CFM, constant value, [s]
    Mdeltat         = 300. #duration of timesteps for this script, choose an arbitrary starting value (<= tstep) [s]
    Mdeltat_new     = 300. #further, we adjust deltat_new and then assign this value to deltat, start with deltat_new == deltat [s]
    deltat_max      = 1*tstep #maximum time step allowed
    deltat_min      = 1e-20 #minimum time step allowed (1e-10 in D'Amboise 2017, not specified in Wever 2014) [s]
    Mtime_counter   = 0 # time for REMF, will be reset at 0 every tstep seconds
    Pdeltat         = 300. #duration of timesteps for flow in PFdom, this is adapted according to max value of Ptheta [s]
    Pdeltat_new     = 300.
    Ptime_counter   = 0 # time for PF, will be reset at 0 every tstep seconds
    timer           = 0 # this tracks the time passing --> remains between 0 and dtCFM [s]

    ### Calculate pore space available in every layer --> water content at saturation 
    porosity           = 1 - rho/RHO_I # Definition of porosity [/]
    porespace_vol      = porosity * dz # Pore space of each layer [m]
    porosity_refr      = porosity*RHO_I/RHO_W_KGM # space available for liq water volume once refrozen, Wever 2014 (9) [/]
    porosity_refr      = np.maximum(porosity_refr,17e-3) # allow space for minimum water content required in both domains for numerical stability, 17e-3 is equivalent to 900.0 density
    porespace_refr_vol = porosity_refr*dz # Available pore space of each layer [m]
    
    theta_sat   = porosity_refr # value of volumetric water content in saturated conditions [/]
    Mtheta_sat  = (1-bigF)*theta_sat
    Ptheta_sat  = bigF*theta_sat
    totrunoff   = 0.

    ### Convergence of the solution between successive iterations ###
    crtn_head   = 1e-3 #Criterion on head pressure, see Wever 2014 Appendix
    crtn_theta  = 1e-5 #Criterion on water content, see Wever 2014 Appendix
    crtn_MB     = 1e-8 # Mass balance criterion, as in SNOWPACK ReSolver line 576
    
    theta_min_fr    = 1e-4 # Threshold value to allow freezing, Wever 2014 Appendix
    lwc_min_fr      = theta_min_fr*dz # corresponding lwc threshold
    Mlwc_min_fr     = 1*lwc_min_fr 
    Plwc_min_fr     = 1*lwc_min_fr
    
    ### Use rain climatic inputy ###
    try:
       raintoadd = self.rainSec[iii] * S_PER_YEAR * RHO_I_MGM # [m]
    except:
        raintoadd = 0.

    melt_vol_a      += raintoadd # we add the rain to the melt volume, everything is in m we (total value for this time step)    
    
    totrefrozen_lwc = 0 # will calculate total lwc refrozen [mWE]
    refrozenlay     = np.zeros_like(dz) # will calculate total lwc refrozen in every layer [mWE]
    
    surfrunoff = 0. # all lateral runoff of this CFMstep  

    if rho[0] >= rhoimp:
        surfrunoff += melt_vol_a
        totinflux  = 0.
    elif rho[0] <= rhoimp:
        totinflux  = melt_vol_a # total surface input for this CFM timestep (should be double checked) [m]
        
    mean_influx  = totinflux/dtCFM # mean value of the liquid water flux at surface [m/s]
        
    liqout = 0. # water flowing out at the bottom of the domain [m/s]
        
     ### To calculate head pressure: van Genuchten parameters and correction factor
    alpha_vG = 4.4e6*(rho/(2*grain))**(-0.98) # Hirashima 2014 (5) ; typical value ~35.0 
    n_vG     = 1 + 2.7e-3*(rho/(2*grain))**(0.61) # Hirashima 2014 (6) ; typical value ~4.0
    m_vG     = 1 - 1/n_vG # Wever 2014 (8) ; typical value ~0.75
    h_e   = 5.8e-3 # air entry pressure for pore size of 5mm at 273K [m], Wever 2014 
    Sc    = (1 + (alpha_vG*h_e)**n_vG)**(-m_vG) #Saturation at cut-off point [/], see Ippisch et al., 2006 eq(11)

    h_we        = 0.0437/(2*grain*1e3) + 0.01074 # water entry suction for snow, Hirashima 2014 (15) [m]
    MSat_we     = (1+(alpha_vG*abs(h_we))**n_vG)**-m_vG / Sc # Effective saturation corresponding to the water entry suction
    MSat_westag = np.append(np.delete(MSat_we,0),1) # MSat_we values staggered, MSat_westag[i+1] == MSat_we[i]
    mu          = 0.001792 # Dynamic viscosity of water [kg m-1 s-1] , can be added to constants list
    Ksat        = RHO_W_KGM*GRAVITY/mu * 3.0*(grain)**2*np.exp(-0.013*rho) # Hydraulic conductivity at saturation (>0) [m s-1], Formula of Calonne et al. 2012, see Wever 2015 (7) and D'Amboise 2017 (10)

    theta_min   = crtn_theta/10*np.ones_like(dz) # minimum initial value of theta used for numerical stability
    ### MFdom variables ###
    Mtheta      = Mlwc/dz # Volumetric liquid water content [/]
    Mprewetting = dz*(theta_min-Mtheta) # amount of water required for prewetting of MFdom [m]
    Mprewetting = np.maximum(0,Mprewetting) # exclude negative values (where Mtheta already higher than theta_min)
    Mtheta      = np.maximum(Mtheta,theta_min) # modify theta
    Mlwc        = Mtheta*dz # modify lwc
    ## Define residual water content as in Wever 2014 ##
    Mthetar     = np.minimum((np.ones_like(Mtheta)*0.02),0.9*Mtheta) # initial residual water content [/], Wever 2014 (10)
    ## Calculate effective saturation ##
    MeffSat     = (Mtheta-Mthetar)/(Mtheta_sat-Mthetar) # effective saturation of the MFdom layers []
    
    ### PFdom variables ###
    Ptheta      = Plwc/dz # Volumetric liquid water content [/]
    Pprewetting = dz*(theta_min-Ptheta) # amount of water required for prewetting of PFdom [m]
    Pprewetting = np.maximum(0,Pprewetting) # exclude negative values (where Ptheta already higher than theta_min)
    Ptheta      = np.maximum(Ptheta,theta_min) # modify theta
    Plwc        = Ptheta*dz # modify lwc
    # residual water content is 0 in PFdom
    ## Calculate effective saturation ##
    PeffSat     = Ptheta/Ptheta_sat # effective saturation of the MFdom layers []

    totprewetting = Mprewetting+Pprewetting # total water added for numerical stability [m]

    totlwc0         = sum(Mlwc) + sum(Plwc)
    totflux         = 0
    massconsalarm   = 0
    provrunoff      = 0.

    ### Move water from layers where effSat is above 1 ###
    # Should be executed before head calculations, otherwise RuntimeWarning
    if np.any(Mtheta > Mtheta_sat): # Evacuate water in layers where we exceed maximal saturation towards other layers
        Mtheta,Mthetar,MeffSat,Mlwc,totrunoff = Msatexcess(dz,rho,Mtheta,Mtheta_sat,crtn_theta,rhoimp,totrunoff)
    
    if np.any(Mtheta[rho>rhoimp]>theta_min[rho>rhoimp]):
        Mtheta,Mthetar,MeffSat,Mlwc,totrunoff = Micedryer(dz,rho,Mtheta,Mtheta_sat,crtn_theta,rhoimp,totrunoff) # Evacuate water from layers of density > rhoimp
      
    if np.any(Ptheta > Ptheta_sat):
        Ptheta,PeffSat,Plwc,totrunoff = Psatexcess(dz,rho,Ptheta,Ptheta_sat,crtn_theta,rhoimp,totrunoff)
        
    ## Spot the aquifer #
    ice1    = np.zeros_like(dz) # bit array for ice layers
    ice1[np.where(rho>=rhoimp)[0]] = 1
    if np.any(ice1==0):
        noflow = np.where(ice1==0)[0][-1] + 1 # this layer and all layers below are saturated or are ice layers
    elif np.all(ice1>0): # If all the domain is ice or saturated
        noflow = 0 # all layers are saturated or are ice layers
     
    indaqM  = 1*noflow
    indaqP  = 1*noflow
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
    if np.any(PeffSat>0.95):
        icebottom = np.zeros_like(dz)
        icebottom[np.where(rho>=rhoimp)[0]] = 1 # ice layers
        Psatbottom = np.zeros_like(dz)
        Psatbottom[np.where(PeffSat>=0.95)[0]] = 1 # saturated layers
        Psaticebottom = Psatbottom+icebottom
        if np.any(Psaticebottom==0):
            indaqP = np.where(Psaticebottom==0)[0][-1]+1 #from this layer, all layers below are saturated or ice
            if indaqP > 0:
                if Ptheta[indaqP-1] > theta_min_fr and Ptheta[indaqP-2] <= theta_min_fr: # case where rest of the storage could not fill entirely a layer -> consider this layer as part of the aquifer
                    indaqP -= 1
        elif np.all(Psaticebottom>0):
            indaqP = 0
        tostore += sum((Ptheta[indaqP:]-theta_min[indaqP:])*dz[indaqP:]) # store all that water, which will be redistributed at end of flow routine
        ### Update all variables
        Ptheta[indaqP:] = 1*theta_min[indaqP:]
        Plwc[indaqP:] = Ptheta[indaqP:]*dz[indaqP:]
        PeffSat[indaqP:] = Ptheta[indaqP:]/Ptheta_sat[indaqP:]
    if tostore > 0:
        ### Consider the top of the aquifer as the lowest layer where we can have incoming flow
        noflow = max(indaqP,indaqM) # top of the aquifer is the highest layer saturated in both PFdom and MFdom    

    ### Head pressure ###
    Mhead  = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
    Phead  = -1*1/alpha_vG * ((Sc * PeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
    
    ## Avoid letting surface influx come in when the entire firn column is saturated
    Mthetar_old = Mthetar
    Pbottom = len(dz)-1
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
                #print('Bottom = 0 case, probably mass cons problem due to calculations, will be fixed shortly')
        elif np.all(icesat>0): # If all the domain is ice or saturated
            bottom = 0 # we will not proceed to RE
            if surfrunoff == 0.: # if the input has not been attributed to surface runoff (because surface layer is not at rhoimp)
                surfrunoff += mean_influx*dtCFM # [m]
                totinflux  = 0.
                mean_influx = 0.
                #print('All melt in runoff because all domain is ice or saturated, totinflux is:',totinflux)         
    
    #Assign old values
    Mtheta_old      = 1*Mtheta
    Mthetar_old     = 1*Mthetar
    MeffSat_old     = 1*MeffSat
    Mlwc_old        = 1*Mlwc
    Mhead_old       = 1*Mhead
    Mhead_old_imin1 = np.append(0,np.delete(Mhead,-1))
    Ptheta_old      = 1*Ptheta
    PeffSat_old     = 1*PeffSat
    Plwc_old        = 1*Plwc
    Phead_old       = 1*Phead
    Phead_old_imin1 = np.append(0,np.delete(Phead,-1))
    
    Mdtheta_dh      = np.zeros_like(dz)
    Mtheta_previous = np.zeros_like(dz)
    Mhead_previous  = np.zeros_like(dz)
    delta_Mtheta    = np.zeros_like(dz)
    delta_Mhead2    = np.zeros_like(dz)
    Pdtheta_dh      = np.zeros_like(dz)
    Ptheta_previous = np.zeros_like(dz)
    Phead_previous  = np.zeros_like(dz)
    delta_Ptheta    = np.zeros_like(dz)
    delta_Phead2    = np.zeros_like(dz)    
    
    ##### Start the Global time loop #####
    while timer < dtCFM:
        
        if timer + tstep >= dtCFM: # make sure we end precisely at the time step of the CFM
            tstep = dtCFM - timer
        
        bottom = max(0,min(bottom,noflow-1)) # don't proceed to RE in the aquifer
        ##### Time loop for MFdom ####
        sat_layers = np.array([])
        while Mtime_counter < tstep:
            Mflag       = 0 # Too many iterations or other cases where finer time step is required
            provrunoff  = 0.
                        
            if np.any(Mtheta<1e-10): # For some unknown reason, this happens on some rare occasions (only in the surface layer when layer below is an ice layer I think)
                instb = np.where(Mtheta<1e-10)[0]
                for ii in instb:
                    Mtheta[ii]  = theta_min[ii] # set theta back to numerical stability threshold and update all MFdom variables
                    Mthetar[ii] = 0.9*Mtheta[ii]
                    MeffSat[ii] = (Mtheta[ii]-Mthetar[ii])/(Mtheta_sat[ii]-Mthetar[ii])
                    Mlwc[ii]    = Mtheta[ii]*dz[ii]
                    Mhead[ii]   = -1*1/alpha_vG[ii] * ((Sc[ii] * MeffSat[ii])**(-1/m_vG[ii])-1)**(1/n_vG[ii]) # [m] Wever 2014 (3)
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
            
            Mdeltat         = 1*Mdeltat_new # Use the dynamically adjusted time step
            ### Assign old values ###
            Mtheta_old      = np.copy(Mtheta)
            Mthetar_old     = np.copy(Mthetar)
            MeffSat_old     = np.copy(MeffSat)
            Mlwc_old        = np.copy(Mlwc)
            Mhead_old       = np.copy(Mhead)
            Mhead_old_imin1 = np.append(0,np.delete(Mhead,-1))
            
            ### Update of Mthetar as in Wever 2014 ###
            Mthetar[0:(bottom+1)] = np.minimum((np.ones_like(Mtheta[0:(bottom+1)])*0.02),0.9*Mtheta[0:(bottom+1)])
            
            ### Update of effSat and head (because Mthetar might have changed) ##
            MeffSat[0:(bottom+1)]   = (Mtheta[0:(bottom+1)]-Mthetar[0:(bottom+1)])/(Mtheta_sat[0:(bottom+1)]-Mthetar[0:(bottom+1)]) # []
            Mhead[0:(bottom+1)]     = -1*1/alpha_vG[0:(bottom+1)] * ((Sc[0:(bottom+1)] * MeffSat[0:(bottom+1)])**(-1/m_vG[0:(bottom+1)])-1)**(1/n_vG[0:(bottom+1)]) # [m] Wever 2014 (3)
            
            
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
                Mlwc[0]     += Mliq_in*Mdeltat # We immediately give all the input of the time step to the surface layer
                Mtheta[0]   = Mlwc[0]/dz[0]
                Mthetar[0]  = np.minimum(0.02,0.9*Mtheta[0])
                MeffSat[0]  = (Mtheta[0]-Mthetar[0])/(Mtheta_sat[0]-Mthetar[0])
                Mhead[0]    = -1*1/alpha_vG[0] * ((Sc[0] * MeffSat[0])**(-1/m_vG[0])-1)**(1/n_vG[0]) # [m] Wever 2014 (3)
                Mliq_in     = 0. # There is no input to give anymore
                if ic == 0: # If all the layers are dry, we don't need to proceed to RE (water only in layer[0] and layer[1] is impermeable)
                    bottom = 0

                Mtheta_old[0]   = 1*Mtheta[0]
                Mthetar_old[0]  = 1*Mthetar[0]
                MeffSat_old[0]  = 1*MeffSat[0]
                Mlwc_old[0]     = 1*Mlwc[0]
                Mhead_old[0]    = 1*Mhead[0]
                
            Miteration   = 0 # Let's count the number of iterations required
    
            delta_Mhead  = np.ones_like(dz) # Difference in head pressure between two iterative steps, calculated by TDMA
            delta_Mhead2 = np.ones_like(dz) # Difference in head pressure between two iterative steps, directly calculated (in case head is artificially adjusted)
            delta_Mtheta = np.ones_like(dz) # Difference in volumetric liquid water content between two iterative steps
            MB_Merror    = np.ones_like(dz) # Mass balance error between discretised flux (with K values of last iteration) and changes calculated by the Picard scheme
            
            ### Spot the impermeable layers ###
            saturated   = np.where(MeffSat >= 0.95)[0]
            sat_layers  = np.unique(np.append(sat_layers,saturated))
            ice_layers  = np.where(rho>=rhoimp) # density threshold should be used here
                
            if impthick == 0.:
                thick_ice = ice_layers[0]
                
            elif impthick > 0.:            
                if len(ice_layers[0]) > 0: # If there are ice layers
                    thick_ice   = np.array([]) # this will contain all the indices of layers that are part of thick ice layers (thickness threshold is to be fixed)
                    start_ice   = np.array([]) # this will contain all the indices where an ice layer starts
                    end_ice     = np.array([]) # this will contain all the indices where an ice layer ends
                    cc          = ice_layers[0][0] # first ice layer

                    while cc <= ice_layers[0][-1]: # loop through all ice layers
                        if cc in ice_layers[0]:
                            start = cc # index of the start of the ice layer
                            while (cc+1 in ice_layers[0]) == True: # goes through all the layers belonging to an individual ice layers
                                cc += 1
                            end         = cc # index of the end of the ice layer
                            start_ice   = np.append(start_ice,start) # add the start of the layer to the start list
                            end_ice     = np.append(end_ice,end) # add the end of the layer to the end list
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
                dzdom           = dz[0:(bottom+1)]
                depthdom        = np.cumsum(dzdom)-dzdom/2 # depth of the center of every layer [m]
                depthdom_imin1  = np.delete(np.append(0,depthdom),-1) # depth staggered 1 level below [m]
                zstepdom        = depthdom-depthdom_imin1 # distance between centers of adjacent layers (between center[i] and center[i-1]), will be used for the PFdom [m]
                zstepdom_iplus1 = np.append(np.delete(zstepdom,0),dzdom[-1]/2) # zstep vector staggered 1 level lower         
                
                Mheaddown       = np.delete(Mhead,0) # Mhead values except last layer -> 1 element shorter than Mhead
                Mheadup         = np.delete(Mhead,-1) # Mhead values except first layer -> 1 element shorter than Mhead
                Mhead_imin1     = np.append(0,Mheadup) # Mhead vector staggered 1 level lower
                Mhead_iplus1    = np.append(Mheaddown,0) # Mhead vector staggered 1 level higher
                
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
                    #print('Mflag is 1!')
                    # time has not been incremented
                                    
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
                       
        
        ##### Time loop for PFdom ####

        Pliq_in = 0.
        Pliq_out = 0.
        
        if np.any(Ptheta>theta_min):
            # We proceed to RE in PFdom only if there is water in it (more than the water added for numerical stability)
            while Ptime_counter < tstep:
                Pflag = 0 # Too many iterations or other cases where finer time step is required
                
                if np.any(Ptheta<1e-10): # For some unknown reason, this happens on some rare occasions (only in the surface layer when layer below is an ice layer I think)
                    instb = np.where(Ptheta<1e-10)[0]
                    for ii in instb:
                        Ptheta[ii] = theta_min[ii] # set theta back to numerical stability threshold and update all MFdom variables
                        PeffSat[ii] = Ptheta[ii]/Ptheta_sat[ii]
                        Plwc[ii] = Ptheta[ii]*dz[ii]
                        Phead[ii] = -1*1/alpha_vG[ii] * ((Sc[ii] * PeffSat[ii])**(-1/m_vG[ii])-1)**(1/n_vG[ii]) # [m] Wever 2014 (3)
                        print('We had to fix an instability in PFdom, layers were:',instb)
                
                ### Assign old values ###
                Ptheta_old = np.copy(Ptheta)
                PeffSat_old = np.copy(PeffSat)
                Plwc_old = np.copy(Plwc)
                Phead_old = np.copy(Phead)
                Phead_old_imin1 = np.append(0,np.delete(Phead,-1))

                ice1 = np.zeros_like(dz) # bit array for ice layers
                satP = np.zeros_like(dz) # bit array for saturated layers in PFdom
                ice1[np.where(rho>=rhoimp)[0]] = 1 # ice layers
                satP[np.where(PeffSat>=0.95)[0]] = 1 # saturated layers in MFdom
                icesatP = ice1+satP # lowest layer with 0 value will be last layer taken into account in RE
                if np.any(icesatP==0):
                    Pbottom = np.where(icesatP==0)[0][-1] # this layer and all layers below are saturated or ice
                elif np.all(icesatP>0):
                    Pbottom = 0

                Pwetlayers = np.where(Ptheta[0:Pbottom] > (theta_min_fr + crtn_theta))[0] # all the layers having a significant theta value
                if np.size(Pwetlayers) == 0.: # in the case all our layers are dry
                    ic = 0 # last wet layer can be considered the surface layer as there will be some meltwater input
                elif np.size(Pwetlayers) > 0.: # in the case we have some layers already wet
                    ic = Pwetlayers[-1] # ic is the index of the deepest wet layer
                lowestdepth = depth[ic] + 0.5 # we take 0.5 meter below the lowest wet layer as the domain to perform calculations
                if lowestdepth < depth[-1]: # we want to define the maximum depth of the calculation domain
                    lc = np.where(depth >= lowestdepth)[0][0] # lc is the index of the last layer of the calculation domain
                elif lowestdepth >= depth[-1]: # the maximum depth of the calculation domain must not be greater than the maximum depth of the real domain
                    lc = len(dz)-1
                Pbottom = min(lc,Pbottom) # replace bottom by lc except in the case where we have ice layers at bottom that are then taken as limit of domain even if the limit just calculated goes further below
                
                Pbottom = max(0,min(Pbottom,noflow-1)) # test
                
                Pdeltat = 1*Pdeltat_new # Use the dynamically adjusted time step
                
                ### RE routine ###
            
                Piteration = 0 # Let's count the number of iterations required
                    
                delta_Phead  = np.ones_like(dz) # Difference in head pressure between two iterative steps, calculated by TDMA
                delta_Phead2 = np.ones_like(dz) # Difference in head pressure between two iterative steps, directly calculated (in case head is artificially adjusted)
                delta_Ptheta = np.ones_like(dz) # Difference in volumetric liquid water content between two iterative steps
                MB_Perror = np.ones_like(dz) # Mass balance error between discretised flux (with K values of last iteration) and changes calculated by the Picard scheme
                
                                
                while (np.max(np.append(0,abs(delta_Ptheta[0:Pbottom+1][PeffSat[0:Pbottom+1]<0.99]))) > crtn_theta or Piteration<2 or np.max(abs(MB_Perror[0:Pbottom+1]))>crtn_MB or np.max(np.append(0,abs(delta_Phead2[0:Pbottom+1][PeffSat[0:Pbottom+1]>=0.99]))) > crtn_head) and Pflag==0 and Pbottom>0:
                # we append a 0 value for the criterion check because this avoids error using max in case no single layer 
                ##### Iterative loop for RE in MFdom #####
                # 3 conditions: difference below convergence criterion, min 2 iterations, special case flag is not on
                    ### "center of layers = nodes" approach, as Wever 2014 and D'Amboise 2017 --> staggering ###
                    dzdom     = dz[0:(Pbottom+1)]
                    depthdom = np.cumsum(dzdom)-dzdom/2 # depth of the center of every layer [m]
                    depthdom_imin1 = np.delete(np.append(0,depthdom),-1) # depth staggered 1 level below [m]
                    zstepdom = depthdom-depthdom_imin1 # distance between centers of adjacent layers (between center[i] and center[i-1]), will be used for the PFdom [m]
                    zstepdom_iplus1 = np.append(np.delete(zstepdom,0),dzdom[-1]/2) # zstep vector staggered 1 level lower         
                    
                    Pheaddown      = np.delete(Phead,0) # Mhead values except last layer -> 1 element shorter than Mhead
                    Pheadup        = np.delete(Phead,-1) # Mhead values except first layer -> 1 element shorter than Mhead
                    Phead_imin1 = np.append(0,Pheadup) # Mhead vector staggered 1 level lower
                    Phead_iplus1   = np.append(Pheaddown,0) # Mhead vector staggered 1 level higher
                    
                    ### Hydraulic conductivity of every layer [Wever 2014 (11)] ###
                    PKdom = Ksat[0:(Pbottom+1)] * PeffSat[0:(Pbottom+1)]**0.5 * (1-(1-PeffSat[0:(Pbottom+1)]**(1/m_vG[0:(Pbottom+1)]))**m_vG[0:(Pbottom+1)])**2 # (>0) [m s-1]
                    #PKdom = 1*PK[0:(Pbottom+1)] # hydraulic saturation only on our calculation domain
                                    
                    ### Hydraulic conductivity at the interfaces ###                                
                    PKup         = np.delete(PKdom,-1) # K values except last layer -> 1 element shorter than K
                    PKdown       = np.delete(PKdom,0) # K values except first layer -> 1 element shorter than K
                    ## Upstream weighted mean for MK at interfaces: Forsyth 1995 (15), Szymkiewicz 2009 (3c)
                    PK_inter     = PKup # Take value of layer above
                    otherPK = np.where((Phead_old[1:(Pbottom+1)]-Phead_old_imin1[1:(Pbottom+1)])/dzdom[1:] -1 >0) # Szymkiewicz 2009 (3c)
                    PK_inter[otherPK[0]] = PKdown[otherPK[0]] # in these case, take value of layer below
                    
                    PKtop        = np.append(0,PK_inter) # Ktop will determine influx in every layer Ktop[0] is at surface layer (set to 0)
                    PKbot        = np.append(PK_inter,0) # Kbot will determine outflux in every layer Kbot[-1] is at bottom layer (set to 0)
                    
                    ### Calculate analytically derivative of theta with respect to head (Ci in Zarba 1988 (3.10) and Celia 1990 (17)) ###
                    Pdtheta_dh[0:(Pbottom+1)] = -1* (Ptheta_sat[0:(Pbottom+1)])/Sc[0:(Pbottom+1)] * (-m_vG[0:(Pbottom+1)]*n_vG[0:(Pbottom+1)]*alpha_vG[0:(Pbottom+1)]**n_vG[0:(Pbottom+1)]) * \
                        (abs(Phead[0:(Pbottom+1)]))**(n_vG[0:(Pbottom+1)]-1) / (1+(alpha_vG[0:(Pbottom+1)]*(abs(Phead[0:(Pbottom+1)])))**n_vG[0:(Pbottom+1)])**(m_vG[0:(Pbottom+1)]+1)
                    ##-1* as we use absolute value of head which is negative
                    
                    ### Flux calculations (part 1) ###
                    Pqbot = PKbot*((Phead[0:(Pbottom+1)]-Phead_iplus1[0:(Pbottom+1)])/zstepdom_iplus1+1) # flux at bottom of every layer [m/s]
                    Pqtop = PKtop*((Phead_imin1[0:(Pbottom+1)]-Phead[0:(Pbottom+1)])/zstepdom+1) # flux at top of every layer (simply staggered with respect to Mqbot) [m/s]
                    Pqtop[0] = Pliq_in # flux boundary condition [m/s]
                    Pqbot[-1] = Pliq_out # flux boundary condition [m/s]
                    
                    ### Parameters for the Tridiagonal Matrix Algorithm, see Zarba (1988) ###
                    #Caution: for Zarba: i=0 at bottom, here: i=0 at surface --> be consistent in signs and with staggered grids
                    #         alpha, beta, gamma not necessarily same (spatially speaking) as in Zarba 1988 but signs compensate
                    #         ex: alpha(CFM) = beta(ZARBA) (spatially) but head_stagdown(CFM) = hi+1(ZARBA)
            
                    # Intermediary calculations #
                    Palpha = PKtop/(dzdom*zstepdom)
                    Pbeta  = PKbot/(dzdom*zstepdom_iplus1)
                    Pgamma = (PKtop-PKbot)/dzdom
                    Pphi   = (Pdtheta_dh[0:(Pbottom+1)])/Pdeltat
                
                    # Diagonals of the matrix: a is sub, b is principal, c is sup #
                    a_diag = -1*Palpha
                    a_diag = np.delete(a_diag,0)
                    b_diag = Palpha+Pbeta+Pphi 
                    c_diag = -1*Pbeta
                    c_diag = np.delete(c_diag,-1)
                    
                    # Calculation of the residual (right-hand side of the equation) #
                    Rd = (Ptheta_old[0:(Pbottom+1)]-Ptheta[0:(Pbottom+1)])/Pdeltat - Palpha*(Phead[0:(Pbottom+1)]-Phead_imin1[0:(Pbottom+1)]) - Pbeta*(Phead[0:(Pbottom+1)]-Phead_iplus1[0:(Pbottom+1)]) + Pgamma 
                    
                    # Surface boundary condition: constant flux [Zarba (1988) Table 3.2]
                    Rd[0] = (Ptheta_old[0]-Ptheta[0])/Pdeltat - Pbeta[0]*(Phead[0]-Phead_iplus1[0]) - PKbot[0]/dz[0] + Pliq_in/dz[0]
                    # Bottom boundary condition: constant flux
                    Rd[Pbottom] = (Ptheta_old[Pbottom]-Ptheta[Pbottom])/Pdeltat - Palpha[Pbottom]*(Phead[Pbottom]-Phead_imin1[Pbottom]) + PKtop[Pbottom]/dzdom[Pbottom] - Pliq_out/dzdom[Pbottom]
                    
                    ### Solve tridiagonal matrix ###
                    #delta_Phead = NPtrid(a_diag,b_diag,c_diag,Rd)
                    delta_Phead = TDMAsolver(a_diag,b_diag,c_diag,Rd)
                    
                    ### Save values of head and theta at the last iterative step ###
                    Ptheta_previous = 1*Ptheta
                    Phead_previous = 1*Phead
                    
                    ### Assign new value to pressure head ###
                    Phead[0:(Pbottom+1)] = Phead[0:(Pbottom+1)] + delta_Phead #[m]
                    Ptheta[0:(Pbottom+1)] = (Ptheta_sat[0:(Pbottom+1)])*(1/Sc[0:(Pbottom+1)])*(1+(alpha_vG[0:(Pbottom+1)]*(abs(Phead[0:(Pbottom+1)])))**n_vG[0:(Pbottom+1)])**(-m_vG[0:(Pbottom+1)]) # []
                    
                    PeffSat[0:(Pbottom+1)] = (Ptheta[0:(Pbottom+1)])/(Ptheta_sat[0:(Pbottom+1)]) # []
                    Plwc[0:(Pbottom+1)]   = Ptheta[0:(Pbottom+1)]*dz[0:(Pbottom+1)] # [m]            
                    
                    delta_Ptheta[0:(Pbottom+1)] = Ptheta[0:(Pbottom+1)]-Ptheta_previous[0:(Pbottom+1)]
                    delta_Phead2[0:(Pbottom+1)] = Phead[0:(Pbottom+1)]-Phead_previous[0:(Pbottom+1)] # Useful if head artificially adjusted (eg to minimum value), normally delta_Mhead2==delta_Mhead
                    
                    ### Flux calculations (part 2) ###
                    change_Plwc = Plwc-Plwc_old # change in lwc over the time step [m]
                    MB_Perror = Pqbot*Pdeltat+change_Plwc[0:(Pbottom+1)]-Pqtop*Pdeltat # change in lwc not explained by flux discretisation [m]
                    totMB_Perror = sum(MB_Perror) # sum of the errors in whole column [m]
            
                    Piteration += 1 # count iterations
                    
                    ### Back step from Wever 2014 (don't use back step for MB criterion, let model iterate further!)
                    if Piteration == 16 or np.any(Phead[0:(Pbottom+1)]>0) or np.any(abs(delta_Phead2[0:(Pbottom+1)])>1000) or np.max(PeffSat[0:(Pbottom+1)])>1.0:#test
                                                
                        # Assign old values and we will start again with a shorter time step #
                        Ptheta[0:(Pbottom+1)] = 1*Ptheta_old[0:(Pbottom+1)]
                        PeffSat[0:(Pbottom+1)] = 1*PeffSat_old[0:(Pbottom+1)]
                        Plwc[0:(Pbottom+1)] = 1*Plwc_old[0:(Pbottom+1)]
                        Phead[0:(Pbottom+1)] = 1*Phead_old[0:(Pbottom+1)]
                        Pdeltat_new = Pdeltat/2
                        Pflag = 1
                        #print('Pflag is 1!')
                        # time has not been incremented
                                        
                ### According to number of iterations, adjust the next time step
                if Piteration <= 5 and Pflag == 0: # Increase next time step but don't exceed max time step value
                    Pdeltat_new = 1.25*Pdeltat
                    Pdeltat_new = min(deltat_max,Pdeltat_new)
                if Piteration > 5 and Piteration <= 10 and Pflag == 0: # Keep same time step
                    Pdeltat_new = 1*Pdeltat
                if Piteration > 10 and Piteration <= 15 and Pflag == 0: # Deccrease next time step but don't exceed min time step value
                    Pdeltat_new = Pdeltat/2    
                    Pdeltat_new = max(deltat_min,Pdeltat_new)            
                if Pflag == 1: #Now we are in the time loop --> don't go until the end of time loop and restart for the process with finer time step
                    Pdeltat_new = Pdeltat/2    
                    Pdeltat_new = max(deltat_min,Pdeltat_new)
                    if Pdeltat_new <= 10*deltat_min:
                        print('Too small Pdeltat, make the simulation crash, self.modeltime[iii] is:',self.modeltime[iii])
                        Ptheta_sat[0:10] = -1e9
                    continue # go back to top of Matrix time loop
                
                Ptime_counter += Pdeltat # increase the Mtime_counter
                
                if Pbottom == noflow-1: # If bottom is at the limit where we allow water to flow
                    if Ptheta[Pbottom] > theta_min_fr: # If water reached that limit -> store it in the aquifer water
                        # store the water to make effSat equal the effSat of overlying layer (no sharp wetting front)
                        Pwaterout = dz[Pbottom]*(Ptheta[Pbottom]-(PeffSat[Pbottom-1]*Ptheta_sat[Pbottom]))
                        Pwaterout = max(Pwaterout,0)
                        Pwaterout = min(Pwaterout,dz[Pbottom]*(Ptheta[Pbottom]-theta_min[Pbottom]))
                        if Pwaterout>0:
                            tostore += Pwaterout
                            Ptheta[Pbottom] -= Pwaterout/dz[Pbottom]
                            ### Update all variables
                            Plwc[Pbottom] = Ptheta[Pbottom]*dz[Pbottom]
                            PeffSat[Pbottom] = Ptheta[Pbottom]/Ptheta_sat[Pbottom]
                            Phead[Pbottom]  = -1*1/alpha_vG[Pbottom] * ((Sc[Pbottom] * PeffSat[Pbottom])**(-1/m_vG[Pbottom])-1)**(1/n_vG[Pbottom])
      
                ### We want to end precisely at the end of tstep ###
                if Ptime_counter + Pdeltat_new >= tstep:
                    Pdeltat_new = tstep - Ptime_counter
                    
        elif np.all(Ptheta<=theta_min):
            # We don't proceed to RE for PFdom
            Pbottom = 0
        
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

        ### Transfers between domains ###
        if np.any(MeffSat[0:noflow]>MSat_westag[0:noflow]):
            Mtheta,Mthetar,MeffSat,Mlwc,Ptheta,PeffSat,Plwc = \
                entrysuction(dz,Mtheta,Mthetar,Mthetar_old,MeffSat,Mtheta_sat,Mlwc,Ptheta,PeffSat,Plwc,Ptheta_sat,crtn_theta,noflow,MSat_westag)
            ## Recalculate pressure head values
            Mhead  = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
            Phead  = -1*1/alpha_vG * ((Sc * PeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
        
        if np.any(Mtheta>0.02):
            Mtheta,Mthetar,MeffSat,Mlwc,Ptheta,PeffSat,Plwc = \
                layerequaliser_eq(dz,Mtheta,Mthetar,Mthetar_old,MeffSat,Mtheta_sat,Mlwc,Ptheta,PeffSat,Plwc,Ptheta_sat,crtn_theta,noflow)
            Mhead  = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
            Phead  = -1*1/alpha_vG * ((Sc * PeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
        ### Issue of SNOWPACK equalising saturation even in dry PFdom layers ###
        # Assume we equalise at layer i -> flow to i+1 at same speed (head depends on effSat and both effSat are equal)
        # But in i+1 Mthetar will be increased -> from i+1 PF will be faster even if originally we only equalised saturations
        # Thus PF water will infiltrate PFdom in layers still dry in MFdom 
        # What is then the point of using the water entry suction criterion ?
        
        if np.any(PeffSat[0:noflow]>PSatlim):
            Mtheta,Mthetar,MeffSat,Mlwc,Ptheta,PeffSat,Plwc = \
                PFleave(dz,rho,Tz,Mtheta,Mthetar,Mthetar_old,MeffSat,Mtheta_sat,Mlwc,Ptheta,PeffSat,Plwc,Ptheta_sat,crtn_theta,rhoimp,noflow,PSatlim)
            ## Recalculate pressure head values
            Mhead  = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
            Phead  = -1*1/alpha_vG * ((Sc * PeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
                
        if (np.any(Tz[Ptheta>(crtn_theta/10+1e-5)][0:noflow]<273.15) and bigN>0): # add the 1e-5 term otherwise imprecision of float causes this to be always called
            deltatime = 1*tstep #just for a trial
            Mtheta,Mthetar,MeffSat,Mlwc,Ptheta,PeffSat,Plwc = \
                PFleaveheat(dz,rho,Tz,Mtheta,Mthetar,Mthetar_old,MeffSat,Mtheta_sat,Mlwc,Ptheta,PeffSat,Plwc,Ptheta_sat,crtn_theta,kth,bigF,bigN,noflow,rhoimp,deltatime)
            ## Recalculate pressure head values
            Mhead  = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
            Phead  = -1*1/alpha_vG * ((Sc * PeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
                
    
        ##### Proceed to refreezing (only in MFdom) #####
        # refreezing process according to cold content of every layer
        if timer == dtCFM: # If we are at the end of the flow routine, we allow for full refreezing
            lwc_min_fr = theta_min*dz # Thus we set minimal water amount that must remain after refreezing to theta_min
   
        rho,Tz,Mhead,Mtheta,Mthetar,MeffSat,Mlwc,Mtheta_sat,Ptheta,PeffSat,Plwc,Ptheta_sat,Ksat,theta_sat,alpha_vG,n_vG,m_vG,Sc,totrefrozen_lwc,refrozenlay,totrunoff = \
            Mrefreezing(dz,zstep,rho,grain,Tz,Mthetar_old,Mlwc,lwc_min_fr,Ptheta,PeffSat,Plwc,h_e,bigF,mu,crtn_theta,rhoimp,totrefrozen_lwc,refrozenlay,totrunoff)
    
        ## For refreezing in PFdom if cold wave propagates from surface  
        if np.any(Plwc>lwc_min_fr):
            if Tz[0] < 273.15:
                if np.any(Tz >= 273.15):
                    dryfront = np.where(Tz >= 273.15)[0][0] - 1
                elif np.all(Tz < 273.15):
                    dryfront = len(dz)-1
                if dryfront > 0:
                    #if dryfront>=np.where(Ptheta>lwc_min_fr)[0][0]:
                    if dryfront>=np.where(Plwc>lwc_min_fr)[0][0]:
                        if timer == dtCFM:
                            # Now we dry out high density (>900/rhoimp TBD) layers of PFdom water because otherwise, this water cannot refreeze due to lack of porespace and we keep forever tiny amounts of water in the firn column (Extradrying impossible)
                            Ptheta,PeffSat,Plwc,totrunoff = Picedryer(dz,rho,Ptheta,Ptheta_sat,crtn_theta,rhoPdr,totrunoff)
                        #print('Sum of Plwc until dryfront before refreezingPFdom():',sum(Plwc[0:dryfront]))
                        #print('Proceed to Prefreezing, dryfront,Plwc[590]:',dryfront,Plwc[590])
                        rho,Tz,Mhead,Mtheta,Mthetar,MeffSat,Mlwc,Mtheta_sat,Ptheta,PeffSat,Plwc,Ptheta_sat,Ksat,theta_sat,alpha_vG,n_vG,m_vG,Sc,totrefrozen_lwc,refrozenlay,totrunoff = \
                            Prefreezing(dz,rho,grain,Tz,Mthetar_old,Mtheta,Mlwc,lwc_min_fr,Ptheta,PeffSat,Plwc,bigF,h_e,mu,crtn_theta,dryfront,totrefrozen_lwc,refrozenlay,rhoimp,totrunoff)                    
                        #print('Sum of Plwc until dryfront after refreezingPFdom():',sum(Plwc[0:dryfront]))

        ## Required updates after refreezing ##
        MSat_we = (1+(alpha_vG*abs(h_we))**n_vG)**-m_vG / Sc # Effective saturation corresponding to the water entry suction
        MSat_westag = np.append(np.delete(MSat_we,0),1) # MSat_we values staggered, MSat_westag[i+1] == MSat_we[i]
        kth = 0.021 + 2.5 * (rho/1000.)**2 # Thermal conductivity [W m-1 K-1] Make sure to use same as in CFM for consistency

        ## Runoff from Zuo and Oerlemans formula ##
        Mtheta,Mthetar,MeffSat,Mlwc,totrunoff = \
            runoff(dz,rho,Mhead,Mtheta,Mthetar,Mthetar_old,MeffSat,Mlwc,Mtheta_sat,theta_min_fr,crtn_theta,slope,rhoimp,noflow,tstep,totrunoff)
        Mhead  = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
    
    
    if tostore > 0: #If there is some water to store in the bottom layers (aquifer)
        Mlwc,Plwc,totrunoff = distribute_tostore(dz,rho,tostore,Mlwc,Plwc,rhoimp,bigF,totrunoff) # redistribute the water stored
        ### Update all variables (required for the refreezing) ###
        Mtheta = Mlwc/dz
        Mthetar = np.minimum((np.ones_like(Mtheta)*0.02),0.9*Mtheta)
        MeffSat = (Mtheta-Mthetar)/(Mtheta_sat-Mthetar)
        Mhead = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)
        Ptheta = Plwc/dz
        PeffSat = Ptheta/Ptheta_sat
        Phead = -1*1/alpha_vG * ((Sc * PeffSat)**(-1/m_vG)-1)**(1/n_vG)
        ### Proceed to refreezing ###
        ##Note that this has to be done after the runoff to avoid emptying aquifer
        rho,Tz,Mhead,Mtheta,Mthetar,MeffSat,Mlwc,Mtheta_sat,Ptheta,PeffSat,Plwc,Ptheta_sat,Ksat,theta_sat,alpha_vG,n_vG,m_vG,Sc,totrefrozen_lwc,refrozenlay,totrunoff = \
            Mrefreezing(dz,zstep,rho,grain,Tz,Mthetar_old,Mlwc,lwc_min_fr,Ptheta,PeffSat,Plwc,h_e,bigF,mu,crtn_theta,rhoimp,totrefrozen_lwc,refrozenlay,totrunoff)
        ## For refreezing in PFdom if cold wave propagates from surface  
        if np.any(Plwc>lwc_min_fr):
            if Tz[0] < 273.15:
                if np.any(Tz >= 273.15):
                    dryfront = np.where(Tz >= 273.15)[0][0] - 1
                elif np.all(Tz < 273.15):
                    dryfront = len(dz)-1
                if dryfront > 0:
                    #if dryfront>=np.where(Ptheta>lwc_min_fr)[0][0]:
                    if dryfront>=np.where(Plwc>lwc_min_fr)[0][0]:
                        if timer == dtCFM:
                            # Now we dry out high density (>900/rhoimp TBD) layers of PFdom water because otherwise, this water cannot refreeze due to lack of porespace and we keep forever tiny amounts of water in the firn column (Extradrying impossible)
                            Ptheta,PeffSat,Plwc,totrunoff = Picedryer(dz,rho,Ptheta,Ptheta_sat,crtn_theta,rhoPdr,totrunoff)
                        rho,Tz,Mhead,Mtheta,Mthetar,MeffSat,Mlwc,Mtheta_sat,Ptheta,PeffSat,Plwc,Ptheta_sat,Ksat,theta_sat,alpha_vG,n_vG,m_vG,Sc,totrefrozen_lwc,refrozenlay,totrunoff = \
                            Prefreezing(dz,rho,grain,Tz,Mthetar_old,Mtheta,Mlwc,lwc_min_fr,Ptheta,PeffSat,Plwc,bigF,h_e,mu,crtn_theta,dryfront,totrefrozen_lwc,refrozenlay,rhoimp,totrunoff)                    
 
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
    Plwc -= np.minimum(Plwc,Pprewetting) # all Pprewetting are at most theta_min but make sure to avoid negative Plwc
        
    r2 = grain**2 # reconversion of grain radius to squared grain radius
    Plwc_mem = 1*Plwc

    lwctot = Mlwc+Plwc

    if ((max(lwctot/dz)<crtn_theta) and (max(lwctot)>0)):
        layers_to_dry = np.where(lwctot>0)
        #print('Extradrying of [mWE]:',sum(lwctot[layers_to_dry[0]]))
        for index in layers_to_dry:
            lwctot[index] = 0.
            Plwc_mem[index] = 0.
        
    ##### Set the variables back on the entire column #####
    rhoC,dzC,TzC,massC,lwcC,Plwc_memC,r2C,refrozenC = combineCFM(split_list,rho,dz,Tz,mass,lwctot,Plwc_mem,r2,refrozenlay)

    self.rho,self.dz,self.Tz,self.mass,self.LWC,self.PLWC_mem,self.r2 = lengthendom(self,rhoC,dzC,TzC,massC,lwcC,Plwc_memC,r2C)
    self.refrozen = np.append(refrozenC,np.zeros(len(self.dz)-len(refrozenC))) # in deep layers not involved in flow routine, no refreezing occured
    
    return self.rho,self.age,self.dz,self.Tz,self.z,self.mass,self.dzn,self.LWC,self.PLWC_mem,self.r2,self.refrozen,totrunoff
    
    
    
    
    
    
    
    


    
