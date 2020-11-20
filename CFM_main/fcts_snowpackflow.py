#!/usr/bin/env python
'''
This script contains all the functions required to make the pereferential flow scheme of snowpack work
https://models.slf.ch/p/snowpack/source/tree/HEAD/branches/dev/snowpack/snowpackCore/ReSolver1d.cc
'''
import numpy as np
from constants import *
from scipy.sparse import diags as diags
from scipy.linalg import solve as solve
from scipy.linalg import solve_banded as solve_banded
import math

from constants import *

def NPtrid(a,b,c,d):
    '''
    function to solve tridiagonal matrix --> find x matrix in Ax=d equation

    a,b,c,d are arrays, should be made of floats, not integer
    a and c are 1 index shorter than b and d
    -   [b0 c0 0. 0. ...]
    -   [a0 b1 c1 0. ...]
    A = [...............]
    -   [. 0. 0. an-1 bn]
    '''

    diagonals  = [a, b, c]
    tridiagmat = diags(diagonals,[-1,0,1]).toarray()
    tridiagsol = solve(tridiagmat,d)
    
    return tridiagsol

def TDMAsolver(a,b,c,d):
    '''
    Way faster than NPtrid!!
    
    TDMA solver --> find x matrix in Ax=d equation
    a,b,c,d are arrays, should be made of floats, not integer
    a and c are 1 index shorter than b and d
    -   [b0 c0 0. 0. ...]
    -   [a0 b1 c1 0. ...]
    A = [...............]
    -   [. 0. 0. an-1 bn]    
    A is nxn tridiagonal matrix with a - b - c as diagonals, x is nx1 matrix, d is nx1 matrix
    Sources:
    https://gist.github.com/cbellei/8ab3ab8551b8dfc8b081c518ccd9ada9
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm#Python
    '''
    
    n = len(d)
    ac,bc,cc,dc = map(np.array, (a,b,c,d))
    flag = 0.0 # In case there would be a division by 0 at a certain point
    if bc[n-1] == 0.: # division by 0
            flag = 1
    for it in range(1,n):
        if bc[it-1] == 0.: # division by 0
            flag = 1
        
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1]
        dc[it] = dc[it] - mc*dc[it-1]
        
    xc = bc
    xc[-1] = dc[-1]/bc[-1]
    
    for il in range(n-2,-1,-1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]
        
    if flag == 0:
        return xc
    else:
        return -1

def splitCFM(rhoC,dzC,TzC,massC,lwcC,Plwc_memC,r2C,vert_res):
    '''
    F for fine grid
    C for coarse grid
    We split the layers of the CFM in layers of a certain maximal thickness.
    Maximal thickness must be specified in vert_res.
    With the implementation of upstream weighted mean for K at interfaces, we can take large vert_res value!
    '''
    
    ### Resolution we want for the layers for RE, should maybe be given as input to funcion ###
    ### Number of sublayers vector ###
    split_list  = np.array([]) # list that contains nb of sublayers of every layer (index corresponds to index of layer in CFM)
    
    ### Vectors of the variables that will be used in Fine grid ###
    dzF         = np.array([]) # thickness vector for the F
    massF       = np.array([]) # mass vector for the F 
    LWCF        = np.array([]) # LWC vector for the F
    rhoF        = np.array([]) # density vector for the F
    TzF         = np.array([]) # temperature vector for the F
    r2F         = np.array([]) # squared radius vector for the F
    PLWC_memF   = np.array([]) # PLWC_mem vector for the flow routine
    # must not be executed for refrozen as refrozen is specific to each CFM step -> starts as zero everywhere at every flow routine
    
    for ii in range(len(dzC)):
        
        if dzC[ii] > vert_res: # Cases where the layer of the CFM is thicker than vert_res
            nb_sublayers = math.ceil(dzC[ii]/vert_res) # nb of layers in which each CFM value is split: we need layers of vert_res thickness at max (round value to upper integer)
            split_list = np.append(split_list,nb_sublayers) # fill in split_list
            
            ### Define new values of the variables: caution to the difference between additive and non-additive ###
            # Additive variables: value of CFM layer repartitioned between sublayers #
            newdz       = dzC[ii]/nb_sublayers*np.ones(nb_sublayers) # thickness value of the layers for the F (will be close to vert_res)
            newmass     = massC[ii]/nb_sublayers*np.ones(nb_sublayers) # mass value of the layers for the F
            newLWC      = lwcC[ii]/nb_sublayers*np.ones(nb_sublayers) # LWC value of the layers for the F
            newPLWC_mem = Plwc_memC[ii]/nb_sublayers*np.ones(nb_sublayers) # PLWC_mem value of the layers for the F

            # Non-additive variables: every sublayer keeps the same value as the CFM layer #
            newrho      = rhoC[ii]*np.ones(nb_sublayers) # density value of the layers for the F
            newTz       = TzC[ii]*np.ones(nb_sublayers) # temperature value of the layers for the F
            newr2       = r2C[ii]*np.ones(nb_sublayers) # grain size value of the layers for the F
            
            ### Fill in the vectors that will be used in the RE routine ###
            dzF         = np.concatenate((dzF,newdz)) # thickness [m]
            massF       = np.concatenate((massF,newmass)) # mass [kg]
            LWCF        = np.concatenate((LWCF,newLWC)) # LWC [m]
            rhoF        = np.concatenate((rhoF,newrho)) # density [kg m-3]
            TzF         = np.concatenate((TzF,newTz)) # temperature [K]
            r2F         = np.concatenate((r2F,newr2)) # squared radius [m2]
            PLWC_memF   = np.concatenate((PLWC_memF,newPLWC_mem)) # LWC [m]
            
        if dzC[ii] <= vert_res: # Cases where the layer of the CFM is not thicker than vert_res -> does not need to be split in sublayers
            nb_sublayers    = 1 # Only one sublayer per layer of CFM
            split_list      = np.append(split_list,nb_sublayers) # fill in split_list with a value of 1 
            
            ### Store immediately the CFM values in the vectors for F ###
            dzF         = np.append(dzF,dzC[ii]) # thickness [m]
            massF       = np.append(massF,massC[ii]) # mass [kg]
            LWCF        = np.append(LWCF,lwcC[ii]) # LWC [m]
            rhoF        = np.append(rhoF,rhoC[ii]) # density [kg m-3]
            TzF         = np.append(TzF,TzC[ii]) # temperature [K]
            r2F         = np.append(r2F,r2C[ii]) # squared radius [m2]
            PLWC_memF   = np.append(PLWC_memF,Plwc_memC[ii]) # LWC [m]
            
    return split_list,rhoF,dzF,TzF,massF,LWCF,PLWC_memF,r2F

def combineCFM(split_list,rhoF,dzF,TzF,massF,lwcF,Plwc_memF,r2F,refrozenF):
    '''
    F for fine grid
    C for coarse grid
    Here, we need the vectors that are outputs of the flow routine:
    - combine them back to reintegrate these to the CFM
    - Note that we also need the split_list ! (nb of sublayers into which every CFM layer was split by the split function)
    - Normally, RE routine (including freezing) should not affect dz and r2 variables but only use them
    -> not necessary to combine them and to give them back to CFM (which can keep its own self.dz and self.r2)
    '''
    
    ### Vectors that will be reattributed to the self. vectors ###
    dzC         = np.array([]) # thickness vector for the C
    massC       = np.array([]) # mass vector for the C 
    lwcC        = np.array([]) # LWC vector for the C
    Plwc_memC   = np.array([]) # PLWC_mem vector for the C
    rhoC        = np.array([]) # density vector for the C
    TzC         = np.array([]) # temperature vector for the C
    r2C         = np.array([]) # squared radius vector for the C
    refrozenC   = np.array([]) # refrozen water amount per layer vector for the C
    
    jj = 0 # index in split_list
    
    for number in split_list: # number stands for the number of sublayers (the individual values of the split_list vector)
        jjlast = int(jj+number) # last index in the RE vectors that corresponds to the same CFM layer as the jjth index in RE vectors
        
        ### Combine values of the variables: caution to the difference between additive and non-additive ###
        # Additive values: sum the values of all the sublayers that we combine
        combdz       = np.sum(dzF[jj:jjlast]) # combining the dz values from same C layer (dzF[jj] == dzF[jj+1] == dzF[jj+2] == ... == dzF[jjlast])
        combmass     = np.sum(massF[jj:jjlast]) # combining the mass values from same C layer
        comblwc      = np.sum(lwcF[jj:jjlast]) # combining the LWC values from same C layer
        combPlwc_mem = np.sum(Plwc_memF[jj:jjlast]) # combining the PLWC_mem values from same C layer
        combrefrozen = np.sum(refrozenF[jj:jjlast]) # combining the refrozen values from same C layer
        # Non-additive variables: take the mean of the values of all the sublayers that we combine
        combrho      = np.mean(rhoF[jj:jjlast])
        combTz       = np.mean(TzF[jj:jjlast])
        combr2       = np.mean(r2F[jj:jjlast])
        
        ### Fill in the vectors that will be given back to the C ###
        dzC         = np.append(dzC,combdz) # thickness vector for dz [m]
        massC       = np.append(massC,combmass) # mass vector for mass [kg]
        lwcC        = np.append(lwcC,comblwc) # LWC vector for LWC [m]
        Plwc_memC   = np.append(Plwc_memC,combPlwc_mem) # PLWC_mem vector for PLWC_mem [m]
        rhoC        = np.append(rhoC,combrho) # density vector for rho [kg m-3]
        TzC         = np.append(TzC,combTz) # temperature vector for Tz [K]
        r2C         = np.append(r2C,combr2) # squared radius vector for r2 [m2]
        refrozenC   = np.append(refrozenC,combrefrozen) # refrozen vector for refrozen [mWE]
        
        jj = jjlast
        
    return rhoC, dzC, TzC, massC, lwcC, Plwc_memC, r2C, refrozenC

def restrictdom(self):
    '''
    Limit the domain to reduce computational time of flow solving
    Spot the deepest layer below pore close-off density (rhoimp). The bottom of the domain will be the next layer.
    Thus, we do not take all the column where rho is constantly above pore close-off density (rhoimp).
    But we make sure not to exclude any layer that has liquid water
    '''
    rhobottom = 830.
    
    iirho = 0
    iilwc = 0

    if np.any(self.rho<rhobottom):        
        iirho = np.where(self.rho<rhobottom)[0][-1] # Last layer below 830
    if np.any(self.LWC>0):
        iilwc = np.where(self.LWC>0)[0][-1] # Last layer with water content
    ii = max(iirho,iilwc) # Last layer below 830 or with water content
    
    # ii is now the last layer where rho is below rhobottom kg/m3 -> limit the domain until ii+1
    rho_short       = self.rho[0:ii+2]
    dz_short        = self.dz[0:ii+2]
    Tz_short        = self.Tz[0:ii+2]
    mass_short      = self.mass[0:ii+2]
    lwc_short       = self.LWC[0:ii+2]
    r2_short        = self.r2[0:ii+2]
    Plwc_mem_short  = self.PLWC_mem[0:ii+2]
    
    return(rho_short,dz_short,Tz_short,mass_short,lwc_short,Plwc_mem_short,r2_short)
    
def lengthendom(self,rho_short,dz_short,Tz_short,mass_short,lwc_short,Plwc_mem_short,r2_short):
    '''
    After having solved the flow, we concatenate back new values with deep values that were not taken into account 
    during the flow routine (because all their rho values was above rhoimp from a certain depth).
    !! This function has to be called before the flow routine but AFTER the melting !!
    '''
    ii              = len(dz_short)-1 # last layer that may have been affected by flow routine
    rho_full        = np.append(rho_short,self.rho[ii+1:]) # modify the old variables that were defined on the entire column
    dz_full         = np.append(dz_short,self.dz[ii+1:])
    Tz_full         = np.append(Tz_short,self.Tz[ii+1:])
    mass_full       = np.append(mass_short,self.mass[ii+1:])
    lwc_full        = np.append(lwc_short,self.LWC[ii+1:])
    r2_full         = np.append(r2_short,self.r2[ii+1:])
    Plwc_mem_full   = np.append(Plwc_mem_short,self.PLWC_mem[ii+1:])
    
    return(rho_full,dz_full,Tz_full,mass_full,lwc_full,Plwc_mem_full,r2_full)


def Msatexcess(dz,rho,Mtheta,Mtheta_sat,crtn_theta,rhoimp,totrunoff):
    '''
    For MFdom
    In case the water content of some layers exceeds the water content at saturation, we move the water to
    the layers below (in priority) and in the layers above (if there is not enough pore space in all the layers below)
    '''
    waterexcess = np.where(Mtheta > Mtheta_sat) # spot layers where we exceed saturation
    ice1        = np.zeros_like(dz)
    sat1        = np.zeros_like(dz)

    ice1[np.where(rho>=rhoimp)[0]] = 1
    sat1[np.where(Mtheta/Mtheta_sat>=0.95)[0]] = 1

    icesat      = ice1+sat1

    if np.any(icesat==0):
        lowest  = np.where(icesat == 0)[0][-1]
    elif np.all(icesat>0):
        lowest  = 0
    
    for index in waterexcess[0]:
        lwcexc  = 1.001*(Mtheta[index]-Mtheta_sat[index]) * dz[index] # move excess of water, with safety margin
        lwcexc  = min((Mtheta[index]-crtn_theta/10)*dz[index],lwcexc) # we still need minimum theta for numerical stability
        Mtheta[index] -= lwcexc/dz[index] # we remove that excess of water but still have to distribute it in the column
        tobelow = 1 # we first try to distribute in layers situated below
        bb      = 1

        if (index+bb)>lowest or (index+bb>len(dz)-1): # if there are no layers below, no below distribution
            tobelow = 0
        while lwcexc > 0. and tobelow == 1: # as long as there is excess to distribute and we did not reach bottom
            if rho[index+bb]<rhoimp: # We only transfer water in layers below pore close-off density
                transf = np.minimum(lwcexc,(0.98*(Mtheta_sat[index+bb]-Mtheta[index+bb])*dz[index+bb])) # do not oversaturate receiving layers, safety margin
                transf = np.maximum(transf,0.) # make sure not to have negative values
            elif rho[index+bb]>=rhoimp: # No transfer of water in layers above pore close-off density
                transf = 0.
            Mtheta[index+bb]    += transf/dz[index+bb] # add the water
            bb                  += 1 # go to layer below
            lwcexc              -= transf # part of lwcexc has been distributed
            if (index+bb)>lowest or (index+bb>len(dz)-1): # if we reach bottom, stop below distribution
                tobelow = 0
        
        toabove = 1 # if there is still some lwcexc to distribute but below distribution not possible anymore, same process in layers above
        #aa = 1
        aa = np.maximum(1,index-lowest) # start to look for space above the aquifer
        if index-aa < 0:
            toabove = 0
        while lwcexc > 0. and toabove == 1:
            if rho[index-aa]<rhoimp: # We only transfer water in layers below pore close-off density
                transf = np.minimum(lwcexc,(0.98*(Mtheta_sat[index-aa]-Mtheta[index-aa])*dz[index-aa]))
                transf = np.maximum(transf,0.)
            elif rho[index-aa]>=rhoimp: # No transfer of water in layers above pore close-off density
                transf = 0.
            Mtheta[index-aa] += transf/dz[index-aa]
            aa += 1
            lwcexc -= transf
            if index-aa < 0:
                toabove = 0
        if lwcexc > 0: # if we could not distribute all the excess
            totrunoff += lwcexc # it is added to runoff
            #print('Problem: some layers exceed saturation but nowhere to store it -> add it to runoff?')
            
    ### Update of thetar ###
    Mthetar = np.minimum((np.ones_like(Mtheta)*0.02),0.9*Mtheta) # initial residual water content [/], Wever 2014 (10)
    
    ### Update of effSat and lwc 
    MeffSat = (Mtheta-Mthetar)/(Mtheta_sat-Mthetar)
    Mlwc = Mtheta*dz
    
    return Mtheta,Mthetar,MeffSat,Mlwc,totrunoff

def Psatexcess(dz,rho,Ptheta,Ptheta_sat,crtn_theta,rhoimp,totrunoff):
    '''
    This works the same as Msatexcess but for PFdom
    In case the water content of some layers exceeds the water content at saturation, we move the water to
    the layers below (in priority) and in the layers above (if there is not enough pore space in all the layers below)
    '''
    waterexcess = np.where(Ptheta > Ptheta_sat) # spot layers where we exceed saturation
    ice1        = np.zeros_like(dz)
    sat1        = np.zeros_like(dz)

    ice1[np.where(rho>=rhoimp)[0]] = 1
    sat1[np.where(Ptheta/Ptheta_sat>=0.95)[0]] = 1

    icesat      = ice1+sat1
    if np.any(icesat==0):
        lowest  = np.where(icesat == 0)[0][-1]
    elif np.all(icesat>0):
        lowest  = 0
        
    for index in waterexcess[0]:
        lwcexc        = 1.001*(Ptheta[index]-Ptheta_sat[index]) * dz[index] # move excess of water, with safety margin
        lwcexc        = min((Ptheta[index]-crtn_theta/10)*dz[index],lwcexc) # we still need minimum theta for numerical stability
        Ptheta[index] -= lwcexc/dz[index] # we remove that excess of water but still have to distribute it in the column
        tobelow       = 1 # we first try to distribute in layers situated below
        bb            = 1
        if (index+bb)>lowest or (index+bb>len(dz)-1): # if there are no layers below, no below distribution
            tobelow = 0
        while lwcexc > 0. and tobelow == 1: # as long as there is excess to distribute and we did not reach bottom
            if rho[index+bb]<rhoimp: # We only transfer water in layers below pore close-off density
                transf = np.minimum(lwcexc,(0.99*(Ptheta_sat[index+bb]-Ptheta[index+bb])*dz[index+bb])) # do not oversaturate receiving layers, safety margin
                transf = np.maximum(transf,0.) # make sure not to have negative values
            elif rho[index+bb]>=rhoimp: # No transfer of water in layers above pore close-off density
                transf = 0.
            Ptheta[index+bb] += transf/dz[index+bb] # add the water
            bb += 1 # go to layer below
            lwcexc -= transf # part of lwcexc has been distributed
            if index+bb > len(dz)-2: # if we reach bottom, stop below distribution + don't allow water transfer in last layer (supposed to be surface of ice sheet)
                tobelow = 0
        
        toabove = 1 # if there is still some lwcexc to distribute but below distribution not possible anymore, same process in layers above
        #aa = 1
        aa = np.maximum(1,index-lowest) # start to look for space above the aquifer
        if index-aa < 0:
            toabove = 0
        while lwcexc > 0. and toabove == 1:
            if rho[index-aa]<rhoimp: # We only transfer water in layers below pore close-off density
                transf = np.minimum(lwcexc,(0.99*(Ptheta_sat[index-aa]-Ptheta[index-aa])*dz[index-aa]))
                transf = np.maximum(transf,0.)
            elif rho[index-aa]>=rhoimp: # No transfer of water in layers above pore close-off density
                transf = 0.
            Ptheta[index-aa] += transf/dz[index-aa]
            aa += 1
            lwcexc -= transf
            if index-aa < 0:
                toabove = 0
        if lwcexc > 0: # if we could not distribute all the excess
            totrunoff += lwcexc # it is added to runoff
            #print('Problem: some layers exceed saturation but nowhere to store it -> add it to runoff?')
            
    ### Update of effSat and lwc 
    PeffSat = Ptheta/Ptheta_sat
    Plwc = Ptheta*dz
    
    return Ptheta,PeffSat,Plwc,totrunoff 

def Micedryer(dz,rho,Mtheta,Mtheta_sat,crtn_theta,rhoimp,totrunoff):
    '''
    This works the same as satexcess but in order to make sure that layers at pore close-off density (rho>rhoimp) are dry in MFdom.
    '''
    icelay      = np.where(rho>rhoimp)[0] # spot ice layer
    wetlay      = np.where(Mtheta>crtn_theta/10)[0] # spot layers above minimum saturation
    ice_to_dry  = np.intersect1d(icelay,wetlay) # spot ice layers above minimum saturation
    
    ice1    = np.zeros_like(dz)
    sat1    = np.zeros_like(dz)
    ice1[np.where(rho>=rhoimp)[0]] = 1
    sat1[np.where(Mtheta/Mtheta_sat>=0.95)[0]] = 1
    icesat  = ice1+sat1

    if np.any(icesat==0):
        lowest = np.where(icesat == 0)[0][-1]
    elif np.all(icesat>0):
        lowest = 0
    
    for index in ice_to_dry:
        lwcexc          = (Mtheta[index]-crtn_theta/10) * dz[index] # move excess of water, with safety margin
        Mtheta[index]   -= lwcexc/dz[index] # we remove that excess of water but still have to distribute it in the column
        tobelow         = 1 # we first try to distribute in layers situated below
        bb              = 1

        if (index+bb)>lowest or (index+bb>len(dz)-1): # if there are no layers below, no below distribution
            tobelow = 0
        while lwcexc > 0. and tobelow == 1: # as long as there is excess to distribute and we did not reach bottom
            if rho[index+bb]<rhoimp:
                transf = np.minimum(lwcexc,(0.99*(Mtheta_sat[index+bb]-Mtheta[index+bb])*dz[index+bb])) # do not oversaturate receiving layers, safety margin
                transf = np.maximum(transf,0.) # make sure not to have negative values
            elif rho[index+bb]>=rhoimp:
                transf = 0.
            Mtheta[index+bb]    += transf/dz[index+bb] # add the water
            bb                  += 1 # go to layer below
            lwcexc              -= transf # part of lwcexc has been distributed
            if index+bb > len(dz)-1: # if we reach bottom, stop below distribution
                tobelow = 0
        
        toabove = 1 # if there is still some lwcexc to distribute but below distribution not possible anymore, same process in layers above
        aa = np.maximum(1,index-lowest) # start to look for space above the aquifer
        if index-aa < 0:
            toabove = 0
        while lwcexc > 0. and toabove == 1:
            if rho[index-aa]<rhoimp:
                transf = np.minimum(lwcexc,(0.99*(Mtheta_sat[index-aa]-Mtheta[index-aa])*dz[index-aa]))
                transf = np.maximum(transf,0.)
            elif rho[index-aa]>=rhoimp:
                transf = 0.
            Mtheta[index-aa] += transf/dz[index-aa]
            aa += 1
            lwcexc -= transf
            if index-aa < 0:
                toabove = 0
        if lwcexc > 0: # if we could not distribute all the excess
            totrunoff += lwcexc # it is added to runoff
    
    ### Update of thetar ###
    Mthetar = np.minimum((np.ones_like(Mtheta)*0.02),0.9*Mtheta) # initial residual water content [/], Wever 2014 (10)
    
    ### Update of effSat and lwc 
    MeffSat = (Mtheta-Mthetar)/(Mtheta_sat-Mthetar)
    Mlwc = Mtheta*dz
    
    return Mtheta,Mthetar,MeffSat,Mlwc,totrunoff

def Picedryer(dz,rho,Ptheta,Ptheta_sat,crtn_theta,rhoPdr,totrunoff):
    '''
    This works the same as satexcess but in order to make sure that layers at pore close-off density (rho>rhoimp) are dry in MFdom.
    '''
    icelay      = np.where(rho>rhoPdr)[0] # spot ice layer
    wetlay      = np.where(Ptheta>crtn_theta/10)[0] # spot layers above minimum saturation
    ice_to_dry  = np.intersect1d(icelay,wetlay) # spot ice layers above minimum saturation
    
    ice1        = np.zeros_like(dz)
    sat1        = np.zeros_like(dz)
    ice1[np.where(rho>=rhoPdr)[0]] = 1
    sat1[np.where(Ptheta/Ptheta_sat>=0.95)[0]] = 1
    icesat      = ice1+sat1
    if np.any(icesat==0):
        lowest = np.where(icesat == 0)[0][-1]
    elif np.all(icesat>0):
        lowest = 0
    
    for index in ice_to_dry:
        lwcexc = (Ptheta[index]-crtn_theta/10) * dz[index] # move excess of water, with safety margin
        Ptheta[index] -= lwcexc/dz[index] # we remove that excess of water but still have to distribute it in the column
        if np.any(rho[index:]<rhoPdr):
            tobelow = 1 # we first try to distribute in layers situated below
            bb = 1
            if (index+bb)>lowest or (index+bb>len(dz)-1): # if there are no layers below, no below distribution
                tobelow = 0
            while lwcexc > 0. and tobelow == 1: # as long as there is excess to distribute and we did not reach bottom
                if rho[index+bb]<rhoPdr:
                    transf = np.minimum(lwcexc,(0.9*(Ptheta_sat[index+bb]-Ptheta[index+bb])*dz[index+bb])) # do not oversaturate receiving layers, safety margin
                    transf = np.maximum(transf,0.) # make sure not to have negative values
                elif rho[index+bb]>=rhoPdr:
                    transf = 0.
                Ptheta[index+bb] += transf/dz[index+bb] # add the water
                bb += 1 # go to layer below
                lwcexc -= transf # part of lwcexc has been distributed
                if index+bb > len(dz)-1: # if we reach bottom, stop below distribution
                    tobelow = 0
        
        toabove = 1 # if there is still some lwcexc to distribute but below distribution not possible anymore, same process in layers above
        aa = np.maximum(1,index-lowest) # start to look for space above the aquifer
        if index-aa < 0:
            toabove = 0
        while lwcexc > 0. and toabove == 1:
            if rho[index-aa]<rhoPdr:
                transf = np.minimum(lwcexc,(0.9*(Ptheta_sat[index-aa]-Ptheta[index-aa])*dz[index-aa]))
                transf = np.maximum(transf,0.)
            elif rho[index-aa]>=rhoPdr:
                transf = 0.
            Ptheta[index-aa] += transf/dz[index-aa]
            aa += 1
            lwcexc -= transf
            if index-aa < 0:
                toabove = 0
        if lwcexc > 0: # if we could not distribute all the excess
            totrunoff += lwcexc # it is added to runoff
    
    ### Update of effSat and lwc 
    PeffSat = Ptheta/Ptheta_sat
    Plwc = Ptheta*dz
    
    return Ptheta,PeffSat,Plwc,totrunoff


def entrysuction(dz,Mtheta,Mthetar,Mthetar_old,MeffSat,Mtheta_sat,Mlwc,Ptheta,PeffSat,Plwc,Ptheta_sat,crtn_theta,aquif,MSat_westag):
    ''' 
    When the pressure exceeds the water entry suction of layer below, water penetrates from MFdom layer to PFdom of layer below
    We convert water entry suction in more intuitive saturation value
    If after transfer, saturation of the layer below in PFdom is still inferior to saturation in layer above in MFdom, we equalise saturations
    '''
    
    entry = np.where(MeffSat[0:aquif-1]>MSat_westag[0:aquif-1])[0] # layers where water entry suction is exceeded, use 
    for ii in entry:        
        transfer = dz[ii]*(MSat_westag[ii]*(Mtheta_sat[ii]-Mthetar[ii])+Mthetar[ii]) #amount of water to transfer so MeffSat[ii]==MSat_we[ii+1]
        transfer = min(transfer,0.95*(Ptheta_sat[ii+1]-Ptheta[ii+1])*dz[ii+1]) # no oversaturation
        transfer = min(transfer,(Mtheta[ii]-crtn_theta/10)*dz[ii]) # preserve numerical stability
        Mtheta[ii] -= transfer/dz[ii] # convert in water content
        ## Update values
        Mthetar[ii] = min((0.02),0.9*Mtheta[ii])
        if Mtheta[ii]<Mthetar[ii]+1e-6:
            if Mtheta[ii]>crtn_theta/10:
                Mthetar[ii] = Mtheta[ii] - crtn_theta/10
            if Mtheta[ii]<=crtn_theta/10:
                Mthetar[ii] = 0
        MeffSat[ii] = (Mtheta[ii]-Mthetar[ii])/(Mtheta_sat[ii]-Mthetar[ii])
        Mlwc[ii] = Mtheta[ii]*(dz[ii]) # update lwc
        
        Ptheta[ii+1] += transfer/dz[ii+1] # transfer towards PFdom
        PeffSat[ii+1] = Ptheta[ii+1] / Ptheta_sat[ii+1]
        Plwc[ii+1] = Ptheta[ii+1]*dz[ii+1]
        
        ## If saturation in MFdom of layer ii is still abobe saturation in PFdom of layer ii+1, we equalise (Wever 2016) ##
        if (PeffSat[ii+1] < MeffSat[ii]) and (Mtheta[ii]>crtn_theta/10):
            # Wever 2016, equations (3),(4),(5) -> corrected versions (there were some errors in (5))
            lwctot = Mtheta[ii]*dz[ii] + Ptheta[ii+1]*dz[ii+1]
            Mtheta[ii] = (dz[ii+1]*(Mthetar[ii]*Ptheta_sat[ii+1])+lwctot*(Mtheta_sat[ii]-Mthetar[ii])) / (dz[ii]*(Mtheta_sat[ii]-Mthetar[ii])+dz[ii+1]*Ptheta_sat[ii+1])
            if Mtheta[ii] < crtn_theta/10:
                Mtheta[ii] = crtn_theta/10 
            Ptheta[ii+1] = (lwctot-Mtheta[ii]*dz[ii])/dz[ii+1]
            
            ## Update all variables ##
            Mthetar[ii] = min((0.02),0.9*Mtheta[ii])
            if Mtheta[ii]<Mthetar[ii]+1e-6:
                if Mtheta[ii]>crtn_theta/10:
                    Mthetar[ii] = Mtheta[ii] - crtn_theta/10
                if Mtheta[ii]<=crtn_theta/10:
                    Mthetar[ii] = 0
            MeffSat[ii] = (Mtheta[ii]-Mthetar[ii])/(Mtheta_sat[ii]-Mthetar[ii])
            Mlwc[ii] = Mtheta[ii]*(dz[ii]) # update Plwc
            PeffSat[ii+1] = Ptheta[ii+1] / Ptheta_sat[ii+1]
            Plwc[ii+1] = Ptheta[ii+1]*dz[ii+1]
            #print('In entrysuction, equalisation of saturation between layers', ii, ii+1, 'MeffSat and PeffSat are:',MeffSat[ii], PeffSat[ii+1])
                
    return Mtheta,Mthetar,MeffSat,Mlwc,Ptheta,PeffSat,Plwc

def layerequaliser_eq(dz,Mtheta,Mthetar,Mthetar_old,MeffSat,Mtheta_sat,Mlwc,Ptheta,PeffSat,Plwc,Ptheta_sat,crtn_theta,aquif):
    '''
    Whenever the effective saturation in the MFdom exceeds the saturation in the PFdom within a same layer, we equalise saturations
    Not sure this is physically realistic but that is what is written in Wever 2016
    https://models.slf.ch/p/snowpack/source/tree/HEAD/branches/dev/snowpack/snowpackCore/ReSolver1d.cc
    + confirmed by Nander Wever, email 27 Sept 2018
    '''
    highsat = np.where(MeffSat[0:aquif]>PeffSat[0:aquif])[0] # Where MeffSat exceeds PeffSat
    oknum   = np.where(Mtheta[0:aquif]>crtn_theta/10)[0] # Don't transfer water where we are at value of numerical stability
    okres   = np.where(Mtheta[0:aquif]>0.02)[0]
    toequalise = np.intersect1d(highsat,okres)
    for ii in toequalise:
        while (abs(PeffSat[ii]-MeffSat[ii])>1e-5): # use of a while loop because Mthetar might change -> MeffSat not as high as what was fixed by first exchange
            # Wever 2016, equations (3),(4),(5) -> corrected versions (there were some errors in (5))
            lwctot = Mlwc[ii] + Plwc[ii]
            Mtheta[ii] = (dz[ii]*(Mthetar[ii]*Ptheta_sat[ii])+lwctot*(Mtheta_sat[ii]-Mthetar[ii])) / (dz[ii]*(Mtheta_sat[ii]-Mthetar[ii])+dz[ii]*Ptheta_sat[ii])
            if Mtheta[ii] < crtn_theta/10:
                Mtheta[ii] = crtn_theta/10         
            Ptheta[ii] = (lwctot-Mtheta[ii]*dz[ii])/dz[ii]
            
            ## Update all variables
            Mthetar[ii] = min((0.02),0.9*Mtheta[ii])
            if Mtheta[ii]<Mthetar[ii]+1e-6:
                if Mtheta[ii]>crtn_theta/10:
                    Mthetar[ii] = Mtheta[ii] - crtn_theta/10
                if Mtheta[ii]<=crtn_theta/10:
                    Mthetar[ii] = 0
            MeffSat[ii] = (Mtheta[ii]-Mthetar[ii])/(Mtheta_sat[ii]-Mthetar[ii])
            Mlwc[ii] = Mtheta[ii]*(dz[ii]) # update Plwc
            PeffSat[ii] = Ptheta[ii] / Ptheta_sat[ii]
            Plwc[ii] = Ptheta[ii]*dz[ii]

    return Mtheta,Mthetar,MeffSat,Mlwc,Ptheta,PeffSat,Plwc   



def PFleave(dz,rho,Tz,Mtheta,Mthetar,Mthetar_old,MeffSat,Mtheta_sat,Mlwc,Ptheta,PeffSat,Plwc,Ptheta_sat,crtn_theta,rhoimp,aquif,PSatlim):
    '''
    When the saturation in the PF domain reaches a threshold value (PSatlim), backflow towards MFdom occurs.
    First we transfer as much water as the cold content of the layer can accomodate
    Second, if PSatlim is still exceeded, we equalise saturations in both domains
    Note that at the end, saturation might not be the same if we equalise because Mthetar can change subsequently
    '''
    
    ### First : calculate cold content of every layer ###
    ## Layers mass ##
    mass = rho*dz
    ## Calculate the refreezing potential in every layer ##
    cold_content            = CP_I * mass * (T_MELT - Tz) # cold content of each box, i.e. how much heat to bring it to 273K [J]
    cold_content_sum        = cold_content.cumsum(axis=0) # cumulative cold content, starting from the surface [J]
    refreeze_mass_pot       = cold_content / LF_I # how much mass of the meltwater could be refrozen due to cold content [kg]
    refreeze_mass_pot       = np.maximum(refreeze_mass_pot,0.)
    refreeze_vol_pot        = refreeze_mass_pot/1000. # how much meters of the meltwater could be refrozen due to cold content [m]
    
    backflowlayers = np.where(PeffSat[0:aquif] > PSatlim)[0]
    for ii in backflowlayers:
        if rho[ii] < rhoimp: # Don't transfer water in MFdom of ice layers
        
            transf = refreeze_vol_pot[ii] # transfer amount that can be accomodated by cold content
            transf = min(transf,(Mtheta_sat[ii]-Mtheta[ii])*dz[ii]) # no oversaturation
            transf = min(transf,(Ptheta[ii]-crtn_theta/10)*dz[ii]) # preserve numerical stability
            transf = max(transf,0.) # make sure no negative values
            Ptheta[ii] -= transf/dz[ii]
            Mtheta[ii] += transf/dz[ii]
            
            ## Update all variables
            Mthetar[ii] = min((0.02),0.9*Mtheta[ii])
            if Mtheta[ii]<Mthetar[ii]+1e-6:
                if Mtheta[ii]>crtn_theta/10:
                    Mthetar[ii] = Mtheta[ii] - crtn_theta/10
                if Mtheta[ii]<=crtn_theta/10:
                    Mthetar[ii] = 0
            MeffSat[ii] = (Mtheta[ii]-Mthetar[ii])/(Mtheta_sat[ii]-Mthetar[ii])
            Mlwc[ii] = Mtheta[ii]*(dz[ii]) # update Plwc
            PeffSat[ii] = Ptheta[ii] / Ptheta_sat[ii]
            Plwc[ii] = Ptheta[ii]*dz[ii]
            if (PeffSat[ii]>PSatlim and PeffSat[ii]>MeffSat[ii]): # If we still exceed the threshold saturation in PFdom, we equalise saturations
#                    if PeffSat[ii] < MeffSat[ii]:
#                        print('PeffSat[ii] < MeffSat[ii] but we equalise in PFleave')
                while (abs(PeffSat[ii]-MeffSat[ii])>1e-5): # use of a while loop because Mthetar might change -> MeffSat not as high as what was fixed by first exchange
                    # Wever 2016, equations (3),(4),(5) -> corrected versions (there were some errors in (5))
                    lwctot = Mlwc[ii] + Plwc[ii]
                    Mtheta[ii] = (dz[ii]*(Mthetar[ii]*Ptheta_sat[ii])+lwctot*(Mtheta_sat[ii]-Mthetar[ii])) / (dz[ii]*(Mtheta_sat[ii]-Mthetar[ii])+dz[ii]*Ptheta_sat[ii])       
                    Ptheta[ii] = (lwctot-Mtheta[ii]*dz[ii])/dz[ii]
                    
                    ## Update all variables
                    Mthetar[ii] = min((0.02),0.9*Mtheta[ii])
                    if Mtheta[ii]<Mthetar[ii]+1e-6:
                        if Mtheta[ii]>crtn_theta/10:
                            Mthetar[ii] = Mtheta[ii] - crtn_theta/10
                        if Mtheta[ii]<=crtn_theta/10:
                            Mthetar[ii] = 0
                    MeffSat[ii] = (Mtheta[ii]-Mthetar[ii])/(Mtheta_sat[ii]-Mthetar[ii])
                    Mlwc[ii] = Mtheta[ii]*(dz[ii]) # update Plwc
                    PeffSat[ii] = Ptheta[ii] / Ptheta_sat[ii]
                    Plwc[ii] = Ptheta[ii]*dz[ii]
            
    return Mtheta,Mthetar,MeffSat,Mlwc,Ptheta,PeffSat,Plwc
    

def PFleaveheat(dz,rho,Tz,Mtheta,Mthetar,Mthetar_old,MeffSat,Mtheta_sat,Mlwc,Ptheta,PeffSat,Plwc,Ptheta_sat,crtn_theta,kth,bigF,bigN,aquif,rhoimp,deltatime):
    '''
    This mimics refreezing in the PF domain by transferring water back to the MF dom depending on how much heat would be lost by water in PFdom
    (see paragraph 2.3 of Wever 2016)
    It depends on the tuning parameter bigN which represents number of preferential flow paths per unit area
    Based on equations (6) and (7) of Wever 2016 but caution, misformulation in (7)
    '''
    coldlayers = np.where(Tz[0:aquif]<273.15)[0]
    oknum   = np.where(Ptheta[0:aquif]>crtn_theta/10)[0]
    layers = np.intersect1d(coldlayers,oknum)
    for ii in layers:
        if rho[ii]<rhoimp:
            bigQ = kth[ii]*abs(Tz[ii]-273.15)/(((1+bigF[ii])/(2*math.pi))**0.5-(bigF[ii]/math.pi)**0.5) # Wever 2016 (6)
            transf = dz[ii] * 2*bigN*((math.pi*bigF[ii])**0.5*bigQ*deltatime)/(LF_I*RHO_W_KGM) # Wever 2016 (7) corrected version + expressed in water amount instead of volumetric water content
            transf = min(transf,(Mtheta_sat[ii]-Mtheta[ii])*dz[ii]) # no oversaturation
            transf = min(transf,(Ptheta[ii]-crtn_theta/10)*dz[ii]) # preserve numerical stability
            transf = max(transf,0.) # make sure no negative values
            Ptheta[ii] -= transf/dz[ii]
            Mtheta[ii] += transf/dz[ii]
            #print('layer and transf are:',ii, transf)
            
            ## Update all variables
            Mthetar[ii] = min((0.02),0.9*Mtheta[ii])
            if Mtheta[ii]<Mthetar[ii]+1e-6:
                if Mtheta[ii]>crtn_theta/10:
                    Mthetar[ii] = Mtheta[ii] - crtn_theta/10
                if Mtheta[ii]<=crtn_theta/10:
                    Mthetar[ii] = 0
            MeffSat[ii] = (Mtheta[ii]-Mthetar[ii])/(Mtheta_sat[ii]-Mthetar[ii])
            Mlwc[ii] = Mtheta[ii]*(dz[ii]) # update Plwc
            PeffSat[ii] = Ptheta[ii] / Ptheta_sat[ii]
            Plwc[ii] = Ptheta[ii]*dz[ii]
        
    return Mtheta,Mthetar,MeffSat,Mlwc,Ptheta,PeffSat,Plwc
        

def Mrefreezing(dz,zstep,rho,grain,Tz,Mthetar_old,Mlwc,lwc_min_fr,Ptheta,PeffSat,Plwc,h_e,bigF,mu,crtn_theta,rhoimp,totrefrozen_lwc,refrozenlay,totrunoff):
    ''' 
    Proceed to refreezing according to cold content of every layer.
    Adjust porosity and hydraulic properties accordingly.
    Verify that we don't oversaturate PFdom. If so, we call for Psatexcess
    '''
    
    ### Layers mass ###
    mass = rho*dz
    
    ### Calculate the refreezing potential in every layer ###
    cold_content            = CP_I * mass * (T_MELT - Tz) # cold content of each box, i.e. how much heat to bring it to 273K [J]
    cold_content_sum        = cold_content.cumsum(axis=0) # cumulative cold content, starting from the surface [J]
    refreeze_mass_pot       = cold_content / LF_I # how much mass of the meltwater could be refrozen due to cold content [kg]
    refreeze_mass_pot       = np.maximum(refreeze_mass_pot,0.)
    refreeze_mass_pot_sum   = refreeze_mass_pot.cumsum(axis=0) # cumulative amount of meltwater that can be refrozen due to cold content [kg]
    rho_pot                 = (mass + refreeze_mass_pot) / dz # density value of the boxes if the refreezemass refroze [kg/m3]
    porosity_pot            = 1 - rho_pot / RHO_I # porosity value of the boxes if the refreezemass refroze [/]
    porespace_vol_pot       = porosity_pot * dz # pore space of the boxes if the refreezemass refroze [m]
    
    refreeze_vol_pot        = refreeze_mass_pot/1000. # how much meters of the meltwater could be refrozen due to cold content [m]
    refreeze_vol_pot_sum    = refreeze_vol_pot.cumsum(axis=0) # cumulative amount of meltwater that can be refrozen due to cold content [m]
    
    ### Refreezing process ###
    refrozen_vol = np.minimum(np.maximum((Mlwc-lwc_min_fr),0),refreeze_vol_pot) # Volume of water refrozen [mWE]
    
    porlimit = np.zeros_like(refreeze_mass_pot) # Now we want to avoid exceeding RHO_I <-> reaching negative porosity
    excrefr = np.where(porosity_pot<0)[0] # Spot where this might happen
    if np.any(refrozen_vol < 0):
        print('Negative refrozen_vol before excrefr')
    for ii in excrefr:
        porlimit[ii] = 1e-3*((RHO_I)*dz[ii]-mass[ii]) # [mWE] limit of possible refreeze due to pore space availability, safety margin provided by min value for porosity_refr
        if porlimit[ii] < refrozen_vol[ii]: # If cold content and available LWC are to make the layer exceed RHO_I-1e-3
            refrozen_vol[ii] = min(refrozen_vol[ii],porlimit[ii]) # we limit the volume we will refreeze
            if refrozen_vol[ii] < 0:
                print('Error due to a negative refrozen_vol')
            #refrozen_vol[ii] = min(refrozen_vol[ii],0.)
    if np.any(refrozen_vol < 0):
        print('Negative refrozen_vol after excrefr')
    
    refrozen_mass = 1000*refrozen_vol # Corresponding mass of water refrozen [kg]
    Mlwc = Mlwc-refrozen_vol # New value of lwc [m]
    refreeze_vol_pot = refreeze_vol_pot-refrozen_vol # what can still be refrozen after the refreezing (!=0 if lwc is limiting factor) [m]
    refreeze_mass_pot = 1000*refreeze_vol_pot # what can still be refrozen after the refreezing (!=0 if lwc is limiting factor) [kg]
    cold_content = refreeze_mass_pot*LF_I # remaining cold content [J]
    
    totrefrozen_lwc += np.sum(refrozen_vol) # [mWE]
    refrozenlay += refrozen_vol # [mWE]
    
    lat_heat = refrozen_mass*LF_I # latent heat released in every layer [J]
    if np.any(refrozen_mass < 0):
                print('Error due to a negative refrozen_mass')
    mass = mass + refrozen_mass # new mass: we added the mass of refrozen water
    Tz = T_MELT - cold_content/(CP_I*mass) # the remaining cold content is equivalent to the energy to raise new mass from new Tz until T_MELT [K]
    rho = mass/dz
    Mtheta = Mlwc/dz
    
    ### Calculate pore space available in every layer --> water content at saturation 
    porosity           = 1 - rho/RHO_I # Definition of porosity [/]
    porespace_vol      = porosity * dz # Pore space of each layer [m]
    porosity_refr      = porosity*RHO_I/RHO_W_KGM # space available for liq water volume once refrozen, Wever 2014 (9) [/]
    #test
    #porosity_refr      = np.maximum(porosity_refr,1e-4) # allow space for minimum water content required in both domains for numerical stability, 1e-4 is equivalent to 916.9 density
    porosity_refr      = np.maximum(porosity_refr,17e-3) # allow space for minimum water content required in both domains for numerical stability, 17e-3 is equivalent to 900.0 density
    porespace_refr_vol = porosity_refr*dz # Available pore space of each layer [m]
    
    ### Re-execute all necessary calculations ###
    theta_sat = porosity_refr # value of volumetric water content in saturated conditions [/]
    alpha_vG = 4.4e6*(rho/(2*grain))**(-0.98) # Hirashima 2014 (5) ; typical value ~35.0 
    n_vG     = 1 + 2.7e-3*(rho/(2*grain))**(0.61) # Hirashima 2014 (6) ; typical value ~4.0
    m_vG     = 1 - 1/n_vG # Wever 2014 (8) ; typical value ~0.75
    Sc    = (1 + (alpha_vG*h_e)**n_vG)**(-m_vG) #Saturation at cut-off point [/], see Ippisch et al., 2006 eq(11)
    Ksat = RHO_W_KGM*GRAVITY/mu * 3.0*(grain)**2*np.exp(-0.013*rho) # Hydraulic conductivity at saturation (>0) [m s-1], Formula of Calonne et al. 2012, see Wever 2015 (7) and D'Amboise 2017 (10)
    
    Mtheta_sat = (1-bigF)*theta_sat
    Ptheta_sat = bigF*theta_sat
    
    if np.any(Mtheta > Mtheta_sat): # If any saturation exceeds maximal saturation, we have to proceed to a redistribution of water
        Mtheta,Mthetar,MeffSat,Mlwc,totrunoff = Msatexcess(dz,rho,Mtheta,Mtheta_sat,crtn_theta,rhoimp,totrunoff) # could happen if water left in a layer that reached RHO_I-1e-3
    
    ## Update of Mthetar as in Wever 2014 ##
    Mthetar = np.minimum((np.ones_like(Mtheta)*0.02),0.9*Mtheta) # residual water content [/]
    ## Update of effSat and head ##
    MeffSat = (Mtheta-Mthetar)/(Mtheta_sat-Mthetar)
    Mhead  = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)


    PeffSat = Ptheta/Ptheta_sat # This might have change as theta_sat has been decreased (but not Mtheta_we) -> Ptheta_sat decreased
    # Ptheta_sat decreased -> PeffSat increased -> We might have oversaturated PF domain ! I think this would be really rare but implemented just in case.
    if np.any(PeffSat>1):
        #print('sum(Plwc) before Psatexcess:',sum(Plwc))
        Ptheta,PeffSat,Plwc,totrunoff = Psatexcess(dz,rho,Ptheta,Ptheta_sat,crtn_theta,rhoimp,totrunoff)
        #print('sum(Plwc) after Psatexcess:',sum(Plwc))
            
    return rho,Tz,Mhead,Mtheta,Mthetar,MeffSat,Mlwc,Mtheta_sat,Ptheta,PeffSat,Plwc,Ptheta_sat,Ksat,theta_sat,alpha_vG,n_vG,m_vG,Sc,totrefrozen_lwc,refrozenlay,totrunoff
    

def runoff(dz,rho,Mhead,Mtheta,Mthetar,Mthetar_old,MeffSat,Mlwc,Mtheta_sat,theta_min_fr,crtn_theta,slope,rhoimp,aquif,deltatime,totrunoff):
    '''
    This is the Greenland specific runoff function of Zuo and Oerlemans 1996 (21)
    We proceed to runoff only for layers of which water content is above theta_min_fr. CAUTION: this might be changed!!
    We don't apply runoff in the aquifer
    '''
    c1Z = 1.5*24*3600 # value in Zuo and Oerlemans 1996, converted in seconds [s]
    c2Z = 25*24*3600 # Zuo and Oerlemans 1996 [s]
    c3Z = 140 # Zuo and Oerlemans 1996 [/]
    tcarac = c1Z + c2Z*np.exp(-1*c3Z*slope) # Zuo and Oerlemans 1996 [s]
    wet_layers = np.where(Mtheta>theta_min_fr)[0] # other possibility: apply runoff to all layers where there is water
    #abice = np.where(rho>=rhoimp)[0] - 1
    #absat = np.where(MeffSat>=0.9)[0] - 1
    #abicesat = np.unique(np.append(abice,absat))
    runoff_layers = np.copy(wet_layers)
    #print('Runoff layers are:',runoff_layers)
    runoff = np.zeros_like(dz) # lwc that will be moved by lateral runoff [m]
    for index in runoff_layers:
        if (index < aquif): # and (rho[index]<830)):
            runoff[index] = deltatime*((Mtheta[index]*dz[index])-(Mthetar[index]*dz[index]))/tcarac # Zuo and Oerlemans 1996 (21), Lefebvre 2003 (1), Langen 2017 (13) [m]
            runoff[index] = min(runoff[index],((Mtheta[index]*dz[index])-(Mthetar[index]*dz[index]))) # Don't let runoff decrease Mtheta below Mthetar
            runoff[index] = min(runoff[index],(Mtheta[index]-crtn_theta/10)*dz[index]) # make sure to preserve numerical stability
            runoff[index] = max(runoff[index],0.) # make sure no negative runoff value
    
            Mtheta[index] -= runoff[index]/dz[index] # remove the runoff from the MFdom
            
            ## Update of all variables
            Mthetar[index] = np.minimum(0.02,0.9*Mtheta[index])
            #if Mtheta[index]<Mthetar[index]+1e-6:
                #if Mtheta[index]>crtn_theta/10:
                    #Mthetar[index] = Mtheta[index] - crtn_theta/10
                #if Mtheta[index]<=crtn_theta/10:
                    #Mthetar[index] = 0
            Mlwc[index] = Mtheta[index]*dz[index]
            MeffSat[index] = (Mtheta[index]-Mthetar[index])/(Mtheta_sat[index]-Mthetar[index])

    totrunoff += sum(runoff) # total amount of runoff
    return Mtheta,Mthetar,MeffSat,Mlwc,totrunoff


def Prefreezing(dz,rho,grain,Tz,Mthetar_old,Mtheta,Mlwc,lwc_min_fr,Ptheta,PeffSat,Plwc,bigF,h_e,mu,crtn_theta,dryfront,totrefrozen_lwc,refrozenlay,rhoimp,totrunoff):
    ''' 
    Not used in preferential flow scheme of Wever 2016 but might be a good idea to use.
    Proceed to refreezing of Plwc until dryfront according to cold content of every layer.
    Adjust porosity and hydraulic properties accordingly.
    Also needs the variables Mtheta and dryfront as input parameters (vs refreezing())
    '''
    
    ### Layers mass ###
    mass = rho[0:dryfront+1]*dz[0:dryfront+1]
    
    ### Calculate the refreezing potential in every layer ###
    cold_content            = CP_I * mass[0:dryfront+1] * (T_MELT - Tz[0:dryfront+1]) # cold content of each box, i.e. how much heat to bring it to 273K [J]
    cold_content_sum        = cold_content.cumsum(axis=0) # cumulative cold content, starting from the surface [J]
    refreeze_mass_pot       = cold_content / LF_I # how much mass of the meltwater could be refrozen due to cold content [kg]
    refreeze_mass_pot       = np.maximum(refreeze_mass_pot,0.)
    refreeze_mass_pot_sum   = refreeze_mass_pot.cumsum(axis=0) # cumulative amount of meltwater that can be refrozen due to cold content [kg]
    rho_pot                 = (mass[0:dryfront+1] + refreeze_mass_pot[0:dryfront+1]) / dz[0:dryfront+1] # density value of the boxes if the refreezemass refroze [kg/m3]
    porosity_pot            = 1 - rho_pot / RHO_I # porosity value of the boxes if the refreezemass refroze [/]
    porespace_vol_pot       = porosity_pot * dz[0:dryfront+1] # pore space of the boxes if the refreezemass refroze [m]
    
    refreeze_vol_pot        = refreeze_mass_pot/1000. # how much meters of the meltwater could be refrozen due to cold content [m]
    refreeze_vol_pot_sum    = refreeze_vol_pot.cumsum(axis=0) # cumulative amount of meltwater that can be refrozen due to cold content [m]
    
    ### Refreezing process ###
    refrozen_vol = np.minimum(np.maximum(Plwc[0:dryfront+1]-lwc_min_fr[0:dryfront+1],0),refreeze_vol_pot) # Volume of water refrozen [mWE]
    
    porlimit = np.zeros_like(refreeze_mass_pot) # Now we want to avoid exceeding RHO_I <-> reaching negative porosity
    excrefr = np.where(porosity_pot<0)[0] # Spot where this might happen
    if np.any(refrozen_vol < 0):
        print('Negative refrozen_vol before excrefr')
    for ii in excrefr:
        #porlimit[ii] = 1e-3*((RHO_I-1e-3)*dz[ii]-mass[ii]) # [mWE] limit of possible refreeze due to pore space availability with little safety margin
        porlimit[ii] = 1e-3*((RHO_I)*dz[ii]-mass[ii]) # [mWE] limit of possible refreeze due to pore space availability, safety margin provided by min value for porosity_refr
        if porlimit[ii] < refrozen_vol[ii]: # If cold content and available LWC are to make the layer exceed RHO_I-1e-3
            refrozen_vol[ii] = min(refrozen_vol[ii],porlimit[ii]) # we limit the volume we will refreeze
            if refrozen_vol[ii] < 0:
                print('Error due to a negative refrozen_vol')
            #refrozen_vol[ii] = min(refrozen_vol[ii],0.)
    if np.any(refrozen_vol < 0):
        print('Negative refrozen_vol after excrefr')
    
    refrozen_mass = 1000*refrozen_vol # Corresponding mass of water refrozen [kg]
    Plwc[0:dryfront+1] = Plwc[0:dryfront+1]-refrozen_vol # New value of lwc [m]
    refreeze_vol_pot = refreeze_vol_pot-refrozen_vol # what can still be refrozen after the refreezing (!=0 if lwc is limiting factor) [m]
    refreeze_mass_pot = 1000*refreeze_vol_pot # what can still be refrozen after the refreezing (!=0 if lwc is limiting factor) [kg]
    cold_content = refreeze_mass_pot*LF_I # remaining cold content [J]
    
#    print('refrozen_vol[0:10] is:',refrozen_vol[0:10])
    totrefrozen_lwc += np.sum(refrozen_vol) # [mWE]
    refrozenlay[0:dryfront+1] += refrozen_vol # [mWE]
    
    lat_heat = refrozen_mass*LF_I # latent heat released in every layer [J]
    if np.any(refrozen_mass < 0):
                print('Error due to a negative refrozen_mass')
    mass = mass + refrozen_mass # new mass: we added the mass of refrozen water
    Tz[0:dryfront+1] = T_MELT - cold_content/(CP_I*mass) # the remaining cold content is equivalent to the energy to raise new mass from new Tz until T_MELT [K]
    rho[0:dryfront+1] = mass/dz[0:dryfront+1]
    Ptheta[0:dryfront+1] = Plwc[0:dryfront+1]/dz[0:dryfront+1]
    
    ### Calculate pore space available in every layer --> water content at saturation 
    porosity           = 1 - rho/RHO_I # Definition of porosity [/]
    porespace_vol      = porosity * dz # Pore space of each layer [m]
    porosity_refr      = porosity*RHO_I/RHO_W_KGM # space available for liq water volume once refrozen, Wever 2014 (9) [/]
    #test
    #porosity_refr      = np.maximum(porosity_refr,1e-4) # allow space for minimum water content required in both domains for numerical stability, 1e-4 is equivalent to 916.9 density
    porosity_refr      = np.maximum(porosity_refr,17e-3) # allow space for minimum water content required in both domains for numerical stability, 17e-3 is equivalent to 900.0 density
    porespace_refr_vol = porosity_refr*dz # Available pore space of each layer [m]
    
    ### Re-execute all necessary calculations ###
    theta_sat = porosity_refr # value of volumetric water content in saturated conditions [/]
    Mtheta_sat = (1-bigF)*theta_sat
    Ptheta_sat = bigF*theta_sat
    alpha_vG = 4.4e6*(rho/(2*grain))**(-0.98) # Hirashima 2014 (5) ; typical value ~35.0 
    n_vG     = 1 + 2.7e-3*(rho/(2*grain))**(0.61) # Hirashima 2014 (6) ; typical value ~4.0
    m_vG     = 1 - 1/n_vG # Wever 2014 (8) ; typical value ~0.75
    Sc    = (1 + (alpha_vG*h_e)**n_vG)**(-m_vG) #Saturation at cut-off point [/], see Ippisch et al., 2006 eq(11)
    Ksat = RHO_W_KGM*GRAVITY/mu * 3.0*(grain)**2*np.exp(-0.013*rho) # Hydraulic conductivity at saturation (>0) [m s-1], Formula of Calonne et al. 2012, see Wever 2015 (7) and D'Amboise 2017 (10)
    
    Mtheta_sat = (1-bigF)*theta_sat
    Ptheta_sat = bigF*theta_sat
    
    if np.any(Mtheta > Mtheta_sat): # If any saturation exceeds maximal saturation, we have to proceed to a redistribution of water
        Mtheta,Mthetar,MeffSat,Mlwc,totrunoff = Msatexcess(dz,rho,Mtheta,theta_sat,crtn_theta,rhoimp,totrunoff) # could happen if water left in a layer that reached RHO_I-1e-3
    
    ## Update of Mthetar as in Wever 2014 ##
    Mthetar = np.minimum((np.ones_like(Mtheta)*0.02),0.9*Mtheta) # initial residual water content [/], Wever 2014 (10)
    #low_Mtheta = np.where(Mtheta<Mthetar+crtn_theta/10) # for low theta values: same procedure as Wever 2014, appendix A3
    #for indices in low_Mtheta[0]:
    #    if Mtheta[indices]>crtn_theta/10:
    #        Mthetar[indices] = Mtheta[indices] - crtn_theta/10
    #    if Mtheta[indices]<=crtn_theta/10:
    #        Mthetar[indices] = 0     
    ## Update of effSat and head ##
    MeffSat = (Mtheta-Mthetar)/(Mtheta_sat-Mthetar)
    Mhead  = -1*1/alpha_vG * ((Sc * MeffSat)**(-1/m_vG)-1)**(1/n_vG) # [m] Wever 2014 (3)

    PeffSat = Ptheta/Ptheta_sat # This might have change as theta_sat has been decreased  -> Ptheta_sat decreased
            
    return rho,Tz,Mhead,Mtheta,Mthetar,MeffSat,Mlwc,Mtheta_sat,Ptheta,PeffSat,Plwc,Ptheta_sat,Ksat,theta_sat,alpha_vG,n_vG,m_vG,Sc,totrefrozen_lwc,refrozenlay,totrunoff


def distribute_tostore(dz,rho,tostore,Mlwc,Plwc,rhoimp,bigF,totrunoff):
    
    ### Calculate pore space available in every layer --> water content at saturation 
    porosity           = 1 - rho/RHO_I # Definition of porosity [/]
    porespace_vol      = porosity * dz # Pore space of each layer [m]
    porosity_refr      = porosity*RHO_I/RHO_W_KGM # space available for liq water volume once refrozen, Wever 2014 (9) [/]
    #porosity_refr      = np.maximum(porosity_refr,1e-4) # allow space for minimum water content required in both domains for numerical stability, 1e-4 is equivalent to 916.9 density
    porosity_refr      = np.maximum(porosity_refr,17e-3) # allow space for minimum water content required in both domains for numerical stability, 17e-3 is equivalent to 900.0 density
    porespace_refr_vol = porosity_refr*dz # Available pore space of each layer [m]
    
    spaceavail = 0.999*porespace_refr_vol # Use a safety margin to avoid calculation problems for the head
    jj = len(dz)-1 # start filling from the bottom layer
    while tostore > 0: # as long as there is water to store, we continue the distribution
        if rho[jj] < rhoimp: # don't put water in ice layers
            toPF = min(bigF[jj]*spaceavail[jj]-Plwc[jj],tostore) # first fill the PFdom part of the porosity
            Plwc[jj] += toPF
            tostore -= toPF # tostore has been (partly) emptied
            toMF = min((1-bigF[jj])*spaceavail[jj]-Mlwc[jj],tostore) # then fill the MFdom part of the porosity
            Mlwc[jj] += toMF
            tostore -= toMF # tostore has been (partly) emptied
        jj -= 1 # go to layer above
        if jj == -1: # if we reach surface
            totrunoff += tostore # put the rest of to store as runoff
            tostore = 0. # nothing to store anymore
    
    return(Mlwc,Plwc,totrunoff)




