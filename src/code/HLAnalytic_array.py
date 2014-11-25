# This particular version of the Herron and Langway code is written to take input matrices/vectors of accumulation rate
# and temperature (e.g. from a grid) and calculate the delta age and bubble close-off depth for each T/A pair. Surface density is fixed. 
# Max Stevens, 7/23/14

import sys
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

 
def rhoHLAnalytic(T,A,rho_0):
  
    h=np.arange(0,maxdepth+dz,dz) # grid
    R = 8.314       #gas constant 
    rho_i = 0.917   #ice density (Mg m^{-3}
    rho_c = 0.550   #critical density (stage 1 - stage 2)
    rho_i_kg = 917.0
    
    dh=np.diff(h)
    dh=np.append(dh,dh[-1])
    
    # site specific rate-constants, eqns 6a and 6b from Herron+Langway
    k_0 = 11  * np.exp(-10160/(R*T))  
    k_1 = 575 * np.exp(-21400/(R*T))
    
    # Given site conditions T, A(ccumulation w.e.) and surface density, calculate the density-depth profile and age.
    
    #Stage 1 Densification
    #depth of critical density, eqn 8 from Herron and Langway
    h0_55 = 1/(rho_i*k_0)*(np.log(rho_c/(rho_i-rho_c))-np.log(rho_0/(rho_i-rho_0)))
    Z_0 = np.exp(rho_i*k_0*h + np.log(rho_0/(rho_i-rho_0))) 
    
    #age of critical density, eq. 9
    t0_55 = 1/(k_0*A)*np.log((rho_i-rho_0)/(rho_i-rho_c))
    rho_h0 = (rho_i * Z_0)/(1+Z_0)
    t_0 =   1/(k_0*A)*np.log((rho_i-rho_0)/(rho_i-rho_h0))
    Z_1 = np.exp(rho_i*k_1*(h-h0_55)/np.sqrt(A) +  np.log(rho_c/(rho_i-rho_c)))
    
    #combine Z for Z_0 less than critical density and Z_1 greater than critical density
    Z = np.append(Z_0[h<h0_55], Z_1[h>h0_55])
    
    # determine rho
    rho_h = (rho_i * Z)/(1+Z)
    t_p = 1/(k_1*np.sqrt(A))*np.log((rho_i-rho_c)/(rho_i-rho_h))+ t0_55 #Equation 11
    age = np.append(t_0[h<h0_55], t_p[h>h0_55]) #Eq. 12, firn age
    
    #Calculate the close-off depth and age
    bcoMartRho = 1/( 1/(917.0) + T*6.95E-7 - 4.3e-5) #Bubble close off densityfrom Martinerie
    
    rho_bco=bcoMartRho #use Martinerie close off (what we should do)
    #rho_bco=815.0 #use close off density of 815 (easier to check answers again matlab to see if I have bugs in this code)
    rho_lid=rho_bco-14.0 #LID density is 14 less than BCO.
    
    rho_hkg=rho_h*1000 #density in kg m^-3
    #HLDIP=np.cumsum(((rho_i_kg-rho_hkg)/rho_i_kg)*dh) #Depth-integrated porosity
    #HLDIPtot=HLDIP[-1] #total DIP
    ind=np.min(np.where(rho_hkg>=rho_bco)) #index of the close off depth
    ind2=np.min(np.where(rho_hkg>=rho_lid)) #index of LID
    BCOHL=h[ind] #bubble close-off depth (where rho=815)
    LIDHL=h[ind2] # lock in depth using the rho_cod minus 14 kg/m^3 method
    LIZ_th=BCOHL-LIDHL #thickness of the lock-in zone
    #bco_vec=815.0*np.ones(np.size(h)) #for plotting
    
    dage=age[ind]

    return BCOHL,LIDHL,dage #we are only saving the close-off depth and age (not the density profiles)

if __name__ == "__main__":
    ## Variables to change
    rho_0=0.36 # Surface density in g cm^-3
    Tc = -30 # Temp in C
    #A = 0.2 # m W.E., e.g. South Pole is 0.08 m WE multiply I.E. by 0.917 to get W.E.
    #T = Tc + 273.15 #Temperature in K
    
    ####
    # This bit just for testing. Load vectors or matrices of accumulation and temperature here.
    AA=np.array([[0.2, 0.15, 0.17],[0.09, 0.25, -9999]]) #AA and TT must be the same size!
    TT=np.array([[-30.0, -35.0, -40.0],[-43.0, -27.0, -30.0]])
    
    #AA=np.loadtxt('Path/to/Data/Data.txt')
    #TT = ...

    ####
    
    dz=0.2 # grid spacing, m
    maxdepth=200 # maximum depth of grid, m
    h=np.arange(0,maxdepth+dz,dz) # grid
    
    rows,cols=np.shape(AA)
    sz=np.size(AA)
    
    AAvec=np.reshape(AA,sz)
    TTvec=np.reshape(TT,sz)
    
    BCOvec=np.ones(sz)
    dagevec=np.ones(sz)
    
    for ii in range(0,sz):
    #vecHLA=np.vectorize(rhoHLAnalytic)
        A=AAvec[ii]
        Tc=TTvec[ii]
        if A==-9999 or Tc==-9999:
            BCOvec[ii]=-9999
            dagevec[ii]=-9999            
        
        else:
            T = Tc + 273.15 #Temperature in K
            BCOHL,LIDHL,dage=rhoHLAnalytic(T,A,rho_0)
            #BCOHL,LIDHL,dage=vecHLA(TT[ii],AA[ii],rho_0)
            BCOvec[ii]=BCOHL
            dagevec[ii]=dage
        
    
    BCOmat=np.reshape(BCOvec,(rows,cols)) 
    dagemat=np.reshape(dagevec,(rows,cols))
    
            
        
        