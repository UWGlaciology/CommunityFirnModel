import sys
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

np.set_printoptions(threshold=np.inf)

## Variables to change
rho_0=0.35 # Surface density in g cm^-3
Tc = -30 # Temp in C
A = 0.08 # m W.E., e.g. South Pole is 0.08 m WE

dz=0.1 # grid spacing, m
maxdepth=200 # maximum depth of grid, m
tracker=0

Avec=np.arange(0.005,0.02,0.005)
Tvec=np.arange(-67,-45,1.)


## Model
# T = 273.15+Tc # T in K

h=np.arange(0,maxdepth+dz,dz) # grid
rho_ice=917.0

R = 8.314       #gas constant
rho_i = 0.917   #ice density (Mg m^{-3}
rho_c = 0.550   #critical density (stage 1 - stage 2)

dh=np.diff(h)
dh=np.append(dh, dh[-1])

dage_hold=np.zeros((np.size(Avec)*np.size(Tvec),4))

for ii in range(0,np.size(Avec)):
    for jj in range(0,np.size(Tvec)):
        
        A=Avec[ii]
        Tc=Tvec[jj]
        T = 273.15+Tc # T in K
        # site specific rate-constants, eqns 6a and 6b from Herron+Langway
        k_0 = 11  * np.exp(-10160/(R*T))
        k_1 = 575 * np.exp(-21400/(R*T))
        
        # Given site conditions T, A(ccumulation w.e.) and surface density,
        # can calculate the density-depth profile and age.
        
        #Stage 1 Densification
        #depth of critical density, eqn 8 from Herron and Langway
        h0_55 = 1/(rho_i*k_0)*(np.log(rho_c/(rho_i-rho_c))-np.log(rho_0/(rho_i-rho_0)))
        
        
        Z_0 = np.exp(rho_i*k_0*h + np.log(rho_0/(rho_i-rho_0)))
        
        #age of critical density, eq. 9
        t0_55 = 1/(k_0*A)*np.log((rho_i-rho_0)/(rho_i-rho_c))
        rho_h0 = (rho_i * Z_0)/(1+Z_0)
        t_0 =   1/(k_0*A)*np.log((rho_i-rho_0)/(rho_i-rho_h0))
        
        
        Z_1 = np.exp(rho_i*k_1*(h-h0_55)/np.sqrt(A) +  np.log(rho_c/(rho_i-rho_c)))
        
        #combine Z for Z_0 less than critical density and Z_1 greater than critical
        #density
        Z = np.concatenate((Z_0[h<h0_55],Z_1[h>h0_55]))
        
        # determine rho
        rho_h = (rho_i * Z)/(1+Z)
        
        t_p = 1/(k_1*np.sqrt(A))*np.log((rho_i-rho_c)/(rho_i-rho_h))+ t0_55 #Equation 11
        
        age = np.concatenate((t_0[h<h0_55], t_p[h>h0_55])) #Eq. 12, firn age
        
        rho_hkg=rho_h*1000 #density in kg m^-3
        HLDIP=np.cumsum(((rho_ice-rho_hkg)/rho_ice)*dh) #Depth-integrated porosity
        HLDIPtot=HLDIP[:] #total DIP
        bcoRho = 1/( 1/(917.0) + T*6.95E-7 - 4.3e-5)
        ind2=np.nonzero(rho_hkg>=bcoRho)
        if np.size(ind2) == 0:
            BCOHL=np.nan
            dage=np.nan
        else:
            ind=np.min(ind2)
            BCOHL=h[ind] #bubble close-off depth (where rho=815)
            dage=age[ind]
        
        bco_vec=815*np.ones(np.size(h)) #for plotting
        asdf=np.array((Avec[ii],Tvec[jj], dage,BCOHL))
        dage_hold[tracker,:]=asdf
        tracker=tracker+1
  
adsf=dage_hold[:,3]
xxx=np.isnan(adsf)
teller=np.where(xxx==False)
dage_hold_new=dage_hold[teller,:]


X,Y=np.meshgrid(dage_hold[:,0],dage_hold[:,1])
ZI=griddata((dage_hold[:,0],dage_hold[:,1]),dage_hold[:,3],(X,Y),method='linear')


#
#
#d2=sort(dage_hold(:,3))
#figure(12)clfplot(d2,'.')
#
plt.figure(14)
plt.contourf(X,Y,ZI,256)
plt.xlabel('Accu')
plt.ylabel('Temp')
# set(gca,'ydir','reverse')
plt.colorbar()
plt.show()
#
#fig123=figure(123)
#clf
#hist(d2,20)

# fig1=figure(1)
# hold on
# box on
# grid on
# plot(rho_hkg,h,'b','linewidth',2)
# plot(bco_vec,h,'k','linewidth',1)
# set(gca,'ydir','reverse','fontsize',16)
# xlabel('Density (kg m^{-3})','fontsize',18)
# ylabel('Depth (m)','fontsize',18)
# title('Herron and Langway (1980) firn depth/density profile','fontsize',18)
#
# fig2=figure(2)
# hold on
# box on
# grid on
# plot(rho_hkg,age,'r','linewidth',2)
# # plot(bco_vec,h,'k','linewidth',1)
# set(gca,'ydir','reverse','fontsize',16)
# xlabel('Density (kg m^{-3})','fontsize',18)
# ylabel('Age (years)','fontsize',18)
# title('Herron and Langway (1980) firn depth/age profile','fontsize',18)
# text(300,1000,sprintf('Delta age = #g',dage),'fontsize',18)
#
# fig3=figure(3)
# hold on
# box on
# grid on
# plot(age,h,'r','linewidth',2)
# # plot(bco_vec,h,'k','linewidth',1)
# set(gca,'ydir','reverse','fontsize',16)
# ylabel('Depth (m)','fontsize',18)
# xlabel('Age (years)','fontsize',18)
# title('Herron and Langway (1980) firn depth/age profile','fontsize',18)
# text(300,110,sprintf('Delta age = #g',dage),'fontsize',18)




