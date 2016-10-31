import numpy as np
import scipy as sp
import netCDF4 as nc
import os
import sys
import shutil
import matplotlib.pyplot as plt

resultsmain='/Users/maxstev/Documents/Grad_School/PIRE/CFM/CommunityFirnModel/code/firnmodel/RACMO/results/KM_spinups/'
figuredir= '/Users/maxstev/Documents/Grad_School/PIRE/CFM/CommunityFirnModel/code/firnmodel/RACMO/figures'

k1='KM_RACMO_summit01spin'
k2='KM_RACMO_summit10spin'
k3='KM_RACMO_summit20spin'
k4='KM_RACMO_summit30spin'
k5='KM_RACMO_summit40spin'
k6='KM_RACMO_summit50spin'
k7='KM_RACMO_summitLOOPspin'

DIP1=np.genfromtxt(os.path.join(resultsmain,k1,'porosity.csv'), delimiter=',')
DIP2=np.genfromtxt(os.path.join(resultsmain,k2,'porosity.csv'), delimiter=',')
DIP3=np.genfromtxt(os.path.join(resultsmain,k3,'porosity.csv'), delimiter=',')
DIP4=np.genfromtxt(os.path.join(resultsmain,k4,'porosity.csv'), delimiter=',')
DIP5=np.genfromtxt(os.path.join(resultsmain,k5,'porosity.csv'), delimiter=',')
DIP6=np.genfromtxt(os.path.join(resultsmain,k6,'porosity.csv'), delimiter=',')
DIP7a=np.genfromtxt(os.path.join(resultsmain,k7,'porosity.csv'), delimiter=',')

BCO1=np.genfromtxt(os.path.join(resultsmain,k1,'BCO.csv'), delimiter=',')
BCO2=np.genfromtxt(os.path.join(resultsmain,k2,'BCO.csv'), delimiter=',')
BCO3=np.genfromtxt(os.path.join(resultsmain,k3,'BCO.csv'), delimiter=',')
BCO4=np.genfromtxt(os.path.join(resultsmain,k4,'BCO.csv'), delimiter=',')
BCO5=np.genfromtxt(os.path.join(resultsmain,k5,'BCO.csv'), delimiter=',')
BCO6=np.genfromtxt(os.path.join(resultsmain,k6,'BCO.csv'), delimiter=',')
BCO7a=np.genfromtxt(os.path.join(resultsmain,k7,'BCO.csv'), delimiter=',')

time7a=DIP7a[0,:]
ind=np.nonzero(time7a>=1958.0)
DIP7=DIP7a[:,ind]
BCO7=BCO7a[:,ind]
time=DIP1[0,:]
time6=DIP6[0,:]
time7=DIP7[0,:]

DIPdev1=DIP1-DIP1[1,0]
DIPdev2=DIP2-DIP2[1,0]
DIPdev3=DIP3-DIP3[1,0]
DIPdev4=DIP4-DIP4[1,0]
DIPdev5=DIP5-DIP5[1,0]
DIPdev6=DIP6-DIP6[1,0]
DIPdev7=DIP7-DIP7[1,0]

plt.figure(11)    
plt.plot(time,DIP1[1,:],'b')
plt.plot(time,DIP2[1,:],'r')
plt.plot(time,DIP3[1,:],'g')
plt.plot(time,DIP4[1,:],'k')
plt.plot(time,DIP5[1,:],'m')
plt.plot(time6,DIP6[1,:],'c')
plt.plot(time7,DIP7[1,:],'y')
plt.xlabel('Year')
plt.ylabel('Depth-Integrated Porosity (m)')
plt.xlim(1958,2014)
savename1='DIP_KMspin.png'
plt.savefig(os.path.join(figuredir,savename1),dpi=200)  

plt.figure(12)    
plt.plot(time,DIPdev1[1,:],'b')
plt.plot(time,DIPdev2[1,:],'r')
plt.plot(time,DIPdev3[1,:],'g')
plt.plot(time,DIPdev4[1,:],'k')
plt.plot(time,DIPdev5[1,:],'m')
plt.plot(time6,DIPdev6[1,:],'c')
plt.plot(time7,DIPdev7[1,:],'y')
plt.xlabel('Year')
plt.ylabel('Change in Depth-Integrated Porosity (m)')
plt.xlim(1958,2014)
savename1='DIPdeviation_KMspin.png'
plt.savefig(os.path.join(figuredir,savename1),dpi=200) 

plt.figure(13)    
plt.plot(time,BCO1[1,:],'b')
plt.plot(time,BCO2[1,:],'r')
plt.plot(time,BCO3[1,:],'g')
plt.plot(time,BCO4[1,:],'k')
plt.plot(time,BCO5[1,:],'m')
plt.plot(time6,BCO6[1,:],'c')
plt.plot(time7,BCO7[1,:],'y')
plt.xlabel('Year')
plt.ylabel('Close-Off Depth (m)')
plt.xlim(1958,2014)
savename1='BCOdepth_KMspin.png'
plt.savefig(os.path.join(figuredir,savename1),dpi=200) 

plt.show()
