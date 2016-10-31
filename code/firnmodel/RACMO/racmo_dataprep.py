import numpy as np
import scipy as sp
import netCDF4 as nc
import os
import sys
import shutil
import matplotlib.pyplot as plt

#remotemain = '/Volumes/Firn/RACMO/Greenland/results_summit'
localmain = '/Users/maxstev/Documents/Grad_School/PIRE/CFM/CommunityFirnModel/code/firnmodel/RACMO/results'
figuredir= '/Users/maxstev/Documents/Grad_School/PIRE/CFM/CommunityFirnModel/code/firnmodel/RACMO/figures'
thenames=thenames=["HLD","HLS","Li2011","helsen","arts","bar","km15","lig","sim","artt"]
leg1=['H & L Dynamic','H & L Sigfus','Li (2011)','Helsen (2008)','Arthern (2010) SS','Barnola (1991)','Kuipers-Munneke (2015)','Ligtenberg (2011)','Simonsen (2013)','Arthern (2010) Transient']
#thenames=thenames=["HLD"]
colos=['b','r','g','c','m','y','k','orange','navy','hotpink']
DIP={}
BCO={}

for nn in range(len(thenames)):
    #for ii in range(11):
        #for jj in range(11):
    ii=5
    jj=5
    numfol=thenames[nn]+'_%s_%s' % (ii, jj)
    
    #remotefolder=os.path.join(remotemain,nn,numfol)
    #localfolder=os.path.join(localmain,nn,numfol)
    #os.makedirs(localfolder)
    #srcDIP=os.path.join(remotefolder,'porosity.csv')
    #dstDIP=os.path.join(localfolder,'porosity.csv')
    #shutil.copyfile(srcDIP, dstDIP)
    #
    #srcBCO=os.path.join(remotefolder,'BCO.csv')
    #dstBCO=os.path.join(localfolder,'BCO.csv')
    #shutil.copyfile(srcBCO, dstBCO)
    
    localfolder=os.path.join(localmain,thenames[nn],numfol)
    dstDIP=os.path.join(localfolder,'porosity.csv')
    dstBCO=os.path.join(localfolder,'BCO.csv')
    key = thenames[nn]
    data_DIP=np.genfromtxt(dstDIP, delimiter=',')
    data_BCO=np.genfromtxt(dstBCO, delimiter=',')
    DIP[key]=data_DIP
    BCO[key]=data_BCO
    
    time1=data_DIP[0,:]
    DIP1=data_DIP[1,:]
    DIP2=DIP1-DIP1[0]
    time2=data_BCO[0,:]
    BCOdep=data_BCO[4,:]
    BCOage=data_BCO[3,:]
    BCOdep2=BCOdep-BCOdep[0]
    
    coluse=colos[nn]
    
    #plt.figure(1)
    #plt.plot(time1,DIP1,coluse)
    #plt.xlabel('Year')
    #plt.ylabel('Depth-Integrated Porosity (m)')
    #plt.xlim(1958,2014)
    #savename1='DIPracmoL.png'
    ##plt.legend(leg1)
    ##plt.savefig(os.path.join(figuredir,savename1),dpi=220)
        
    #plt.figure(2)
    #plt.plot(time2,BCOdep,coluse)#,label=leg1[nn])
    #plt.xlabel('Year')
    #plt.ylabel('Firn thickness (m)')
    #plt.xlim(1958,2014)
    ###plt.legend(loc='upper left')
    #savename1='BCOdepthracmo.png'
    #plt.savefig(os.path.join(figuredir,savename1),dpi=220)
    #
    
    plt.figure(10)
    plt.plot(time2,BCOdep,coluse)
    plt.xlabel('Year')
    plt.ylabel('Firn thickness (m)')
    plt.xlim(1958,2014)
    savename1='BCOdepthracmo2.png'
    plt.savefig(os.path.join(figuredir,savename1),dpi=220)
    
    plt.figure(4)
    plt.plot(time2,BCOdep2,coluse)
    plt.xlabel('Year')
    plt.ylabel('Firn-thickness change (m)')
    plt.xlim(1958,2014)
    #plt.legend(leg1)
    savename1='BCOdepthchangeL.png'
    plt.savefig(os.path.join(figuredir,savename1),dpi=220)
        
    plt.figure(3)
    plt.plot(time1,DIP2,coluse)
    plt.xlabel('Year')
    plt.ylabel('Depth-Integrated-Porosity change (m)')
    plt.xlim(1958,2014)
    savename1='DIPdeviation.png'
    plt.savefig(os.path.join(figuredir,savename1),dpi=220)
    
#plt.legend(loc='upper left')
#savename1='BCOdepthracmoL.png'
#plt.savefig(os.path.join(figuredir,savename1),dpi=220)
plt.show()
    
    
        

