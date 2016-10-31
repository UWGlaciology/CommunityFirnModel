import csv
import numpy as np
import netCDF4 as nc
import os
import sys
import math
from shutil import rmtree

spot = os.path.dirname(sys.argv[0])
os.chdir(spot)
smb_st= '/Users/maxstev/Documents/Grad_School/PIRE/CFM/CommunityFirnModel/code/firnmodel/RACMO/GREENLAND/smb_Summit.RACMO2.3_1958-2014.nc'
temp_st = '/Users/maxstev/Documents/Grad_School/PIRE/CFM/CommunityFirnModel/code/firnmodel/RACMO/GREENLAND/t2m_Summit.RACMO2.3_1958-2014.nc'

saveloc=os.path.join(spot,'RACMO_Summit_csv')

if os.path.exists(saveloc):
    rmtree(saveloc)
os.makedirs(saveloc)

days=30

f=nc.Dataset(smb_st,'r')
g=nc.Dataset(temp_st,'r')

for mm in range(11):
    for nn in range(11):
        tt=f.variables['time'][:]
        tt2=tt.astype('float64')
        smb=f.variables['smb'][:,mm,nn]
        smb2=smb.astype('float64')
        
        tt3=g.variables['time'][:]
        tt4=tt3.astype('float64')
        t2m=g.variables['t2m'][:,mm,nn]
        t2m2=t2m.astype('float64')
        
        tt_out=tt2[0::days]
        tt_out2=tt4[0::days]
        smb_out=np.zeros_like(tt_out)
        t2m_out=np.zeros_like(tt_out)

        for ii in range(len(smb_out)):
            ind=ii*days
            smb_out[ii]=np.sum(smb2[ind:ind+days])*365./days/1.e3/0.917 #this is converting the N-day total accumulation to m ICE per year rate
            
        for jj in range(len(t2m_out)):
            ind=jj*days
            t2m_out[jj]=np.mean(t2m2[ind:ind+days])
            
        smbfile='smb_%s_%s.csv' % (mm,nn)
        tempfile='t2m_%s_%s.csv' % (mm,nn)    
        
        outpath1=os.path.join(saveloc,smbfile)
        outpath2=os.path.join(saveloc,tempfile)
        
        with open(outpath1, "wb") as ff:
            writer=csv.writer(ff)
            writer.writerow(tt_out)
            writer.writerow(smb_out)
            
        with open(outpath2, "wb") as ff2:
            writer=csv.writer(ff2)
            writer.writerow(tt_out2)
            writer.writerow(t2m_out)
            
            
#np.savetxt('smb2.csv',tt,fmt='delimiter=',')




