import matplotlib.pyplot as plt 
import h5py as h5
import os
import sys
import numpy as np
import scipy as sp
import pickle

rfolder = '/Users/maxstev/Documents/Grad_School/Research/FIRN/CFM/CommunityFirnModel/firnmodel/CFM_main/CFM_current_170208/melt_test_DYE2_KP'
rfile = 'CFMresults.hdf5'

fn = os.path.join(rfolder,rfile)

f = h5.File(fn,'r')

LWC = f['LWC'][10000:,1:3000]
lbin = np.copy(LWC)
lbin[lbin>0.0]=1
lbin[lbin<=0.0]  =0
T = f['temperature'][10000:,1:3000]
T[T<230] = 230
depth = f['depth'][10000:,1:3000]
time = f['depth'][10000:,0]

n_grid = np.arange(0,np.ceil(depth[0,-1]),0.1)
ro,co = np.shape(LWC)
lwc_interp = np.zeros((ro,len(n_grid)))
t_interp = np.zeros_like(lwc_interp)
# f_i = sp.interpolate.interp1d()
for jj in range(ro):
	yy=np.interp(n_grid,depth[jj,:],LWC[jj,:])
	zz=np.interp(n_grid,depth[jj,:],T[jj,:])
	lwc_interp[jj,:]=yy
	t_interp[jj,:]=zz

lbin_interp=np.copy(lwc_interp)
lbin_interp[lbin_interp>0.0] = 1
lbin_interp[lbin_interp<=0.0]  = 0


plt.ion()
f3, ax3 = plt.subplots()
cax3=ax3.contourf(time,n_grid,t_interp.T,256)#,extend='both')
# cax3.cmap.set_under('k')
# cax3.cmap.set_over('r')
ax3.invert_yaxis()
# cax3.set_clim(240,273.0)
f3.colorbar(cax3)

# f2, ax2 = plt.subplots()
# cax2 = ax2.contourf(LWC)
# cax2.set_clim(1.0e-5,0.015)
# f2.colorbar(cax2)

# ax3.set_ylim([0,20])
# ax3.invert_yaxis()
# bounds = [240,273.15]
# f3.colorbar(cax,ticks=v)
# ax3.contourf(t_interp[0:1000])

f1, ax1 = plt.subplots()
cax1 = ax1.contourf(time,n_grid,lbin_interp.T)
ax1.invert_yaxis()
# cax1.set_clim(1.0e-5,0.015)
f1.colorbar(cax1)

plt.show()