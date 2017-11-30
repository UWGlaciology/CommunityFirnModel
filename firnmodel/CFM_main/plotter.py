import matplotlib.pyplot as plt 
import h5py as h5
import os
import sys
import numpy as np
import scipy as sp
import pickle

def plotter(rfolder,rfile):

	'''
	Plot results from the CFM
	'''
	fn = os.path.join(rfolder,rfile)

	f = h5.File(fn,'r')

	timesteps = f['depth'][1:,0]

	# with open('/Users/maxstev/Documents/Grad_School/Research/FIRN/GREENLAND_CVN/Data/CVN_DATA/core_data_dict.pkl','rb') as ff:
	# 	d=pickle.load(ff)

	# core = '2016_2'
	# cdepth = d[core]['depth']/100
	# cdensity = d[core]['density']

	depth = f['depth'][1:,1:]
	density = f['density'][1:,1:]
	temperature = f['temperature'][1:,1:]
	# air = f['gasses'][:,1:]

	jj=np.where(timesteps>=200.0)[0][0]

	f1=plt.figure(1)
	ax1 = f1.add_subplot(111)
	plt.plot(density[jj,1:],depth[jj,1:])
	plt.plot(cdensity,cdepth,'r')
	# plt.ylim([0,50])
	plt.gca().invert_yaxis()
	stri = 'date = %.6f' % timesteps[jj]
	plt.gca().text(0.2, 0.1,stri,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax1.transAxes)
	# plt.text(0.2,0.1,stri)
	# plt.savefig('melt.eps')

	plt.show()



	# niter = np.shape(depth)[0]

	# f2=plt.figure(2)
	# ax2 = f2.add_subplot(111)
	# plt.ion()
	# for ii in range(niter):
	# 	stri = '%.2f' % timesteps[ii]
	# 	plt.clf()
	# 	plt.plot(density[ii,:],depth[ii,:])
	# 	# plt.plot(depth[ii,:],(air[ii,:]-1)*1000)
		
	# 	# plt.ylim([0,50])
	# 	# plt.gca().invert_yaxis()
	# 	plt.gca().text(0.2, 0.1,stri, horizontalalignment='center', verticalalignment='center', transform = ax2.transAxes)
	# 	plt.pause(0.05)

	# 	plt.show()
	# while True:
	# 	plt.pause(0.05)

		
	# ap = (air[jj,1:]-1)*1000

	n_grid = np.arange(0,np.ceil(depth[0,-1]),0.1)
	ro,co = np.shape(temperature)
	t_interp = np.zeros((ro,len(n_grid)))
	# f_i = sp.interpolate.interp1d()
	for jj in range(ro):
		yy=np.interp(n_grid,depth[jj,:],temperature[jj,:])
		t_interp[jj,:]=yy

	# t_plot1 = t_interp[:,0:1000]
	# t_plot = t_plot1.T
	# n_plot = n_grid[0:1000]


	# f3, ax3 = plt.subplots()
	# v=np.linspace(240,273.15,10,endpoint=True)
	# cax=ax3.contourf(timesteps,n_plot,t_plot,v)
	# ax3.set_ylim([0,20])
	# ax3.invert_yaxis()
	# bounds = [240,273.15]
	# f3.colorbar(cax,ticks=v)
	# # ax3.contourf(t_interp[0:1000])

	# plt.show()








if __name__ == '__main__':

	rfolder = '/Users/maxstev/Documents/Grad_School/Research/FIRN/CFM/CommunityFirnModel/firnmodel/CFM_main/CFM_current_170208/melt_test_kanu_KP' #path to the results folder

	# rfolder = '/Volumes/FirnSSD/CFMresults/melt_test_kanu_KP'

	rfile = 'CFMresults.hdf5'

	plotter(rfolder,rfile)