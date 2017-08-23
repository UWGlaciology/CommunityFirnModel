import matplotlib.pyplot as plt 
import h5py as h5
import os
import sys
import numpy as np

def plotter(rfolder,rfile):

	'''
	Plot results from the CFM
	'''
	fn = os.path.join(rfolder,rfile)

	f = h5.File(fn,'r')

	timesteps = f['depth'][1:,0]



	depth = f['depth'][:,1:]
	density = f['density'][:,1:]
	air = f['gasses'][:,1:]

	jj=-1

	# f1=plt.figure(1)
	# ax1 = f1.add_subplot(111)
	# plt.plot(density[jj,1:],depth[jj,1:])
	# # plt.ylim([0,50])
	# plt.gca().invert_yaxis()
	# stri = 'date = %.6f' % timesteps[jj]
	# plt.gca().text(0.2, 0.1,stri,
 #     horizontalalignment='center',
 #     verticalalignment='center',
 #     transform = ax1.transAxes)
	# # plt.text(0.2,0.1,stri)
	# plt.savefig('melt.eps')

	# plt.show()

	niter = np.shape(depth)[0]

	# f2=plt.figure(2)
	# ax2 = f2.add_subplot(111)
	# plt.ion()
	# for ii in range(niter):
	# 	stri = '%.2f' % timesteps[ii]
	# 	plt.clf()
	# 	# plt.plot(density[ii,:],depth[ii,:])
	# 	plt.plot(depth[ii,:],(air[ii,:]-1)*1000)
		
	# 	# plt.ylim([0,50])
	# 	# plt.gca().invert_yaxis()
	# 	plt.gca().text(0.2, 0.1,stri, horizontalalignment='center', verticalalignment='center', transform = ax2.transAxes)
	# 	plt.pause(0.05)

	# 	plt.show()
	# while True:
	# 	plt.pause(0.05)

		
	ap = (air[jj,1:]-1)*1000

	f3, ax3 = plt.subplots()
	ax3.plot(depth[jj,1:],ap)

	plt.show()







if __name__ == '__main__':

	rfolder = 'air_test'
	rfile = 'air_test.hdf5'

	plotter(rfolder,rfile)