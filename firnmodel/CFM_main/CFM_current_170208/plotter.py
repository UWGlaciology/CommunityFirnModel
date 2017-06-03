import matplotlib.pyplot as plt 
import h5py as h5
import os
import sys

def plotter(rfolder,rfile):

	'''
	Plot results from the CFM
	'''
	fn = os.path.join(rfolder,rfile)

	f = h5.File(fn,'r')

	timesteps = f['depth'][1:,0]



	depth = f['depth'][-1,1:]
	density = f['density'][-1,1:]

	plt.figure(1)
	plt.plot(density,depth)
	plt.ylim([0,50])
	plt.gca().invert_yaxis()

	plt.show()







if __name__ == '__main__':

	rfolder = 'melt_test'
	rfile = 'melt.hdf5'

	plotter(rfolder,rfile)