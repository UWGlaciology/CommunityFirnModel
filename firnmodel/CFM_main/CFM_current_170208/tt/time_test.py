import numpy as np
import scipy as sp
import os
import csv
import time
import h5py
import tables

def time_test1(stps):

	vecpath = 'vec1.csv'

	vec=np.arange(350.0,500.0,1.0)

	with open(vecpath, 'w') as f:
		writer = csv.writer(f)
		writer.writerow(vec)

	stp = stps

	vec_out = np.zeros((stp,len(vec)))

	for ii in xrange(stp):

		dr = vec*0.01
		vec = vec + dr
		vec = np.concatenate(([350.0],vec[0:-1]))
		vec_out[ii-1,:] = vec

		with open(vecpath, 'a') as f:
			writer = csv.writer(f)
			writer.writerow(vec)

def time_test2(stps):
	vecpath = 'vec2.csv'

	vec=np.arange(350.0,500.0,1.0)

	with open(vecpath, 'w') as f1:
		writer = csv.writer(f1)
		writer.writerow(vec)

		stp = stps

		vec_out = np.zeros((stp,len(vec)))

		for ii in xrange(stp):

			dr = vec*0.01
			vec = vec + dr
			vec = np.concatenate(([350.0],vec[0:-1]))
			vec_out[ii-1,:] = vec

			writer.writerow(vec)

def time_test3(stps):

	vecpath = 'vec3.csv'

	vec=np.arange(350.0,500.0,1.0)

	stp = stps

	vec_out = np.zeros((stp,len(vec)))

	with open(vecpath, 'w') as f2:
		writer = csv.writer(f2)
		writer.writerow(vec)

	for ii in xrange(stp):

		dr = vec*0.01
		vec = vec + dr
		vec = np.concatenate(([350.0],vec[0:-1]))

		vec_out[ii-1,:] = vec

	with open(vecpath, 'a') as f2:
		writer = csv.writer(f2)
		writer.writerow(vec_out)

def time_test4(stps):

	f4 = h5py.File('vec4.hdf5','w')

	vec=np.arange(350.0,500.0,1.0)

	stp = stps

	vec_out = np.zeros((stp+1,len(vec)))
	vec_out[0,:]=vec

	for ii in xrange(stp):

		dr = vec*0.01
		vec = vec + dr
		vec = np.concatenate(([350.0],vec[0:-1]))

		vec_out[ii+1,:] = vec

	f4.create_dataset('vec4', shape=np.shape(vec_out), data=vec_out, dtype='f8')

	f4.close()

def time_test5(stps):

	f5 = h5py.File('vec5.hdf5','w')


	vec=np.arange(350.0,500.0,1.0)

	stp = stps

	vec_out = np.zeros((stp+1,len(vec)))
	dset = f5.create_dataset('vec5', shape=np.shape(vec_out), dtype='f8')
	dset[0,:]=vec

	for ii in xrange(stp):

		dr = vec*0.01
		vec = vec + dr
		vec = np.concatenate(([350.0],vec[0:-1]))

		# vec_out[ii-1,:] = vec
		dset[ii+1,:] = vec
		

	f5.close()

def time_test6(stps):

	f6 = h5py.File('vec6.hdf5','w')

	vec=np.arange(350.0,500.0,1.0)

	stp = stps

	lv=len(vec)

	vw = np.zeros((100,lv))

	dset = f6.create_dataset('vec6', shape=np.shape(stp+1,len(vec)), dtype='f8')
	vw[0,:] = vec
	counter = 1

	for ii in xrange(stp):

		dr = vec*0.01
		vec = vec + dr
		vec = np.concatenate(([350.0],vec[0:-1]))

		vw[counter,:] = vec
		# dset[ii+1,:] = vec
		if counter ==

		else:
			counter += 1 

		dset = f6.create_dataset('hh%d' %w, data=vec)
		
	# dall = f6.create_dataset('vec6',shape=(stp+1,len(vec)))
	# for jj in xrange(stp+1):
	# 	dall[jj,:]=f6['hh%d' %jj][:]

	f6.close()



if __name__ == '__main__':

	stps = 10000
	t1a = time.time()
	time_test1(stps)
	t1b = time.time()
	ttot1 = t1b - t1a
	print 'ttot1 =', ttot1

	t2a = time.time()
	time_test2(stps)
	t2b = time.time()
	ttot2 = t2b - t2a
	print 'ttot2 =', ttot2

	# t3a = time.time()
	# time_test3(stps)
	# t3b = time.time()
	# ttot3 = t3b - t3a
	# print 'ttot3 =', ttot3

	t4a = time.time()
	time_test4(stps)
	t4b = time.time()
	ttot4 = t4b - t4a
	print 'ttot4 =', ttot4

	t5a = time.time()
	time_test5(stps)
	t5b = time.time()
	ttot5 = t5b - t5a
	print 'ttot5 =', ttot5

	t6a = time.time()
	time_test6(stps)
	t6b = time.time()
	ttot6 = t6b - t6a
	print 'ttot6 =', ttot6




