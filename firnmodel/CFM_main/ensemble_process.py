'''
9/8/17
use this script to process the ensemble runs - create mean and standard deviation fields in a new file so all files don't have to be moved around. 
'''

import matplotlib.pyplot as plt 
import h5py as h5
import os
import sys
import numpy as np
from shutil import rmtree

# thefiles2=['r0','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10','r11','r12','r13','r14','r15','r16','r17','r18','r19','r20','r21','r22','r23','r24','r25','r26','r27','r28','r29','r30','r31','r32','r33','r34','r35','r36','r37','r38','r39']
# thefiles2 = ['r0', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r9', 'r10', 'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18', 'r19', 'r20', 'r21', 'r22', 'r23', 'r24', 'r25', 'r26', 'r27', 'r28', 'r29', 'r30', 'r31', 'r32', 'r33', 'r34', 'r35', 'r36', 'r37', 'r38', 'r39', 'r40', 'r41', 'r42', 'r43', 'r44', 'r45', 'r46', 'r47', 'r48', 'r49', 'r50', 'r51', 'r52', 'r53', 'r54', 'r55', 'r56', 'r57', 'r58', 'r59', 'r60', 'r61', 'r62', 'r63', 'r64', 'r65', 'r66', 'r67', 'r68', 'r69', 'r70', 'r71', 'r72', 'r73', 'r74', 'r75', 'r76', 'r77', 'r78', 'r79', 'r80', 'r81', 'r82', 'r83', 'r84', 'r85', 'r86', 'r87', 'r88', 'r89', 'r90', 'r91', 'r92', 'r93', 'r94', 'r95', 'r96', 'r97', 'r98', 'r99']

thefiles2 = ['r0', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r9', 'r10', 'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18', 'r19', 'r20', 'r21', 'r22', 'r23', 'r24', 'r25', 'r26', 'r27', 'r28', 'r29', 'r30', 'r31', 'r32', 'r33', 'r34', 'r35', 'r36', 'r37', 'r38', 'r39', 'r40', 'r41', 'r42', 'r43', 'r44', 'r45', 'r46', 'r47', 'r48', 'r49']

# thefiles2 = ['r50', 'r51', 'r52', 'r53', 'r54', 'r55', 'r56', 'r57', 'r58', 'r59', 'r60', 'r61', 'r62', 'r63', 'r64', 'r65', 'r66', 'r67', 'r68', 'r69', 'r70', 'r71', 'r72', 'r73', 'r74', 'r75', 'r76', 'r77', 'r78', 'r79', 'r80', 'r81', 'r82', 'r83', 'r84', 'r85', 'r86', 'r87', 'r88', 'r89', 'r90', 'r91', 'r92', 'r93', 'r94', 'r95', 'r96', 'r97', 'r98', 'r99']

thefilesr=list(reversed(thefiles2))

thenames=['HLdynamic','Helsen2008','Arthern2010S','Arthern2010T','Simonsen2013','Ligtenberg2011','Barnola1991','KuipersMunneke2015','Li2011','Goujon2003','Crocus']

# thenames=['HLSigfus']

# thesites= ['CRAWFORD','DYE2','EGRIP','EKT','KANU','NASASE','SADDLE','Summit']
thesites = ['Summit']
# thesites = ['WAISDvarden']

cont = 'Greenland'

datasource = 'RACMO'

rootfolder = '/wd1/wrk/maxstev/CommunityFirnModel/firnmodel/CFM_main'

writer = True

def stdadd(n,std_old,x_new,mean_new,mean_old):
	s2 = ((n-1)*(std_old)**2 + (x_new - mean_new)*(x_new - mean_old)) / n
	s=np.sqrt(s2)

	return s

with open("missing_files_%s.txt" %datasource, "w") as tf:
	tf.write('missing files are' + '\n')

for idx0, site in enumerate(thesites):
	for idx1, name in enumerate(thenames):
		# nomean = False
		counter = 0
		for idx2, fil in enumerate(thefilesr):
			try:
				fn = rootfolder + '/%sresults_ens_all/%s/%s/%s/%s/CFMresults.hdf5' %(datasource, cont, site, fil, name)
				# fn = '%sresults_ens_all/%s/%s/%s/CFMresults.hdf5' %(datasource, site, fil, name)

				f = h5.File(fn,'r')
				print('file is %s' %fil)

			except:
				# nomean = True
				print('no file ', fn)
				with open("missing_files_%s.txt" %datasource, "a+") as tf:
					tf.write('################' + '\n')
					# tf.write('no file')
					tf.write('site= %s' %site + '\n')
					tf.write('name= %s' %name + '\n')
					tf.write('file number= %s' %fil + '\n')
					tf.write(fn + '\n')
					tf.write('################' + '\n')
				continue

			time = f['depth'][1:,0]
			depth = f['depth'][1:,1:]
			density = f['density'][1:,1:]
			temperature = f['temperature'][1:,1:]
			compaction_rate = f['compaction_rate'][1:,1:]
			try:
				age = f['age'][1:,1:]
			except:
				pass
				# print('no age file')

			DIP = f['DIP'][1:,1:]
			BCO = f['BCO'][1:,1:]

			if counter==0:
				print('initializing arrays at ', fil)
				depth_all = depth
				density_mean = density
				density_std = np.zeros_like(density_mean)
				temperature_mean = temperature
				temperature_std = np.zeros_like(temperature_mean)
				compaction_rate_mean = np.zeros_like(depth)
				rcr,ccr=np.shape(compaction_rate)
				compaction_rate_mean[0:rcr,0:ccr] = compaction_rate
				compaction_rate_std = np.zeros_like(compaction_rate_mean)
				try:
					age_mean = age
					age_std = zeros_like(age_mean)
				except:
					pass
				DIP_mean = DIP
				DIP_std = np.zeros_like(DIP_mean)
				BCO_mean = BCO
				BCO_std = np.zeros_like(BCO_mean)
				counter += 1

				# rD,cD = np.shape(DIP_mean)
				# DIP_all = np.zeros((rD,cD,len(thefiles2)))
				# DIP_all[:,:,idx2] = DIP

			else:
				# n = idx2 + 1
				n = counter + 1

				row,col = np.shape(depth_all)
				rcr,ccr = np.shape(compaction_rate)
				# print(r)
				density_interp = np.zeros_like(depth_all)
				temperature_interp = np.zeros_like(depth_all)
				compaction_rate_interp = np.zeros_like(depth_all)
				try:
					age_interp = np.zeros_like(depth_all)
				except:
					pass

				for r in range(row):
					density_interp[r,:] = np.interp(depth_all[r,:],depth[r,:],density[r,:])
					temperature_interp[r,:] = np.interp(depth_all[r,:],depth[r,:],temperature[r,:])
					crm = np.zeros_like(depth[r,:])
					crm[0:ccr] = compaction_rate[r,:]
					compaction_rate_interp[r,:] = np.interp(depth_all[r,:],depth[r,:],crm)				
					try:
						age_interp[r,:] = np.interp(depth_all[r,:],depth[r,:],age[r,:])
					except:
						pass

				d_mean_old = density_mean
				density_mean = density_mean + (density_interp - density_mean)/(n)
				density_std = stdadd(n, density_std, density_interp, density_mean, d_mean_old)

				t_mean_old = temperature_mean
				temperature_mean = temperature_mean + (temperature_interp - temperature_mean)/(n)
				temperature_std = stdadd(n, temperature_std, temperature_interp, temperature_mean, t_mean_old)

				c_mean_old = compaction_rate_mean
				compaction_rate_mean = compaction_rate_mean + (compaction_rate_interp - compaction_rate_mean)/(n)
				compaction_rate_std = stdadd(n, compaction_rate_std, compaction_rate_interp, compaction_rate_mean, c_mean_old)

				try:
					a_mean_old = age_mean
					age_mean = age_mean + (age_interp - density_mean) / (n)
					age_std = stdadd(n, age_std, age_interp, age_mean, a_mean_old)

				except:
					pass

				BCO_mean_old = BCO_mean
				BCO_mean = BCO_mean + (BCO-BCO_mean)/(n)
				BCO_std = stdadd(n, BCO_std, BCO, BCO_mean, BCO_mean_old)

				DIP_mean_old = DIP_mean
				DIP_mean = DIP_mean + (DIP-DIP_mean)/(n)
				DIP_std = stdadd(n, DIP_std, DIP, DIP_mean, DIP_mean_old)

				# DIP_all[:,:,idx2] = DIP
				counter += 1

				f.close()

		if writer:
			depth_out = np.c_[time,depth_all]
			density_out = np.c_[time,density_mean]
			density_std_out = np.c_[time,density_std]
			temperature_out = np.c_[time,temperature_mean]
			temperature_std_out = np.c_[time,temperature_std]			
			compaction_rate_out = np.c_[time,compaction_rate_mean]
			compaction_rate_std_out = np.c_[time,compaction_rate_std]
			try:
				age_out = np.c_[time,age_mean]
				age_std_out = np.c_[time,age_std]
			except:
				pass

			BCO_out = np.c_[time,BCO_mean]
			BCO_std_out = np.c_[time,BCO_std]
			DIP_out = np.c_[time,DIP_mean]
			DIP_std_out = np.c_[time,DIP_std]

			# f4 = h5.File('%sresults_ens_all/%s/ensmean/%s/CFMresults_ens_mean.hdf5' %(datasource, site, name))
			rfolder = rootfolder + '/%sresults_ens_all/%s/%s/ensmean049/%s' %(datasource, cont, site, name)
			if os.path.exists(rfolder):
				rmtree(rfolder)
			os.makedirs(rfolder)

			f4 = h5.File(rootfolder + '/%sresults_ens_all/%s/%s/ensmean049/%s/CFMresults.hdf5' %(datasource, cont, site, name),'w')

			f4.create_dataset('depth', data=depth_out)
			f4.create_dataset('density',data=density_out)
			f4.create_dataset('density_std',data=density_std_out)
			f4.create_dataset('temperature',data=temperature_out)
			f4.create_dataset('temperature_std',data=temperature_std_out)			
			f4.create_dataset('compaction_rate',data=compaction_rate_out)
			f4.create_dataset('compaction_rate_std',data=compaction_rate_std_out)
			try:
				f4.create_dataset('age',data=age_out)
				f4.create_dataset('age_std',data=age_std_out)
			except:
				pass
			f4.create_dataset('BCO',data=BCO_out)
			f4.create_dataset('BCO_std',data=BCO_std_out)
			f4.create_dataset('DIP',data=DIP_out)
			f4.create_dataset('DIP_std',data=DIP_std_out)
			f4.create_dataset('counter',data=counter)

			f4.close()





