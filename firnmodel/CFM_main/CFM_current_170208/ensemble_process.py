import matplotlib.pyplot as plt 
import h5py as h5
import os
import sys
import numpy as np
from shutil import rmtree

thefiles2=['r0','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10','r11','r12','r13','r14','r15','r16','r17','r18','r19','r20','r21','r22','r23','r24','r25','r26','r27','r28','r29','r30','r31','r32','r33','r34','r35','r36','r37','r38','r39']

thefilesr=list(reversed(thefiles2))

thenames=['HLdynamic','Helsen2008','Arthern2010S','Arthern2010T','Simonsen2013','Ligtenberg2011','Barnola1991','KuipersMunneke2015','Li2011','Goujon2003','Crocus']

# thenames=['HLdynamic']

thesites= ['CRAWFORD','DYE2','EGRIP','EKT','KANU','NASASE','SADDLE','Summit']

datasource = 'MAR'

def stdadd(n,std_old,x_new,mean_new,mean_old):
	s2 = ((n-1)*(std_old)**2 + (x_new - mean_new)*(x_new - mean_old)) / n
	s=np.sqrt(s2)
	return s

# with open("missing_files.txt", "w") as tf:
# 	tf.write('missing files are')

for idx0, site in enumerate(thesites):
	for idx1, name in enumerate(thenames):
		nomean = False
		for idx2, fil in enumerate(thefilesr):
			try:
				fn = '/Volumes/FirnSSD/CFMresults/%sresults_ens_all/%s/%s/%s/CFMresults.hdf5' %(datasource, site, fil, name)
				# fn = '%sresults_ens_all/%s/%s/%s/CFMresults.hdf5' %(datasource, site, fil, name)

				f = h5.File(fn,'r')

			except:
				nomean = True
				with open("missing_files.txt", "a+") as tf:
					tf.write('################')
					# tf.write('no file')
					tf.write('site= %s' %site)
					tf.write('name= %s' %name)
					tf.write('file number= %s' %fil)
					tf.write(fn)
					tf.write('################')
				break

			time = f['depth'][1:,0]
			depth = f['depth'][1:,1:]
			density = f['density'][1:,1:]
			crate = f['compaction_rate'][1:,1:]
			try:
				age = f['age'][1:,1:]
			except:
				pass
				# print('no age file')

			DIP = f['DIP'][1:,1:]
			BCO = f['BCO'][1:,1:]

			if idx2==0:
				depth_all = depth
				density_mean = density
				density_std = np.zeros_like(density_mean)
				crate_mean = crate
				crate_std = np.zeros_like(crate_mean)
				try:
					age_mean = age
					age_std = zeros_like(age_mean)
				except:
					pass
				DIP_mean = DIP
				DIP_std = np.zeros_like(DIP_mean)
				BCO_mean = BCO
				BCO_std = np.zeros_like(BCO_mean)

				rD,cD = np.shape(DIP_mean)
				DIP_all = np.zeros((rD,cD,len(thefiles2)))
				DIP_all[:,:,idx2] = DIP

			else:
				n = idx2 + 1

				row,col = np.shape(depth_all)
				# print(r)
				density_interp = np.zeros_like(depth_all)
				crate_interp = np.zeros_like(depth_all)
				try:
					age_interp = np.zeros_like(depth_all)
				except:
					pass

				for r in range(row):
					density_interp[r,:] = np.interp(depth_all[r,:],depth[r,:],density[r,:])
					crate_interp[r,:] = np.interp(depth_all[r,:],depth[r,:],crate[r,:])					
					try:
						age_interp[r,:] = np.interp(depth_all[r,:],depth[r,:],age[r,:])
					except:
						pass

				d_mean_old = density_mean
				density_mean = density_mean + (density_interp - density_mean)/(idx2+1)
				density_std = stdadd(n, density_std, density_interp, density_mean, d_mean_old)

				c_mean_old = crate_mean
				crate_mean = crate_mean + (crate_interp - crate_mean)/(idx2+1)
				crate_std = stdadd(n, crate_std, crate_interp, crate_mean, c_mean_old)

				try:
					a_mean_old = age_mean
					age_mean = age_mean + (age_interp - density_mean) / (idx2+1)
					age_std = stdadd(n, age_std, age_interp, age_mean, a_mean_old)

				except:
					pass

				BCO_mean_old = BCO_mean
				BCO_mean = BCO_mean + (BCO-BCO_mean)/(idx2+1)
				BCO_std = stdadd(n, BCO_std, BCO, BCO_mean, BCO_mean_old)

				DIP_mean_old = DIP_mean
				DIP_mean = DIP_mean + (DIP-DIP_mean)/(idx2+1)
				DIP_std = stdadd(n, DIP_std, DIP, DIP_mean, DIP_mean_old)

				DIP_all[:,:,idx2] = DIP

				f.close()

		if not nomean:
			density_out = np.c_[time,density_mean]
			density_std_out = np.c_[time,density_std]
			crate_out = np.c_[time,crate_mean]
			crate_std_out = np.c_[time,crate_std]
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
			rfolder = '/Volumes/Samsung_T1/CFMresults/MAR_Summit/%sresults_ens_all/%s/ensmean/%s' %(datasource, site, name)
			if os.path.exists(rfolder):
				rmtree(rfolder)
			os.makedirs(rfolder)

			f4 = h5.File('/Volumes/Samsung_T1/CFMresults/MAR_Summit/%sresults_ens_all/%s/ensmean/%s/CFMresults_ens_mean.hdf5' %(datasource, site, name),'w')

			f4.create_dataset('density',data=density_out)
			f4.create_dataset('density_std',data=density_std_out)
			f4.create_dataset('crate',data=crate_out)
			f4.create_dataset('crate_std',data=crate_std_out)
			try:
				f4.create_dataset('age',data=age_out)
				f4.create_dataset('age_std',data=age_std_out)
			except:
				pass
			f4.create_dataset('BCO',data=BCO_out)
			f4.create_dataset('BCO_std',data=BCO_std_out)
			f4.create_dataset('DIP',data=DIP_out)
			f4.create_dataset('DIP_std',data=DIP_std_out)

			f4.close()





