import matplotlib.pyplot as plt 
import h5py as h5
import os
import sys
import numpy as np
from shutil import rmtree

thefiles2=['r0','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10','r11','r12','r13','r14','r15','r16','r17','r18','r19','r20','r21','r22','r23','r24','r25','r26','r27','r28','r29','r30','r31','r32','r33','r34','r35','r36','r37','r38','r39']

thefilesr=list(reversed(thefiles2))

thenames=['HLdynamic','Helsen2008','Arthern2010S','Arthern2010T','Simonsen2013','Ligtenberg2011','Barnola1991','KuipersMunneke2015','Li2011','Goujon2003','Crocus']

thesites= ['CRAWFORD','DYE2','EGRIP','EKT','KANU','NASASE','SADDLE','Summit']

datasource = 'MAR'

for idx0, site in enumerate(thesites):
	for idx1, name in enumerate(thenames):
		noavg = False
		for idx2, fil in enumerate(thefilesr):
			try:
				fn = '/Volumes/Samsung_T1/CFMresults/MAR_Summit/%sresults_ens_all/%s/%s/%s/CFMresults.hdf5' %(datasource, site, fil, name)
				# fn = '%sresults_ens_all/%s/%s/%s/CFMresults.hdf5' %(datasource, site, fil, name)

				f = h5.File(fn,'r')

			except:
				noavg = True
				print('################')
				print('no file')
				print('site', site)
				print('name', name)
				print('file number', fil)
				print(fn)
				print('################')
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
				density_avg = density
				crate_avg = crate
				try:
					age_avg = age
				except:
					pass
				DIP_avg = DIP
				BCO_avg = BCO

			else:
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

				density_avg = density_avg + (density_interp - density_avg)/(idx2+1)
				crate_avg = crate_avg + (crate_interp - crate_avg)/(idx2+1)
				try:
					age_avg = age_avg + (age_interp - density_avg) / (idx2+1)
				except:
					pass

				BCO_avg = BCO_avg + (BCO-BCO_avg)/(idx2+1)
				DIP_avg = DIP_avg + (DIP-DIP_avg)/(idx2+1)

				f.close()





		if not noavg:
			density_out = np.c_[time,density_avg]
			crate_out = np.c_[time,crate_avg]
			try:
				age_out = np.c_[time,age_avg]
			except:
				pass

			BCO_out = np.c_[time,BCO_avg]
			DIP_out = np.c_[time,DIP_avg]

			# f4 = h5.File('%sresults_ens_all/%s/ensavg/%s/CFMresults_ens_avg.hdf5' %(datasource, site, name))
			rfolder = '/Volumes/Samsung_T1/CFMresults/MAR_Summit/%sresults_ens_all/%s/ensavg/%s' %(datasource, site, name)
			if os.path.exists(rfolder):
				rmtree(rfolder)
			os.makedirs(rfolder)

			f4 = h5.File('/Volumes/Samsung_T1/CFMresults/MAR_Summit/%sresults_ens_all/%s/ensavg/%s/CFMresults_ens_avg.hdf5' %(datasource, site, name),'w')

			f4.create_dataset('density',data=density_out)
			f4.create_dataset('crate',data=crate_out)
			try:
				f4.create_dataset('age',data=age_out)
			except:
				pass
			f4.create_dataset('BCO',data=BCO_out)
			f4.create_dataset('DIP',data=DIP_out)

			f4.close()





