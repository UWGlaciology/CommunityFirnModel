from constants import *
import numpy as np

def simple_melt(self, iii):

	
    if self.bdotSec[iii]<=self.snowmeltSec[iii]: # more melt than accumulation; no new accumulation; runoff from old box

        meltdiff = self.snowmeltSec[iii] - self.bdotSec[iii] # units m ice.
        # melt_mass = meltdiff * RHO_I *S_PER_YEAR

        self.age = self.age + self.dt
        self.dz_old = self.dz
        self.sdz_old = np.sum(self.dz) # old total column thickness
        self.z_old = self.z
        self.dz = self.mass / self.rho * self.dx
        self.sdz_new = np.sum(self.dz) #total column thickness after densification, before new snow addedv
        self.z = self.dz.cumsum(axis = 0)
        self.z = self.z - self.z[0]
        self.dzNew = 0
        # if self.rho[0]>=800: # old surface is ice, so just runoff
        #     self.dz = self.mass / self.rho * self.dx
        #     self.sdz_new = np.sum(self.dz) #total column thickness after densification, before new snow added


        # else:
        
    elif self.bdotSec[iii]>self.snowmeltSec[iii]: # more accumulation than melt; all melt is from that year's accumulation; remaining accumulation is added
	    self.age = np.concatenate(([0,0],self.age[:-2]))+self.dt
	    self.dz_old = self.dz
	    self.sdz_old = np.sum(self.dz) # old total column thickness
	    self.z_old = self.z
	    self.dzNewSnow = (self.bdotSec[iii]-self.snowmeltSec[iii]) * RHO_I / self.rhos0[iii] * S_PER_YEAR
	    self.dzNewIce = self.snowmeltSec[iii] * S_PER_YEAR
	    self.dz = self.mass / self.rho * self.dx
	    self.sdz_new = np.sum(self.dz) #total column thickness after densification, before new snow added
	    self.dz = np.concatenate(([self.dzNewSnow,self.dzNewIce], self.dz[:-2]))
	    self.z = self.dz.cumsum(axis = 0)
	    self.z = np.concatenate(([0], self.z[:-1]))
	    self.rhos0[iii]=400.0
	    self.rho  = np.concatenate(([self.rhos0[iii],RHO_I], self.rho[:-2])) 
	    massNewSnow = self.bdotSec[iii] * S_PER_YEAR * RHO_I
	    massNewIce = self.snowmeltSec[iii] * S_PER_YEAR *RHO_I
	    self.mass = np.concatenate(([massNewSnow,massNewIce], self.mass[:-2]))
	    self.dzNew = self.dzNewIce + self.dzNewSnow

def percolation(self, iii):
	
	# print('!!! iii ', iii)
	porosity = 1 - self.rho/RHO_I #porosity
	dz_old = self.dz
	# print 'porosity', porosity

	porespace_vol 			= porosity * self.dz #porosity in meters of each box


	melt_volume_IE  	= self.snowmeltSec[iii] * S_PER_YEAR #meters
	melt_volume_WE		= melt_volume_IE * 0.917 #meters
	melt_mass			= melt_volume_WE * 1000. #kg

	heat_to_freeze 		= melt_mass * LF_I #amount of heat needed to refreeze the melt (J)

	# print('melt_mass (orig)', melt_mass)
	if (self.mass_sum==melt_mass).any():
		exactmass = True
	else:
		exactmass = False

	ind1a = np.where(self.mass_sum<=melt_mass)[0] # indicies of boxes that will be melted away
	# if iii>1020:
	# 	print(iii)
	# 	print(self.T_mean)
	# 	print('melt_mass', melt_mass)
	# 	print('mass_sum', self.mass_sum[0:20])
	num_boxes_melted = len(ind1a)+1 #number of boxes that melt away, include the box that is partially melted
	# print('num_boxes_melted', num_boxes_melted)
	ind1 = np.where(self.mass_sum>melt_mass)[0][0] # index which will become the new surface
	# print('ind1a', ind1a)
	# print('num_boxes_melted',num_boxes_melted)
	# print('ind1', ind1)
	# print('self.mass_sum[ind1-1]',self.mass_sum[ind1])
	# print('melt_mass', melt_mass)

	# pm is the partial melt (the box/volume that has a portion melted away)
	pm_mass = self.mass_sum[ind1] - melt_mass # the remaining mass of the PM box
	pm_dz = pm_mass / self.rho[ind1] #remaining thickness
	pm_porespace = (1 - self.rho[ind1]/RHO_I) * pm_dz #porespace in the PM box
	pm_rho = self.rho[ind1] #density of the PM box

	# print('pm_mass', pm_mass)
	# print('pm_dz', pm_dz)
	# print('pm_porespace', pm_porespace)
	# print('pm_rho', pm_rho)

	cold_content_0 		= CP_I * self.mass * (T_MELT - self.Tz) #cold content of each box, i.e. how much heat to bring it to 273K
	cold_content_0_sum = cold_content_0.cumsum(axis=0)
	cold_content = cold_content_0[ind1:] #just the boxes that don't melt away
	cold_content[0]  = CP_I * pm_mass * (T_MELT - self.Tz[ind1]) # the partial melt box has its cold content reassigned.
	cold_content_sum 	= cold_content.cumsum(axis=0)
	# print('cold_content_sum', cold_content_sum)
	# print('heat_to_freeze', heat_to_freeze)

	ind2_rel = np.where(cold_content_sum>heat_to_freeze)[0][0] #freeze horizon index (where the wetting front freezes), index relative to ind1
	ind2 = ind2_rel + ind1 #absolute index on real grid (we have not removed the melted boxes yet)

	
	# print('rho!', self.rho[0:ind2])
	# print('rhomin', np.min(self.rho)

	if (self.rho[ind1:ind2+1]>830.0).any(): #if there is an ice lens somewhere between the new surface and where freezing should occur
		# print('!!!!!  ',iii)
		# print(self.rho[ind1:ind2+1])
		# print("blocking ice lens")
		# print('ind2 (old)', ind2)

		ind2_rel = np.where(self.rho[ind1:]>=830.0)[0][0]
		ind2 = ind2_rel + ind1 #the recalculated freezing front (water can not go past this level)
		# print('rhoind2',self.rho[ind2])

		# print('ind2 (new)', ind2)
		# print('ind2_rel (new)', ind2_rel)
		# cold_content_lens = cold_content_0[ind1:ind2+1].sum() #cold content that is is available in space between the surface and the lens 
		cold_content_lens = cold_content[0:ind2_rel].sum() #cold content that is is available in space between the surface and the lens 
		hh = heat_to_freeze - cold_content_lens 
		# hh = cold_content_lens # how much of the liquid can be refrozen with the available cold content
		
		refreeze_mass = hh / LF_I # this is how much should be able to refreeze in the available pore space above the ice lens.

		# print('melt_mass', melt_mass)
		# print('refreeze_mass', refreeze_mass)

		melt_volume_WE = refreeze_mass / 1000.
		melt_volume_IE = melt_volume_WE / 0.917 #the volume of melt that is not runoff
		runoff_volume_duetolens = melt_mass - refreeze_mass

	else:
		runoff_volume_duetolens	 = 0.0
		# print('ind2 (no change)', ind2)
		# print('ind2_rel (no change)', ind2_rel)

	# print('ind2',ind2)

	pore_indices = np.arange(ind1,ind2+1) # indicies of the boxes that are available to fill with water
	pore_indices_flip = np.flipud(pore_indices)

	# print('pore_indices', pore_indices)
	# print('dz', self.dz[ind1:ind2+1])

	porespace_vol[ind1] = pm_porespace
	porespace_0_sum		= porespace_vol.cumsum(axis=0)
	# print('porespace_0_sum', porespace_0_sum)
	# print('porespace_vol',porespace_vol)
	# print('ind1+1',ind1+1)
	# print('ind2+1',ind2+1)
	porespace = porespace_vol[ind1+1:ind2+1] #space available for the water

	# print 	'porespace1', porespace
	# porespace[0]=pm_porespace
	# print 'porespace2', porespace
	porespace_sum = porespace.cumsum(axis=0)
	# print('pss', porespace_sum)
	porespace_sum_flip = (np.flipud(porespace)).cumsum(axis=0)

	# if self.rho[ind1:ind2].any()>= 830.0:
	# 	''
	# ind2a = ind2 - ind1
	# print 'ind2a', ind2a

	# print 'porespace_sum', porespace_0_sum[ind2a]-porespace_0_sum[ind1]
	# print 'melt_volume', melt_volume_IE
	available_space = porespace_0_sum[ind2]-porespace_0_sum[ind1]

	# print('available_space', available_space)
	# print('melt_volume_IE',melt_volume_IE)

	# if porespace_sum[ind2a]<melt_volume_IE:
	if available_space < melt_volume_IE: # melt volume has already been recalculated based on how much can freeze with the cold content

		# if iii>1440:
		# 	 # or self.bdotSec[iii-1])<0:
		# 	print('iii',iii)
			
		# 	print('not enough pore space')
		# 	# if iii>1440:
		# 	print(self.bdotSec[iii]*S_PER_YEAR*12)
		# 	if self.bdotSec[iii]<0:
		# 		print('negative bdot')

		# 	print('num_boxes_melted', num_boxes_melted)
		# 	print('iii',iii)
		runoff_volume_duetolimitedporespace = (melt_volume_IE - porespace_sum) * 0.917

		self.rho[ind1:ind2+1] = 870.0 #fill all of the boxes with water.
		self.Tz[ind1:ind2+1] = T_MELT

		# split up last box into several
		divider = num_boxes_melted

		self.rho = np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted)))
		self.age = np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_melted)))
		self.dz  = np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted)))
		self.dz[0] = pm_dz
		self.dzn = np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
		self.dzn = self.dzn[0:self.compboxes]
		# if (self.dz<=0).any():
		# 	print self.dz[0:ind2]
		# 	raw_input('negative dz')
		self.Tz  = np.concatenate((self.Tz[ind1:-1],self.Tz[-1]*np.ones(num_boxes_melted)))
		self.bdot_mean = np.concatenate((self.bdot_mean[ind1:-1],self.bdot_mean[-1]*np.ones(num_boxes_melted)))
		self.z = self.dz.cumsum(axis = 0)
		self.z = np.concatenate(([0], self.z[:-1]))
		self.mass = self.rho*self.dz
		if len(self.mass)!=self.gridLen:
			print('num_boxes_melted', num_boxes_melted)
			print('divider',divider)
			print('ind1',ind1)
			print('ind1a',ind1a)
			print('line202')

	
	elif available_space == 0.0: #the top layer is an ice lens, so the melt runs off
		# print('ice lens on top')

		# print('rho',self.rho[0:ind1+4])

				# split up last box into several
		divider = num_boxes_melted
		# ind3 should be removed and replaced with 2 new boxes.
		self.rho = np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted)))
		self.age = np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_melted)))
		self.dz  = np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted)))
		self.dz[0] = pm_dz
		self.dzn = np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
		self.dzn = self.dzn[0:self.compboxes]
		# if (self.dz<=0).any():
		# 	print self.dz[0:ind2]
		# 	raw_input('negative dz')
		self.Tz  = np.concatenate((self.Tz[ind1:-1],self.Tz[-1]*np.ones(num_boxes_melted)))
		self.bdot_mean = np.concatenate((self.bdot_mean[ind1:-1],self.bdot_mean[-1]*np.ones(num_boxes_melted)))
		self.z = self.dz.cumsum(axis = 0)
		self.z = np.concatenate(([0], self.z[:-1]))
		self.mass = self.rho*self.dz
		if len(self.mass)!=self.gridLen:
			print('line224')



	else:
		# if self.bdotSec[iii]<0:
		# print('iii',iii)
		# print('enough pore space')
		# print('iii',iii)
		runoff_volume_duetolimitedporespace = 0
		# print('porespace_sum_flip',porespace_sum_flip[0:5])
		# print('melt_volume_IE',melt_volume_IE)
		ind3a = np.where(porespace_sum_flip>melt_volume_IE)[0][0]
		# print(ind3a)
		ind3 = ind2 - ind3a #the index of the node that is partially filled with water

		# print 'ind3',ind3
		# print 'ind3a',ind3a
		# print 'rho3', self.rho[ind3]

		partial_volume = melt_volume_IE - np.sum(porespace_vol[ind3+1:ind2+1]) # pore space filled in the box that is partially filled

		leftover_porespace = porespace_vol[ind3]-partial_volume #open pore space in the the partially-filled box

		new_node_1_rho = self.rho[ind3] #split up the partial box into 2 parts
		new_node_2_rho = 870.0

		new_node_1_dz = leftover_porespace / (1 - self.rho[ind3]/RHO_I)
		new_node_2_dz = self.dz[ind3] - new_node_1_dz

		# if new_node_1_dz<0 or new_node_2_dz<0:

		# 	print 'new_node_1_dz', new_node_1_dz
		# 	print 'new_node_2_dz', new_node_2_dz
		# 	print 'leftover_porespace', leftover_porespace
		# 	print 'np.sum(porespace[ind3+1:ind2+1])', np.sum(porespace_vol[ind3+1:ind2+1])
		# 	print 'dz_ind3', self.dz[ind3]
		# 	print 'porespace_ind3', porespace_vol[ind3]
		# 	raw_input('negative dz new')

		self.rho[ind3+1:ind2+1] = 870.0
		self.Tz[ind1:ind2+1] = T_MELT

		# split up last box into several
		divider = num_boxes_melted
		# ind3 should be removed and replaced with 2 new boxes.
		self.rho = np.concatenate((self.rho[ind1:ind3] , [new_node_1_rho,new_node_2_rho] , self.rho[ind3+1:-1] , self.rho[-1]*np.ones(num_boxes_melted-1)))
				
		# if (iii>7000 and iii<7435):
		# 	print('num_boxes_melted',num_boxes_melted)
		# 	print('ind1',ind1)
		# 	print('ind2',ind2)
		# 	print('ind3',ind3)
		# 	print(available_space)
		# 	print(melt_volume_IE)
		# 	print('top',len(self.rho[ind1:ind3]))
		# 	print('new',len([new_node_1_rho,new_node_2_rho]))
		# 	print('old',len(self.rho[ind3+1:-1]))
		# 	print('bottom',len(self.rho[-1]*np.ones(num_boxes_melted-1)))
		# 	print('pm_mass',pm_mass)
		# 	print('self.mass',self.mass[0:ind2])
		# 	print('melt_mass',melt_mass)

		# if (iii>7429 and iii<7435):
			# print('rholen',len(self.rho))
		self.age = np.concatenate((self.age[ind1:ind3] , [self.age[ind3],self.age[ind3]] , self.age[ind3+1:-1] , self.age[-1]*np.ones(num_boxes_melted-1)))
		dzhold = self.dz[ind1+1:ind3]
		dzhold2 = self.dz[ind3+1:-1]
		self.dz = np.concatenate((self.dz[ind1:ind3] , [new_node_1_dz,new_node_2_dz] , self.dz[ind3+1:-1] ,self.dz[-1]/divider*np.ones(num_boxes_melted-1)))
		self.dzn = np.concatenate((np.zeros(num_boxes_melted), np.append(dzhold, new_node_1_dz+new_node_2_dz), dzhold2))
		self.dzn = self.dzn[0:self.compboxes]

		
		# print(np.min(self.dz)
		# if (self.dz<=0).any():
		# 	print self.dz[0:ind2]
		# 	raw_input('negative dz')
		self.dz[0] = pm_dz
		self.Tz = np.concatenate((self.Tz[ind1:ind3] , [self.Tz[ind3],self.Tz[ind3]] , self.Tz[ind3+1:-1] , self.Tz[-1]*np.ones(num_boxes_melted-1)))
		self.bdot_mean = np.concatenate((self.bdot_mean[ind1:ind3] , [self.bdot_mean[ind3],self.bdot_mean[ind3]] , self.bdot_mean[ind3+1:-1] , self.bdot_mean[-1]*np.ones(num_boxes_melted-1)))
		self.z = self.dz.cumsum(axis = 0)
		self.z = np.concatenate(([0], self.z[:-1]))
		self.mass = self.rho*self.dz
		# if len(self.mass)!=self.gridLen:
		# 	print('line283')

		# print('self.mass', self.mass)

		# if there is an ice lens:
			# how much water can the porous firn refreeze? the rest is runoff.


	# print('###rho###', self.rho[0:6])
	# runoff_volume = runoff_volume_duetolimitedporespace + runoff_volume_duetolens
	# print 'runoff_volume = ', runoff_volume

	return self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn

def percolation2(self, iii):
	maxpore = 0.6

	porosity 	= 1 - self.rho/RHO_I #porosity (unitless)
	dz_old 		= self.dz

	porespace_vol 		= porosity * self.dz 		# pore space volume (meters) of each box - volume of air + water
	porespace_air		= porespace_vol - self.LWC 	# pore space that is filled with air (meters)
	maxLWC1				= porespace_vol * maxpore 	# maximum volume of water that can be stored in each node (meters)
	maxLWC2				= ((900 * self.dz) - self.mass)/RHO_W_KGM
	maxLWC 				= np.maximum(maxLWC1,maxLWC2)
	maxLWC_mass 		= maxLWC * RHO_W_KGM		# mass of the maximum volume of water
	
	# pot_mass = (self.mass+self.LWC*RHO_W_KGM)/self.dz # potential mass; checking if mass is too high
	# if np.any(pot_mass>917):
	# 	aa = np.where(pot_mass>917)[0][0]
	# 	print('problem at 359')
	# 	print('iii', iii)
	# 	print('rho',self.rho[aa])
	# 	print('dz',self.dz[aa])
	# 	print('pot_mass',pot_mass[aa])
	# 	print('LWC',self.LWC[aa])
	# 	print('maxmass',maxLWC_mass[aa])
	# 	print('maxLWC',maxLWC[aa])
	# 	print('porosity',porosity[aa])
	# 	print('porespace_air',porespace_air[aa])
	# 	input('366')

	rho_old = np.copy(self.rho)

	melt_volume_IE  	= self.snowmeltSec[iii] * S_PER_YEAR # meters
	melt_volume_WE		= melt_volume_IE * 0.917 # meters
	melt_mass			= melt_volume_WE * 1000. # kg

	# melt_mass_a 		= melt_mass # this will get iterated and subtracted from
	# melt_vol_a			= melt_mass_a / RHO_W_KGM

	heat_to_freeze 		= melt_mass * LF_I #amount of heat needed to refreeze the melt (J)

	ind1a 				= np.where(self.mass_sum<=melt_mass)[0] # indicies of boxes that will be melted away
	num_boxes_melted 	= len(ind1a)+1 #number of boxes that melt away, include the box that is partially melted
	ind1 				= np.where(self.mass_sum>melt_mass)[0][0] # index which will become the new surface
	# print('ind1',ind1)

	# pm is the partial melt (the box/volume that has a portion melted away)
	pm_mass 			= self.mass_sum[ind1] - melt_mass # the remaining mass of the PM box
	pm_dz 				= pm_mass / self.rho[ind1] #remaining thickness
	pm_porespace 		= (1 - self.rho[ind1]/RHO_I) * pm_dz #porespace in the PM box
	pm_rho 				= self.rho[ind1] #density of the PM box

	melt_boxes_LWC_vol  = np.sum(self.LWC[0:ind1]) #include the water mass from the boxes that melt (currently does not include from the partial melt box)
	melt_boxes_LWC_mass = melt_boxes_LWC_vol * RHO_W_KGM

	melt_mass_a = melt_mass + melt_boxes_LWC_mass
	melt_vol_a = melt_mass_a / RHO_W_KGM

	self.LWC[ind1] = 0

	cold_content_0 		= CP_I * self.mass * (T_MELT - self.Tz) 	# cold content of each box, i.e. how much heat to bring it to 273K (kJ)
	cold_content_0_sum 	= cold_content_0.cumsum(axis=0)
	cold_content 		= cold_content_0[ind1:] 					# just the boxes that don't melt away
	cold_content[0]  	= CP_I * pm_mass * (T_MELT - self.Tz[ind1]) # the partial melt box has its cold content reassigned.
	cold_content_sum 	= cold_content.cumsum(axis=0)

	for ii in range(len(self.dz)): # cycle through each node until all of the meltwater is refrozen or retained.

		ii = ind1 + ii #index on the real grid

		rhoperm = 800. # impermeable boundary density.

		if self.rho[ii]>rhoperm: # this layer is impermeable, so fill it (and the layers above it) up with water

			for kk in range(ii): # start filling element ii-1, ii-2, ii-3, etc. - bottom up.

				available_volume = maxLWC[ii-kk-1] - self.LWC[ii-kk-1] #how much porespace can be filled up.

				if ii-kk-1<0: #can't fill above the surface
					print('Saturated to the surface')
					runoff = melt_vol_a
					melt_vol_a = 0
					melt_mass_a = 0
					mark = 'one'
					break

				if available_volume >= melt_vol_a: # if there is more pore space than melt, the melt just fills that node
					self.LWC[ii-kk-1] 	= self.LWC[ii-kk-1] + melt_vol_a
					if self.LWC[ii-kk-1] > maxLWC[ii-kk-1]:
						print('problem at 418')
					self.Tz[ii-kk-1]	= T_MELT
					melt_vol_a = 0
					melt_mass_a = 0
					mark = 'two'
					break

				elif available_volume<0: # should never happen!
					print('ii',ii)
					print('kk',kk)
					print(self.rho[0:6])
					print(self.LWC[0:6])
					print(porespace_vol[0:6])
					print('negative available_volume! (406)')
					print('pot_mass',pot_mass[ii-kk-1])
					# input()
					available_volume = 0
					mark = 'three'
					continue

				elif available_volume==0: # no space available in this node.
					# self.LWC[ii-kk-1] 	= self.LWC[ii-kk-1]
					# self.LWC[ii-kk-1] 	= maxLWC[ii-kk-1]
					# print('hi')

					continue

				else:					
					self.LWC[ii-kk-1] 	= self.LWC[ii-kk-1] + available_volume
					if self.LWC[ii-kk-1] > maxLWC[ii-kk-1]:
						self.LWC[ii-kk-1] 	= maxLWC[ii-kk-1]

						# print('439')
						# print(self.LWC[ii-kk-1] - maxLWC[ii-kk-1]) 
						# # print(maxLWC[ii-kk-1])
						# input('i439')					
					self.Tz[ii-kk-1]	= T_MELT
					melt_vol_a = melt_vol_a - available_volume

					if self.LWC[ii-kk-1]<0:
						print('negative LWC at 461', iii)

					mark = 'three'

					# rr = (self.mass[ii-kk-1] + self.LWC[ii-kk-1]*RHO_W_KGM) / self.dz[ii-kk-1]
					# if rr>917.0:
					# 	print('too much water')
					# 	input()

			porespace_air		= porespace_vol - self.LWC 	# pore space that is filled with air (meters)
			maxLWC				= porespace_vol * maxpore 	# maximum volume of water that can be stored in each node (meters)
			maxLWC_mass 		= maxLWC * RHO_W_KGM		# mass of the maximum volume of water

			break # break the ii loop after water has saturated the porosity

		#### refreezing:
		if self.LWC[ii]>0: # no refreezing (of new meltwater) if there is any liquid water present
			refreeze_mass = 0
			refreeze_volume = 0

		else: # there is some refreezing

			max_refreeze_volume = maxLWC[ii] - self.LWC[ii] # how much water volume can be accomodated in the pore space
			max_refreeze_mass = max_refreeze_volume * RHO_W_KGM  # how much mass can be accomodated in the pore space

			refreeze_potential_mass = cold_content[ii] / LF_I # how much mass can refreeze in this node due to cold content(kg)
			refreeze_potential_volume = refreeze_potential_mass / RHO_W_KGM # water volume that can refreeze

			# if refreeze_potential_mass > maxLWC_mass[ii]:
				# refreeze_potential_mass = maxLWC_mass[ii]

			if ((melt_mass_a < refreeze_potential_mass) and (refreeze_potential_mass < max_refreeze_mass)): # if the volume of liquid is less than the cold content is capable of freezing and the liquid volume is less than the porespace can accomodate, all of the water refreezes
				
				refreeze_mass 		= melt_mass_a # how much mass refreezes
				cc_n 				= refreeze_mass * LF_I # heat that is used for refreezing, will raise temperature
				self.Tz[ii] 		= self.Tz[ii] + cc_n / (CP_I * self.mass[ii]) # new temperature
				self.mass[ii]		= self.mass[ii] + refreeze_mass				# new mass				
				self.rho[ii] 		=  self.mass[ii] / self.dz[ii] 		# new density after refreezing
				if self.rho[ii]>917:
					print('refreeze_mass 470',refreeze_mass)
					print('rho!',self.rho[ii])
					print('rho_old',rho_old[ii])
				self.LWC[ii]		= 0									# all is frozen, so no LWC
				porosity[ii] 		= 1 - self.rho[ii] / RHO_I 			# new porosity
				porespace_vol[ii] 	= porosity[ii] * self.dz[ii] 		# new porespace, m
				porespace_air[ii]	= porespace_vol[ii] - self.LWC[ii]	# new air porespace
				maxLWC[ii]			= porespace_vol[ii] * maxpore 		# maximum volume that can be stored in each node
				melt_mass_a			= 0 
				mark = 'four'

				break #quit the ii loop (no more melt mass available)

			elif ((melt_mass_a < refreeze_potential_mass) and (refreeze_potential_mass >= max_refreeze_mass)): # enough cold content to freeze, but not enough space to accomodate that water
				refreeze_mass 		= max_refreeze_mass  # how much mass refreezes
				cc_n 				= refreeze_mass * LF_I # heat that is used for refreezing, will raise temperature
				self.Tz[ii] 		= self.Tz[ii] + cc_n / (CP_I * self.mass[ii]) # new temperature
				self.mass[ii]		= self.mass[ii] + refreeze_mass				# new mass
				self.rho[ii] 		=  self.mass[ii] / self.dz[ii] 		# new density after refreezing
				# if self.rho[ii]>917:
				print('refreeze_mass 495',refreeze_mass)
				print('rho!',self.rho[ii])
				self.LWC[ii]		= 0									# all is frozen, so no LWC
				porosity[ii] 		= 1 - self.rho[ii] / RHO_I 			# new porosity
				porespace_vol[ii] 	= porosity[ii] * self.dz[ii] 		# new porespace, m
				porespace_air[ii]	= porespace_vol[ii] - self.LWC[ii]	# new air porespace
				maxLWC[ii]			= porespace_vol[ii] * maxpore 			# maximum volume that can be stored in each node
				melt_mass_a 		= melt_mass_a - refreeze_mass
				mark = 'five'

				continue # water should move to the ii-1 cell and sit there.

			elif melt_mass_a >= refreeze_potential_mass: # there is more melt than can be refrozen, so refreeze whatever is possible
 
				refreeze_mass = refreeze_potential_mass # how much the cold content can refreeze

				if refreeze_mass > max_refreeze_mass: # if the cold content refreeze is more than the pore space allows, use the amout pore space allows
					# extra = refreeze_mass - max_refreeze_mass
					refreeze_mass = max_refreeze_mass
					# melt_mass_a = melt_mass_a + extra
					
				melt_mass_a 		= melt_mass_a - refreeze_mass
				cc_n 				= refreeze_mass * LF_I # heat that is used for refreezing, will raise temperature
				self.Tz[ii] 		= self.Tz[ii] + cc_n / (CP_I * self.mass[ii]) # new temperature
				if self.Tz[ii] > T_MELT:
					self.Tz[ii] = T_MELT
				m0 = np.copy(self.mass)
				self.mass[ii]		= self.mass[ii] + refreeze_mass				# new mass
				# rho_old = self.rho 
				self.rho[ii] 		= self.mass[ii] / self.dz[ii] 		# new density after refreezing
				if self.rho[ii]>917:
					print('refreeze_mass 540',refreeze_mass)
					print('ii',ii)
					print('ind1',ind1)
					print('mass',self.mass[ii])
					print('dz', self.dz[ii])
					print('rho!',self.rho[ii])
					print('rho_old',rho_old[ii])				
				self.LWC[ii]		= 0
				porosity[ii] 		= 1 - self.rho[ii] / RHO_I 					# new porosity
				porespace_vol[ii] 	= porosity[ii] * self.dz[ii] 				# new porespace, m
				porespace_air[ii]	= porespace_vol[ii] - self.LWC[ii]	# new air porespace
				maxLWC[ii]			= porespace_vol[ii] * maxpore 			# maximum volume that can be stored in each node
				mark = 'six'				

			if melt_mass_a <= 0: # stop the loop if there is no more melt mass
				break
		######## end refreezing

		###############################
		#### now deal with liquid water
		if self.rho[ii] > rhoperm:
			Wmi = 0
			Swi = 0
		else:
			Wmi 		= 0.057*(RHO_I - self.rho[ii])/self.rho[ii] + 0.017 # water per snow-plus- water mass irreducible liquid water content, Langen eqn 3 unitless)
			Swi			= Wmi / (1 - Wmi) * (self.rho[ii]*RHO_I)/(1000*(RHO_I - self.rho[ii])) 	#irreducible water saturation, volume of water per porespace volume (unitless)

		# print('Swi', Swi)
		# irreducible_mass = 0.02 * porespace_vol[ii] * RHO_W_KGM 
		# irreducible_mass 	= Swi * porespace_vol[ii] * RHO_W_KGM 		# mass of irreducible water for each volume (potential - does not separate how much is already there)
		irreducible_mass 	= Swi * porespace_air[ii] * RHO_W_KGM

	

	# if np.any(maxlwc>self.LWC):
	# 	input('here')
		maxLWC_mass 		= maxLWC * RHO_W_KGM					# maximum mass that will fit  with porespace available
		irreducible_volume 	= irreducible_mass / RHO_W_KGM
		maxLWC				= maxLWC_mass / RHO_W_KGM

		# retained_water_mass	= maxLWC_mass - self.LWC[ii]*RHO_W_KGM 	# how much mass of water will be retained from the melt volume (kg)
		retained_water_mass = irreducible_mass - self.LWC[ii]*RHO_W_KGM

		if retained_water_mass < 0:
			retained_water_mass = 0
		
		if (retained_water_mass + self.LWC[ii]*RHO_W_KGM) > maxLWC_mass[ii]: # if the sum of the retained meltwater and current LWC is greater than the maximum possible, make ml_new the maximum value.
			# retained_water_mass = maxLWC_mass - self.LWC[ii]*RHO_W_KGM
			retained_water_mass = 0


		retained_water_vol 	= retained_water_mass / RHO_W_KGM 	# volume retained (m, because 1D)

		self.LWC[ii] 		= self.LWC[ii] + retained_water_vol # new LWC volume in this node (meters)
		if self.LWC[ii]<0:
			print('negative LWC at 596', iii)

		if self.LWC[ii]> maxLWC[ii]:
			add_melt_mass = (self.LWC[ii]-maxLWC[ii]) * RHO_W_KGM
			self.LWC[ii] = maxLWC[ii]
			melt_mass_a = melt_mass_a + add_melt_mass
			mark = 'seven'

		if self.LWC[ii]<0:
			print('negative LWC at 604', iii)
			# print('568')
			# input()
		self.Tz[ii]			= T_MELT

		if (self.LWC[ii]*RHO_W_KGM) > irreducible_mass: # if the LWC exceeds the mimimum, add some liquid to the percolating mass

			excess = self.LWC[ii]*RHO_W_KGM - irreducible_mass
			self.LWC[ii] = irreducible_mass / RHO_W_KGM
			if self.LWC[ii]<0:
				self.LWC[ii]=0
				excess = 0
			mark = 'eight'
				# print('negative LWC at 614', iii)

			melt_mass_a = melt_mass_a + excess

		porespace_air[ii]	= porespace_vol[ii] - self.LWC[ii]	# new air porespace
		maxLWC[ii]			= porespace_vol[ii] * maxpore 			# maximum volume that can be stored in each node

		melt_mass_a			= melt_mass_a - retained_water_mass



		if melt_mass_a<=0:
			break

	

	# split up last box into several
	divider = num_boxes_melted
	print(self.modeltime[iii])
	self.rho = np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted)))
	self.LWC = np.concatenate((self.LWC[ind1:-1] , self.LWC[-1]*np.ones(num_boxes_melted)))
	self.age = np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_melted)))
	self.dz  = np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted)))
	self.dz[0] = pm_dz
	self.dzn = np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
	self.dzn = self.dzn[0:self.compboxes]

	self.Tz  = np.concatenate((self.Tz[ind1:-1],self.Tz[-1]*np.ones(num_boxes_melted)))
	self.bdot_mean = np.concatenate((self.bdot_mean[ind1:-1],self.bdot_mean[-1]*np.ones(num_boxes_melted)))
	self.z = self.dz.cumsum(axis = 0)
	self.z = np.concatenate(([0], self.z[:-1]))
	self.mass = self.rho*self.dz
	self.Tz[self.LWC>0] = T_MELT

	pot_mass = (self.mass+self.LWC*RHO_W_KGM)/self.dz #checking if mass is too high
	if np.any(pot_mass>917):
		aa = np.where(pot_mass>917)[0][0]

		print('problem at 615')
		print('iii', iii)
		print('rho',self.rho[aa])
		print('pot_mass',pot_mass[aa])
		print('LWC',self.LWC[aa])
		# print()
		print()

		input()

	maxlwc_mass = (917*self.dz - self.mass)

	if np.any(maxlwc_mass<(self.LWC*RHO_W_KGM)):
		print(mark)
		input('here')

	# self.LWC[self.LWC<0]=0.0

	return self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn, self.LWC

def percolation3(self, iii):
	maxpore = 0.4

	melt_volume_IE  	= self.snowmeltSec[iii] * S_PER_YEAR # meters
	melt_volume_WE		= melt_volume_IE * 0.917 # meters
	melt_mass			= melt_volume_WE * 1000. # kg

	heat_to_freeze 		= melt_mass * LF_I #amount of heat needed to refreeze the melt (J)

	ind1a 				= np.where(self.mass_sum<=melt_mass)[0] # indicies of boxes that will be melted away
	num_boxes_melted 	= len(ind1a)+1 #number of boxes that melt away, include the box that is partially melted
	ind1 				= np.where(self.mass_sum>melt_mass)[0][0] # index which will become the new surface

	# pm is the partial melt (the box/volume that has a portion melted away)
	pm_mass 			= self.mass_sum[ind1] - melt_mass # the remaining mass of the PM box
	pm_dz 				= pm_mass / self.rho[ind1] #remaining thickness
	pm_porespace 		= (1 - self.rho[ind1]/RHO_I) * pm_dz #porespace in the PM box
	pm_rho 				= self.rho[ind1] #density of the PM box
	pm_lwc				= self.LWC[ind1]/self.dz[ind1] * pm_dz

	melt_boxes_LWC_vol  = np.sum(self.LWC[0:ind1]) - pm_lwc #include the water mass from the boxes that melt (currently does not include from the partial melt box)
	melt_boxes_LWC_mass = melt_boxes_LWC_vol * RHO_W_KGM

	melt_mass_a = melt_mass + melt_boxes_LWC_mass
	melt_vol_a = melt_mass_a / RHO_W_KGM

	#regrid now (after melt)
	divider = num_boxes_melted
	self.rho = np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted)))
	self.LWC = np.concatenate((self.LWC[ind1:-1] , self.LWC[-1]*np.ones(num_boxes_melted)))
	self.LWC[0] = pm_lwc
	self.age = np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_melted)))
	self.dz  = np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted)))
	self.dz[0] = pm_dz
	self.dzn = np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
	self.dzn = self.dzn[0:self.compboxes]
	self.Tz  = np.concatenate((self.Tz[ind1:-1],self.Tz[-1]*np.ones(num_boxes_melted)))
	self.bdot_mean = np.concatenate((self.bdot_mean[ind1:-1],self.bdot_mean[-1]*np.ones(num_boxes_melted)))
	self.z = self.dz.cumsum(axis = 0)
	self.z = np.concatenate(([0], self.z[:-1]))
	self.mass = self.rho*self.dz

	# now working all with the new grid
	porosity 	= 1 - self.rho/RHO_I #porosity (unitless)

	porespace_vol 		= porosity * self.dz 		# pore space volume (meters) of each box - volume of air + water
	porespace_air		= porespace_vol - self.LWC 	# pore space that is filled with air (meters)
	maxLWC1				= porespace_vol * maxpore 	# maximum volume of water that can be stored in each node (meters)
	maxLWC2				= ((900 * self.dz) - self.mass)/RHO_W_KGM
	maxLWC 				= np.minimum(maxLWC1,maxLWC2)
	maxLWC_mass 		= maxLWC * RHO_W_KGM		# mass of the maximum volume of water

	cold_content		= CP_I * self.mass * (T_MELT - self.Tz) 	# cold content of each box, i.e. how much heat to bring it to 273K (kJ)
	cold_content_sum 	= cold_content.cumsum(axis=0)

	rho_old = np.copy(self.rho)

	for ii in range(len(self.dz)): # cycle through each node until all of the meltwater is refrozen or retained.

		rhoperm = 600. # impermeable boundary density.

		if self.rho[ii]>rhoperm: # this layer is impermeable, so fill it (and the layers above it) up with water

			for kk in range(ii): # start filling element ii-1, ii-2, ii-3, etc. - bottom up.

				available_volume = maxLWC[ii-kk-1] - self.LWC[ii-kk-1] #how much porespace can be filled up.

				if ii-kk-1<0: #can't fill above the surface
					print('Saturated to the surface')
					runoff = melt_vol_a
					melt_vol_a = 0
					melt_mass_a = 0
					break

				if available_volume >= melt_vol_a: # if there is more pore space than melt, the melt will fill up as much pore space as it can.
					self.LWC[ii-kk-1] 	= self.LWC[ii-kk-1] + melt_vol_a
					self.Tz[ii-kk-1]	= T_MELT
					melt_vol_a = 0
					melt_mass_a = 0
					break

				elif available_volume<0: # should never happen!
					available_volume = 0

					continue

				elif available_volume==0: 
					continue

				else:					
					self.LWC[ii-kk-1] 	= self.LWC[ii-kk-1] + available_volume # accomodate whatever is possible.
					if self.LWC[ii-kk-1] > maxLWC[ii-kk-1]:
						self.LWC[ii-kk-1] 	= maxLWC[ii-kk-1]
				
					self.Tz[ii-kk-1]	= T_MELT
					melt_vol_a = melt_vol_a - available_volume


			# porespace_air		= porespace_vol - self.LWC 	# pore space that is filled with air (meters)
			# maxLWC				= porespace_vol * maxpore 	# maximum volume of water that can be stored in each node (meters)
			# maxLWC_mass 		= maxLWC * RHO_W_KGM		# mass of the maximum volume of water

			break # break the ii loop after water has saturated the porosity

		#### refreezing:
		if self.LWC[ii]>0: # no refreezing (of new meltwater) if there is any liquid water present
			refreeze_mass = 0
			refreeze_volume = 0

		else: # there is some refreezing

			max_refreeze_volume = maxLWC[ii] - self.LWC[ii] # how much water volume can be accomodated in the pore space
			max_refreeze_mass = max_refreeze_volume * RHO_W_KGM  # how much mass can be accomodated in the pore space

			refreeze_potential_mass = cold_content[ii] / LF_I # how much mass can refreeze in this node due to cold content(kg)
			refreeze_potential_volume = refreeze_potential_mass / RHO_W_KGM # water volume that can refreeze


			if ((melt_mass_a < refreeze_potential_mass) and (refreeze_potential_mass < max_refreeze_mass)): # if the volume of liquid is less than the cold content is capable of freezing and the liquid volume is less than the porespace can accomodate, all of the water refreezes
				
				refreeze_mass 		= melt_mass_a # how much mass refreezes
				cc_n 				= refreeze_mass * LF_I # heat that is used for refreezing, will raise temperature
				self.Tz[ii] 		= self.Tz[ii] + cc_n / (CP_I * self.mass[ii]) # new temperature
				self.mass[ii]		= self.mass[ii] + refreeze_mass				# new mass				
				self.rho[ii] 		=  self.mass[ii] / self.dz[ii] 		# new density after refreezing
				self.LWC[ii]		= 0									# all is frozen, so no LWC

				porosity[ii] 		= 1 - self.rho[ii] / RHO_I 			# new porosity
				porespace_vol[ii] 	= porosity[ii] * self.dz[ii] 		# new porespace, m
				porespace_air[ii]	= porespace_vol[ii] - self.LWC[ii]	# new air porespace
				maxLWC[ii]			= porespace_vol[ii] * maxpore 		# maximum volume that can be stored in each node
				melt_mass_a			= 0 

				break #quit the ii loop (no more melt mass available)

			elif ((melt_mass_a < refreeze_potential_mass) and (refreeze_potential_mass >= max_refreeze_mass)): # enough cold content to freeze, but not enough space to accomodate that water
				refreeze_mass 		= max_refreeze_mass  # how much mass refreezes
				cc_n 				= refreeze_mass * LF_I # heat that is used for refreezing, will raise temperature
				self.Tz[ii] 		= self.Tz[ii] + cc_n / (CP_I * self.mass[ii]) # new temperature
				self.mass[ii]		= self.mass[ii] + refreeze_mass				# new mass
				self.rho[ii] 		=  self.mass[ii] / self.dz[ii] 		# new density after refreezing
				self.LWC[ii]		= 0									# all is frozen, so no LWC
				porosity[ii] 		= 1 - self.rho[ii] / RHO_I 			# new porosity
				porespace_vol[ii] 	= porosity[ii] * self.dz[ii] 		# new porespace, m
				porespace_air[ii]	= porespace_vol[ii] - self.LWC[ii]	# new air porespace
				maxLWC[ii]			= porespace_vol[ii] * maxpore 			# maximum volume that can be stored in each node
				melt_mass_a 		= melt_mass_a - refreeze_mass

				continue # water should move to the ii-1 cell and sit there.

			elif melt_mass_a >= refreeze_potential_mass: # there is more melt than can be refrozen, so refreeze whatever is possible
 
				refreeze_mass = refreeze_potential_mass # how much the cold content can refreeze

				if refreeze_mass > max_refreeze_mass: # if the cold content refreeze is more than the pore space allows, use the amout pore space allows					
					refreeze_mass = max_refreeze_mass
					
					
				melt_mass_a 		= melt_mass_a - refreeze_mass
				cc_n 				= refreeze_mass * LF_I # heat that is used for refreezing, will raise temperature
				self.Tz[ii] 		= self.Tz[ii] + cc_n / (CP_I * self.mass[ii]) # new temperature
				if self.Tz[ii] > T_MELT:
					self.Tz[ii] = T_MELT
				self.mass[ii]		= self.mass[ii] + refreeze_mass				# new mass
				self.rho[ii] 		= self.mass[ii] / self.dz[ii] 		# new density after refreezing				
				self.LWC[ii]		= 0
				porosity[ii] 		= 1 - self.rho[ii] / RHO_I 					# new porosity
				porespace_vol[ii] 	= porosity[ii] * self.dz[ii] 				# new porespace, m
				porespace_air[ii]	= porespace_vol[ii] - self.LWC[ii]	# new air porespace
				maxLWC[ii]			= porespace_vol[ii] * maxpore 			# maximum volume that can be stored in each node				

			if melt_mass_a <= 0: # stop the loop if there is no more melt mass
				break
		######## end refreezing

		###############################
		#### now deal with liquid water
		if self.rho[ii] > rhoperm:
			Wmi = 0
			Swi = 0
		else:
			Wmi 		= 0.057*(RHO_I - self.rho[ii])/self.rho[ii] + 0.017 # water per snow-plus- water mass irreducible liquid water content, Langen eqn 3 unitless)
			Swi			= Wmi / (1 - Wmi) * (self.rho[ii]*RHO_I)/(1000*(RHO_I - self.rho[ii])) 	#irreducible water saturation, volume of water per porespace volume (unitless)

		irreducible_mass = 0.02 * porespace_vol[ii] * RHO_W_KGM  		
		# irreducible_mass 	= Swi * porespace_air[ii] * RHO_W_KGM # mass of irreducible water for each volume (potential - does not separate how much is already there)

		maxLWC_mass 		= maxLWC * RHO_W_KGM					# maximum mass that will fit  with porespace available
		irreducible_volume 	= irreducible_mass / RHO_W_KGM
		maxLWC				= maxLWC_mass / RHO_W_KGM

		# retained_water_mass	= maxLWC_mass - self.LWC[ii]*RHO_W_KGM 	# how much mass of water will be retained from the melt volume (kg)
		retained_water_mass = irreducible_mass - self.LWC[ii]*RHO_W_KGM

		if retained_water_mass < 0:
			retained_water_mass = 0
		
		if (retained_water_mass + self.LWC[ii]*RHO_W_KGM) > maxLWC_mass[ii]: # if the sum of the retained meltwater and current LWC is greater than the maximum possible, make ml_new the maximum value.
			# retained_water_mass = maxLWC_mass - self.LWC[ii]*RHO_W_KGM
			retained_water_mass = 0

		retained_water_vol 	= retained_water_mass / RHO_W_KGM 	# volume retained (m, because 1D)

		self.LWC[ii] 		= self.LWC[ii] + retained_water_vol # new LWC volume in this node (meters)

		if self.LWC[ii] > maxLWC[ii]:
			add_melt_mass = (self.LWC[ii]-maxLWC[ii]) * RHO_W_KGM
			self.LWC[ii] = maxLWC[ii]
			melt_mass_a = melt_mass_a + add_melt_mass

		self.Tz[ii]			= T_MELT

		if (self.LWC[ii]*RHO_W_KGM) > irreducible_mass: # if the LWC exceeds the mimimum, add some liquid to the percolating mass

			excess = self.LWC[ii]*RHO_W_KGM - irreducible_mass
			self.LWC[ii] = irreducible_mass / RHO_W_KGM
			if self.LWC[ii]<0:
				self.LWC[ii]=0
				excess = 0
				# print('negative LWC at 614', iii)

			melt_mass_a = melt_mass_a + excess

		porespace_air[ii]	= porespace_vol[ii] - self.LWC[ii]	# new air porespace
		maxLWC[ii]			= porespace_vol[ii] * maxpore 			# maximum volume that can be stored in each node
		melt_mass_a			= melt_mass_a - retained_water_mass

		if melt_mass_a<=0:
			break



	return self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn, self.LWC

def percolation4(self, iii):
	# maxpore 				= 0.2 #upper limit on what percentage of the porosity can be filled with water.
	impermeable_rho			= 725. 

	if np.any(self.LWC<0):
		print('negative lwc')

	melt_volume_IE  		= self.snowmeltSec[iii] * S_PER_YEAR # meters
	melt_volume_WE			= melt_volume_IE * 0.917 	# meters
	melt_mass				= melt_volume_WE * 1000. 	# kg
	heat_to_freeze 			= melt_mass * LF_I 			#amount of heat needed to refreeze the melt (J)

	ind1a 					= np.where(self.mass_sum <= melt_mass)[0] # indicies of boxes that will be melted away
	num_boxes_melted 		= len(ind1a)+1 							# number of boxes that melt away, include the box that is partially melted
	ind1 					= np.where(self.mass_sum > melt_mass)[0][0] # index which will become the new surface

	# pm is the partial melt (the box/volume that has a portion melted away)
	pm_mass 				= self.mass_sum[ind1] - melt_mass # the remaining mass of the PM box
	pm_dz 					= pm_mass / self.rho[ind1] #remaining thickness
	pm_porespace 			= (1 - self.rho[ind1]/RHO_I) * pm_dz #porespace in the PM box
	pm_rho 					= self.rho[ind1] #density of the PM box
	pm_lwc					= self.LWC[ind1]/self.dz[ind1] * pm_dz

	melt_boxes_LWC_vol  	= np.sum(self.LWC[0:ind1+1]) - pm_lwc #include the water mass from the boxes that melt (currently does not include from the partial melt box)
	melt_boxes_LWC_mass 	= melt_boxes_LWC_vol * RHO_W_KGM
	melt_mass_a 			= melt_mass + melt_boxes_LWC_mass
	if melt_mass_a<0:
		print('########')
		print('melt_mass',melt_mass)
		print('melt_boxes_LWC_mass',melt_boxes_LWC_mass)
		print('lwcsum',np.sum(self.LWC[0:ind1]))
		print('pm_lwc',pm_lwc)
		print('########')


	melt_vol_a 				= melt_mass_a / RHO_W_KGM

	### regrid now (after melt)
	divider 				= num_boxes_melted
	self.rho 				= np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted)))
	self.LWC 				= np.concatenate((self.LWC[ind1:-1] , self.LWC[-1]*np.ones(num_boxes_melted)))
	self.LWC[0] 			= pm_lwc
	self.age 				= np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_melted)))
	# self.dz  				= np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted)))
	self.dz  				= np.concatenate((self.dz[ind1:-1] , self.dz[-1]*np.ones(num_boxes_melted)))
	self.dz[0] 				= pm_dz
	self.Dcon 				= np.concatenate((self.Dcon[ind1:-1] , self.Dcon[-1]*np.ones(num_boxes_melted)))
	self.dzn 				= np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
	self.dzn 				= self.dzn[0:self.compboxes]
	self.Tz  				= np.concatenate((self.Tz[ind1:-1] , self.Tz[-1]*np.ones(num_boxes_melted)))
	self.bdot_mean 			= np.concatenate((self.bdot_mean[ind1:-1] , self.bdot_mean[-1]*np.ones(num_boxes_melted)))
	self.z 					= self.dz.cumsum(axis = 0)
	self.z 					= np.concatenate(([0] , self.z[:-1]))
	uni, uni_i, uni_n, uni_c = np.unique(self.z,return_index=True,return_inverse=True,return_counts=True)
	if len(uni)!=len(self.z):
		print('uni')
		print('end',self.z[-10:])
		print('type',self.z.dtype)
		print('dz',self.dz[-10:])
		input('enter to continue')
	self.mass 				= self.rho * self.dz

	##########################################
	### now working all with the new grid ####
	##########################################
	porosity 				= 1 - self.rho / RHO_I #porosity (unitless)
	porespace_vol 			= porosity * self.dz 		# pore space volume (meters) of each box - volume of air + water
	porespace_air			= porespace_vol - self.LWC 	# pore space that is filled with air (meters)

	cold_content			= CP_I * self.mass * (T_MELT - self.Tz) 	# cold content of each box, i.e. how much heat to bring it to 273K (kJ)
	cold_content_sum 		= cold_content.cumsum(axis=0)
	refreeze_mass_pot 		= cold_content / LF_I # how much mass of the meltwater could be refrozen due to cold content
	refreeze_mass_pot_sum 	= refreeze_mass_pot.cumsum(axis=0) 

	# calculate what the values will be after refreeze happens (does make a small difference)
	rho_pot					= (self.mass + refreeze_mass_pot) / self.dz # what the mass of the boxes would be if the refreezemass refroze
	porosity_pot			= 1 - rho_pot / RHO_I
	porespace_vol_pot		= porosity_pot * self.dz
	porespace_air_pot		= porespace_vol_pot - self.LWC

	Wmi 					= 0.057 * (RHO_I - rho_pot) / rho_pot + 0.017 # water per snow-plus- water mass irreducible liquid water content, Langen eqn 3 unitless)
	Swi						= Wmi / (1 - Wmi) * (rho_pot * RHO_I) / (1000 * (RHO_I - rho_pot)) 	#irreducible water saturation, volume of water per porespace volume (unitless), Colbeck 1972
	# Swi[Swi>impermeable_rho] = 0
	maxpore = Swi * 2.0

	maxLWC1					= porespace_vol * maxpore 	# maximum volume of water that can be stored in each node (meters)
	maxLWC2					= ((917.0 * self.dz) - self.mass) / RHO_W_KGM # double check that the LWC does not get too large. 
	maxLWC 					= np.minimum(maxLWC1 , maxLWC2)
	maxLWC[self.rho>impermeable_rho] = 0
	maxLWC_mass 			= maxLWC * RHO_W_KGM		# mass of the maximum volume of water
	maxLWC1_pot				= porespace_vol_pot * maxpore 	# maximum volume of water that can be stored in each node (meters)
	maxLWC2_pot				= ((917.0 * self.dz) - (self.mass + refreeze_mass_pot)) / RHO_W_KGM # double check that the LWC does not get too large. 
	maxLWC_pot 				= np.minimum(maxLWC1_pot , maxLWC2_pot)
	# maxLWC_pot[rho_pot>impermeable_rho] = 0
	maxLWC_mass_pot 		= maxLWC_pot * RHO_W_KGM		# mass of the maximum volume of water



	irreducible_mass_pot 	= Swi * porespace_vol_pot * RHO_W_KGM # mass of irreducible water for each volume (potential - does not separate how much is already there)
	irreducible_vol_pot		= irreducible_mass_pot / RHO_W_KGM
	liquid_storage_vol_pot	= irreducible_vol_pot - self.LWC
	liquid_storage_mass_pot = liquid_storage_vol_pot * RHO_W_KGM

	extra_liquid_mass		= np.sum(self.LWC[self.LWC > irreducible_vol_pot] * RHO_W_KGM - irreducible_mass_pot[self.LWC > irreducible_vol_pot])
	storage_mass_pot		= liquid_storage_mass_pot + refreeze_mass_pot #how much can be refrozen plus how much will stick around due to capillary
	storage_mass_pot_sum	= storage_mass_pot.cumsum(axis=0)
	total_liquid_mass 		= melt_mass_a + extra_liquid_mass
	
	trac = False
	if total_liquid_mass<0:
		print('negative total_liquid_mass 1')
		print('iii',iii)
		print('melt_mass_a',melt_mass_a)
		print('extra_liquid_mass',extra_liquid_mass)
		trac = True


	try:
		ind_p 					= np.where(storage_mass_pot_sum >= total_liquid_mass)[0][0] #the layer that water will percolate to
	except: # all of the liquid is runoff.
		ind_p=0
	###################################

	### if there is an impermeable layer, block water from getting through
	if np.any(self.rho[0:ind_p+1] >= impermeable_rho):

		if trac:
			print('trac')
			print('iii',iii)

		ind_p 					= np.where(self.rho >= impermeable_rho)[0][0] #- 1 # the index of the node that has density greater than the impermeable density
		# print('blocking lens')
		# print('depth',self.z[ind_p])
		id1 					= np.where(self.LWC >  irreducible_vol_pot)[0] # indices where the liquid water content is greater than the irreducible
		# id2 					= np.where(id1 < ind_p)[0]
		id2 = id1[id1<ind_p] 

		extra_liquid_mass		= np.sum(self.LWC[id2] * RHO_W_KGM) - np.sum(irreducible_mass_pot[id2])
		storage_mass_pot		= liquid_storage_mass_pot[0:ind_p] + refreeze_mass_pot[0:ind_p] #how much can be refrozen plus how much will stick around due to capillary
		storage_mass_pot_sum	= storage_mass_pot.cumsum(axis=0)
		total_liquid_mass 		= melt_mass_a + extra_liquid_mass
		if total_liquid_mass<0:
			print('!!!!!!!!!')
			print('negative total_liquid_mass 2')
			print('melt_mass_a',melt_mass_a)
			print('extra_liquid_mass',extra_liquid_mass)
			print('ind_p',ind_p)
			print('id1',id1)
			print('id2',id2)
			print('sumlwc', np.sum(self.LWC[id2] * RHO_W_KGM))
			print('sumirr', np.sum(irreducible_mass_pot[id2]))
			print('!!!!!!!!!')
			input('press enter to continue')

		#first, refreeze where possible
		self.mass[0:ind_p] 		= self.mass[0:ind_p] + refreeze_mass_pot[0:ind_p]
		self.rho[0:ind_p] 		= self.mass[0:ind_p] / self.dz[0:ind_p]
		self.Tz[0:ind_p] 		= T_MELT

		mass_frozen 			= np.sum(refreeze_mass_pot[0:ind_p])
		if mass_frozen >= total_liquid_mass:
			total_liquid_mass 	= 0
		else:
			total_liquid_mass 	= total_liquid_mass - mass_frozen

		# then, fill up the nodes above the ice slab
		maxLWC_mass_pot_f 		= np.flipud(maxLWC_mass_pot[0:ind_p])
		maxLWC_mass_pot_f_sum 	= maxLWC_mass_pot_f.cumsum(axis=0)

		if total_liquid_mass >= np.sum(maxLWC_mass_pot_f): # all porespace gets filled and there is runoff		
			self.LWC[0:ind_p] 		= maxLWC_pot[0:ind_p] # each node gets the maximum allowed
			# stored_water_vol 		= np.sum(self.LWC[0:ind_p]) #can calculate how much runoff there is, need to consider how much LWC there was previously
		
		else: # fill up however much porespace is needed to accomodate the meltwater

			# print('mass pot sum',np.sum(maxLWC_mass_pot_f))
			# print('liquid mass',total_liquid_mass)
			# try:
			# 	print('ml_end',maxLWC_mass_pot_f_sum[-1])
			# except:
			# 	print('ml_all',maxLWC_mass_pot_f_sum)	

			ind_f 					= np.where(maxLWC_mass_pot_f_sum > total_liquid_mass)[0][0] #index on the flipped grid
			ind_g 					= ind_p - 1 - ind_f #index on the real grid.
			self.LWC[ind_g+1:ind_p] = maxLWC_mass_pot[ind_g + 1:ind_p] / RHO_W_KGM # fill the indices up with the maximum allowed water
			lv_mass 				= total_liquid_mass - np.sum(maxLWC_mass_pot[ind_g + 1:ind_p]) 	# leftover volume
			self.LWC[ind_g] 		= lv_mass / RHO_W_KGM 						# put that into the ind_g node
		###################################

	### there is not an impermeable layer, water goes to layer ind_p
	elif ind_p>0: 

		# first, up to ind_p (not inclusive)
		self.mass[0:ind_p] 		= self.mass[0:ind_p] + refreeze_mass_pot[0:ind_p]
		self.rho[0:ind_p] 		= self.mass[0:ind_p] / self.dz[0:ind_p]
		lwc_old 				= np.copy(self.LWC)
		self.LWC[0:ind_p] 		= irreducible_mass_pot[0:ind_p] / RHO_W_KGM
		self.Tz[0:ind_p] 		= T_MELT
		lw_mass_retained 		= np.sum(refreeze_mass_pot[0:ind_p]) + np.sum(irreducible_mass_pot[0:ind_p]) - np.sum(lwc_old[0:ind_p] * RHO_W_KGM)
		lw_mass_remaining 		= total_liquid_mass - lw_mass_retained # mass left that will go into the ind_p node

		# now deal with the very last node where there may be just freezing or both freezing and some amount of retention
		if lw_mass_remaining <= refreeze_mass_pot[ind_p]: # all remaining water freezes
			latent_heat_released 	= lw_mass_remaining * LF_I
			self.Tz[ind_p] 			= self.Tz[ind_p] + latent_heat_released / (CP_I * self.mass[ind_p])
			self.mass[ind_p] 		= self.mass[ind_p] + lw_mass_remaining
			self.rho[ind_p] 		= self.mass[ind_p] / self.dz[ind_p]
			self.LWC[ind_p] 		= 0
			
		else: 	#some refreeze, some sticks around 
			self.mass[ind_p] 		= self.mass[ind_p] + refreeze_mass_pot[ind_p]
			self.rho[ind_p] 		= self.mass[ind_p] / self.dz[ind_p]
			self.LWC[ind_p] 		= (lw_mass_remaining - refreeze_mass_pot[ind_p]) / RHO_W_KGM
			self.Tz[ind_p] 			= T_MELT
		###################################

	self.LWC[self.LWC<0] = 0
	# else:
		# print('lots of runoff')

	return self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn, self.LWC