from constants import *
import numpy as np

def percolation_noliquid(self, iii):

	'''
	This is a simple (and very wrong, probably) percolation scheme that does not
	allow any liquid water to exist. Any water that can be refrozed by the 
	cold content of the firn will be. There is no limit of how much of the
	porosity can be filled with that water (which will refreeze). Any water that
	is not refrozed is considered runoff.
	'''
	
	porosity 			= 1 - self.rho/RHO_I # porosity
	dz_old 				= self.dz
	porespace_vol 		= porosity * self.dz #porosity in meters of each box
	melt_volume_IE  	= self.snowmeltSec[iii] * S_PER_YEAR #meters
	melt_volume_WE		= melt_volume_IE * RHO_I_MGM #meters
	melt_mass			= melt_volume_WE * 1000. #kg
	heat_to_freeze 		= melt_mass * LF_I #amount of heat needed to refreeze the melt (J)

	if (self.mass_sum==melt_mass).any():
		exactmass = True
	else:
		exactmass = False

	ind1a 				= np.where(self.mass_sum<=melt_mass)[0] # indicies of boxes that will be melted away
	num_boxes_melted 	= len(ind1a)+1 #number of boxes that melt away, include the box that is partially melted
	ind1 				= np.where(self.mass_sum>melt_mass)[0][0] 	# index which will become the new surface

	### pm is the partial melt (the box/volume that has a portion melted away)
	pm_mass				= self.mass_sum[ind1] - melt_mass 	# the remaining mass of the PM box
	pm_dz 				= pm_mass / self.rho[ind1] #remaining thickness
	pm_porespace 		= (1 - self.rho[ind1]/RHO_I) * pm_dz #porespace in the PM box
	pm_rho 				= self.rho[ind1] #density of the PM box

	cold_content_0 		= CP_I * self.mass * (T_MELT - self.Tz) #cold content of each box, i.e. how much heat to bring it to 273K
	cold_content_0_sum 	= cold_content_0.cumsum(axis=0)
	cold_content 		= cold_content_0[ind1:] #just the boxes that don't melt away
	cold_content[0]  	= CP_I * pm_mass * (T_MELT - self.Tz[ind1]) # the partial melt box has its cold content reassigned.
	cold_content_sum 	= cold_content.cumsum(axis=0)

	ind2_rel 			= np.where(cold_content_sum>heat_to_freeze)[0][0] #freeze horizon index (where the wetting front freezes), index relative to ind1
	ind2 				= ind2_rel + ind1 #absolute index on real grid (we have not removed the melted boxes yet)

	if (self.rho[ind1:ind2+1]>830.0).any(): #if there is an ice lens somewhere between the new surface and where freezing should occur
		ind2_rel 				= np.where(self.rho[ind1:]>=830.0)[0][0]
		ind2 					= ind2_rel + ind1 #the recalculated freezing front (water can not go past this level)
		cold_content_lens 		= cold_content[0:ind2_rel].sum() #cold content that is is available in space between the surface and the lens 	
		refreeze_mass 			= cold_content_lens / LF_I # this is how much should be able to refreeze in the available pore space above the ice lens.
		melt_volume_WE 			= refreeze_mass / RHO_W_KGM
		melt_volume_IE 			= melt_volume_WE / RHO_I_MGM #the volume of melt that is not runoff
		runoff_volume_duetolens = (melt_mass - refreeze_mass) * RHO_I_KGM

	else:
		runoff_volume_duetolens	= 0.0

	pore_indices 		= np.arange(ind1,ind2+1) # indicies of the boxes that are available to fill with water
	pore_indices_flip 	= np.flipud(pore_indices)
	porespace_vol[ind1] = pm_porespace
	porespace_0_sum		= porespace_vol.cumsum(axis=0)
	porespace 			= porespace_vol[ind1+1:ind2+1] #space available for the water
	porespace_sum 		= porespace.cumsum(axis=0)
	porespace_sum_flip 	= (np.flipud(porespace)).cumsum(axis=0)
	available_space 	= porespace_0_sum[ind2]-porespace_0_sum[ind1]

	if available_space < melt_volume_IE: # melt volume has already been recalculated based on how much can freeze with the cold content

		runoff_volume_duetolimitedporespace = (melt_volume_IE - porespace_sum) * RHO_I_MGM

		self.rho[ind1:ind2+1] 	= 870.0 #fill all of the boxes with water.
		self.Tz[ind1:ind2+1] 	= T_MELT

		# split up last box into several
		divider 		= num_boxes_melted
		self.rho 		= np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted)))
		self.age 		= np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_melted)))
		self.dz  		= np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted)))
		self.dz[0] 		= pm_dz
		self.dzn 		= np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
		self.dzn 		= self.dzn[0:self.compboxes]
		self.Tz  		= np.concatenate((self.Tz[ind1:-1],self.Tz[-1]*np.ones(num_boxes_melted)))
		self.bdot_mean 	= np.concatenate((self.bdot_mean[ind1:-1],self.bdot_mean[-1]*np.ones(num_boxes_melted)))
		self.z 			= self.dz.cumsum(axis = 0)
		self.z 			= np.concatenate(([0], self.z[:-1]))
		self.mass 		= self.rho*self.dz
	
	elif available_space == 0.0: #the top layer is an ice lens, so the melt runs off

		# split up last box into several
		divider 		= num_boxes_melted # ind3 should be removed and replaced with 2 new boxes.
		self.rho 		= np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted)))
		self.age 		= np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_melted)))
		self.dz  		= np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted)))
		self.dz[0] 		= pm_dz
		self.dzn 		= np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
		self.dzn 		= self.dzn[0:self.compboxes]
		self.Tz  		= np.concatenate((self.Tz[ind1:-1],self.Tz[-1]*np.ones(num_boxes_melted)))
		self.bdot_mean 	= np.concatenate((self.bdot_mean[ind1:-1],self.bdot_mean[-1]*np.ones(num_boxes_melted)))
		self.z 			= self.dz.cumsum(axis = 0)
		self.z 			= np.concatenate(([0], self.z[:-1]))
		self.mass 		= self.rho*self.dz

	else:
		runoff_volume_duetolimitedporespace = 0
		ind3a 				= np.where(porespace_sum_flip>melt_volume_IE)[0][0]
		ind3 				= ind2 - ind3a #the index of the node that is partially filled with water
		partial_volume 		= melt_volume_IE - np.sum(porespace_vol[ind3+1:ind2+1]) # pore space filled in the box that is partially filled
		leftover_porespace 	= porespace_vol[ind3]-partial_volume #open pore space in the the partially-filled box

		new_node_1_rho 		= self.rho[ind3] #split up the partial box into 2 parts
		new_node_2_rho 		= 870.0
		new_node_1_dz 		= leftover_porespace / (1 - self.rho[ind3]/RHO_I)
		new_node_2_dz 		= self.dz[ind3] - new_node_1_dz

		self.rho[ind3+1:ind2+1] = 870.0
		self.Tz[ind1:ind2+1] 	= T_MELT

		# split up last box into several
		divider 		= num_boxes_melted # ind3 should be removed and replaced with 2 new boxes.
		self.rho 		= np.concatenate((self.rho[ind1:ind3] , [new_node_1_rho,new_node_2_rho] , self.rho[ind3+1:-1] , self.rho[-1]*np.ones(num_boxes_melted-1)))
		self.age 		= np.concatenate((self.age[ind1:ind3] , [self.age[ind3],self.age[ind3]] , self.age[ind3+1:-1] , self.age[-1]*np.ones(num_boxes_melted-1)))
		dzhold 			= self.dz[ind1+1:ind3]
		dzhold2 		= self.dz[ind3+1:-1]
		self.dz 		= np.concatenate((self.dz[ind1:ind3] , [new_node_1_dz,new_node_2_dz] , self.dz[ind3+1:-1] ,self.dz[-1]/divider*np.ones(num_boxes_melted-1)))
		self.dzn 		= np.concatenate((np.zeros(num_boxes_melted), np.append(dzhold, new_node_1_dz+new_node_2_dz), dzhold2))
		self.dzn 		= self.dzn[0:self.compboxes]
		self.dz[0] 		= pm_dz
		self.Tz 		= np.concatenate((self.Tz[ind1:ind3] , [self.Tz[ind3],self.Tz[ind3]] , self.Tz[ind3+1:-1] , self.Tz[-1]*np.ones(num_boxes_melted-1)))
		self.bdot_mean 	= np.concatenate((self.bdot_mean[ind1:ind3] , [self.bdot_mean[ind3],self.bdot_mean[ind3]] , self.bdot_mean[ind3+1:-1] , self.bdot_mean[-1]*np.ones(num_boxes_melted-1)))
		self.z 			= self.dz.cumsum(axis = 0)
		self.z 			= np.concatenate(([0], self.z[:-1]))
		self.mass 		= self.rho*self.dz

	return self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn

def percolation_bucket(self, iii):

	'''
	This is the bucket scheme that allows liquid water to persist in the firn.
	It includes consideration of irreducible liquid water content (LWC) an maximum
	LWC. Water that encounters a slab of a certain density (impermeable_rho) will
	not percolate through.

	LWC is in volume (m^3), and since we are working in one dimension we assume
	that 
	'''

	maxpore_f 				= 2.0 	# factor by which the maximum filled porespace can exceed the irreducible saturation.
	impermeable_rho			= 725. 	# impermeable lens density.

	if np.any(self.LWC<0):
		print('ERROR: negative LWC')
		print('(model will continue to run)')

	melt_volume_IE  		= self.snowmeltSec[iii] * S_PER_YEAR 	# meters
	melt_volume_WE			= melt_volume_IE * RHO_I_MGM 			# meters
	melt_mass				= melt_volume_WE * 1000. 				# kg
	heat_to_freeze 			= melt_mass * LF_I 						# amount of heat needed to refreeze the melt (J)
	ind1a 					= np.where(self.mass_sum <= melt_mass)[0] 	# indicies of boxes that will be melted away
	num_boxes_melted 		= len(ind1a)+1 								# number of boxes that melt away, include the box that is partially melted
	ind1 					= np.where(self.mass_sum > melt_mass)[0][0] # index which will become the new surface

	### pm is the partial melt (the model volume that has a portion melted away)
	pm_mass 				= self.mass_sum[ind1] - melt_mass 		# the remaining mass of the PM box
	pm_dz 					= pm_mass / self.rho[ind1] 				# remaining thickness
	pm_porespace 			= (1 - self.rho[ind1]/RHO_I) * pm_dz 	# porespace in the PM box
	pm_rho 					= self.rho[ind1] 						# density of the PM box
	pm_lwc					= self.LWC[ind1]/self.dz[ind1] * pm_dz	# LWC of the PM box

	melt_boxes_LWC_vol  	= np.sum(self.LWC[0:ind1+1]) - pm_lwc #include the water mass from the boxes that melt (currently does not include from the partial melt box)
	melt_boxes_LWC_mass 	= melt_boxes_LWC_vol * RHO_W_KGM
	melt_mass_a 			= melt_mass + melt_boxes_LWC_mass
	melt_vol_a 				= melt_mass_a / RHO_W_KGM

	###################################
	### Regrid after melt
	### Melted boxes are accomodated by just adding more (new) boxes at the bottom of the column
	### Beware of this if you are not modeling to firn-ice transition depth.
	divider 				= num_boxes_melted
	self.rho 				= np.concatenate((self.rho[ind1:-1] , self.rho[-1]*np.ones(num_boxes_melted)))
	self.LWC 				= np.concatenate((self.LWC[ind1:-1] , self.LWC[-1]*np.ones(num_boxes_melted)))
	self.LWC[0] 			= pm_lwc
	self.age 				= np.concatenate((self.age[ind1:-1] , self.age[-1]*np.ones(num_boxes_melted)))
	# self.dz  				= np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted))) # this splits the last box into many.
	self.dz  				= np.concatenate((self.dz[ind1:-1] , self.dz[-1]*np.ones(num_boxes_melted))) # this adds new boxes at the bottom.
	self.dz[0] 				= pm_dz
	self.Dcon 				= np.concatenate((self.Dcon[ind1:-1] , self.Dcon[-1]*np.ones(num_boxes_melted)))
	self.dzn 				= np.concatenate((np.zeros(num_boxes_melted), self.dz[1:])) #this is not quite right because is assumes compaction for the pm box is zero.
	self.dzn 				= self.dzn[0:self.compboxes]
	self.Tz  				= np.concatenate((self.Tz[ind1:-1] , self.Tz[-1]*np.ones(num_boxes_melted)))
	self.bdot_mean 			= np.concatenate((self.bdot_mean[ind1:-1] , self.bdot_mean[-1]*np.ones(num_boxes_melted)))
	self.z 					= self.dz.cumsum(axis = 0)
	self.z 					= np.concatenate(([0] , self.z[:-1]))
	self.mass 				= self.rho * self.dz
	###################################

	##########################################
	### now working all with the new grid ####
	##########################################
	porosity 				= 1 - self.rho / RHO_I 		# porosity (unitless)
	porespace_vol 			= porosity * self.dz 		# pore space volume (meters) of each box - volume of air + water
	porespace_air			= porespace_vol - self.LWC 	# pore space that is filled with air (meters)

	cold_content			= CP_I * self.mass * (T_MELT - self.Tz) # cold content of each box, i.e. how much heat to bring it to 273K (kJ)
	cold_content_sum 		= cold_content.cumsum(axis=0)
	refreeze_mass_pot 		= cold_content / LF_I 					# how much mass of the meltwater could be refrozen due to cold content
	refreeze_mass_pot_sum 	= refreeze_mass_pot.cumsum(axis=0) 

	### calculate what the values will be after refreeze happens (pot stands for potential)
	rho_pot					= (self.mass + refreeze_mass_pot) / self.dz # what the mass of the boxes would be if the refreezemass refroze
	porosity_pot			= 1 - rho_pot / RHO_I
	porespace_vol_pot		= porosity_pot * self.dz
	porespace_air_pot		= porespace_vol_pot - self.LWC

	Wmi 					= 0.057 * (RHO_I - rho_pot) / rho_pot + 0.017 # water per snow-plus- water mass irreducible liquid water content, Langen eqn 3 unitless)
	Swi						= Wmi / (1 - Wmi) * (rho_pot * RHO_I) / (1000 * (RHO_I - rho_pot)) 	#irreducible water saturation, volume of water per porespace volume (unitless), Colbeck 1972

	maxpore 				= Swi * 2.0 # upper limit on what percentage of the porosity can be filled with water.

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
	
	try:
		ind_p 	= np.where(storage_mass_pot_sum >= total_liquid_mass)[0][0] # the layer that water will percolate to
	except: # all of the liquid is runoff.
		ind_p 	= 0
	###################################

	### if there is an impermeable layer, block water from getting through
	if np.any(self.rho[0:ind_p+1] >= impermeable_rho):

		ind_p 					= np.where(self.rho >= impermeable_rho)[0][0] #- 1 # the index of the node that has density greater than the impermeable density
		id1 					= np.where(self.LWC >  irreducible_vol_pot)[0] # indices where the liquid water content is greater than the irreducible
		id2 					= id1[id1<ind_p] 

		extra_liquid_mass		= np.sum(self.LWC[id2] * RHO_W_KGM) - np.sum(irreducible_mass_pot[id2])
		storage_mass_pot		= liquid_storage_mass_pot[0:ind_p] + refreeze_mass_pot[0:ind_p] #how much can be refrozen plus how much will stick around due to capillary
		storage_mass_pot_sum	= storage_mass_pot.cumsum(axis=0)
		total_liquid_mass 		= melt_mass_a + extra_liquid_mass

		### first, refreeze where possible
		self.mass[0:ind_p] 		= self.mass[0:ind_p] + refreeze_mass_pot[0:ind_p]
		self.rho[0:ind_p] 		= self.mass[0:ind_p] / self.dz[0:ind_p]
		self.Tz[0:ind_p] 		= T_MELT

		mass_frozen 			= np.sum(refreeze_mass_pot[0:ind_p])
		if mass_frozen >= total_liquid_mass:
			total_liquid_mass 	= 0
		else:
			total_liquid_mass 	= total_liquid_mass - mass_frozen

		### then, fill up the nodes above the ice slab
		maxLWC_mass_pot_f 		= np.flipud(maxLWC_mass_pot[0:ind_p])
		maxLWC_mass_pot_f_sum 	= maxLWC_mass_pot_f.cumsum(axis=0)

		if total_liquid_mass >= np.sum(maxLWC_mass_pot_f): # all porespace gets filled and there is runoff		
			self.LWC[0:ind_p] 		= maxLWC_pot[0:ind_p] # each node gets the maximum allowed
			# stored_water_vol 		= np.sum(self.LWC[0:ind_p]) # can calculate how much runoff there is, need to consider how much LWC there was previously
		
		else: # fill up however much porespace is needed to accomodate the meltwater

			ind_f 					= np.where(maxLWC_mass_pot_f_sum > total_liquid_mass)[0][0] #index on the flipped grid
			ind_g 					= ind_p - 1 - ind_f #index on the real grid.
			self.LWC[ind_g+1:ind_p] = maxLWC_mass_pot[ind_g + 1:ind_p] / RHO_W_KGM # fill the indices up with the maximum allowed water
			lv_mass 				= total_liquid_mass - np.sum(maxLWC_mass_pot[ind_g + 1:ind_p]) 	# leftover volume
			self.LWC[ind_g] 		= lv_mass / RHO_W_KGM 						# put that into the ind_g node
	###################################

	
	### there is not an impermeable layer, water goes to layer ind_p
	elif ind_p>0: 

		### first, up to ind_p (not inclusive)
		self.mass[0:ind_p] 		= self.mass[0:ind_p] + refreeze_mass_pot[0:ind_p]
		self.rho[0:ind_p] 		= self.mass[0:ind_p] / self.dz[0:ind_p]
		lwc_old 				= np.copy(self.LWC)
		self.LWC[0:ind_p] 		= irreducible_mass_pot[0:ind_p] / RHO_W_KGM
		self.Tz[0:ind_p] 		= T_MELT
		lw_mass_retained 		= np.sum(refreeze_mass_pot[0:ind_p]) + np.sum(irreducible_mass_pot[0:ind_p]) - np.sum(lwc_old[0:ind_p] * RHO_W_KGM)
		lw_mass_remaining 		= total_liquid_mass - lw_mass_retained # mass left that will go into the ind_p node

		### now deal with the very last node where there may be just freezing or both freezing and some amount of retention
		if lw_mass_remaining <= refreeze_mass_pot[ind_p]: # all remaining water freezes
			latent_heat_released 	= lw_mass_remaining * LF_I
			self.Tz[ind_p] 			= self.Tz[ind_p] + latent_heat_released / (CP_I * self.mass[ind_p])
			self.mass[ind_p] 		= self.mass[ind_p] + lw_mass_remaining
			self.rho[ind_p] 		= self.mass[ind_p] / self.dz[ind_p]
			self.LWC[ind_p] 		= 0
			
		else: 	# some refreeze, some sticks around 
			self.mass[ind_p] 		= self.mass[ind_p] + refreeze_mass_pot[ind_p]
			self.rho[ind_p] 		= self.mass[ind_p] / self.dz[ind_p]
			self.LWC[ind_p] 		= (lw_mass_remaining - refreeze_mass_pot[ind_p]) / RHO_W_KGM
			self.Tz[ind_p] 			= T_MELT
	###################################

	self.LWC[self.LWC<0] = 0

	return self.rho, self.age, self.dz, self.Tz, self.z, self.mass, self.dzn, self.LWC