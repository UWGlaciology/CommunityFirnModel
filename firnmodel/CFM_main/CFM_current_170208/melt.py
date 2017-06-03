from constants import *
import numpy as np

def simple_melt(self, iii):

	
    if self.bdotSec[iii]<=self.snowmeltSec[iii]: # more melt than accumulation; no new accumulation; runoff from old box

        meltdiff = self.snowmeltSec[iii] - self.bdotSec[iii] # units m ice.
        # meltmass = meltdiff * RHO_I *S_PER_YEAR

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
	

	porosity = 1 - self.rho/RHO_I #porosity
	porespace_0 			= porosity * self.dz #porosity in meters of each box

	melt_volume_IE  	= self.snowmeltSec[iii] * S_PER_YEAR #meters
	melt_volume_WE		= melt_volume_IE * 0.917 #meters
	melt_mass			= melt_volume_WE * 1000. #kg
	heat_to_freeze 		= melt_mass * LF_I #amount of heat needed to refreeze the melt (J)

	ind1a = np.where(self.mass_sum<melt_mass)[0] # indicies of boxes that will be melted away
	num_boxes_melted = len(ind1a)+1 #number of boxes that melt away
	print 'num_boxes_melted', num_boxes_melted
	ind1 = np.where(self.mass_sum>melt_mass)[0][0] # index which will be the new surface

	pm_mass = melt_mass - self.mass_sum[ind1-1] # pm = partial melt (the box/volume that has a portion melted away)
	pm_dz = pm_mass / self.rho[ind1]
	pm_porespace = (1 - self.rho[ind1]/RHO_I) * pm_dz

	cold_content_0 		= CP_I * self.mass * (K_TO_C - self.Tz) #cold content of each box, i.e. how much heat to bring it to 273K
	cold_content = cold_content_0[ind1:] #just the boxes that don't melt away
	cold_content[0]  = CP_I * pm_mass * (K_TO_C - self.Tz[ind1]) #cold content of each box, i.e. how much heat
	cold_content_sum 	= cold_content.cumsum(axis=0)

	ind2 = np.where(cold_content_sum>heat_to_freeze)[0][0] #freeze horizon index (where the wetting front freezes)

	pore_indices = np.arange(ind1,ind2+1) # indicies of the boxes that are available to fill with water
	pore_indices_flip = np.flipud(pore_indices)

	porespace = porespace_0[ind1:ind2+1] #space available for the water
	porespace[0]=pm_porespace
	porespace_sum = porespace.cumsum(axis=0)
	porespace_sum_flip = (np.flipud(porespace)).cumsum(axis=0)

	# if self.rho[ind1:ind2].any()>= 830.0:
	# 	''

	if porespace_sum[ind2] < melt_volume_IE: #there is more melt than available pore space

		print 'not enough pore space'
		runoff_volume = (melt_volume_IE - porespace_sum) *0.917

		self.rho[ind1:ind2+1] = RHO_I
		self.Tz[ind1:ind2+1] = 273.

		# split up last box into several
		divider = num_boxes_melted

		self.rho = np.concatenate((self.rho[ind1:-1] , self.rho[-1]/divider*np.ones(num_boxes_melted)))
		self.age = np.concatenate((self.age[ind1:-1] , self.age[-1]/divider*np.ones(num_boxes_melted)))
		self.dz  = np.concatenate((self.dz[ind1:-1] , self.dz[-1]/divider*np.ones(num_boxes_melted)))
		self.dz[0] = pm_dz
		self.Tz  = np.concatenate((self.Tz[ind1:-1],self.rho[-1]/divider*np.ones(num_boxes_melted)))
		self.z = self.dz.cumsum(axis = 0)
		self.z = np.concatenate(([0], self.z[:-1]))
		self.mass = self.rho*self.dz

	else: #enough pore space to accomodate the water. 
		
		ind3a = np.where(porespace_sum_flip>melt_volume_IE)
		ind3 = ind2 - ind3a #the index of the node that is partially filled with water

		partial_volume = melt_volume_IE - np.sum(porespace[ind3+1:ind2+1]) # pore space filled in the box that is partially filled

		leftover_porespace = porespace[ind3]-partial_volume #open pore space in the the partial

		new_node_1_rho = self.rho[ind3] #split up the partial voume
		new_node_2_rho = RHO_I

		new_node_1_dz = leftover_porespace / (1 - self.rho[ind3]/RHO_I)
		new_node_2_dz = self.dz[ind3] - new_node_1_dz

		self.rho[ind3+1:ind2+1] = RHO_I
		self.Tz[ind1:ind2+1] = 273.

		# split up last box into several
		divider = num_boxes_melted+1
		# ind3 should be removed and replaced with 2 new boxes.
		self.rho = np.concatenate((self.rho[ind1:ind3] , [new_node_1_rho,new_node_2_rho] , self.rho[ind3+1:] , self.rho[-1]*np.ones(num_boxes_melted-1)))
		self.age = np.concatenate((self.age[ind1:ind3] , [self.age[ind3],self.age[ind3]] , self.age[ind3+1:] , self.age[-1]*np.ones(num_boxes_melted-1)))
		self.dz = np.concatenate((self.dz[ind1:ind3] , [new_node_1_dz,new_node_2_dz] , self.dz[ind3+1:-1] ,self.dz[-1]/divider*np.ones(num_boxes_melted-1)))
		self.dz[0] = pm_dz
		self.Tz = np.concatenate((self.Tz[ind1:ind3] , [self.Tz[ind3],self.Tz[ind3]] , self.Tz[ind3+1:] , self.Tz[-1]*np.ones(num_boxes_melted-1)))
		self.z = self.dz.cumsum(axis = 0)
		self.z = np.concatenate(([0], self.z[:-1]))
		self.mass = self.rho*self.dz

		# if there is an ice lens:
			# how much water can the porous firn refreeze? the rest is runoff.

			# remember to add the accumulation

	return self.rho, self.age, self.dz, self.Tz, self.z, self.mass




