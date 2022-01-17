#!/usr/bin/env python
'''
Surface Energy balance class
'''

import numpy as np 
# from solver import solver
from solver import transient_solve_TR
# from Gasses import gasses 
# from Sites import sites 
from reader import read_input, read_init
import json
import scipy.interpolate as interpolate
from constants import *
import os
import sys

class SEB:
	'''
	Class to handle surface energy balance in the CFM
	'''

	def __init__(self,spin,config):

		self.c = config
		
		localpath = os.path.join(self.c['InputFileFolder'],self.c['InputFileNameTemp'])
		seb_in = pd.read_csv(os.path.join(os.getcwd(),localpath))

		# Would be good to build functionality to take pandas inputs directly? Problem with datetime long ago?

		self.SWdown = seb_in['SWdown'].values

		# self.D_sh = 15 # Sensible heat flux coefficient, Born et al. 2019 [W m^-2 K^-1] 


	def ()