'''
Max Stevens, 3/6/17
Take files with .txt and output as csv
'''

import csv
import numpy as np
import os

dd = os.path.dirname(os.path.realpath(__file__))
print dd
# dd='/Users/maxstev/Documents/Grad_School/Research/FIRN/CFM/CommunityFirnModel/firnmodel/CFM_main/CFM_current_170208/tyler_CU'

b_in=np.loadtxt(dd+'/CFM_age_accum.txt')
b_in[0,:] = b_in[0,:]-b_in[0,0]
b_in[1,:] = 0.18

t_in=np.loadtxt(dd+'/CFM_age_temp.txt')
t_in[0,:] = t_in[0,:]-t_in[0,0]
t_in[1,:] = -34.6

# i_in=np.loadtxt(dd+'/CFM_age_isotope.txt')
# i_in[0,:] = i_in[0,:]-i_in[0,0]

np.savetxt('CFM_age_accum_constant.csv',b_in,delimiter=',',fmt='%1.4f')
np.savetxt('CFM_age_temp_constant.csv',t_in,delimiter=',',fmt='%1.4f')
# np.savetxt('CFM_age_isotope.csv',i_in,delimiter=',',fmt='%1.4f')
