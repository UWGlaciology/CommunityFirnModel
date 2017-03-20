'''
Max Stevens, 3/6/17
Take files with .txt and output as csv
'''

import csv
import numpy as np 

dd='/Users/maxstevens//Documents/Grad_School/Research/FIRN/CFM/CommunityFirnModel/firnmodel/CFM_main/CFM_current_170208/tyler_CU'

b_in=np.loadtxt(dd+'/CFM_age_accum.txt')
# print b_in[0,0]
# b_out=
b_in[0,:] = b_in[0,:]-b_in[0,0]
t_in=np.loadtxt(dd+'/CFM_age_temp.txt')
# print t_in[0,0]
t_in[0,:] = t_in[0,:]-t_in[0,0]
i_in=np.loadtxt(dd+'/CFM_age_isotope.txt')
# print i_in[0,0]
i_in[0,:] = i_in[0,:]-i_in[0,0]

np.savetxt('CFM_age_accum.csv',b_in,delimiter=',',fmt='%1.4f')
np.savetxt('CFM_age_temp.csv',t_in,delimiter=',',fmt='%1.4f')
np.savetxt('CFM_age_isotope.csv',i_in,delimiter=',',fmt='%1.4f')
