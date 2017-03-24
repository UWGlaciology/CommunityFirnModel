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
print 'b_in', b_in[:,0:4]
bnc = np.copy(b_in[:,0][...,None])
bnc[0,0] = 0.0
b_new = np.hstack((bnc,b_in))
# b_in[0,:] = b_in[0,:]-b_in[0,0]
b_con = np.copy(b_new)
b_con[1,:] = 0.18

t_in=np.loadtxt(dd+'/CFM_age_temp.txt')
tnc = np.copy(t_in[:,0][...,None])
tnc[0,0] = 0.0
t_new = np.hstack((tnc,t_in))

# t_in[0,:] = t_in[0,:]-t_in[0,0]
t_con = np.copy(t_new)
t_con[1,:] = -34.6

i_in=np.loadtxt(dd+'/CFM_age_isotope.txt')
inc = np.copy(i_in[:,0][...,None])
inc[0,0] = 0.0
i_new = np.hstack((inc,i_in))
# i_in[0,:] = i_in[0,:]-i_in[0,0]


np.savetxt('CFM_age_accum.csv',b_new,delimiter=',',fmt='%1.4f')
np.savetxt('CFM_age_temp.csv',t_new,delimiter=',',fmt='%1.4f')
np.savetxt('CFM_age_accum_constant.csv',b_con,delimiter=',',fmt='%1.4f')
np.savetxt('CFM_age_temp_constant.csv',t_con,delimiter=',',fmt='%1.4f')
np.savetxt('CFM_age_isotope.csv',i_in,delimiter=',',fmt='%1.4f')
