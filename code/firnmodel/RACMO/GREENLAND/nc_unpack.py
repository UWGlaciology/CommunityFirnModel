import netCDF4 as nc
import numpy as np
import scipy.io
import csv
import math
import sys
import decimal
import os

spot = os.path.dirname(sys.argv[0]) #Add Folder
os.chdir(spot)

#rootgrp = nc.Dataset('smb_Summit.RACMO2.3_1958-2014.nc', 'r')
def write_interp(name_t, years, intp_temp):
	with open(name_t, 'w') as csvfile:
		#year_temp = []
		#for i in xrange(len(years)):
		#	year_temp.append(decimal.Decimal(str(years[i])))
		interp_writer = csv.writer(csvfile, delimiter=',')
		interp_writer.writerow(years)
		interp_writer.writerow(intp_temp)


if len(sys.argv) != 1:
	w_lat = sys.argv[1]
	w_lon = sys.argv[2]


nc_fn_smb = '/Users/maxstev/Documents/Grad_School/PIRE/CFM/CommunityFirnModel/code/firnmodel/RACMO/GREENLAND/smb_Summit.RACMO2.3_1958-2014.nc'
nc_fn_temp = '/Users/maxstev/Documents/Grad_School/PIRE/CFM/CommunityFirnModel/code/firnmodel/RACMO/GREENLAND/t2m_Summit.RACMO2.3_1958-2014.nc'
#nc_s = scipy.io.netcdf.netcdf_file(nc_fn_smb, 'r')
#nc_t = scipy.io.netcdf.netcdf_file(nc_fn_temp, 'r')

nc_s = nc.Dataset(nc_fn_smb, 'r')
nc_t = nc.Dataset(nc_fn_temp, 'r')


#print(nc.variables)

# for the SMB file:
# 3D
smb = nc_s.variables['smb'][:] 
# 2D 
lat_smb = nc_s.variables['lat'][:]
# 2D 
lon_smb = nc_s.variables['lon'][:]
# 1D - no need to alter
time_smb = nc_s.variables['time'][:]
smb = np.array(smb)
summ_smb = []

# t2m files:
t2m = nc_t.variables['t2m'][:]
lat_t2m = nc_t.variables['lat'][:]
lon_t2m = nc_t.variables['lon'][:]
time_t2m = nc_t.variables['time'][:]

#for i in xrange(len(nc_t.variables['time'][:])):
#	time_t2m.append(decimal.Decimal(str(nc_t.variables['time'][i])))
summ_t2m = []
print time_t2m[0:4]


#print summ_t2m[0]
#print nc_t.variables['time'][:]
#time_t2m = np.around(nc_t.variables['time'][:], 8)
#print time_t2m[:]


# loop that takes in all data at certain point.
#for stps in range(len(smb)):
#	summ_smb.append(smb[stps][5][5])

for stps in range(len(t2m)):
	summ_t2m.append((t2m[stps][5][5]))

# Check the time steps

#for i in range(len(time_t2m)):
#	if time_t2m[i] != time_smb[i]:
#	print decimal.Decimal(time_t2m[i])

# TESTING Purposes. Can be deleted
#count = 0
#for i in range(len(summ_smb)):
#	print summ_smb[i]
#	count+=1

#print 'Total number of entries: ' + str(count)


#smb = np.array(smb)

# Wrties the csv file
#write_interp('smb_test2.csv', time_smb, summ_smb)
#write_interp('t2m_test.csv', time_t2m, summ_t2m)


