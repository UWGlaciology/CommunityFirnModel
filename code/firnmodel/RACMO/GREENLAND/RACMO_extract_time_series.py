import netCDF4 as nc
import numpy as np
import scipy.io
import csv
import math
import sys
import decimal
import os

spot = os.path.dirname(os.path.abspath(sys.argv[0])) #Add Folder

os.chdir(spot)

if len(sys.argv) != 1:
	w_lat = sys.argv[1]
	w_lon = sys.argv[2]


nc_fn_smb = 'smb_Summit.RACMO2.3_1958-2014.nc'
nc_fn_temp = 't2m_Summit.RACMO2.3_1958-2014.nc'
#nc_s = scipy.io.netcdf.netcdf_file(nc_fn_smb, 'r')
#nc_t = scipy.io.netcdf.netcdf_file(nc_fn_temp, 'r')

nc_s = nc.Dataset(nc_fn_smb, 'r')
nc_t = nc.Dataset(nc_fn_temp, 'r')

node1=5
node2=5

# for the SMB file:
# 3D
smb = nc_s.variables['smb'][:,node1,node2] 
# 2D 
lat_smb = nc_s.variables['lat'][node1,node2]
# 2D 
lon_smb = nc_s.variables['lon'][node1,node2]
# 1D - no need to alter
time_smb = nc_s.variables['time'][:]
smb = np.array(smb)
tt=np.array(time_smb)
# summ_smb = []

# t2m files:
t2m = nc_t.variables['t2m'][:,node1,node2]
lat_t2m = nc_t.variables['lat'][node1,node2]
lon_t2m = nc_t.variables['lon'][node1,node2]
time_t2m = nc_t.variables['time'][:]

array_out=np.row_stack((time_smb,smb))
np.savetxt('Summit_SMB_dailyRACMO.csv',array_out,delimiter=',')


