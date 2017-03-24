### Max Stevens
### 3/23/17

import netCDF4 as nc
import numpy as np
import scipy.io
import csv
import math
import sys
import decimal
import os
import matplotlib.pyplot as plt
from dateutil import rrule
from datetime import datetime, timedelta
import pandas as pd

spot = os.path.dirname(os.path.realpath(__file__)) #Add Folder
print spot
# os.chdir(spot)

#def write_interp(name_t, years, intp_temp):
#	with open(name_t, 'w') as csvfile:
#		#year_temp = []
#		#for i in xrange(len(years)):
#		#	year_temp.append(decimal.Decimal(str(years[i])))
#		interp_writer = csv.writer(csvfile, delimiter=',')
#		interp_writer.writerow(years)
#		interp_writer.writerow(intp_temp)



ddir = '/Users/maxstev/Documents/Grad_School/Research/FIRN/GREENLAND_CVN/Kristin-RACMOv2.3-1958-2013/'
ddir = '/Users/maxstev/Documents/Grad_School/Research/FIRN/GREENLAND_CVN/Kristin-RACMOv2.3-1958-2013/'

nc_fn_smb = spot + '/ZGRN11_smb_monthly_1958-2013.nc'
nc_fn_temp = spot + '/ZGRN11_tskin_monthly_1958-2013.nc'
nc_s = nc.Dataset(nc_fn_smb, 'r')
nc_t = nc.Dataset(nc_fn_temp, 'r')

#smb has dimensions of (time,lat,lon)
#temperature has dimensions of (time,lat,lon)

# SMB file:
smb = nc_s.variables['smb'][:] 
lat_smb = nc_s.variables['LAT'][:]
lon_smb = nc_s.variables['LON'][:]
time_smb = nc_s.variables['time'][:]

# tskin files:
tskin = nc_t.variables['tskin'][:]
lat_tskin = nc_t.variables['LAT'][:]
lon_tskin = nc_t.variables['LON'][:]
time_tskin = nc_t.variables['time'][:]

### for a specific point of interest
lat_int=72.57972
lon_int=-38.50454

dist_lat_mat=(lat_smb-lat_int)**2.0
dist_lon_mat=(lon_smb-lon_int)**2.0

dist=(dist_lat_mat+dist_lon_mat)**(0.5)
ii,jj=np.unravel_index(dist.argmin(),dist.shape)
###

# t_start = datetime(1958,1,15,12,00,00)

# tend1 = t_start + timedelta(days=365.*(2014-1958))

# timevec = np.zeros(length(time_tskin))

# for idx, dt in enumerate(rrule.rrule(rrule.MONTHLY, dtstart=t_start, until=tend1)):
# 	pass
dates = pd.date_range('1958-01-01',periods=len(time_smb),freq='MS')+pd.DateOffset(days=14)
# dates = pd.date_range('1958-01',periods=len(time_smb),freq='MS')

s1=smb[:,ii,jj]
t1=tskin[:,ii,jj]

smb_df=pd.DataFrame(s1,index=dates,columns=['smb'])
tskin_df=pd.DataFrame(t1,index=dates)

smb_df['smb'].groupby([smb_df.index.year, smb_df.index.month])#.unstack()
# racmo=pd.DataFrame({'SMB':smb_ind,'TSKIN':tskin_ind})