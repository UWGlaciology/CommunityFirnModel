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
from datetime import datetime, timedelta, date
import pandas as pd
# import datetime

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
# dates = pd.date_range('1958-01-01',periods=len(time_smb),freq='MS')+pd.DateOffset(days=14)
dates = pd.date_range('1958-01-01','1977-12-31',freq='MS')+pd.DateOffset(days=14)
# dates = pd.date_range('1958-01',periods=len(time_smb),freq='MS')

s1=smb[:,ii,jj]
t1=tskin[:,ii,jj]

smbdata = {'date':dates,'smb':s1[0:len(dates)]}
tskindata = {'date':dates,'tskin':t1[0:len(dates)]}

smb_df=pd.DataFrame(smbdata,columns=['date','smb'])
# s2 = smb_df.copy()
tskin_df=pd.DataFrame(tskindata,columns=['date','tskin'])
# tskin_df=pd.DataFrame(t1,index=dates)

# smb_df = smb_df.set_index([smb_df.date.dt.year, smb_df.date.dt.month]).smb.unstack()
# tskin_df = tskin_df.set_index([tskin_df.date.dt.year, tskin_df.date.dt.month]).tskin.unstack()


smb_df = smb_df.set_index([smb_df.date.dt.month, smb_df.date.dt.year]).smb.unstack()
# smb_df = smb_df.set_index([smb_df.date.dt.strftime('%b'), smb_df.date.dt.year]).smb.unstack()
tskin_df = tskin_df.set_index([tskin_df.date.dt.month, tskin_df.date.dt.year]).tskin.unstack()

# smb_df_stat = smb_df.copy()
smb_df['average'] = smb_df.mean(numeric_only=True, axis=1)
smb_df['std'] = smb_df.std(numeric_only=True, axis=1)

tskin_df['average'] = tskin_df.mean(numeric_only=True, axis=1)
tskin_df['std'] = tskin_df.std(numeric_only=True, axis=1)

years = 1000

st_year = dates.year[0] - years
# st_date = date(st_year,1,1)
# en_date = date(dates.year[0]-1,12,31)

# sdates = pd.date_range(st_date,en_date,freq='MS')+pd.DateOffset(days=14)

filler_smb = np.zeros([12,years])

for jj in xrange(12):
	# print jj
	randfill_smb = np.random.normal(smb_df.loc[jj+1,'average'],smb_df.loc[jj+1,'std'],years)
	filler_smb[jj,:] = randfill_smb

smb_spindata = np.ndarray.flatten(filler_smb,'F')

months = dates[0:12].to_pydatetime()
helper = np.vectorize(lambda x: x.timetuple().tm_yday)

decis = (helper(months)-0.5)/365

yrvec = np.arange(st_year,st_year+years)
yrtile = np.tile(yrvec[...,None],[1,12])

allspintime = yrtile + decis
allspintime_vec = np.ndarray.flatten(allspintime)



# s_dt_smb = pd.DataFrame

# df['smb'].groupby([smb_df.index.year, smb_df.index.month])#.unstack()
# racmo=pd.DataFrame({'SMB':smb_ind,'TSKIN':tskin_ind})

# df = pd.DataFrame(
#          dict(date=pd.date_range('2013-01-01', periods=42, freq='M'),
#               pb=np.random.rand(42)))

# df.set_index([df.date.dt.month, df.date.dt.year]).pb.unstack()