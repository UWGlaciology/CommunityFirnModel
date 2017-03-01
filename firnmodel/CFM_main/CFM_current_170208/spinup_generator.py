import netCDF4 as nc
import numpy as np
import scipy.io
import csv
import math
import sys
import decimal
import os
import csv
# import matplotlib.pyplot as plt
from dateutil import rrule
from datetime import datetime, timedelta
import time
import re
import calendar

def toYearFraction(date):
    # def sinceEpoch(date): # returns seconds since epoch
    #     return time.mktime(date.timetuple())
    # s = sinceEpoch

    year = date.year
    # startOfThisYear = datetime(year=year, month=1, day=1)
    # startOfNextYear = datetime(year=year+1, month=1, day=1)

    # daysinyear = (startOfNextYear - startOfThisYear).days

    if calendar.isleap(year):
    	daysinyear = 366.
    else:
    	daysinyear = 365.
    
    fraction = date.timetuple().tm_yday/daysinyear

    return date.year + fraction

def spinup_generator(ddir,latind,lonind, spintype, ss_climate_length, spin_length):

	nc_s = nc.Dataset(os.path.join(ddir,'ZGRN11_smb_monthly_1958-2013.nc'),'r')
	nc_t = nc.Dataset(os.path.join(ddir,'ZGRN11_tskin_monthly_1958-2013.nc'),'r')

	# SMB file:
	smb = nc_s.variables['smb'][:,latind,lonind] 
	lat_smb = nc_s.variables['LAT'][:]
	lon_smb = nc_s.variables['LON'][:]
	time_smb = nc_s.variables['time'][:]

	# tskin files:
	tskin = nc_t.variables['tskin'][:,latind,lonind] 
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

	tunits = nc_s.variables['time'].units
	tu = map(int, re.findall('\d+',tunits))
	tst = datetime(tu[0],tu[1],tu[2])
	# print 'tst', tst
	# tst = datetime(1958,1,15)
	# tend1 = tst + timedelta(days=365.25*(2013-1958))
	decdates_main=np.empty(len(time_smb))

	for idx, dt in enumerate(rrule.rrule(rrule.MONTHLY, dtstart=tst, count=len(time_smb))):
		if idx<5:
			print dt
		dd=toYearFraction(dt)
		decdates_main[idx]=dd

	spinstart = datetime(tst.year-spin_length, tst.month, tst.day)
	decdates_spin=np.empty(spin_length*12+1)

	for idx, dt in enumerate(rrule.rrule(rrule.MONTHLY, dtstart=spinstart, until=tst)):
		# print dt
		dd=toYearFraction(dt)
		decdates_spin[idx]=dd
	print decdates_spin[-3:]

	with open('temperature_test.csv','w') as f:
		writer = csv.writer(f)
		writer.writerow(decdates_main)
		writer.writerow(tskin)

	with open('bdot_test.csv','w') as f:
		writer = csv.writer(f)
		writer.writerow(decdates_main)
		writer.writerow(tskin)
 


if __name__ == '__main__':

	ddir='/Users/maxstev/Documents/Grad_School/Research/FIRN/GREENLAND_CVN/Data/Kristin-RACMOv2.3-1958-2013'

	latind = 150

	lonind = 150

	spintype = 'random' # 'random' or 'loop'

	ss_climate_length = 20 # years

	spin_length = 600

	spinup_generator(ddir, latind, lonind, spintype, ss_climate_length, spin_length)


