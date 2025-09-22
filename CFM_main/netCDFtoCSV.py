#!/usr/bin/env python
'''
C. Max Stevens, 3/26/20

This script is an example of how to take climate data from
Regional Climate Models (RCMs) and convert it to CSV files
that can be used to force the CFM.

I wrote it to work with MAR or RACMO daily data. You may need to 
modify it a bit depending on the data you have. For example, if you
have monthly data the units of the SMB outputs may be different.

The important thing is that any SMB outputs (e.g. accumulation, 
melt) need to have units of meters ice equivalent per year.

The gist is that the files are loaded using the netCDF4 package.
The time from those files is converted to an array of python 
datetimes. All of the climate fields are put into a pandas 
dataframe with the datetime as the index. An additional column is
added to the dataframe that is the decimal date.

The various fields are then saved to .csv files.
'''

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
import fnmatch
import time
from scipy.spatial import cKDTree
######################################

datatype = 'MAR' # RACMO or MAR
lat_int = 72.57972 #location that you want the climate fields for
lon_int = -38.50454
######################################

def find_indices(points,lon,lat,tree=None):
    if tree is None:
        # lon,lat = lon.T,lat.T
        lonlat = np.column_stack((lon.ravel(),lat.ravel()))
        tree = cKDTree(lonlat,balanced_tree=False,compact_nodes=False)
    dist,idx = tree.query(points,k=1)
    ind = np.column_stack(np.unravel_index(idx,lon.shape))
    # print(ind)
    for i,j in ind:
        ii=i
        jj=j

    return ii,jj

def toYearFraction(date):
    '''
    Function to convert a datetime to decimal date
    '''
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

######################################

if datatype == 'RACMO':

    path_to_files = '/Volumes/Samsung_T1/RACMO/Greenland/'

    yrs = ['1957','1961','1971','1981','1991','2001','2011']

    for kk, yr in enumerate(yrs):

        fn_smb = 'smb.KNMI-{}.FGRN11.BN_RACMO2.4_FGRN11.DD.nc'.format(yr)
        fn_tskin = 'tskin.KNMI-{}.FGRN11.BN_RACMO2.4_FGRN11.DD.nc'.format(yr)
        fn_melt = '/snowmelt.KNMI-{}.FGRN11.BN_RACMO2.4_FGRN11.DD.nc'.format(yr)
        nc_s = nc.Dataset(path_to_files+fn_smb, 'r')
        nc_t = nc.Dataset(path_to_files+fn_tskin, 'r')
        nc_m = nc.Dataset(path_to_files+fn_melt, 'r')

        ### SMB file:
        smb_in = nc_s.variables['smb'][:] # units for racmo are kg/m^2/s
        lat = nc_s.variables['lat'][:]
        lon = nc_s.variables['lon'][:]
        times = nc_s.variables['time'][:]
        ### tskin files:
        tskin_in = nc_t.variables['tskin'][:]
        smelt_in = nc_m.variables['snowmelt'][:]
        if kk==0:
            ii,jj = find_indices((lon_int,lat_int),lon,lat) # get the indices for the site of interest
            ### If you know the indices, you can put them in and comemnt out the above.
            ### e.g.:
            # ii = 164 # Summit, Greenland
            # jj = 170

        s1=smb_in[:,0,ii,jj]*(365.25*24*3600)/1000/0.917 #put into units of m IE/year
        t1=tskin_in[:,0,ii,jj]
        m1 = smelt_in[:,0,ii,jj]*(365.25*24*3600)/1000/0.917
        m1[m1<0]=0

        dti = nc.num2date(times,nc_s['time'].units)

        df_y = pd.DataFrame({'smb':s1,'tskin':t1,'smelt':m1},index=dti)

        df_y['decdate'] = [toYearFraction(jj) for jj in dti]

        if kk == 0:
            df = df_y.copy()
        else:
            df = pd.concat([df,df_y])

        nc_s.close()
        nc_t.close()

        df.to_csv('RACMO_smb_example.csv',columns=['decdate','smb'],header=False,index=False)
        df.to_csv('RACMO_tskin_example.csv',columns=['decdate','tskin'],header=False,index=False)
        df.to_csv('RACMO_melt_example.csv',columns=['decdate','smelt'],header=False,index=False)
######################################

######################################
elif datatype=='MAR':
    mar_dir = '/Volumes/Samsung_T1/MAR310/Greenland/ERA_1958-2019-15km/daily_15km/'
    years_mar = np.arange(1960,2019)
    for kk, year in enumerate(years_mar):
        for file1 in os.listdir(mar_dir):             
            if fnmatch.fnmatch(file1, 'MAR*%s.nc' %year):
                fn = file1

        MARdata = nc.Dataset(mar_dir+fn,'r')

        lat = MARdata['LAT'][:,:]
        lon = MARdata['LON'][:,:]          
        if kk==0:
            ii,jj = find_indices((lon_int,lat_int),lon,lat) # get 

        t1 = MARdata['ST2'][:,0,ii,jj]
        s1 = MARdata['SMB'][:,0,ii,jj]*365.25/1000./0.917 #converted to m IE per year
        m1 = MARdata['ME'][:,0,ii,jj]*365.25/1000./0.917 #converted to m IE per year
        m1[m1<0]=0

        dti = nc.num2date(MARdata['TIME'][:],units=MARdata['TIME'].units)

        df_y = pd.DataFrame({'smb':s1,'tskin':t1,'smelt':m1},index=dti)

        df_y['decdate'] = [toYearFraction(jj) for jj in dti]

        if kk == 0:
            df = df_y.copy()
        else:
            df = pd.concat([df,df_y])

        MARdata.close()

        df.to_csv('MAR_smb_example.csv',columns=['decdate','smb'],header=False,index=False)
        df.to_csv('MAR_tskin_example.csv',columns=['decdate','tskin'],header=False,index=False)
        df.to_csv('MAR_melt_example.csv',columns=['decdate','smelt'],header=False,index=False)

##########################################
##########################################
##########################################
'''
Here is some code that will create .csv files that have 1000 years of spin up
before the main model run. In this example, we take the first 20 years RCM fields
and repeat them over and over again (50 times). That repeating series is concatenated
to the beginning of the RCM time series. 

This repeating series technique is what the RACMO (Ligtenberg and others, 2011;
Kuipers Munneke and others, 2015) team uses, to the best of my knowledge.

A more concrete example: if your RCM data goes from 1960 to 2019, and you want
to spin up for 1000 years, you need a time series going from year 960 to 1959. 
We take the 1960 to 1979 SMB, and that becomes the 960 to 979 SMB, and the 980 to 
999 SMB, and so on. 

If you use these files with the CFM, you can set 'yearsSpin' in the .json file to
some small number. The model will initialize, run the 'spin' module quickly, but
then the effective spin up occurs in the main model run.
'''

d_startyear = np.modf(df.decdate.iloc[0])[1] # Start year of the time series from RCM
yearstorepeat = 20.0 # How many of the years of the RCM time series you want to repeat
rep_inds = np.where(np.modf(df.decdate)[1].values < d_startyear+yearstorepeat)[0] #indices of those years
dec_reps = np.modf(df['decdate'].iloc[rep_inds])[0].values #decimal parts of those years
yr_reps =  np.modf(df['decdate'].iloc[rep_inds])[1].values # year portion of hte decimal year

spin_years = 1000 # How many years the spinup should be

if np.modf(spin_years/yearstorepeat)[0]!=0: # make sure that spin years is an even multiple of yearstorepeat
    spin_years = np.modf(spin_years/yearstorepeat)[1]*yearstorepeat

reps = spin_years/yearstorepeat #how many times the time series will be repeated

spin_y = (yr_reps + np.flipud(-1.0*yearstorepeat*np.arange(reps+1))[:,None]).ravel()[:int(-1*len(yr_reps))]

spin_time = spin_y + np.tile(dec_reps,int(reps))
all_time = np.concatenate((spin_time,df['decdate'].values)) # this is the vector of decimal times

# Repeat the climate data 
smb_spin = np.tile(df['smb'].iloc[rep_inds].values,int(reps))
tskin_spin = np.tile(df['tskin'].iloc[rep_inds].values,int(reps))
smelt_spin = np.tile(df['smelt'].iloc[rep_inds].values,int(reps))

# Create arrays of the time and climate fields
smb_out = np.array([all_time,np.concatenate((smb_spin,df['smb'].values))])
tskin_out = np.array([all_time,np.concatenate((tskin_spin,df['tskin'].values))])
smelt_out = np.array([all_time,np.concatenate((smelt_spin,df['smelt'].values))])

# Save them to .csv
np.savetxt('smb_example_withspin.csv',smb_out,delimiter=',',fmt='%1.4f')
np.savetxt('tskin_example_withspin.csv',tskin_out,delimiter=',',fmt='%1.4f')
np.savetxt('smelt_example_withspin.csv',smelt_out,delimiter=',',fmt='%1.4f')



