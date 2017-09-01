'''
Max Stevens
8/23/17
Use this file to create spin-up files for the CFM
This is the right one to use!
Melt for spin up is calculated by performing an rbf regression on the temerature data, using the actual data as the training set.
'''


import netCDF4 as nc
import numpy as np
import scipy.io
import csv
import math
import sys
import decimal
import os
import sys
import matplotlib.pyplot as plt
from dateutil import rrule
from datetime import datetime, timedelta, date
import pandas as pd
import fnmatch
from scipy.spatial import cKDTree
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.svm import SVR

# import datetime
def find_indices(points,lon,lat,tree=None):
    if tree is None:
        # lon,lat = lon.T,lat.T
        lonlat = np.column_stack((lon.ravel(),lat.ravel()))
        tree = cKDTree(lonlat)
    dist,idx = tree.query(points,k=[1])
    ind = np.column_stack(np.unravel_index(idx,lon.shape))
    print(ind)
    for i,j in ind:
    	ii=i
    	jj=j
    # return [(i,j) for i,j in ind]
    return ii,jj #, [(i,j) for i,j in ind]

writer=True
spot = os.path.dirname(os.path.realpath(__file__)) #Add Folder
print(spot)
datatype = 'MAR'
# datatype = 'RACMO'
print('datatype is ', datatype)
# sites=['Summit','DYE2','KANU','EKT','NASASE','SADDLE','CRAWFORD','EGRIP']
sites=['DYE2']
SPY = 365.25*24*3600

for site in sites:

	# site = 'EGRIP'
	print("site = ", site)
	resultsdir = 'inputdata/'+ datatype + 'input/' + site
	try: 
		os.makedirs(resultsdir)
	except:
		pass
		# print(resultsdir' already exists. Do you want to continue and overwrite? (y/n)')
		# ip = input(resultsdir+' already exists. Do you want to continue and overwrite? (y/n)')
		# if ip == 'y':
		# 	pass
		# else:
		# 	sys.exit()

	if site == 'Summit':
		lat_int = 72.57972 #summit
		lon_int = -38.50454	
	elif site == 'DYE2':
		lat_int = 66.4806 #this is DYE2
		lon_int = -46.2831
	elif site == 'KANU':
		lat_int = 67.000383
		lon_int = -47.02615
	elif site == 'EKT':
		lat_int = 66.9854
		lon_int = -44.39465
	elif site == 'NASASE':
		lat_int = 66.47768
		lon_int = -42.49635
	elif site == 'SADDLE':
		lat_int = 65.9994
		lon_int = -44.50248
	elif site == 'CRAWFORD':
		lat_int = 69.87615
		lon_int = -47.03112
	elif site == 'EGRIP':
		lat_int = 75.62556
		lon_int = -35.97803


	# os.chdir(spot)

	#def write_interp(name_t, years, intp_temp):
	#	with open(name_t, 'w') as csvfile:
	#		#year_temp = []
	#		#for i in xrange(len(years)):
	#		#	year_temp.append(decimal.Decimal(str(years[i])))
	#		interp_writer = csv.writer(csvfile, delimiter=',')
	#		interp_writer.writerow(years)
	#		interp_writer.writerow(intp_temp)


	#### Where to find the data
	# ddir = '/Users/maxstev/Documents/Grad_School/Research/FIRN/GREENLAND_CVN/Kristin-RACMOv2.3-1958-2013/'

	if datatype=='RACMO':
		yrs = ['1957','1961','1971','1981','1991','2001','2011']
		# ddir = '/Users/maxstev/Documents/Grad_School/Research/FIRN/CFM/CommunityFirnModel/firnmodel/CFM_main/CFM_current_170208/RACMO'
		ddir = '/Volumes/Samsung_T1/RACMO/'
		sd={}
		td={}
		md={}
		tid={}


		# evap.KNMI-1957.FGRN11.BN_RACMO2.4_FGRN11.DD.nc
		####
		for idx, kk in enumerate(yrs):
			nc_fn_smb = ddir + '/smb.KNMI-%s.FGRN11.BN_RACMO2.4_FGRN11.DD.nc' %kk
			# nc_fn_smb = spot + '/smb_Summit.RACMO2.3_1958-2014.nc'
			nc_fn_temp = ddir + '/tskin.KNMI-%s.FGRN11.BN_RACMO2.4_FGRN11.DD.nc' %kk
			nc_fn_melt = ddir + '/snowmelt.KNMI-%s.FGRN11.BN_RACMO2.4_FGRN11.DD.nc' %kk

			# nc_fn_smb2 = ddir + '/smb.KNMI-1961.FGRN11.BN_RACMO2.4_FGRN11.DD.nc'
			# # nc_fn_smb = spot + '/smb_Summit.RACMO2.3_1958-2014.nc'
			# nc_fn_temp2 = ddir + '/tskin.KNMI-1961.FGRN11.BN_RACMO2.4_FGRN11.DD.nc'
			# nc_fn_melt2 = ddir + '/snowmelt.KNMI-1961.FGRN11.BN_RACMO2.4_FGRN11.DD.nc'

			# nc_fn_temp = spot + '/t2m_Summit.RACMO2.3_1958-2014.nc'
			nc_s = nc.Dataset(nc_fn_smb, 'r')
			nc_t = nc.Dataset(nc_fn_temp, 'r')
			nc_m = nc.Dataset(nc_fn_melt, 'r')

			# nc_s2 = nc.Dataset(nc_fn_smb2, 'r')

			#smb has dimensions of (time,lat,lon)
			#temperature has dimensions of (time,lat,lon)

			# SMB file:
			smb_in = nc_s.variables['smb'][:] 
			lat_smb = nc_s.variables['lat'][:]
			lon_smb = nc_s.variables['lon'][:]
			time_smb = nc_s.variables['time'][:]

			# tskin files:
			tskin_in = nc_t.variables['tskin'][:]
			# lat_tskin = nc_t.variables['lat'][:]
			# lon_tskin = nc_t.variables['lon'][:]
			# time_tskin = nc_t.variables['time'][:]

			# melt files:
			smelt_in = nc_m.variables['snowmelt'][:]
			# lat_smelt = nc_m.variables['lat'][:]
			# lon_smelt = nc_m.variables['lon'][:]
			# time_smelt = nc_m.variables['time'][:]

			lat = lat_smb
			lon = lon_smb

			if idx==0:

				ii,jj = find_indices((lon_int,lat_int),lon,lat)

				print(ii,jj)
				print(lat[ii,jj])
				print(lon[ii,jj])


			### for a specific point of interest
			# lat_int=72.57972 #this is summit
			# lon_int=-38.50454

			# lat_int=66.4806 #this is DYE2
			# lon_int=-46.2831

			# dist_lat_mat=(lat_smb-lat_int)**2.0
			# dist_lon_mat=(lon_smb-lon_int)**2.0

			# dist=(dist_lat_mat+dist_lon_mat)**(0.5)
			# ii,jj=np.unravel_index(dist.argmin(),dist.shape) # indices closest to specified point
			# print(ii,jj)
			# tt = nc_s['time']
			s1=smb_in[:,0,ii,jj]*24*3600 #put into units of kg/m^2/day
			t1=tskin_in[:,0,ii,jj]
			m1 = smelt_in[:,0,ii,jj]*24*3600 #put into units of m IE/day
			m1[m1<0]=0 # current version of RACMO goes through end of 2016

			sd['smb%s' %kk]=s1
			td['tskin%s' %kk]=t1
			md['smelt%s' %kk]=m1
			tid['time%s' %kk]=time_smb

		smb_all = np.concatenate([v for k,v in sorted(sd.items())], 0)
		tskin_all = np.concatenate([v for k,v in sorted(td.items())], 0)
		smelt_all = np.concatenate([v for k,v in sorted(md.items())], 0)
		time_all = np.concatenate([v for k,v in sorted(tid.items())], 0)

		date_start = pd.to_datetime('19500101',format='%Y%m%d')
		temp = pd.to_timedelta(time_all,unit='D')
		dates_daily = date_start + temp

		date_end = 2016

		dates = pd.date_range('1958-01-01','1977-12-31',freq='MS')+pd.DateOffset(days=14)
		dates_all = pd.date_range('1958-01-01','%s-12-31' %date_end,freq='MS')+pd.DateOffset(days=14)
		dates_years = pd.date_range('1958-01-01','%s-12-31' %date_end,freq='AS')

		st_date=pd.to_datetime('19580101',format='%Y%m%d')
		trun_en_date=pd.to_datetime('19771231',format='%Y%m%d')

		#### smb ##############
		smb_data_daily = {'date':dates_daily,'smb':smb_all}
		smb_df_daily = pd.DataFrame(smb_data_daily,columns=['smb'],index=smb_data_daily['date'])
		# smb_df_daily = pd.DataFrame(smb_data,columns=['date','smb'])
		smb_df_daily = smb_df_daily[smb_df_daily.index>=st_date]

		smb_df_monthly_all = smb_df_daily.resample('MS').sum() # this is the smb for each month, in kg/m^2 (i.e. the value for December of 2016 is the amount of snow that fell that month)
		smb_df_monthly_all.index = smb_df_monthly_all.index+pd.DateOffset(days=14)
		smb_df_monthly_all['smb']= smb_df_monthly_all['smb']/917*12 # convert to m I.E. per year 

		smb_df = smb_df_monthly_all[st_date:trun_en_date] #1958-1977
		smbdata={'date':smb_df.index,'smb':smb_df['smb'].values}
		smb_df = smb_df.set_index([smb_df.index.month,smb_df.index.year]).unstack()

		# smb_df_monthly_all.set_index([smb_df_monthly_all.index.month,smb_df_monthly_all.index.year]).unstack()

		s1=smb_df_monthly_all['smb'][smb_df_monthly_all.index>st_date].values
		#######################

		###### tskin ##########
		tskin_data_daily = {'date':dates_daily,'tskin':tskin_all}
		tskin_df_daily = pd.DataFrame(tskin_data_daily,columns=['tskin'],index=tskin_data_daily['date'])
		# tskin_df_daily = pd.DataFrame(tskin_data,columns=['date','tskin'])
		tskin_df_daily = tskin_df_daily[tskin_df_daily.index>=st_date]

		tskin_df_monthly_all = tskin_df_daily.resample('MS').mean()
		tskin_df_monthly_all.index = tskin_df_monthly_all.index+pd.DateOffset(days=14)

		tskin_df = tskin_df_monthly_all[st_date:trun_en_date]
		tskindata={'date':tskin_df.index,'tskin':tskin_df['tskin'].values}
		tskin_df = tskin_df.set_index([tskin_df.index.month,tskin_df.index.year]).unstack()

		# tskin_df_monthly_all.set_index([tskin_df_monthly_all.index.month,tskin_df_monthly_all.index.year]).unstack()

		t1=tskin_df_monthly_all['tskin'][tskin_df_monthly_all.index>st_date].values
		########################

		####### snowmelt ##########
		smelt_data_daily = {'date':dates_daily,'smelt':smelt_all}
		smelt_df_daily = pd.DataFrame(smelt_data_daily,columns=['smelt'],index=smelt_data_daily['date'])
		# smelt_df_daily = pd.DataFrame(smelt_data,columns=['date','smelt'])
		smelt_df_daily = smelt_df_daily[smelt_df_daily.index>=st_date]

		smelt_df_monthly_all = smelt_df_daily.resample('MS').sum()
		smelt_df_monthly_all.index = smelt_df_monthly_all.index+pd.DateOffset(days=14)
		smelt_df_monthly_all['smelt'] = smelt_df_monthly_all['smelt'] / 917 * 12 # convert to m I.E. per year
		smelt_df_monthly_all['smelt'][smelt_df_monthly_all['smelt']<=1.0e-3]=0.0

		smelt_df = smelt_df_monthly_all[st_date:trun_en_date]
		smeltdata={'date':smelt_df.index,'smelt':smelt_df['smelt'].values}
		smelt_df = smelt_df.set_index([smelt_df.index.month,smelt_df.index.year]).unstack()

		smelt_df_monthly_all.set_index([smelt_df_monthly_all.index.month,smelt_df_monthly_all.index.year]).unstack()

		m1=smelt_df_monthly_all['smelt'][smelt_df_monthly_all.index>st_date].values
		##################

		# ['2012-07-01':'2012-08-01']

		##### end RACMO #########

	elif datatype == 'MAR': # current version of MAR goes through end of 2015
		'''
		s1,t1,m1 are the full time series of the data (1958-2016, monthly), values only
		dates is the monthly datetimeindex, 1958-1977
		dates all is the monthly datetimeindex, 1958-2016
		dates_years is the annual datetimeindex, 1958-2016
		smbdata is a dict, keys 'date' and 'smb', 1958-1977
		smb_df is dataframe, monthly, 1958 - 1977
		smb_df_monthly_all is a DataFrame, monthly, 1958-2016
		'''

		mar_dir = '/Volumes/Samsung_T1/mar2'
		years = np.arange(1958,2016)

		t1=np.zeros(len(years)*12)
		s1=np.zeros(len(years)*12)
		m1=np.zeros(len(years)*12)

		for kk in range(len(years)):
			year = years[kk]
			for file in os.listdir(mar_dir):
			    if fnmatch.fnmatch(file, 'MAR*%s.nc' %year):
			        fn = file
			        # print type(fn)
			        # print file
			fp = os.path.join(mar_dir,fn)
			# print fp
			rgr=nc.Dataset(fp,'r')

			lat = rgr['LAT'][:,:]
			lon = rgr['LON'][:,:]

			if kk==0:
				ii,jj = find_indices((lon_int,lat_int),lon,lat)
				print(lat[ii,jj])
				print(lon[ii,jj])
				print('x=',rgr['x'][jj])
				print('y=',rgr['y'][ii])

			t1[12*kk:12*kk+12] = rgr['STcorr'][:,ii,jj]
			s1[12*kk:12*kk+12] = rgr['SMBcorr'][:,ii,jj]*12/1000./0.917 #converted to m IE per year
			m1[12*kk:12*kk+12] = rgr['MEcorr'][:,ii,jj]*12/1000./0.917 #converted to m IE per year
			m1[m1<0]=0

			date_end=2015
	###

			dates = pd.date_range('1958-01-01','1977-12-31',freq='MS')+pd.DateOffset(days=14)
			dates_all = pd.date_range('1958-01-01','%s-12-31' %date_end,freq='MS')+pd.DateOffset(days=14)
			dates_years = pd.date_range('1958-01-01','%s-12-31' %date_end,freq='AS')

			smbdata = {'date':dates,'smb':s1[0:len(dates)]}
			tskindata = {'date':dates,'tskin':t1[0:len(dates)]}
			smeltdata = {'date':dates,'smelt':m1[0:len(dates)]}

			smb_df=pd.DataFrame(smbdata,columns=['date','smb'])
			tskin_df=pd.DataFrame(tskindata,columns=['date','tskin'])
			smelt_df=pd.DataFrame(smeltdata,columns=['date','smelt'])

			smb_df_monthly_all=pd.DataFrame(s1,index=dates_all,columns=['smb'])
			tskin_df_monthly_all=pd.DataFrame(t1,index=dates_all,columns=['tskin'])
			melt_df_monthly_all=pd.DataFrame(m1,index=dates_all,columns=['smelt'])

			adf = pd.concat([smb_df_monthly_all.resample('AS').sum(),tskin_df_monthly_all.resample('AS').mean(),melt_df_monthly_all.resample('AS').sum()],axis=1)

			smb_year = np.array([dates_years.year,adf['smb']/1000])
			tskin_year = np.array([dates_years.year,adf['tskin']])
			smelt_year = np.array([dates_years.year,adf['smelt']/1000])

			# np.savetxt('smb_year.csv',smb_year,delimiter=',',fmt='%1.4f')
			# np.savetxt('tskin_year.csv',tskin_year,delimiter=',',fmt='%1.4f')
			# np.savetxt('smelt_year.csv',smelt_year,delimiter=',',fmt='%1.4f')

			# smb_df = smb_df.set_index([smb_df.date.dt.year, smb_df.date.dt.month]).smb.unstack()
			# tskin_df = tskin_df.set_index([tskin_df.date.dt.year, tskin_df.date.dt.month]).tskin.unstack()


			smb_df = smb_df.set_index([smb_df.date.dt.month, smb_df.date.dt.year]).smb.unstack()
			# smb_df = smb_df.set_index([smb_df.date.dt.year]).smb.unstack()
			# smb_df = smb_df.set_index([smb_df.date.dt.strftime('%b'), smb_df.date.dt.year]).smb.unstack()
			tskin_df = tskin_df.set_index([tskin_df.date.dt.month, tskin_df.date.dt.year]).tskin.unstack()
			smelt_df = smelt_df.set_index([smelt_df.date.dt.month, smelt_df.date.dt.year]).smelt.unstack()

			###### END MAR ##########




	
	s2 = smb_df.copy()
	t2 = tskin_df.copy()
	m2 = smelt_df.copy()

	s_flip = np.fliplr(s2.as_matrix())
	t_flip = np.fliplr(t2.as_matrix())
	m_flip = np.fliplr(m2.as_matrix())

	tskin_df_values = tskin_df.copy()
	smelt_df_values = smelt_df.copy()

	# smb_df_stat = smb_df.copy()
	smb_df['average'] = smb_df.mean(numeric_only=True, axis=1)
	smb_df['std'] = smb_df.std(numeric_only=True, axis=1)

	tskin_df['average'] = tskin_df.mean(numeric_only=True, axis=1)
	tskin_df['std'] = tskin_df.std(numeric_only=True, axis=1)

	smelt_df['average'] = smelt_df.mean(numeric_only=True, axis=1)
	smelt_df['std'] = smelt_df.std(numeric_only=True, axis=1)
	smelt_df.loc[smelt_df['average']<1.0e-3]=0

	### time vector that data will be written
	# years = 0
	years = 1000

	st_year = dates.year[0] - years
	# st_date = date(st_year,1,1)
	# en_date = date(dates.year[0]-1,12,31)

	# sdates = pd.date_range(st_date,en_date,freq='MS')+pd.DateOffset(days=14)
	months = dates[0:12].to_pydatetime()
	helper = np.vectorize(lambda x: x.timetuple().tm_yday)
	decis = (helper(months)-0.5)/365

	yrvec = np.arange(st_year,st_year+years)
	yrtile = np.tile(yrvec[...,None],[1,12])
	allspintime = yrtile + decis
	allspintime_vec = np.ndarray.flatten(allspintime)

	yrvec2 = np.arange(dates.year[0],date_end+1)
	yrtile2 = np.tile(yrvec2[...,None],[1,12])
	allspintime2 = yrtile2 + decis
	allspintime2_vec = np.ndarray.flatten(allspintime2)

	time_out = np.concatenate((allspintime_vec,allspintime2_vec))
	#####

	##### make the random time series #####
	sno = 40 # number of random time series to make

	for ii in range(sno):

		filler_smb = np.zeros([12,years])
		filler_tskin = np.zeros([12,years])
		filler_smelt = np.zeros([12,years])

		for jj in range(12):
			# print jj
			randfill_smb = np.random.normal(smb_df.loc[jj+1,'average'],smb_df.loc[jj+1,'std'],years)
			filler_smb[jj,:] = randfill_smb
			randfill_tskin = np.random.normal(tskin_df.loc[jj+1,'average'],tskin_df.loc[jj+1,'std'],years)
			filler_tskin[jj,:] = randfill_tskin
			try:
				randfill_smelt = np.random.normal(smelt_df.loc[jj+1,'average'],smelt_df.loc[jj+1,'std'],years)
			except:
				randfill_smelt = np.zeros(years)

			# [randfill_tskin,randfill_smelt] = np.random.multivariate_normal( [tskin_df.loc[jj+1,'average'] , smelt_df.loc[jj+1,'average']],np.cov())
			
			randfill_smelt[randfill_smelt<0] = 0.0
			filler_smelt[jj,:] = randfill_smelt 

		smb_spindata = np.ndarray.flatten(filler_smb,'F')
		tskin_spindata = np.ndarray.flatten(filler_tskin,'F')
		smelt_spindata = np.ndarray.flatten(filler_smelt,'F')

		# create model for melt/temperature relationship
		# clf = SVR(kernel='rbf', C=1, epsilon=0.2)
		# clf.fit(t1.reshape(-1,1),m1)

		# smelt_spindata = clf.predict(tskin_spindata.reshape(-1,1))
		# smelt_spindata[smelt_spindata<1.0e-3]=0.0

		# # if clftest:
		if datatype=='MAR':
			thresh=-20.0
		else:
			thresh=273-20.0
		mask = t1>thresh
		mm = m1[mask]
		tt = t1[mask]

		mask2=mm>0

		p=np.poly1d(np.polyfit(tt[mask2],np.log(mm[mask2]),2))
		wm = tskin_spindata>thresh
		cm = tskin_spindata<=thresh
		smelt_spindata[wm]=np.exp(p(tskin_spindata[wm]))
		smelt_spindata[cm]=0.0

		if site==('EGRIP' or 'Summit'):
			smelt_spindata[:]=0.0

		# plt.scatter(t1, m1,  color='black')
		# plt.scatter(tt,mm, color='red')
		# plt.scatter(tskin_spindata, smelt_spindata, color='blue',marker='+')
		# plt.show()

		# plt.scatter(tt[mask2],np.log(mm[mask2]))
		# plt.scatter(tt[mask2],p(tt[mask2]))
		# plt.show()

		smb_d = np.concatenate((smb_spindata,s1))
		tskin_d = np.concatenate((tskin_spindata,t1))
		smelt_d = np.concatenate((smelt_spindata,m1))

		smb_out = np.array([time_out,smb_d])
		tskin_out = np.array([time_out,tskin_d])
		smelt_out = np.array([time_out,smelt_d])

		sfn = resultsdir + '/' + site + '_smb_%s_' %datatype+str(ii)+'.csv'
		tfn = resultsdir + '/' + site + '_tskin_%s_' %datatype+str(ii)+'.csv'
		mfn = resultsdir + '/' + site + '_melt_%s_' %datatype+str(ii)+'.csv'

		# sfn = 'melt_test_smb_%s' %datatype +'.csv'
		# tfn = 'melt_test_tskin_%s' %datatype +'.csv'
		# mfn = 'melt_test_melt_%s' %datatype +'.csv'

		# print(sfn)

		##########
		if writer:
			np.savetxt(sfn,smb_out,delimiter=',',fmt='%1.4f')
			np.savetxt(tfn,tskin_out,delimiter=',',fmt='%1.4f')
			np.savetxt(mfn,smelt_out,delimiter=',',fmt='%1.4f')
		##########

		# rho_vec = np.random.normal(329.4,53.0,len(smb_d))
		# rho_out = np.array([time_out,rho_vec])
	
	# dd=[smelt_df['average'],tskin_df['average']]
	# ddd = pd.concat(dd,axis=1,join_axes=[tskin_df.index])
	# print(ddd)
	# input('press enter to continue')


	##### Make the looping time series
	s_rev = np.ndarray.flatten(s_flip,'F') # this is the reversed time series
	t_rev = np.ndarray.flatten(t_flip,'F')
	m_rev = np.ndarray.flatten(m_flip,'F')
	s_loop = np.concatenate((smbdata['smb'],s_rev)) # this is the 40 - year forward-backward series
	t_loop = np.concatenate((tskindata['tskin'],t_rev))
	m_loop = np.concatenate((smeltdata['smelt'],m_rev))

	s_loop_spin = np.tile(s_loop, len(allspintime_vec)//len(s_loop))
	t_loop_spin = np.tile(t_loop, len(allspintime_vec)//len(t_loop))
	m_loop_spin = np.tile(m_loop, len(allspintime_vec)//len(m_loop))

	s_loop_out_d = np.concatenate((s_loop_spin,s1))
	t_loop_out_d = np.concatenate((t_loop_spin,t1))
	m_loop_out_d = np.concatenate((m_loop_spin,m1)) 

	smb_loop_out = np.array([time_out,s_loop_out_d])
	tskin_loop_out = np.array([time_out,t_loop_out_d])
	melt_loop_out = np.array([time_out,m_loop_out_d])

	if writer:
		np.savetxt(resultsdir + '/' + site+'_smb_%s_loop.csv' %datatype,smb_loop_out,delimiter=',',fmt='%1.4f')
		np.savetxt(resultsdir + '/' + site+'_tskin_%s_loop.csv' %datatype,tskin_loop_out,delimiter=',',fmt='%1.4f')
		np.savetxt(resultsdir + '/' + site+'_melt_%s_loop.csv' %datatype,melt_loop_out,delimiter=',',fmt='%1.4f')
	######

	###### make constant time series for spin up
	s_out_con = np.concatenate((np.mean(smbdata['smb'])*np.ones_like(allspintime_vec),s1))
	t_out_con = np.concatenate((np.mean(tskindata['tskin'])*np.ones_like(allspintime_vec),t1))
	m_out_con = np.concatenate((np.mean(smeltdata['smelt'])*np.ones_like(allspintime_vec),m1))

	# s_out_con_d = np.concatenate((s_out_con,s1))
	# t_out_con_d = np.concatenate((t_out_con,t1))

	smb_con_out = np.array([time_out,s_out_con])
	tskin_con_out = np.array([time_out,t_out_con])
	melt_con_out = np.array([time_out,m_out_con]) 

	if writer:
		np.savetxt(resultsdir + '/' + site + '_smb_%s_con.csv' %datatype,smb_con_out,delimiter=',',fmt='%1.4f')
		np.savetxt(resultsdir + '/' + site + '_tskin_%s_con.csv' %datatype,tskin_con_out,delimiter=',',fmt='%1.4f')
		np.savetxt(resultsdir + '/' + site + '_melt_%s_con.csv' %datatype,melt_con_out,delimiter=',',fmt='%1.4f')

	#####

	# np.savetxt('Summit_spin_smb_base.csv',smb_out,delimiter=',',fmt='%1.4f')
	# np.savetxt('Summit_spin_temp_base.csv',tskin_out,delimiter=',',fmt='%1.4f')
	# np.savetxt('Summit_spin_rho_base.csv',rho_out,delimiter=',',fmt='%1.4f')

	# stind=12240
	# reclen = 12*36

	# ss=smb_d[stind:]
	# tt=tskin_d[stind:]

	##### Mocon: monthly constant (each month has a constant value through the spin up)
	smb_mocon = np.tile(smb_df['average'],years)
	tskin_mocon = np.tile(tskin_df['average'],years)
	melt_mocon = np.tile(smelt_df['average'],years)

	melt_mocon[melt_mocon<1.0e-3]=0.0

	s_mocon_out_d = np.concatenate((smb_mocon,s1))
	t_mocon_out_d = np.concatenate((tskin_mocon,t1))
	m_mocon_out_d = np.concatenate((melt_mocon,m1)) 

	smb_mocon_out = np.array([time_out,s_mocon_out_d])
	tskin_mocon_out = np.array([time_out,t_mocon_out_d])
	melt_mocon_out = np.array([time_out,m_mocon_out_d])

	if writer:
		np.savetxt(resultsdir + '/' + site+'_smb_%s_mocon.csv' %datatype,smb_mocon_out,delimiter=',',fmt='%1.4f')
		np.savetxt(resultsdir + '/' + site+'_tskin_%s_mocon.csv' %datatype,tskin_mocon_out,delimiter=',',fmt='%1.4f')
		np.savetxt(resultsdir + '/' + site+'_melt_%s_mocon.csv' %datatype,melt_mocon_out,delimiter=',',fmt='%1.4f')




	### Below was making file for Lora Koenig, Summit meeting (email exchanges late March 2017)
	# tnoise05=np.random.normal(0,0.5,reclen)
	# tskin_out_05 = tskin_out.copy()
	# tskin_out_05[1,stind:]=tskin_out_05[1,stind:]+tnoise05
	# np.savetxt('Summit_spin_temp_05.csv',tskin_out_05,delimiter=',',fmt='%1.4f')


	# tnoise10=np.random.normal(0,1.0,reclen)
	# tskin_out_10 = tskin_out.copy()
	# tskin_out_10[1,stind:]=tskin_out_10[1,stind:]+tnoise10
	# np.savetxt('Summit_spin_temp_10.csv',tskin_out_10,delimiter=',',fmt='%1.4f')

	# tnoise15=np.random.normal(0,1.5,reclen)
	# tskin_out_15 = tskin_out.copy()
	# tskin_out_15[1,stind:]=tskin_out_15[1,stind:]+tnoise15
	# np.savetxt('Summit_spin_temp_15.csv',tskin_out_15,delimiter=',',fmt='%1.4f')

	# smb_rm = np.mean(ss)
	# # smb_std = np.std(ss)

	# snoise05=np.random.normal(0,0.05*smb_rm,reclen)
	# smb_out_05 = smb_out.copy()
	# smb_out_05[1,stind:]=smb_out_05[1,stind:]+snoise05
	# np.savetxt('Summit_spin_smb_05.csv',smb_out_05,delimiter=',',fmt='%1.4f')


	# snoise10=np.random.normal(0,0.1*smb_rm,reclen)
	# smb_out_10 = smb_out.copy()
	# smb_out_10[1,stind:]=smb_out_10[1,stind:]+snoise10
	# np.savetxt('Summit_spin_smb_10.csv',smb_out_10,delimiter=',',fmt='%1.4f')

	# snoise20=np.random.normal(0,0.2*smb_rm,reclen)
	# smb_out_20 = smb_out.copy()
	# smb_out_20[1,stind:]=smb_out_20[1,stind:]+snoise20
	# np.savetxt('Summit_spin_smb_20.csv',smb_out_20,delimiter=',',fmt='%1.4f')

	# rho_rm = 330.

	# rnoise05=np.random.normal(0,0.05*rho_rm,reclen)
	# rho_out_05 = rho_out.copy()
	# rho_out_05[1,stind:]=rho_out_05[1,stind:]+rnoise05
	# np.savetxt('Summit_spin_rho_05.csv',rho_out_05,delimiter=',',fmt='%1.4f')

	# rnoise10=npy.random.normal(0,0.1*rho_rm,reclen)
	# rho_out_10 = rho_out.copy()
	# rho_out_10[1,stind:]=rho_out_10[1,stind:]+rnoise10
	# np.savetxt('Summit_spin_rho_10.csv',rho_out_10,delimiter=',',fmt='%1.4f')

	# rnoise20=np.random.normal(0,0.2*rho_rm,reclen)
	# rho_out_20 = rho_out.copy()
	# rho_out_20[1,stind:]=rho_out_20[1,stind:]+rnoise20
	# np.savetxt('Summit_spin_rho_20.csv',rho_out_20,delimiter=',',fmt='%1.4f')


	# s_dt_smb = pd.DataFrame

	# df['smb'].groupby([smb_df.index.year, smb_df.index.month])#.unstack()
	# racmo=pd.DataFrame({'SMB':smb_ind,'TSKIN':tskin_ind})

	# df = pd.DataFrame(
	#          dict(date=pd.date_range('2013-01-01', periods=42, freq='M'),
	#               pb=np.random.rand(42)))

	# df.set_index([df.date.dt.month, df.date.dt.year]).pb.unstack()