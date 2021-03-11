#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
2/24/2021

This script takes outputs from a regional climate model (RCM) - e.g. MERRA, 
MAR - for a particular site and puts that data into a pandas dataframe. 

The output can be fed to RCMpkl_to_spin.py to generate a time series to force the CFM

YOU MAY HAVE TO EDIT THIS SCRIPT A LOT TO MAKE IT WORK WITH YOUR FILE STRUCTURE 
AND WHAT CLIMATE FILES YOU HAVE.

And, for now there are little things you need to search out and change manually,
like the reference climate interval. Sorry!

@author: maxstev
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
import time
import xarray as xr
import glob
import hl_analytic as hla


def find_indices(points,lon,lat,tree=None):
    '''
    find the grid point nearest a given coordinate.
    '''
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

    return ii,jj #, [(i,j) for i,j in ind]


def read_netcdfs_merra(files, dim, ii, jj, vv, transform_func=None):
    '''
    Read merra files and concatenate into a pandas dataframe
    '''
    def process_one_path(path):
        with xr.open_dataset(path) as ds:
            # transform_func should do some sort of selection or
            # aggregation
            if transform_func is not None:
                ds = transform_func(ds)
            # load all data from the transformed dataset, to ensure we can
            # use it after closing each original file
            ds = ds[vv][:,ii,jj]
            ds.load()
            return ds
    datasets = [process_one_path(p) for p in files]
    combined = xr.concat(datasets, dim)
    df1 = combined.to_dataframe()
    return (df1.drop(labels=['lon','lat'],axis=1)).sort_index()

def read_netcdfs_mar(files, dim, ii, jj, vv):
    '''
    Read mar files and concatenate into a pandas dataframe
    '''
    def process_one_path(path):
        with xr.open_dataset(path) as ds:
            dsd = {}
            for v in vv:
                # print(v)
                if len(ds[v].dims)==4:
                    dsd[v] = ds[v][:,0,ii,jj].to_dataframe()
                else:
                    dsd[v] = ds[v][:,ii,jj].to_dataframe()
            df_list = [v for k,v in dsd.items()]
            df1 = pd.concat(df_list, axis=1)
            return df1[df1.columns.intersection(vv)]
    datasets = [process_one_path(p) for p in files]
    return (pd.concat(datasets)).sort_index()

def effectiveT(T):
    '''
    The Arrhenius mean temperature.
    '''
    Q   = -1 * 60.0e3
    R   = 8.314
    k   = np.exp(Q/(R*T))
    km  = np.mean(k)
    return Q/(R*np.log(km))

def getClimate(lat_int,lon_int,writer=True,datatype='MERRA',timeres='1D',melt=False,runtype='local',dsource = None):
    '''
    Load data from MERRA or MAR or whatever.
    Put it into a pandas dataframe, called df_CLIM. index must be datetimeindex for 
    resampling.
    df_CLIM can have any number of columns: BDOT, TSKIN, SMELT, RAIN, 
    SUBLIMATION (use capital letters. We use SMELT because melt is a pandas function)
    Hopefully this makes it easy to adapt for the different climate products.
    write df_CLIM into a pickle for future use.

    Reference for Summit, Greenland (my favorite test site):
    lat = 72.57972
    lon = -38.50454

    DYE-2 (my favorite wet test site):
    lat = 66.5
    lon = -46.2

    UNITS FOR MASS FLUXES IN THE DATAFRAMES ARE kg/m^2 PER TIME STEP SIZE IN
    THE DATA FRAME. e.g. if you have hourly data in the dataframe, the units
    for accumulation are kg/m^2/hour - the mass of precip that fell during that 
    time interval.

    Parameters
    ----------
    lat_int: float
        the latitude of the site you want to build a climate history for
    lon_int: float
        the longitude of the site you want to build a climate history for
    writer: boolean
        Whether or not you want to write the pandas dataframe to a pickle
    datatype: string
        The type of RCM data you are using 'MERRA' or 'MAR' for now.
    melt: boolean
        Whether or not to put melt into the pandas dataframe
    Tinterp: 'mean', 'effective', or 'weighted'
        how to resample the temperature; mean is regular mean, 'effective' is 
        Arrhenius mean; 'weighted' is accumulation-weighted mean
    runtype: 'local' or 'remote'
        Allows you easily switch between directory structures if you are testing
        code locally and running on a remote server
    dsource: 'ERA10k', 'ERA6k', or 'NCEP20k'
        MAR has several flavors; choose which one.

    Returns
    -------
    df_CLIM: pandas dataframe
        Dataframe containing the time series of each pertinent variable for 
        the site, pulled from the RCM data. Index is a datetimeindex.


    '''

    if not writer:
        print('Files will not be written!')
    SPY = 365.25*24*3600

    todaystring = date.today().strftime("%Y%m%d")
    # write_out_dir = 'inputdata{}/'.format(todaystring) + datatype + 'input'
    write_out_dir = 'pickle'
    if writer:
        try: 
            os.makedirs(write_out_dir)
        except:
            pass

    if datatype == 'MERRA':
        '''
        smb has dimensions of (time,lat,lon)
        smb has units of kg m^-2 s^-1 per day (because I sum the hourly values to get a value for each day, but do not divde by 24 after that) (pretty sure, at least!)
        temperature has dimensions of (time,lat,lon)
        temperature has units K

        '''

        ### Set directory to find climate files.
        if lat_int < 0: # Antarctica
            if runtype=='local':
                ddir = 'PATH/TO/LOCAL/DATA/MERRA/Antarctica/Hourly'
            elif runtype=='remote':
                ddir = 'PATH/TO/REMOTE/DATA/MERRA/Antarctica/Hourly'
            elif runtype=='differentremote':
                ddir = 'PATH/TO/OTHER/REMOTE/DATA/CFM/MERRA/Antarctica/Hourly'
            
            # Adjust these as you see fit to set the Reference Climate Interval (RCI)
            spin_date_st = 1980 
            spin_date_end = 2019

        else: # Greenland
            if runtype=='local':
                # ddir = 'PATH/TO/LOCAL/DATA/MERRA/Greenland/Hourly'
                ddir = '/Volumes/Samsung_T1/MERRA/Greenland/Hourly'
            elif runtype=='remote':
                ddir = 'PATH/TO/REMOTE/DATA/MERRA/Greenland/Hourly'
            elif runtype == 'differentremote':
                ddir = 'PATH/TO/OTHER/REMOTE/DATA/MERRA/Greenland/Hourly'
            
            # Adjust these as you see fit to set the Reference Climate Interval (RCI)
            spin_date_st = 1980
            spin_date_end = 1995

        # input_datetimes = [dparser.parse((re.search(r'\d{8}',xx)).group()) for xx in ff] # this will extract the dates for each file
        # yy = np.array([float((re.search(r'\d{8}',xx)).group()[0:4]) for xx in glob.glob(ddir+'/TS/*.nc*')])
        # yrs = np.arange(min(yy),max(yy)+1)

        fn_ll = glob.glob(ddir + '/SMB/*.nc*')[0]
        nc_ll = nc.Dataset(fn_ll,'r')
        lat_ll = nc_ll.variables['lat'][:]
        lon_ll = nc_ll.variables['lon'][:]
        ii, lat_val = min(enumerate(lat_ll), key=lambda x: abs(x[1]-lat_int))
        jj, lon_val = min(enumerate(lon_ll), key=lambda x: abs(x[1]-lon_int))
        nc_ll.close()       
        print('lat_val: ', lat_val)
        print('lon_val: ', lon_val)

        if runtype=='local':
            # pickle_folder = '/PUT/PICKLES/HERE/MERRA/IDSpickle/pickle/'
            pickle_folder = 'example_pickle/'
        else:
            pickle_folder = 'pickle/'
        pickle_name = pickle_folder + 'MERRA2_TS_PREC_df_{}_{}.pkl'.format(lat_val,lon_val)
        if not os.path.exists(pickle_folder):
            os.makedirs(pickle_folder)

        if os.path.isfile(pickle_name):
            print('pickle found')
            writer = False
            loadnetcdf = False
            df_CLIM = pd.read_pickle(pickle_name)
            try:
                df_BDOT = pd.DataFrame(df_CLIM['PRECTOT'])
                df_TS = pd.DataFrame(df_CLIM['TS'])
                df_CLIM.rename(columns={'PRECTOT':'BDOT','TS':'TSKIN'},inplace=True)
            except Exception:
                df_BDOT = pd.DataFrame(xx['BDOT'])
                df_TS = pd.DataFrame(xx['TSKIN'])

            if df_CLIM.BDOT.resample('1A').sum().mean()<1:
                df_CLIM.BDOT = df_CLIM.BDOT *3600 #get rid of seconds dimension - MERRA is hourly, so this gives precip per hour.

        else:
            flist_TS = glob.glob(ddir+'/TS/*.nc*')
            df_TS = read_netcdfs_merra(flist_TS, dim='time',ii=ii,jj=jj,vv='TS')
            df_TS.rename(columns={'TS':'TSKIN'},inplace=True)

            flist_SMB = glob.glob(ddir+'/SMB/*.nc*')
            df_BDOT = read_netcdfs_merra(flist_SMB, dim='time',ii=ii,jj=jj,vv='PRECTOT') # [kg m^-2 s^-1]
            df_BDOT = (df_BDOT.rename(columns={'PRECTOT':'BDOT'}))*3600 # [kg m^-2 hour^-1] (this is amount of precip per MERRA time interval)

            df_CLIM = df_TS.join(df_BDOT)
        # ACCVAR = 'PRECTOT'
        # TVAR = 'TS'
             
        # df_MELT = None
        # df_RAIN = None
        ####################
        #### end MERRA #####

    elif datatype == 'MAR':
        spin_date_st = 1980
        spin_date_end = 1995
        print('Using MAR')
        if lat_int < 0:
            print('no Antarctic MAR data')
            sys.exit()            
        else:
            if runtype=='local':
                ddir = '/Volumes/Samsung_T1/MAR311/Greenland/Daily'

        if not dsource:
            dsource = 'ERA10k'
            print('using MAR ', dsource)

        if dsource == 'ERA10k':
            d2 = '/ERA_1958-2019-10km/'
            vv = ['ME','SF','ST2','RF','SU']
        elif dsource == 'ERA6k':
            d2 = '/ERA_1979-2020-6km/'
            vv = ['ME','SF','ST2','RF']
        elif dsource == 'NCEP20k':
            d2 = '/NCEP1_1948-2020_20km/'
            vv = ['ME','SF','ST2','RF','SU']

        pickle_folder = ddir + '/pickles' + d2
        print(pickle_folder)
        if not os.path.exists(pickle_folder):
            os.makedirs(pickle_folder)
        # searchdir = ddir + d2 + '/*.nc'
        flist = glob.glob(ddir + d2 + '*.nc')
        rgr = nc.Dataset(flist[0],'r')
        lat = rgr['LAT'][:,:]
        lon = rgr['LON'][:,:]
        ii,jj = find_indices((lon_int,lat_int),lon,lat)
        lat_val = lat[ii,jj]
        lon_val = lon[ii,jj]
        print('lat_val: ', lat_val)
        print('lon_val: ', lon_val)
        rgr.close()

        PN = pickle_folder + 'MAR_{}_CLIM_df_{}_{}.pkl'.format(dsource,lat_val,lon_val)
        if os.path.isfile(PN):
            df_CLIM = pd.read_pickle(PN)
            print('Pickle found!')
            df_BDOT = pd.DataFrame(df_CLIM.BDOT)
            df_TS = pd.DataFrame(df_CLIM.TSKIN)
        
        # vv = ['ST2','SMB']
        else:
            df_CLIM = (read_netcdfs_mar(flist,'TIME',ii=ii,jj=jj,vv=vv))[str(spin_date_st):]

            if 'SMB' in df_CLIM.columns:
                df_BDOT = pd.DataFrame(df_CLIM['SMB']/1000*917).rename(columns = ['BDOT']) #put into units kg/m^2/day (i.e. per time resolution in the files))
                df_MELT = None
                df_RAIN = None
            else:
                if 'SU' in df_CLIM.columns:
                    df_BDOT = pd.DataFrame(((df_CLIM['SF']-df_CLIM['SU'])/1000*917),columns=['BDOT']) #put into units kg/m^2/day (i.e. per time resolution in the files))
                    df_CLIM['BDOT'] = df_BDOT.BDOT.values
                    df_CLIM.drop(['SF','SU'],axis=1,inplace=True)
                else:
                    df_BDOT = pd.DataFrame((df_CLIM['SF'])/1000*917).rename(columns={'SF':'BDOT'}) #put into units kg/m^2/day (i.e. per time resolution in the files))
                    df_CLIM['BDOT'] = df_BDOT.BDOT.values
                    df_CLIM.drop(['SF'],axis=1,inplace=True)
                df_CLIM['ME'] = df_CLIM['ME']/1000*917 #put into units kg/m^2/day (i.e. per time resolution in the files))
                df_CLIM['RF'] = df_CLIM['RF']/1000*917 #put into units kg/m^2/day (i.e. per time resolution in the files))
                # df_MELT = pd.DataFrame(df_CLIM['ME']/1000*917/3600).rename(columns={'ME':'MELT'}) #put into equivalent units to the merra data (kg/m^2/s)
                # df_RAIN = pd.DataFrame(df_CLIM['RF']/1000*917/3600).rename(columns={'RF':'RAIN'}) #put into equivalent units to the merra data (kg/m^2/s)
            df_TS = pd.DataFrame(df_CLIM['ST2']).rename(columns = {'ST2':'TSKIN'}) + 273.15

            drn = {'ME':'SMELT','SU':'SUBLIMATION','SF':'BDOT','RF':'RAIN','ST2':'TSKIN','SMB':'BDOT'}
            df_CLIM.rename(mapper=drn,axis=1,inplace=True)
            df_CLIM.TSKIN = df_CLIM.TSKIN + 273.15
    ###############
    ### end MAR ###
    ###############


    if writer:
        if datatype =='MERRA':
            df_CLIM.to_pickle(pickle_folder + 'MERRA2_CLIM_df_{}_{}.pkl'.format(lat_val,lon_val))
        elif datatype == 'MAR':
            df_CLIM.to_pickle(pickle_folder + 'MAR_{}_CLIM_df_{}_{}.pkl'.format(dsource,lat_val,lon_val))

    return df_CLIM
    # return CD, stepsperyear, depth_S1, depth_S2, desired_depth


if __name__ == '__main__':

    LLpair = sys.argv[1]
    nn = np.fromstring(LLpair,dtype =float, sep=' ')
    lat_int = nn[0]
    lon_int = nn[1]
    writer=True
    datatype='MERRA'
    runtype = 'local'

    df_CLIM = getClimate(lat_int,lon_int,writer = True, runtype = runtype)





