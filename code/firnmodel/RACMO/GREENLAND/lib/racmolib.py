#! /usr/bin/env python

import sys
import os
import glob
from datetime import datetime, timedelta

import numpy as np
import gdal

import malib
import warplib
import timelib
import moslib

#This is the correct way to determine surface elevation due to smb
def racmo_surface(smb):
    rho_s = 0.350
    rho_i = 0.917
    surface = 0
    #Note: first value is Jan 1, need to assign a lower minima, account for snowfall from previous Sept to Dec
    #Should interpolate what we have from Jan to minima_idx[1]
    ice_surface = -60 / rho_s
    z = [0]
    for i in smb:
        if i > 0:
            surface += i / rho_s
        else:
            if surface == ice_surface:
                surface += i / rho_i
            else:
                surface += i / rho_s
        z.append(surface)
        if surface <= ice_surface:
            ice_surface = surface
    #Note: the smb is technically the "end of the day" value, so remove starting point
    return np.array(z[1:])

def racmo_lininterp_ts(dt, ds=None):
    dt_list = get_racmo_dt_list()
    closest_ts_idx = timelib.get_closest_dt_idx(dt, dt_list)
    if dt_list[closest_ts_idx] > dt:
        closest_ts_idx -= 1
    dt1 = dt_list[closest_ts_idx]
    dt2 = dt_list[closest_ts_idx+1]
    fn1 = fn_list[closest_ts_idx]
    fn2 = fn_list[closest_ts_idx+1]
    ds1 = gdal.Open(fn1)
    ds2 = gdal.Open(fn1)
    if ds is not None:
        ds1, ds2 = warplib.memwarp_multi([ds1, ds2], extent=ds, res=ds, t_srs=ds)
        ds_a = malib.ds_getma(ds)
    a1 = malib.ds_getma(ds1)
    a2 = malib.ds_getma(ds1)
    m = (a2 - a1)/abs(dt2 - dt1).total_seconds()
    a = a1 + m*abs(dt - dt1).total_seconds()
    #drv = ds1.GetDriver()
    #ds1.GetRasterBand(1).WriteArray(a)
    #return ds1
    out_a = np.ma.array(a, mask=ds_a.mask)
    return out_a 

#Given an input filename, find the corresponding RACMO data
def get_racmo_fn(fn, fn_list):
    dt = timelib.fn_getdatetime(fn)
    ds = gdal.Open(fn)
    dt_list = [timelib.fn_getdatetime(os.path.split(fn)[-1]) for fn in fn_list]
    max_diff = timedelta(days=3.0) 
    out_fn = None
    if dt < (min(dt_list) - max_diff) or dt > (max(dt_list) + max_diff):
        print "Unable to find suitable data for %s" % fn 
    else:
        idx = timelib.get_closest_dt_idx(dt, dt_list)
        out_fn = fn_list[idx]
    return out_fn 

#Preparation for RACMO quieries
def get_racmo_fn_list(var='zs'):
    #topdir = '/Volumes/insar5/dshean/Greenland/RACMO/smb_reproj'
    topdir = '/Volumes/insar5/dshean/Antarctica/RACMO2/DATA/FDM_RACMO2.3'
    if not os.path.exists(topdir):
        topdir = '/nobackup/deshean/amund_ext_dem_final/racmo'
    datadir = os.path.join(topdir, var+'_tif')
    racmo_fn_list = glob.glob(datadir+'/*_[0-9][0-9][0-9].tif')
    racmo_fn_list.sort()
    return racmo_fn_list

#Given an input filename, create RACMO tif for same masked extent
def dem_fn_get_racmo(dem_fn, var='zs', racmo_fn_list=None):
    print dem_fn
    ext = var+'.tif'
    out_fn = moslib.find_existing(dem_fn, ext)
    if not out_fn:
        if not racmo_fn_list:
            racmo_fn_list = get_racmo_fn_list(var)
        dem_ma = malib.fn_getma(dem_fn)
        dem_racmo_fn = get_racmo_fn(dem_fn, racmo_fn_list)
        if dem_racmo_fn:
            out_fn = os.path.splitext(dem_fn)[0]+'_racmo_%s.tif' % var
            #Might want to use intermediate resolution here to avoid artifacts
            racmo_ds = warplib.memwarp_multi_fn([dem_racmo_fn,], res=dem_fn, extent=dem_fn, t_srs=dem_fn, r='cubicspline')[0]
            racmo_ma = np.ma.array(malib.ds_getma(racmo_ds), mask=dem_ma.mask, fill_value=dem_ma.fill_value)
            malib.writeGTiff(racmo_ma.astype(float), out_fn, racmo_ds)
    return out_fn
