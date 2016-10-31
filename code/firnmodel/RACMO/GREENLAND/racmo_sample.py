#! /usr/bin/env python

import sys
import os
import glob

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from lib import malib
from lib import geolib

def find_breaks(smb):
    from scipy.signal import argrelextrema
    import cookb_signalsmooth as sm
    smb_sm = sm.smooth(smb, window_len=7) 
    #Find contiguous locations where smb_sm is below zero
    maxima_idx = argrelextrema(smb_sm, np.greater)[0]
    minima_idx = argrelextrema(smb_sm, np.less)[0]
    return minima_idx

#def density_corr(smb, minima_idx, maxima_idx):
def density_corr(smb, minima_idx):
    rho_s = 0.350
    rho_i = 0.917
    smb_cum = smb.cumsum(axis=0)
    #import ipdb; ipdb.set_trace()
    #maxima = smb_cum[maxima_idx]
    minima = smb_cum[minima_idx]
    #Note: first value is Jan 1, need to assign a lower minima, account for snowfall from previous Sept to Dec
    #Should interpolate what we have from Jan to minima_idx[1]
    minima[0] -= 60
    print minima_idx
    print minima
    snow_idx = []
    ice_idx = []
    for m in range(len(minima_idx)-1): 
        #Note: np.where returns tuple of ndarray
        #Evaluate only for the region from one minima to the next
        snow_idx.append(np.where(smb_cum[minima_idx[m]:minima_idx[m+1]] >= minima[m])[0] + minima_idx[m])
        ice_idx.append(np.where(smb_cum[minima_idx[m]:minima_idx[m+1]] < minima[m])[0] + minima_idx[m])
        #smb_cum_scale[smb_cum[minima_idx[m]:minima_idx[m+1]] > minima[m]] /= rho_s
        #smb_cum_scale[smb_cum[minima_idx[m]:minima_idx[m+1]] <= minima[m]] /= rho_i
    #The rest should be snow
    m += 1
    snow_idx.append(np.where(smb_cum[minima_idx[m]:] >= minima[m])[0] + minima_idx[m]) 
    snow_idx = np.concatenate(snow_idx)
    ice_idx = np.concatenate(ice_idx)
    print snow_idx.size
    print ice_idx.size
    print smb.size
    smb_scale = np.ma.copy(smb)
    smb_scale[snow_idx] /= rho_s
    smb_scale[ice_idx] /= rho_i
    #Could also do this
    #smb_scale[~snow_idx] /= rho_i
    return smb_scale

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

racmo_dir = '/Volumes/insar5/dshean/Greenland/RACMO/smb_reproj'
fn_list = glob.glob('smb*_[0-9][0-9][0-9].tif')
fn_list.sort()

#Load the stack
stack = malib.DEMStack(fn_list)

coord = [-49.54335,69.12397,-181064.724,-2278605.582]
#This is mapx, mapy
#These can be arrays
map_coord = [coord[2], coord[3]]
px_coord = geolib.mapToPixel(map_coord[0], map_coord[1], stack.gt)
px_coord_int = np.rint(px_coord).astype(int)
#This needs to be int idxy, idxx
smb = stack.ma_stack[:,px_coord_int[1], px_coord_int[0]]
smb_cum = smb.cumsum(axis=0) / 1000.

#map_coords = np.broadcast_arrays(np.arange(stack.ma_stack.shape[0])[:, None], y, x)
#samp = ndimage.map_coordinates(stack.ma_stack, map_coords, order=1)

x = stack.date_list

plt.figure()
plt.plot_date(x, smb, linestyle='-', marker=None)


#Could probably do this with dates
#minima_idx = find_breaks(smb)
#This should be the date of the first snowfall each season
#minima_idx = [252, 650, 975, 1345]
#minima_idx.insert(0, 0)
#smb_scale = density_corr(smb, minima_idx) 
#smb_scale_cum = smb_scale.cumsum(axis=0) / 1000.

fig = plt.figure()
#plt.plot_date(x, smb_cum, color='b', linestyle='-', marker=None)
surf_scale_cum = racmo_surface(smb) / 1000.
ax = fig.add_subplot(111)
plt.plot_date(x, surf_scale_cum, color='b', linestyle='-', marker=None)

plt.ylabel('SMB Relative Elevation (m)')
#plt.hlines([0, 3.07, 6.74, 9.05, 12.21], colors='k', linestyles='dashed')

import matplotlib.dates as mdates
#This is used for interactive display of x-value in plot window 
date_str = '%Y-%m-%d %H:%M'
date_fmt = mdates.DateFormatter(date_str)
months = mdates.MonthLocator() 
months_int = mdates.MonthLocator((4,7,10))  # every n months 
years = mdates.YearLocator()   # every year
years_int = mdates.MonthLocator((1, 7))   # every year
fig.autofmt_xdate()
ax.fmt_xdata = mdates.DateFormatter(date_fmt)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_major_locator(years_int)
#ax.xaxis.set_major_formatter(date_fmt)
ax.fmt_xdata = date_fmt
ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

ax.yaxis.grid(True)

plt.autoscale()
#plt.ylim(-14, 2)

plt.title('RACMO/2 Surface Mass Balance Elevation Change')
fig_fn = os.path.join(racmo_dir, 'smb_rel_elevation_2009-2012.pdf')
plt.savefig(fig_fn)

