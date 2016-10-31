#! /usr/bin/env python

#Utility for reading RACMO netcdf files

import sys, os
import numpy as np
import scipy.io 
import netCDF4
import matplotlib.pyplot as plt
from lib import timelib

#See the following for KDTree approach to search arrays:
#http://stackoverflow.com/questions/10818546/finding-index-of-nearest-point-in-numpy-arrays-of-x-and-y-coordinates
def ll2px(p_lat, p_lon):
    #distance = [(lon - x)**2 + (lat - y)**2 for x,y in zip(p_lon, p_lat)]
    x = p_lon
    y = p_lat
    #distance = [(lon - x)**2 + (lat - y)**2 for x,y in (p_lon, p_lat)]
    distance = [(lon - x)**2 + (lat - y)**2]
    idx_a = [np.unravel_index(a.argmin(), lon.shape) for a in distance]
    #return np.array(idx_a).T
    return zip(*idx_a)

#Extract values at indices
def getvals(var, iy, ix):
    return var[:, iy, ix].T

nc_fn = sys.argv[1]
csv_fn = sys.argv[2]

nc = scipy.io.netcdf.netcdf_file(nc_fn)

#keys = nc.variables
smb = nc.variables['smb'][:]
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
smb_time = netCDF4.num2date(nc.variables['time'][:], nc.variables['time'].units)
smb_doy = [t.strftime('%j') for t in smb_time]
#time = nc.variables['time'][:]

#ndv = -9999
#smb_ma = numpy.ma.masked_equal(smb.data, ndv)

#Load csv
#Need the dtype=None here to automatically handle site name string
pts = np.genfromtxt(csv_fn, dtype=None, missing_values=-32767.0, usemask=True, delimiter=',', names=True)
header = list(pts.dtype.names)
sites = [pts[header[0]].tolist()]
#Extract lat lon
lat_idx = 1
lon_idx = 2
p_lat = [pts[header[lat_idx]].tolist()]
p_lon = [pts[header[lon_idx]].tolist()]
#Loop through each station, or do everything with vectors
smb_pts = getvals(smb, *ll2px(p_lat, p_lon))
smb_pts_cum = smb_pts.cumsum(axis=1)

dem_dt = [timelib.fn_getdatetime(s) for s in header[lon_idx+1:]]
dem_doy = np.array([t.strftime('%j') for t in dem_dt])

plt.figure(0, figsize=(10,7.5))
plt.suptitle('RACMO2/GL 2011 Daily Net Surface Mass Balance (SMB)')
plt.title('WHOI/UW West Greenland Lakes GPS sites')
plt.ylabel('SMB (mm w.e.)')
plt.xlabel('DOY 2011')
plt.xlim(0,365)

plt.figure(1, figsize=(10,7.5))
plt.suptitle('RACMO2/GL 2011 Cumulative Surface Mass Balance (SMB)')
plt.title('WHOI/UW West Greenland Lakes GPS sites')
plt.ylabel('SMB (mm w.e.)')
plt.xlabel('DOY 2011')
plt.xlim(0,365)

plt.figure(2, figsize=(10,7.5))
plt.suptitle('WorldView DEM Surface Elevation')
plt.title('WHOI/UW West Greenland Lakes GPS sites')
plt.ylabel('Elevation (m, WGS84)')
plt.xlabel('DOY 2011')
plt.xlim(0,365)

for site, mysmb, mysmb_cum, mylat, mylon in zip(sites, smb_pts, smb_pts_cum, p_lat, p_lon):
    plt.figure(0)
    plt.plot(smb_doy, mysmb, label='%s (%0.5f, %0.5f)' % (site, mylat, mylon))
    plt.figure(1)
    plt.plot(smb_doy, mysmb_cum, label='%s (%0.5f, %0.5f)' % (site, mylat, mylon))

plt.figure(2)
for n, site in enumerate(sites[0:2]):
    dem_pts = np.ma.array(list(pts[n])[lon_idx+1:])
    print dem_pts
    #This scales to same value
    #dem_pts = (dem_pts - dem_pts[0])*1000
    plt.plot(dem_doy[~dem_pts.mask], dem_pts.compressed(), label='%s DEM elevation' % site)

plt.figure(0)
plt.legend(loc=3, prop={'size':12})
#plt.legend(loc=3, prop={'size':10}, ncol=2)
#plt.savefig('racmo_2011_dailysmb_allsites.pdf')

plt.figure(1)
plt.legend(loc=3, prop={'size':12})
#plt.legend(loc=3, prop={'size':10}, ncol=2)
#plt.savefig('racmo_2011_cumsmb_allsites.pdf')

plt.figure(2)
plt.legend(loc=3, prop={'size':12})
plt.savefig('racmo_2011_DEMelev.pdf')

plt.show()
