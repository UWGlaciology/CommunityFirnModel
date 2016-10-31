#! /usr/bin/env python

#Utility for reading RACMO netcdf files

#Need to mask after interpolation
#Look up alpha shape, concave hull
#CGAL has implementation with python bindings

import sys, os
import numpy as np
import scipy.io 
import netCDF4
import matplotlib.pyplot as plt
from datetime import datetime
import osr

from lib import geolib
from lib import malib
from lib import timelib

#Write out vrt
def writevrt(out_csv):
    out_vrt = os.path.splitext(out_csv)[0]+'.vrt'
    #srs=proj
    #This is ECEF, the coord sys of ASP PC
    #srs='EPSG:4978'
    #WGS84
    srs='EPSG:4326'
    #srs='EPSG:3413'
    f = open(out_vrt, 'w')
    f.write('<OGRVRTDataSource>\n')
    f.write('   <OGRVRTLayer name="%s">\n' % os.path.splitext(out_csv)[0])
    f.write('       <SrcDataSource>%s</SrcDataSource>\n' % out_csv)
    f.write('       <GeometryType>wkbPoint</GeometryType>\n')
    f.write('       <LayerSRS>%s</LayerSRS>\n' % srs)
    f.write('       <GeometryField encoding="PointFromColumns" x="lon" y="lat"/>\n')
    #f.write('       <GeometryField encoding="PointFromColumns" x="x" y="y"/>\n')
    f.write('   </OGRVRTLayer>\n')
    f.write('</OGRVRTDataSource>\n')
    f.close()

def writeout(x, y, var_ma, fn, t_srs=None):
    val = var_ma.compressed()
    #Now interpolate in projected coords
    out, gt = interpgrid(x, y, val)
    out_fn = os.path.splitext(fn)[0]+'.tif'
    malib.writeGTiff(out, out_fn, gt=gt, proj=t_srs.ExportToWkt())

#Convert lat/lon to desired projected coordinate sys
def ll2proj(lon, lat, t_srs):
    s_srs = geolib.wgs_srs 
    #t_srs = osr.SpatialReference()
    #t_srs.ImportFromProj4(proj)
    ct = osr.CoordinateTransformation(s_srs, t_srs)
    coord = np.array(ct.TransformPoints(zip(lon, lat)))
    return coord[:,0], coord[:,1] 

#res was 11000.0 for RACMO, maybe 5500.0 to avoid aliasing
def interpgrid(x, y, z, res=100.0):
    import scipy.interpolate
    #y *= -1
    #Input resolution is roughly 11 km, sample output at 5.5 km
    dx = dy = res 
    x_out = np.arange(x.min(), x.max()+dx, dx)
    y_out = np.arange(y.min(), y.max()+dy, dy)
    #m_x, m_y = np.meshgrid(x_out, y_out)
    #out = scipy.interpolate.griddata(x,y,z,x_out,y_out,method='linear')
    out = scipy.interpolate.griddata((x,y),z,(x_out[None,:],y_out[:,None]), method='cubic', fill_value=np.nan)
    #The -1 here flips the data, gdal data model is positive y down
    out = np.ma.fix_invalid(out)[::-1]
    #dilate
    #gt = [x_out[0], dx, 0, y_out[0], dy, 0]
    #Need the 0.5 here b/c GDAL coords have origin at UL pixel corner
    gt = [x_out[0]-0.5*dx, dx, 0, y_out[-1]+0.5*dy, 0, -dy]
    return out, gt

nc_fn_list = sys.argv[1:]
#Define output projection, EPSG:3413 here
#proj = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '
#This is alternative projection for Greenland
#proj = '+proj=stere +lat_0=90 +lat_ts=72 +lon_0=-37.5 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '
#EPSG:3031 for Antarctica
proj = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '
t_srs = osr.SpatialReference()
t_srs.ImportFromProj4(proj)

for i,nc_fn in enumerate(nc_fn_list):
    print "Processing %i of %i: %s" % (i+1, len(nc_fn_list), nc_fn) 
    nc = scipy.io.netcdf.netcdf_file(nc_fn)
    keys = nc.variables.keys()
    latkey = 'lat'
    lonkey = 'lon'
    latkey = 'latgrid'
    lonkey = 'longrid'
    timekey = 'time'
    #key = 'SMB'
    #latkey = 'lat2d'
    #lonkey = 'lon2d'
    lat = malib.checkma(nc.variables[latkey][:]).compressed()
    lon = malib.checkma(nc.variables[lonkey][:]).compressed()
    keys.remove(latkey)
    keys.remove(lonkey)
    #Compute projected coords for all lat/lon points
    x, y = ll2proj(lon, lat, t_srs)
    timeflag = False
    if timekey in keys:
        timeflag = True
        keys.remove(timekey)
    for key in keys:
        print key
        #key = 'smb'
        #key = 'FirnAir'
        #key = 'zs'
        var = nc.variables[key][:]
        ndv = -999
        if timeflag: 
            #This works for older RACMO2 output - daily, monthly SMB
            #nc_time = netCDF4.num2date(nc.variables['time'][:], nc.variables['time'].units)
            #nc_time = nc_time[-74:-1]
            #This is for decimal years
            print nc.variables[timekey][:]
            nc_time = timelib.np_decyear2dt(nc.variables[timekey][:])
            #Just take last ~5 years
            #nc_time = timelib.np_decyear2dt(nc.variables['time'][-795:-1])
            #Go back to early 2008
            #nc_time = timelib.np_decyear2dt(nc.variables['time'][-1120:-1])
            #nc_time = timelib.np_decyear2dt(nc.variables['time'][-72:-1])
            nc_doy = [int(t.strftime('%j')) for t in nc_time]
            #The 2012 smb netcdf includes surrounding coastlines
            #Output size is different
            for n,t in enumerate(nc_time):
                print "Processing %s" % t
                var_ma = np.ma.masked_equal(var, ndv)
                #lat = np.ma.array(lat, mask=var_ma[0].mask).compressed()
                #lon = np.ma.array(lon, mask=var_ma[0].mask).compressed()
                #csv_fn = '%s_%s_%03i.csv' % (key, datetime.strftime(t, '%Y%m%d'), nc_doy[n])
                csv_fn = '%s_%s_%s_%03i.csv' % (nc_fn, key, datetime.strftime(t, '%Y%m%d'), nc_doy[n])
                #Write out original lat/lon/var points
                np.savetxt(csv_fn, zip(lon, lat, val), fmt='%0.6f,%0.6f,%0.3f', delimiter=',', header='lon,lat,%s' % key,comments='')
                #np.savetxt(csv_fn, zip(x, y, val), fmt='%0.3f,%0.3f,%0.3f', delimiter=',', header='x,y,%s' % key,comments='')
                writevrt(csv_fn)
                writeout(x, y, var_ma[n], csv_fn, t_srs)
        else:
            var_ma = np.ma.masked_equal(var, ndv)
            fn = '%s_%s.tif' % (os.path.splitext(nc_fn)[0], key)
            writeout(x, y, var_ma, fn, t_srs)
