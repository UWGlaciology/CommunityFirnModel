from scipy.io import netcdf
import numpy as np

f=netcdf.netcdf_file('smb_Summit.RACMO2.3_1958-2014.nc','r')
# lat=f.variables['lat'][:]
# lon=f.variables['lon'][:]
# smb=f.variables['smb'][:]
# time=f.variables['time'][:]

f.variables['lat'][:]
f.variables['lon'][:]
f.variables['smb'][:]
f.variables['time'][:]



