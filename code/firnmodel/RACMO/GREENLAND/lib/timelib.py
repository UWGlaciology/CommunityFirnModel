#! /usr/bin/env python

#David Shean
#dshean@gmail.com

#Library with time conversion utilities

import os
import numpy as np
from datetime import datetime, timedelta
import time
import dateutil.parser
import matplotlib.dates

#Get timezone for a given lat/lon
#lon,lat = geolib.get_center(ds, t_srs=geolib.wgs_srs)
def getTimeZone(lat, lon):
    import urllib2
    import xml.etree.ElementTree as ET
    #http://api.askgeo.com/v1/918/aa8292ec06199d1207ccc15be3180213c984832707f0cbf3d3859db279b4b324/query.xml?points=37.78%2C-122.42%3B40.71%2C-74.01&databases=Point%2CTimeZone%2CAstronomy%2CNaturalEarthCountry%2CUsState2010%2CUsCounty2010%2CUsCountySubdivision2010%2CUsTract2010%2CUsBlockGroup2010%2CUsPlace2010%2CUsZcta2010
    req = "http://api.askgeo.com/v1/918/aa8292ec06199d1207ccc15be3180213c984832707f0cbf3d3859db279b4b324/query.xml?points="+str(lat)+"%2C"+str(lon)+"&databases=TimeZone"
    opener = urllib2.build_opener()
    f = opener.open(req)
    tree = ET.parse(f)
    root = tree.getroot()
    #Check response
    tzid = None
    if root.attrib['code'] == '0':
        tz = list(root.iter('TimeZone'))[0]
        shortname = tz.attrib['ShortName']
        tzid = tz.attrib['TimeZoneId']
    return tzid 

#Return local timezone time
def getLocalTime(utc_dt, tz):
    import pytz
    local_tz = pytz.timezone(tz)
    local_dt = utc_dt.replace(tzinfo=pytz.utc).astimezone(local_tz)
    return local_dt

#Compute local time based on longitude
def ul_time(utc_dt, lon):
    #return utc_dt + timedelta(hours=lon / np.pi * 12)
    offset = timedelta(hours=(lon*(24.0/360)))
    return utc_dt + offset

#Compute local solar time
#This should be relevant for illumination 
def solarTime(utc_dt, lat, lon):
    import ephem
    o = ephem.Observer()
    o.date = utc_dt
    o.lat = str(lat)
    o.lon = str(lon)
    sun = ephem.Sun()
    sun.compute(o)
    hour_angle = o.sidereal_time() - sun.ra
    rad = str(ephem.hours(hour_angle + ephem.hours('12:00')).norm)
    t = datetime.strptime(rad, '%H:%M:%S.%f')
    solar_dt = datetime.combine(utc_dt.date(), t.time()) 
    return solar_dt 

#Note: this returns current date if not found 
#If only year is provided, will return current month, day
def strptime_fuzzy(s):
    dt = dateutil.parser.parse(str(s), fuzzy=True) 
    return dt 

def fn_getdatetime(fn):
    dt_list = fn_getdatetime_list(fn)
    if dt_list:
        return dt_list[0]
    else:
        return None

#Return datetime object extracted from arbitrary filename
def fn_getdatetime_list(fn):
    #Want to split last component
    fn = os.path.split(os.path.splitext(fn)[0])[-1]
    import re
    #WV01_12JUN152223255-P1BS_R1C1-102001001B3B9800__WV01_12JUN152224050-P1BS_R1C1-102001001C555C00-DEM_4x.tif
    #Need to parse above with month name 
    #Note: made this more restrictive to avoid false matches:
    #'20130304_1510_1030010020770600_1030010020CEAB00-DEM_4x'
    #This is a problem, b/c 2015/17/00:
    #WV02_20130315_10300100207D5600_1030010020151700
    #This code should be obsolete before 2019 
    #Assume new filenames
    #fn = fn[0:13]
    #Use cascading re find to pull out timestamps
    #Note: Want to be less restrictive here - could have a mix of YYYYMMDD_HHMM, YYYYMMDD and YYYY in filename
    #Should probably search for all possibilities, then prune
    dstr = None
    dstr = re.findall(r'(?:^|_)(?:19|20)[0-9][0-9](?:0[1-9]|1[012])(?:0[1-9]|[12][0-9]|3[01])[_T](?:0[0-9]|1[0-9]|2[0-3])[0-5][0-9]', fn)
    if not dstr:
        dstr = re.findall(r'(?:^|_)(?:19|20)[0-9][0-9](?:0[1-9]|1[012])(?:0[1-9]|[12][0-9]|3[01])_', fn)
    if not dstr:
        dstr = re.findall(r'(?:^|_)(?:19|20)[0-9][0-9]_', fn)
    #This is a hack to remove leading _
    dstr = [d.lstrip('_').rstrip('_') for d in dstr]
    #This returns an empty list of nothing is found
    return [strptime_fuzzy(s) for s in dstr]

def get_closest_dt_fn(fn, fn_list):
    dt = fn_getdatetime(fn)
    dt_list = np.array([fn_getdatetime(fn) for fn in fn_list])
    idx = get_closest_dt(dt, dt_list) 
    return fn_list[idx]

def get_closest_dt_idx(dt, dt_list):
    import malib
    dt_list = malib.checkma(dt_list, fix=False)
    dt_diff = np.abs(dt - dt_list)
    return dt_diff.argmin()

#Return center date between two datetime
#Useful for velocity maps
def center_date(dt1, dt2):
    return dt1 + (dt2 - dt1)/2

#Seconds since epoch
def sinceEpoch(dt):
    return time.mktime(dt.timetuple())

def dt2decyear(dt):
    year = dt.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)
    yearElapsed = sinceEpoch(dt) - sinceEpoch(startOfThisYear)
    yearDuration = sinceEpoch(startOfNextYear) - sinceEpoch(startOfThisYear)
    fraction = yearElapsed/yearDuration
    return year + fraction 

def decyear2dt(t):
    year = int(t)
    rem = t - year 
    base = datetime(year, 1, 1)
    dt = base + timedelta(seconds=(base.replace(year=base.year+1) - base).total_seconds() * rem)
    #This works for np array input
    #year = t.astype(int)
    #rem = t - year 
    #base = np.array([datetime(y, 1, 1) for y in year])
    return dt

#Better to use astro libe or jdcal for julian to gregorian conversions 
#Source: http://code.activestate.com/recipes/117215/
def dt2jd(dt):
    """Returns the Julian day number of a date."""
    a = (14 - dt.month)//12
    y = dt.year + 4800 - a
    m = dt.month + 12*a - 3
    return dt.day + ((153*m + 2)//5) + 365*y + y//4 - y//100 + y//400 - 32045

def jd2dt(n):
    """Returns a date corresponding to the given Julian day number."""
    n = int(round(float(n)))
    a = n + 32044
    b = (4*a + 3)//146097
    c = a - (146097*b)//4
    d = (4*c + 3)//1461
    e = c - (1461*d)//4
    m = (5*e + 2)//153
    day = e + 1 - (153*m + 2)//5
    month = m + 3 - 12*(m//10)
    year = 100*b + d - 4800 + m/10
    return datetime(year, month, day)

#This has not been tested 
def gps2dt(gps_week, gps_ms):
    gps_epoch = datetime(1980,1,6,0,0,0)
    gps_week_s = timedelta(seconds=gps_week*7*24*60*60)
    gps_ms_s = timedelta(milliseconds=gps_ms) 
    return gps_epoch + gps_week_s + gps_ms_s

#Matlab ordinal to Python datetime
#Need to account for AD 0 and AD 1 discrepancy between the two
def mat2dt(o):
    return o2dt(o) - timedelta(days=366)

#Python datetime to matlab ordinal
def dt2mat(dt):
    return dt2o(dt + timedelta(days=366))

def dt2o(dt):
    #return datetime.toordinal(dt)
    #This works for arrays of dt
    return matplotlib.dates.date2num(dt)

#Need to split ordinal into integer and decimal parts
def o2dt(o):
    #omod = np.modf(o)
    #return datetime.fromordinal(int(omod[1])) + timedelta(days=omod[0])
    return matplotlib.dates.num2date(o)

#Return integer DOY (julian)
def dt2j(dt):
    return int(dt.strftime('%j'))

#Year and day of year to datetime
#Add comment to http://stackoverflow.com/questions/2427555/python-question-year-and-day-of-year-to-date
#ordinal allows for days>365 and decimal days
def j2dt(yr, j):
    return o2dt(dt2o(datetime(int(yr), 1, 1))+j-1)
    #The solution below can't deal with jd>365
    #jmod = np.modf(j)
    #return datetime.strptime(str(yr)+str(int(jmod[1])), '%Y%j') + timedelta(days=jmod[0])

def print_dt(dt):
    return dt.strftime('%Y%m%d_%H%M')

#Generate a new files with time ordinal written to every pixel
#If dt_ref is provided, return time interval in decimal days
#Should add functionality to do relative doy
def gen_ts_fn(fn, dt_ref=None, ma=False):
    import gdal
    import malib
    print "Generating timestamp for: %s" % fn
    fn_ts = os.path.splitext(fn)[0]+'_ts.tif' 
    if not os.path.exists(fn_ts) or dt_ref is not None:
        ds = gdal.Open(fn)
        #Should be ok with float ordinals here
        a = malib.ds_getma(ds)
        ts = fn_getdatetime(fn)
        #Want to check that dt_ref is valid datetime object
        if dt_ref is not None:
            t = ts - dt_ref 
            t = t.total_seconds()/86400.
            fn_ts = os.path.splitext(fn)[0]+'_ts_rel.tif'
        else:
            t = dt2o(ts)
        a[~np.ma.getmaskarray(a)] = t
        #Probably want to be careful about ndv here - could be 0 for rel
        #ndv = 1E20
        ndv = -9999.0 
        a.set_fill_value(ndv)
        malib.writeGTiff(a, fn_ts, ds) 
    if ma:
        return a
    else:
        return fn_ts

#np vectorize form of functions
np_mat2dt = np.vectorize(mat2dt)
np_dt2mat = np.vectorize(dt2mat)
np_dt2o = np.vectorize(dt2o)
np_o2dt = np.vectorize(o2dt)
np_j2dt = np.vectorize(j2dt)
np_decyear2dt = np.vectorize(decyear2dt)
np_dt2decyear = np.vectorize(dt2decyear)
