#! /usr/bin/env python

"""
David Shean
dshean@gmail.com
10/29/12

Library of various useful raster geospatial functions
"""

#Need to make sure all geom have spatial reference included

import numpy as np
import gdal
import ogr
import osr

#Define WGS84 srs
wgs_srs = osr.SpatialReference()
wgs_srs.SetWellKnownGeogCS('WGS84')

#Define ECEF srs
ecef_srs=osr.SpatialReference()
ecef_srs.ImportFromEPSG(4978)

#Define ITRF2008 srs
itrf_srs=osr.SpatialReference()
itrf_srs.ImportFromEPSG(5332)

#Define EGM96 srs
#Note: must have gtx grid files in /usr/local/share/proj
#cd /usr/local/share/proj
#wget http://download.osgeo.org/proj/vdatum/egm96_15/egm96_15.gtx
#wget http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx
#See: http://lists.osgeo.org/pipermail/gdal-dev/2011-August/029856.html
egm96_srs=osr.SpatialReference()
egm96_srs.ImportFromProj4("+proj=longlat +datum=WGS84 +no_defs +geoidgrids=egm96_15.gtx")

#Define EGM2008 srs
egm08_srs=osr.SpatialReference()
egm08_srs.ImportFromProj4("+proj=longlat +datum=WGS84 +no_defs +geoidgrids=egm08_25.gtx")

#Define N Polar Stereographic srs
nps_srs=osr.SpatialReference()
nps_srs.ImportFromEPSG(3413)

#Define N Polar Stereographic srs
sps_srs=osr.SpatialReference()
sps_srs.ImportFromEPSG(3031)

#To do for transformations below:
#Check input order of lon, lat
#Need to broadcast z=0.0 if z is not specified
#Check that all inputs have same length

def cT_helper(x, y, z, in_srs, out_srs):
    x, y, z = np.atleast_1d(x), np.atleast_1d(y), np.atleast_1d(z)
    #Handle cases where z is 0 - probably a better way to use broadcasting for this
    if x.shape[0] != z.shape[0]:
        #Watch out for masked array input here
        orig_z = z[0]
        z = np.zeros_like(x)
        z[:] = orig_z
    orig_shape = x.shape
    cT = osr.CoordinateTransformation(in_srs, out_srs)
    #x2, y2, z2 = zip(*[cT.TransformPoint(*xyz) for xyz in zip(x, y, z)])
    x2, y2, z2 = zip(*[cT.TransformPoint(*xyz) for xyz in zip(np.ravel(x),np.ravel(y),np.ravel(z))])
    if len(x2) == 1:
        x2, y2, z2 = x2[0], y2[0], z2[0] 
    else:
        x2 = np.array(x2).reshape(orig_shape)
        y2 = np.array(y2).reshape(orig_shape)
        z2 = np.array(z2).reshape(orig_shape)
    return x2, y2, z2

def ll2ecef(lon, lat, z=0.0):
    return cT_helper(lon, lat, z, wgs_srs, ecef_srs)
    
def ecef2ll(x, y, z):
    return cT_helper(x, y, z, ecef_srs, wgs_srs)

def ll2itrf(lon, lat, z=0.0):
    return cT_helper(lon, lat, z, wgs_srs, itrf_srs)

def itrf2ll(x, y, z):
    return cT_helper(x, y, z, itrf_srs, wgs_srs)

#Note: the lat/lon values returned here might be susceptible to rounding errors 
#Or are these true offsets due to dz?
#120.0 -> 119.99999999999999
#def geoid2ell(lon, lat, z=0.0, geoid=egm08_srs):
def geoid2ell(lon, lat, z=0.0, geoid=egm96_srs):
    llz = cT_helper(lon, lat, z, geoid, wgs_srs)
    return lon, lat, llz[2]

def ell2geoid(lon, lat, z=0.0, geoid=egm96_srs):
    llz = cT_helper(lon, lat, z, wgs_srs, geoid)
    return lon, lat, llz[2]

def ll2nps(lon, lat, z=0.0):
    #Should throw error here
    if np.any(lat < 0.0):
        print "Warning: latitude out of range for output projection"
    return cT_helper(lon, lat, z, wgs_srs, nps_srs)

def nps2ll(x, y, z=0.0):
    return cT_helper(x, y, z, nps_srs, wgs_srs)

def ll2sps(lon, lat, z=0.0):
    if np.any(lat > 0.0):
        print "Warning: latitude out of range for output projection"
    return cT_helper(lon, lat, z, wgs_srs, sps_srs)

def sps2ll(x, y, z=0.0):
    return cT_helper(x, y, z, sps_srs, wgs_srs)

def scale_ps_ds(ds):
    clat, clon = get_center(ds)
    return scale_ps(clat)

def scale_ps(lat):
    """
    From Ben Smith email on 7/12/12: PS scale m file
    
    This function calculates the scaling factor for a polar stereographic
    projection (ie. SSM/I grid) to correct area calculations. The scaling
    factor is defined (from Snyder, 1982, Map Projections used by the U.S.
    Geological Survey) as:

    k = (mc/m)*(t/tc), where:

    m = cos(lat)/sqrt(1 - e2*sin(lat)^2)
    t = tan(Pi/4 - lat/2)/((1 - e*sin(lat))/(1 + e*sin(lat)))^(e/2)
    e2 = 0.006693883 is the earth eccentricity (Hughes ellipsoid)
    e = sqrt(e2)
    mc = m at the reference latitude (70 degrees)
    tc = t at the reference latitude (70 degrees)

    The ratio mc/tc is precalculated and stored in the variable m70_t70.

    """
    #Note: want to adapt for np arrays
    if lat > 0:
        m70_t70 = 1.9332279 
    else:
        # for 71 deg, southern PS  -- checked BS 5/2012
        m70_t70 = 1.93903005  

    #for WGS84, a=6378137, 1/f = 298.257223563 -> 1-sqrt(1-e^2) = f
    #-> 1-(1-f)^2 = e2 =    0.006694379990141
    #e2 = 0.006693883
    e2 = 0.006694379990141  # BS calculated from WGS84 parameters 5/2012
    e = np.sqrt(e2)

    #Hack to deal with pole
    if abs(lat) >= 90.0:
        lat = 89.999999999

    lat = abs(np.deg2rad(lat))
    slat = np.sin(lat)
    clat = np.cos(lat)

    m = clat/np.sqrt(1. - e2*slat**2)
    t = np.tan(np.pi/4 - lat/2)/((1. - e*slat)/(1. + e*slat))**(e/2)
    k = m70_t70*t/m

    scale=(1./k)
    return scale

def lon360to180(lon):
    lon[lon > 180.0] -= 360.0
    return lon

def lon180to360(lon):
    lon[lon < 0.0] += 360.0
    return lon

#Want to accept np arrays for these
def dd2dms(dd):
    n = dd < 0
    dd = abs(dd)
    m,s = divmod(dd*3600,60)
    d,m = divmod(m,60)
    if n:
        d = -d
    return d,m,s

def dms2dd(d,m,s):
    if d < 0:
        sign = -1
    else:
        sign = 1
    dd = sign * (int(abs(d)) + float(m) / 60 + float(s) / 3600)
    return dd

#Note: this needs some work, not sure what input str format was supposed to be
def dms2dd_str(dms_str):
    import re
    dms_str = re.sub(r'\s', '', dms_str)
    if re.match('[swSW]', dms_str):
        sign = -1
    else:
        sign = 1
    (degree, minute, second, frac_seconds, junk) = re.split('\D+', dms_str, maxsplit=4)
    #dd = sign * (int(degree) + float(minute) / 60 + float(second) / 3600 + float(frac_seconds) / 36000)
    dd = dms2dd(degree*sign, minute, second+frac_seconds) 
    return dd

#Note: These should work with input np arrays
#Note: these functions are likely in osr/pyproj
#GDAL model used here - upper left corner of upper left pixel for mX, mY (and in GeoTransform)
def mapToPixel(mX, mY, geoTransform):
    mX = np.asarray(mX)
    mY = np.asarray(mY)
    if geoTransform[2] + geoTransform[4] == 0:
        pX = ((mX - geoTransform[0]) / geoTransform[1]) - 0.5
        pY = ((mY - geoTransform[3]) / geoTransform[5]) - 0.5
        #pX = (mX - geoTransform[0]) / geoTransform[1]
        #pY = (mY - geoTransform[3]) / geoTransform[5]
    else:
        pX, pY = applyGeoTransform(mX, mY, invertGeoTransform(geoTransform))
    #return int(pX), int(pY)
    return pX, pY

def pixelToMap(pX, pY, geoTransform):
    pX = np.asarray(pX, dtype=float)
    pY = np.asarray(pY, dtype=float)
    pX += 0.5
    pY += 0.5
    mX, mY = applyGeoTransform(pX, pY, geoTransform)
    return mX, mY

#Should probably keep this clean and deal with 0.5 px offsets in pixelToMap
def applyGeoTransform(inX, inY, geoTransform):
    inX = np.asarray(inX)
    inY = np.asarray(inY)
    #inX = np.array(inX) + 0.5
    #inY = np.array(inY) + 0.5
    outX = geoTransform[0] + inX * geoTransform[1] + inY * geoTransform[2]
    outY = geoTransform[3] + inX * geoTransform[4] + inY * geoTransform[5]
    return outX, outY

def invertGeoTransform(geoTransform):
    # we assume a 3rd row that is [1 0 0]
    # compute determinate
    det = geoTransform[1] * geoTransform[5] - geoTransform[2] * geoTransform[4]
    if abs(det) < 0.000000000000001:
        return
    invDet = 1.0 / det
    # compute adjoint and divide by determinate
    outGeoTransform = [0, 0, 0, 0, 0, 0]
    outGeoTransform[1] = geoTransform[5] * invDet
    outGeoTransform[4] = -geoTransform[4] * invDet
    outGeoTransform[2] = -geoTransform[2] * invDet
    outGeoTransfrom[5] = geoTransform[1] * invDet
    outGeoTransform[0] = (geoTransform[2] * geoTransform[3] - geoTransform[0] * geoTransform[5]) * invDet
    outGeoTransform[3] = (-geoTransform[1] * geoTransform[3] + geoTransform[0] * geoTransform[4]) * invDet
    return outGeoTransform

def block_stats(x,y,z,ds,stat='median'):
    import scipy.stats as stats
    extent = get_extent(ds)
    #[[xmin, xmax], [ymin, ymax]]
    range = [[extent[0], extent[2]], [extent[1], extent[3]]]
    #bins = (ns, nl)
    bins = (ds.RasterXSize, ds.RasterYSize)
    if stat == 'max':
        stat = np.max
    elif stat == 'min':
        stat = np.min
    #block_count, xedges, yedges, bin = stats.binned_statistic_2d(x,y,z,'count',bins,range)
    block_stat, xedges, yedges, bin = stats.binned_statistic_2d(x,y,z,stat,bins,range)
    #Get valid blocks
    #if (stat == 'median') or (stat == 'mean'):
    if stat in ('median', 'mean', np.max, np.min):
        idx = ~np.isnan(block_stat)
    else:
        idx = (block_stat != 0)
    idx_idx = idx.nonzero()
    #Cell centers
    res = [(xedges[1] - xedges[0]), (yedges[1] - yedges[0])]
    out_x = xedges[:-1]+res[0]/2.0
    out_y = yedges[:-1]+res[1]/2.0
    out_x = out_x[idx_idx[0]]
    out_y = out_y[idx_idx[1]]
    out_z = block_stat[idx]
    return out_x, out_y, out_z

#Note: the above method returns block_stat, which is already a continuous grid
#Just need to account for ndv and the upper left x_edge and y_edge

def block_stats_grid(x,y,z,ds,stat='median'):
    mx, my, mz = block_stats(x,y,z,ds,stat)
    gt = ds.GetGeoTransform()
    pX, pY = mapToPixel(mx, my, gt)
    shape = (ds.RasterYSize, ds.RasterXSize)
    ndv = -9999.0
    a = np.full(shape, ndv)
    a[pY.astype('int'), pX.astype('int')] = mz
    return np.ma.masked_equal(a, ndv) 

def block_stats_gen(x,y,z,stat='median',bins=None,res=None):
    #res should be (x_res, y_res)
    #bins should be (x_bins, y_bins)
    range = np.array([[x.min(), x.max()],[y.min(), y.max()]])
    #dim = np.array([int(round(range[0][1] - range[0][0])), int(round(range[1][1] - range[1][0]))])
    dim = np.array([int(np.ceil(range[0][1] - range[0][0])), int(np.ceil(range[1][1] - range[1][0]))])
    if (res is not None) and (bins is None):
        res = np.array(res)
        bins = (np.ceil(dim/res.astype('float'))).astype(int)
    if (bins is not None) and (res is None):
        bins = np.array(bins)
        res = dim/(bins).astype('float')
    #Note: bins is (nx, ny)
    block_stat, xedges, yedges, bin = stats.binned_statistic_2d(x,y,z,stat,bins,range)
    #gt should be upper left pixel coords
    gt = (range[0][0], res[1], 0.0, range[1][1], 0.0, -res[0])

#Modify proj/gt of dst_fn in place
def copyproj(src_fn, dst_fn, gt=True):
    src_ds = gdal.Open(src_fn, gdal.GA_ReadOnly)
    dst_ds = gdal.Open(dst_fn, gdal.GA_Update)
    dst_ds.SetProjection(src_ds.GetProjection())
    if gt:
        dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    src_ds = None
    dst_ds = None

#This should be a function of a new geom class
#Assumes geom has srs defined
#Modifies geom in place
def geom_transform(geom, t_srs):
    s_srs = geom.GetSpatialReference()
    if not s_srs.IsSame(t_srs):
        ct = osr.CoordinateTransformation(s_srs, t_srs)
        geom.Transform(ct)
        geom.AssignSpatialReference(t_srs)

#Return dataset bbox envelope as geom
def get_geom(ds, t_srs=None):
    gt = ds.GetGeoTransform()
    ds_srs = get_ds_srs(ds)
    if t_srs is None:
        t_srs = ds_srs
    ns = ds.RasterXSize
    nl = ds.RasterYSize
    x = np.array([0, ns, ns, 0, 0])
    y = np.array([0, 0, nl, nl, 0])
    mx, my = pixelToMap(x, y, gt)
    geom_wkt = 'POLYGON(({0}))'.format(', '.join(['{0} {1}'.format(*a) for a in zip(mx,my)]))
    geom = ogr.CreateGeometryFromWkt(geom_wkt)
    geom.AssignSpatialReference(ds_srs)
    if not ds_srs.IsSame(t_srs):
        geom_transform(geom, t_srs)
    return geom

#Get geom from shapefile
#Need to handle multi-part geom
#http://osgeo-org.1560.x6.nabble.com/Multipart-to-singlepart-td3746767.html
def shp2geom(shp_fn):
    ds = ogr.Open(shp_fn)
    lyr = ds.GetLayer()
    srs = lyr.GetSpatialRef()
    lyr.ResetReading()
    geom_list = []
    for feat in lyr:
        geom = feat.GetGeometryRef()
        #Duplicate the geometry, or segfault
        #See: http://trac.osgeo.org/gdal/wiki/PythonGotchas
        geom_dup = ogr.CreateGeometryFromWkt(geom.ExportToWkt())
        geom_dup.AssignSpatialReference(srs)
        geom_list.append(geom_dup)
    #geom = ogr.ForceToPolygon(' '.join(geom_list))    
    #Dissolve should convert multipolygon to single polygon 
    #return geom_list[0]
    ds = None
    return geom_list

#Write out a shapefile for input geometry
#Useful for debugging
def geom2shp(geom, out_fn, fields=False):
    import os
    import timelib
    driverName = "ESRI Shapefile"
    drv = ogr.GetDriverByName(driverName)
    if os.path.exists(out_fn):
        drv.DeleteDataSource(out_fn)
    out_ds = drv.CreateDataSource(out_fn)
    out_lyrname = os.path.splitext(os.path.split(out_fn)[1])[0]
    geom_srs = geom.GetSpatialReference()
    geom_type = geom.GetGeometryType()
    out_lyr = out_ds.CreateLayer(out_lyrname, geom_srs, geom_type)
    if fields:
        field_defn = ogr.FieldDefn("name", ogr.OFTString)
        field_defn.SetWidth(128)
        out_lyr.CreateField(field_defn)
        #field_defn = ogr.FieldDefn("date", ogr.OFTString)
        field_defn = ogr.FieldDefn("date", ogr.OFTInteger)
        field_defn.SetWidth(32)
        out_lyr.CreateField(field_defn)
    out_feat = ogr.Feature(out_lyr.GetLayerDefn())
    out_feat.SetGeometry(geom)
    if fields:
        out_feat_name = os.path.splitext(out_fn)[0]
        #out_feat_date = str(timelib.fn_getdatetime(out_fn))
        #out_feat_date = int(timelib.fn_getdatetime(out_fn).strftime('%Y%m%d%H%M'))
        out_feat_date = int(timelib.fn_getdatetime(out_fn).strftime('%Y%m%d'))
        out_feat.SetField("name", out_feat_name)
        out_feat.SetField("date", out_feat_date)
    out_lyr.CreateFeature(out_feat)
    out_ds = None
    #return status?

#get_outline is an attempt to reproduce the PostGIS Raster ST_MinConvexHull function
#Could potentially do the following:
#Extract random pts from unmasked elements, get indices
#Run scipy convex hull
#Convert hull indices to mapped coords

#See this:
#http://stackoverflow.com/questions/3654289/scipy-create-2d-polygon-mask

#This generates a wkt polygon outline of valid data for the input raster
#Need to implement geoma, or take ma as optional argument don't want to load again
def get_outline(ds, t_srs=None, scale=1.0, simplify=False, convex=False):
    import malib
    gt = np.array(ds.GetGeoTransform())
    
    #Want to limit the dimensions of a, as notmasked_edges is slow
    a = malib.gdal_getma_sub(ds, scale=scale)

    #Create empty geometry
    geom = ogr.Geometry(ogr.wkbPolygon)
    #Check to make sure we have unmasked data
    if a.count() != 0:
        #Scale the gt for reduced resolution
        #The UL coords should remain the same, as any rounding will trim LR
        if (scale != 1.0):
            gt[1] *= scale
            gt[5] *= scale
        #Get srs
        ds_srs = get_ds_srs(ds)
        if t_srs is None:
            t_srs = ds_srs
        #Find the unmasked edges
        #Note: using only axis=0 from notmasked_edges will miss undercuts - see malib.get_edgemask
        #Better ways to do this - binary mask, sum (see numpy2stl)
        #edges0, edges1, edges = malib.get_edges(a)
        px = np.ma.notmasked_edges(a, axis=0)
        coord = []
        #Combine edge arrays, reversing order and adding first point to complete polygon
        x = np.concatenate((px[0][1][::1], px[1][1][::-1], [px[0][1][0]]))
        #x = np.concatenate((edges[0][1][::1], edges[1][1][::-1], [edges[0][1][0]]))
        y = np.concatenate((px[0][0][::1], px[1][0][::-1], [px[0][0][0]]))
        #y = np.concatenate((edges[0][0][::1], edges[1][0][::-1], [edges[0][0][0]]))
        #Use np arrays for computing mapped coords
        mx, my = pixelToMap(x, y, gt)
        #Create wkt string
        geom_wkt = 'POLYGON(({0}))'.format(', '.join(['{0} {1}'.format(*a) for a in zip(mx,my)]))
        geom = ogr.CreateGeometryFromWkt(geom_wkt)
        if not ds_srs.IsSame(t_srs):
            ct = osr.CoordinateTransformation(ds_srs, t_srs)
            geom.Transform(ct)
        #Make sure geometry has correct srs assigned
        geom.AssignSpatialReference(t_srs)
        if not geom.IsValid():
            tol = gt[1] * 0.1
            geom = geom.Simplify(tol)
        #Need to get output units and extent for tolerance specification
        if simplify:
            #2 pixel tolerance
            tol = gt[1] * 2
            geom = geom.Simplify(tol)
        if convex:
            geom = geom.ConvexHull()
    else:
        print "No unmasked values found"
    return geom

#Given an input line geom, generate points at fixed interval
def line2pts(geom, dl=None):
    #Extract list of (x,y) tuples at nodes
    nodes = geom.GetPoints()
    #print "%i nodes" % len(nodes)
   
    #Point spacing in map units
    if dl is None:
        nsteps=1000
        dl = geom.Length()/nsteps

    #This only works for equidistant projection!
    #l = np.arange(0, geom.Length(), dl)

    #Initialize empty lists
    l = []
    mX = []
    mY = []

    #Add first point to output lists
    l += [0]
    x = nodes[0][0]
    y = nodes[0][1]
    mX += [x]
    mY += [y]

    #Remainder
    rem_l = 0
    #Previous length (initially 0)
    last_l = l[-1]
    
    #Loop through each line segment in the feature
    for i in range(0,len(nodes)-1):
        x1, y1 = nodes[i]
        x2, y2 = nodes[i+1]
      
        #Total length of segment
        tl = np.sqrt((x2-x1)**2 + (y2-y1)**2)

        #Number of dl steps we can fit in this segment
        #This returns floor 
        steps = int((tl+rem_l)/dl)

        if steps > 0:
            dx = ((x2-x1)/tl)*dl
            dy = ((y2-y1)/tl)*dl
            rem_x = rem_l*(dx/dl)
            rem_y = rem_l*(dy/dl)
            
            #Loop through each step and append to lists
            for n in range(1, steps+1):
                l += [last_l + (dl*n)]
                #Remove the existing remainder
                x = x1 + (dx*n) - rem_x
                y = y1 + (dy*n) - rem_y
                mX += [x]
                mY += [y]

            #Note: could just build up arrays of pX, pY for entire line, then do single z extraction
            #Update the remainder
            rem_l += tl - (steps * dl)
            last_l = l[-1]
        else:
            rem_l += tl 

    return l, mX, mY 
    
#Return resolution stats for an input dataset list
def get_res_stats(ds_list, t_srs=None):
    if t_srs is None:
        t_srs = get_ds_srs(ds_list[0]) 
    res = np.array([get_res(ds, t_srs=t_srs) for ds in ds_list])
    #Check that all projections are identical
    #gt_array = np.array([ds.GetGeoTransform() for ds in args])
    #xres = gt_array[:,1]
    #yres = -gt_array[:,5]
    #if xres == yres:
    #res = np.concatenate((xres, yres))
    min = np.min(res)
    max = np.max(res)
    mean = np.mean(res)
    med = np.median(res)
    return (min, max, mean, med)

#Get resolution of dataset in specified coordinate system
#mpd = 111319.9
def get_res(ds, t_srs=None, square=False):
    gt = ds.GetGeoTransform()
    ds_srs = get_ds_srs(ds)
    #This is Xres, Yres
    res = [gt[1], np.abs(gt[5])]
    if square:
        res = [np.mean(res), np.mean(res)]
    if t_srs is not None and not ds_srs.IsSame(t_srs):
        if True:
            #This diagonal approach is similar to the approach in gdaltransformer.cpp
            #Bad news for large extents near the poles
            ullr = get_ullr(ds, t_srs)
            diag = np.sqrt((ullr[0]-ullr[2])**2 + (ullr[1]-ullr[3])**2)
            #extent = get_extent(ds, t_srs)
            #diag = np.sqrt((extent[2]-extent[0])**2 + (extent[3]-extent[1])**2)
            res = diag / np.sqrt(ds.RasterXSize**2 + ds.RasterYSize**2)
            res = [res, res]
        else:        
            #Compute from center pixel
            ct = osr.CoordinateTransformation(ds_srs, t_srs)
            pt = get_center(ds)
            #Transform center coordinates
            pt_ct = ct.TransformPoint(*pt)
            #Transform center + single pixel offset coordinates
            pt_ct_plus = ct.TransformPoint(pt[0] + gt[1], pt[1] + gt[5])
            #Compute resolution in new units 
            res = [pt_ct_plus[0] - pt_ct[0], np.abs(pt_ct_plus[1] - pt_ct[1])]
    return res

#Return center coordinates
def get_center(ds, t_srs=None):
    gt = ds.GetGeoTransform()
    ds_srs = get_ds_srs(ds)
    center = [gt[0] + (gt[1] * ds.RasterXSize/2.0), gt[3] + (gt[5] * ds.RasterYSize/2.0)]
    #include t_srs.Validate() and t_srs.Fixup()
    if t_srs is not None and not ds_srs.IsSame(t_srs):
        ct = osr.CoordinateTransformation(ds_srs, t_srs)
        center = list(ct.TransformPoint(*center)[0:2])
    return center

#Return ullr extent of dataset
#ul_x, ul_y, lr_x, lr_y
#If t_srs is specified, output will be converted to specified srs
def get_ullr(ds, t_srs=None):
    gt = ds.GetGeoTransform()
    ds_srs = get_ds_srs(ds) 
    ul = [gt[0], gt[3]]
    lr = [gt[0] + (gt[1] * ds.RasterXSize), gt[3] + (gt[5] * ds.RasterYSize)]
    if t_srs is not None and not ds_srs.IsSame(t_srs):
        ct = osr.CoordinateTransformation(ds_srs, t_srs)
        #Check to see if ct creation failed
        #if ct == NULL:
        #Check to see if transform failed
        #if not ct.TransformPoint(ullr[0], ullr[1]):
        #Need to check that transformed coordinates fall within appropriate bounds
        #Note: egm96-5 geoid bounds are 
        #[-0.4583333333333333, 90.45833333333333]
        #[362.0416666666667, -90.45833333333333]
        #When transforming to EPSG:3413, get error about lat<-90
        ul = ct.TransformPoint(*ul)[0:2]
        lr = ct.TransformPoint(*lr)[0:2]
    return list(ul + lr)

#Get srs object for input dataset
def get_ds_srs(ds):
    ds_srs = osr.SpatialReference()
    ds_srs.ImportFromWkt(ds.GetProjectionRef())
    return ds_srs

#Return min/max extent of dataset
#xmin, xmax, ymin, ymax
#If t_srs is specified, output will be converted to specified srs
def get_extent(ds, t_srs=None):
    gt = ds.GetGeoTransform()
    ds_srs = get_ds_srs(ds) 
    ul = [gt[0], gt[3]]
    ll = [gt[0], gt[3] + (gt[5] * ds.RasterYSize)]
    ur = [gt[0] + (gt[1] * ds.RasterXSize), gt[3]]
    lr = [gt[0] + (gt[1] * ds.RasterXSize), gt[3] + (gt[5] * ds.RasterYSize)]

    if t_srs is not None and not ds_srs.IsSame(t_srs):
        ct = osr.CoordinateTransformation(ds_srs, t_srs)
        #Check to see if ct creation failed
        #if ct == NULL:
        #Check to see if transform failed
        #if not ct.TransformPoint(extent[0], extent[1]):
        #Need to check that transformed coordinates fall within appropriate bounds
        ul = ct.TransformPoint(*ul)
        ll = ct.TransformPoint(*ll)
        ur = ct.TransformPoint(*ur)
        lr = ct.TransformPoint(*lr)
    
    xmin = min(ul[0], ll[0], ur[0], lr[0])
    xmax = max(ul[0], ll[0], ur[0], lr[0])
    ymin = min(ul[1], ll[1], ur[1], lr[1])
    ymax = max(ul[1], ll[1], ur[1], lr[1])
    return [xmin, ymin, xmax, ymax] 

def get_extent_gt(gt, nx, ny):
    ul = [gt[0], gt[3]]
    ll = [gt[0], gt[3] + (gt[5] * ny)]
    ur = [gt[0] + (gt[1] * nx), gt[3]]
    lr = [gt[0] + (gt[1] * nx), gt[3] + (gt[5] * ny)]
    
    xmin = min(ul[0], ll[0], ur[0], lr[0])
    xmax = max(ul[0], ll[0], ur[0], lr[0])
    ymin = min(ul[1], ll[1], ur[1], lr[1])
    ymax = max(ul[1], ll[1], ur[1], lr[1])

    return [xmin, ymin, xmax, ymax] 

def get_extent_geom(geom):
    #Envelope is ul_x, ur_x, lr_y, ll_y (?)
    env = geom.GetEnvelope()
    #return xmin, ymin, xmax, ymax 
    return [env[0], env[2], env[1], env[3]]

def get_extent_ds(ds, t_srs=None):
    ds_srs = get_ds_srs(ds)
    geom = get_geom(ds, t_srs)
    return get_extent_geom(geom)

#Quick and dirty filter to check for points inside bbox
def pt_within_extent(x, y, extent):
    xmin, ymin, xmax, ymax = extent
    idx = ((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)) 
    return x[idx], y[idx]

#Pad extent
#Want to rewrite to allow for user-specified map units in addition to percentage
def pad_extent(extent, perc=0.1, uniform=False):
    e = np.array(extent)
    dx = e[2] - e[0]
    dy = e[3] - e[1]
    if uniform:
        dx = dy = np.mean([dx, dy])
    return e + (perc * np.array([-dx, -dy, dx, dy]))

"""
Extent notes:
gdalwarp uses '-te xmin ymin xmax ymax'
gdalbuildvrt uses '-te xmin ymin xmax ymax'
gdal_translate uses '-projwin ulx uly lrx lry' or '-projwin xmin ymax xmax ymin'
"""

#Need to migrate to use geom
#Compute union for arbitrary number of input datsets
#Output is [xmin, ymin, xmax, ymax]
def compute_union(*args):
    ref_ds = args[0]
    ref_srs = get_ds_srs(ref_ds) 
    union = get_extent(ref_ds) 
    for ds in args[1:]:
        #Want to check to make sure projections are identical
        #Note, the gt values are for top left corner of top left pixel
        #Center of pixel is actually(gt[0] + 0.5*gt[1], gt[3] + 0.5*gt[5])
        #The returned union is inclusive of both top left and bottom right pixels
        r = get_extent(ds, ref_srs)
        #This was for ul_lr
        #union = [min(r[0], union[0]), max(r[1], union[1]), max(r[2], union[2]), min(r[3], union[3])]
        union = [min(r[0], union[0]), min(r[1], union[1]), max(r[2], union[2]), max(r[3], union[3])]
    return union 

#Compute intersection for arbitrary number of input datsets
def compute_intersection(*args):
    ref_ds = args[0]
    ref_srs = get_ds_srs(ref_ds) 
    intsect = get_extent(ref_ds) 
    for ds in args[1:]:
        #Note, the gt values are for top left corner of top left pixel
        #Center of pixel is actually(gt[0] + 0.5*gt[1], gt[3] + 0.5*gt[5])
        #The returned intsect is inclusive of both top left and bottom right pixels
        r = get_extent(ds, ref_srs)
        #This was for ul_lr
        #intsect = [max(intsect[0], r[0]), min(intsect[1], r[1]), min(intsect[2], r[2]), max(intsect[3], r[3])]
        intsect = [max(intsect[0], r[0]), max(intsect[1], r[1]), min(intsect[2], r[2]), min(intsect[3], r[3])]
    #Check to make sure images actually intersect
    #if (intsect[2] <= intsect[0]) or (intsect[1] <= intsect[3]):
    if (intsect[2] <= intsect[0]) or (intsect[3] <= intsect[1]):
        #Throw error?
        intsect = None
    return intsect

#Temporary function to compute union from input list of geometries
#Switch to using lists of geom or lists of ds instead of expanded args
#Create additional method to do this from a list of ds
#Add option to return envelope, don't need additional functions to do this
#Note: this can return multipolygon geometry!
def union_geom(*args, **kwargs):
    convex=False
    union = args[0]
    for geom in args[1:]:
        union = union.Union(geom)
    if convex:
        union = union.ConvexHull()
    return union

#Check to make sure we have at least 2 input ds
def compute_union_geom(*args, **kwargs):
    ref_ds = args[0]
    ref_srs = get_ds_srs(ref_ds) 
    if 't_srs' in kwargs: 
        if kwargs['t_srs'] is not None:
            if not ref_srs.IsSame(kwargs['t_srs']):
                ref_srs = kwargs['t_srs']
    geom0 = get_geom(ref_ds, t_srs=ref_srs)
    union = geom0
    for ds in args[1:]:
        r = get_geom(ds, t_srs=ref_srs)
        union = union.Union(r)
    #Envelope is ul_x, ur_x, lr_y, lr_x
    #Define new geom class with better Envelope options?
    env = union.GetEnvelope()
    return [env[0], env[2], env[1], env[3]]

#Check to make sure we have at least 2 input ds
def compute_intersection_geom(*args, **kwargs):
    ref_ds = args[0]
    ref_srs = get_ds_srs(ref_ds) 
    if 't_srs' in kwargs: 
        if kwargs['t_srs'] is not None:
            if not ref_srs.IsSame(kwargs['t_srs']):
                ref_srs = kwargs['t_srs']
    geom0 = get_geom(ref_ds, t_srs=ref_srs)
    intsect = geom0
    valid = False
    for ds in args[1:]:
        r = get_geom(ds, t_srs=ref_srs)
        if intsect.Intersects(r):
            valid = True
            intsect = intsect.Intersection(r)
    #Check to make sure intersection is valid
    #***
    #This means checking that stereographic projections don't extend beyond equator
    #***
    #If it is the same as the original, then nothing intersects
    #if intsect.Equals(geom0):
    if not valid:
        return None 
    else:
        #Envelope is ul_x, ur_x, lr_y, lr_x
        #Define new geom class with better Envelope options?
        env = intsect.GetEnvelope()
        return [env[0], env[2], env[1], env[3]]

#Write for arbitrary number of inputs
def clip_to_intersection(ds1, ds2):
    intersection = compute_intersection(ds1, ds2)
    #if intersection is None:
    #convert that rectangle into pixels for each image by subtracting the top and left coordinates and dividing by the pixel size, rounding up.
    #band.ReadRaster(px1[0], px1[1], px1[2] - px1[0], px1[3] - px1[1], px1[2] - px1[0], px1[3] - px1[1], <band's datatype here>)

#This is necessary because extent precision is different
def extent_compare(e1, e2):
    e1_f = '%0.6f %0.6f %0.6f %0.6f' % tuple(e1)
    e2_f = '%0.6f %0.6f %0.6f %0.6f' % tuple(e2)
    return e1_f == e2_f

#This is necessary because extent precision is different
def res_compare(r1, r2):
    r1_f = '%0.6f' % r1
    r2_f = '%0.6f' % r2
    return r1_f == r2_f

#Clip raster by shape
#Note, this is a hack that uses gdalwarp command line util
#It is possible to do this with GDAL/OGR python API, but this works for now
#See: http://stackoverflow.com/questions/2220749/rasterizing-a-gdal-layer
def clip_raster_by_shp(dem_fn, shp_fn):
    import os, subprocess
    #This is ok when writing to outdir, but clip_raster_by_shp.sh writes to raster dir
    #try:
    #    with open(dem_fn) as f: pass
    #except IOError as e:
    cmd = 'clip_raster_by_shp.sh '+dem_fn+' '+shp_fn
    print cmd
    subprocess.call(cmd, shell=True)
    print
    dem_clip_fn = os.path.splitext(dem_fn)[0]+'_shpclip.tif'
    dem_clip_ds = gdal.Open(dem_clip_fn, gdal.GA_ReadOnly)
    return dem_clip_ds

#This will rasterize a geom for a given ma and geotransform
#Proper way would be to take in ds, transform geom to ds_srs, then convert to pixel coord
#Another need for Geoma
#See ogr_explode.py
#def rasterize_geom(ma, geom, gt=[0,1,0,0,1,0]):
def geom2mask(geom, ds):
    from PIL import Image, ImageDraw
    #width = ma.shape[1]
    #height = ma.shape[0]
    width = ds.RasterXSize
    height = ds.RasterYSize
    gt = ds.GetGeoTransform()

    img = Image.new('L', (width, height), 0)
    draw = ImageDraw.Draw(img)
    #Check to make sure we have polygon
    #3 is polygon, 6 is multipolygon
    #At present, this doesn't handle internal polygons
    #Want to set these to 0
    if (geom.GetGeometryType() == 3): 
        for ring in geom:
            pts = np.array(ring.GetPoints())
            px = np.array(mapToPixel(pts[:,0], pts[:,1], gt))
            px_poly = px.T.astype(int).ravel().tolist()
            draw.polygon(px_poly, outline=1, fill=1)
    elif (geom.GetGeometryType() == 6):
        for poly in geom:
            for ring in poly:
                pts = np.array(ring.GetPoints())
                px = np.array(mapToPixel(pts[:,0], pts[:,1], gt))
                px_poly = px.T.astype(int).ravel().tolist()
                draw.polygon(px_poly, outline=1, fill=1)
    # polygon = [(x1,y1),(x2,y2),...] or [x1,y1,x2,y2,...]
    mask = np.array(img).astype(bool)
    return ~mask

#These gdaldem functions should be able to ingest masked array
#Just write out temporary file, or maybe mem vrt?
#NOTE: probably want to smooth input DEM here!
def gdaldem_wrapper(fn, product='hs'):
    import os, subprocess
    import malib
    out_fn = os.path.splitext(fn)[0]+'_%s.tif' % product
    try:
        with open(fn) as f: pass
    except IOError as e:
        print "Unable to open %s" %fn

    valid_opt = ['hillshade', 'hs', 'slope', 'aspect', 'color-relief', 'TRI', 'TPI', 'roughness']
    bma = None
    if product in valid_opt: 
        if product == 'hs':
            product = 'hillshade'
        cmd = ['gdaldem', product, fn, out_fn]
        print ' '.join(cmd)
        subprocess.call(cmd, shell=False)
        ds = gdal.Open(out_fn, gdal.GA_ReadOnly)
        bma = malib.ds_getma(ds, 1)
    else:
        print "Invalid gdaldem option specified"
    return bma 

def gdaldem_slope(fn):
    return gdaldem_wrapper(fn, 'slope')

def gdaldem_aspect(fn):
    return gdaldem_wrapper(fn, 'aspect')

#Perhaps this should be generalized, and moved to malib
def bilinear(px, py, band_array, gt):
    '''Bilinear interpolated point at(px, py) on band_array
    example: bilinear(2790501.920, 6338905.159)'''
    #import malib
    #band_array = malib.checkma(band_array)
    ndv = band_array.fill_value
    ny, nx = band_array.shape
    # Half raster cell widths
    hx = gt[1]/2.0
    hy = gt[5]/2.0
    # Calculate raster lower bound indices from point
    fx =(px -(gt[0] + hx))/gt[1]
    fy =(py -(gt[3] + hy))/gt[5]
    ix1 = int(np.floor(fx))
    iy1 = int(np.floor(fy))
    # Special case where point is on upper bounds
    if fx == float(nx - 1):
        ix1 -= 1
    if fy == float(ny - 1):
        iy1 -= 1
    # Upper bound indices on raster
    ix2 = ix1 + 1
    iy2 = iy1 + 1
    # Test array bounds to ensure point is within raster midpoints
    if(ix1 < 0) or(iy1 < 0) or(ix2 > nx - 1) or(iy2 > ny - 1):
        return ndv
    # Calculate differences from point to bounding raster midpoints
    dx1 = px -(gt[0] + ix1*gt[1] + hx)
    dy1 = py -(gt[3] + iy1*gt[5] + hy)
    dx2 =(gt[0] + ix2*gt[1] + hx) - px
    dy2 =(gt[3] + iy2*gt[5] + hy) - py
    # Use the differences to weigh the four raster values
    div = gt[1]*gt[5]
    return(band_array[iy1,ix1]*dx2*dy2/div +
            band_array[iy1,ix2]*dx1*dy2/div +
            band_array[iy2,ix1]*dx2*dy1/div +
            band_array[iy2,ix2]*dx1*dy1/div)

#Note: the existing egm96-5 dataset has problematic extent
#warplib writes out correct res/extent, but egm96 is empty
#Eventually accept geoma
def wgs84_to_egm96(dem_ds, geoid_dir=None):
    import os
    import malib
    import warplib

    #Check input dem_ds - make sure WGS84

    geoid_dir = os.getenv('ASP_DATA')
    if geoid_dir is None:
        print "No geoid directory available"
        print "Set ASP_DATA or specify"
    
    egm96_fn = geoid_dir+'/geoids-1.1/egm96-5.tif' 
    try:
        with open(egm96_fn) as f: pass
    except IOError as e:
        sys.exit("Unable to find "+egm96_fn)
    egm96_ds = gdal.Open(egm96_fn)

    #Warp egm96 to match the input ma
    ds_list = warplib.memwarp_multi([dem_ds, egm96_ds], res='first', extent='first', t_srs='first') 

    #Try two-step with extent/res in wgs84
    #ds_list = warplib.memwarp_multi([dem_ds, egm96_ds], res='first', extent='intersection', t_srs='last') 
    #ds_list = warplib.memwarp_multi([dem_ds, ds_list[1]], res='first', extent='first', t_srs='first')

    print "Extracting ma from dem and egm96 ds"
    dem = malib.ds_getma(ds_list[0])
    egm96=malib.ds_getma(ds_list[1])

    print "Removing offset"
    dem_egm96 = dem - egm96
   
    return dem_egm96

#Run ASP dem_geoid adjustment utility
#Note: this is multithreaded
def dem_geoid(dem_fn):
    import os
    import malib
    out_prefix = os.path.splitext(dem_fn)[0]
    adj_fn = out_prefix +'-adj.tif'
    if not os.path.exists(adj_fn):
        import subprocess
        cmd_args = ["-o", out_prefix, dem_fn]
        cmd = ['dem_geoid'] + cmd_args
        #cmd = 'dem_geoid -o %s %s' % (out_prefix, dem_fn)
        print ' '.join(cmd)
        subprocess.call(cmd, shell=False)
    adj_ds = gdal.Open(adj_fn, gdal.GA_ReadOnly)
    return adj_ds
    #return malib.ds_getma(adj_ds, 1)

def get_xy_grids(ds, stride=1, getval=False):
    gt = ds.GetGeoTransform()
    #stride = stride_m/gt[1]
    pX = np.arange(0, ds.RasterXSize, stride)
    pY = np.arange(0, ds.RasterYSize, stride)
    psamp = np.meshgrid(pX, pY)
    mX, mY = pixelToMap(psamp[0], psamp[1], gt)
    return mX, mY
