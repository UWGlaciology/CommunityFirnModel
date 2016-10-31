#! /usr/bin/env python

"""
David Shean
dshean@gmail.com
2/15/13

To do:
Much better type checking
Implement multiprocessing!
Proper vdatum tags in geotiff header
filename preservation
Instead of doing the 'first' stuff, check actual values before writing out original dataset
"""

import sys
import os 
import argparse
import math

import gdal
import osr

import geolib
import malib

#Note: can run into filesystem limits for number of open files
#http://superuser.com/questions/433746/is-there-a-fix-for-the-too-many-open-files-in-system-error-on-os-x-10-7-1
gdal.SetConfigOption('GDAL_MAX_DATASET_POOL_SIZE', '2048')
#Need to set in the shell
#ulimit -S -n 2048
#cmd = ['ulimit', '-S', '-n', 2048]
#import subprocess
#subprocess.call(cmd)

mem_drv = gdal.GetDriverByName('MEM')
gtif_drv = gdal.GetDriverByName ("GTiff")
gdal_opt = ['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=IF_SAFER']

#Use this to warp to file - no need to write to memory then write to file 
#gdal is much better about memory management
#def diskwarp(src_ds, dst_fn=None, res=None, extent=None, t_srs=None, r='cubic', driver=gtif_drv):
#    if dst_fn is None:
#        dst_fn = os.path.splitext(src_ds.GetFileList()[0])[0]+'_warp.tif'
#    dst_ds = warp(src_ds, res, extent, t_srs, r, driver, dst_fn)
#    #return dst_ds
#    dst_ds = None
#    return dst_ds

#Preserve these input options
def memwarp(src_ds, res=None, extent=None, t_srs=None, r=None, driver=mem_drv):
    return warp(src_ds, res, extent, t_srs, r, driver)

def warp(src_ds, res=None, extent=None, t_srs=None, r='cubic', driver=mem_drv, dst_fn=None):
    src_srs = geolib.get_ds_srs(src_ds)
    
    if t_srs is None:
        t_srs = geolib.get_ds_srs(src_ds)
    
    src_gt = src_ds.GetGeoTransform()
    #Note: get_res returns [x_res, y_res]
    #Could just use gt here and average x_res and y_res
    src_res = geolib.get_res(src_ds, t_srs=t_srs, square=True)[0]

    if res is None:
        res = src_res

    if extent is None:
        extent = geolib.get_extent_ds(src_ds, t_srs=t_srs)
    
    #Note: GDAL Lanczos creates block artifacts
    #Wait for gdalwarp to support gaussian resampling
    #Want to use Lanczos for downsampling
    #if src_res < res:
    #    gra = gdal.GRA_Lanczos
    #See http://blog.codinghorror.com/better-image-resizing/
    # Suggests cubic for downsampling, bilinear for upsampling
    #    gra = gdal.GRA_Cubic
    #Cubic for upsampling
    #elif src_res >= res:
    #    gra = gdal.GRA_Bilinear

    #At this point, the resolution and extent values must be float
    #Extent must be list
    res = float(res)
    extent = [float(i) for i in extent]

    #Might want to move this to memwarp_multi, keep memwarp basic w/ gdal.GRA types

    #Note:GRA_CubicSpline created huge block artifacts for the St. Helen's compute_dh WV cases
    #Stick with CubicSpline for both upsampling/downsampling for now
    if r == 'near':
        #Note: Nearest respects nodata when downsampling
        gra = gdal.GRA_NearestNeighbour
    elif r == 'bilinear':
        gra = gdal.GRA_Bilinear
    elif r == 'cubic':
        gra = gdal.GRA_Cubic
    elif r == 'cubicspline':
        gra = gdal.GRA_CubicSpline
    elif r == 'average':
        gra = gdal.GRA_Average
    elif r == 'lanczos':
        gra = gdal.GRA_Lanczos
    elif r == 'mode':
        #Note: Mode respects nodata when downsampling, but very slow
        gra = gdal.GRA_Mode
    else:
        sys.exit("Invalid resampling method")

    #Create progress function
    prog_func = gdal.TermProgress
    
    if dst_fn is None:
        #This is a dummy fn if only in mem, but can be accessed later via GetFileList()
        #Actually, no, doesn't look like the filename survivies
        dst_fn = ''
    
    #Compute output image dimensions
    dst_nl = int(round((extent[3] - extent[1])/res))
    dst_ns = int(round((extent[2] - extent[0])/res))
    #dst_nl = int(math.ceil((extent[3] - extent[1])/res))
    #dst_ns = int(math.ceil((extent[2] - extent[0])/res))
    #dst_nl = int(math.floor((extent[3] - extent[1])/res))
    #dst_ns = int(math.floor((extent[2] - extent[0])/res))
    print 'nl: %i ns: %i res: %0.3f' % (dst_nl, dst_ns, res)
    #Create output dataset
    src_b = src_ds.GetRasterBand(1)
    src_dt = src_b.DataType

    dst_ds = driver.Create(dst_fn, dst_ns, dst_nl, src_ds.RasterCount, src_dt) 

    dst_ds.SetProjection(t_srs.ExportToWkt())
    #Might be an issue to use src_gt rotation terms here with arbitrary extent/res
    dst_gt = [extent[0], res, src_gt[2], extent[3], src_gt[4], -res]
    dst_ds.SetGeoTransform(dst_gt)
    
    for n in range(1, src_ds.RasterCount+1):
        src_b = src_ds.GetRasterBand(n)
        src_ndv = malib.get_ndv_b(src_b)
        if src_ndv is not None:
            b = dst_ds.GetRasterBand(n)
            b.SetNoDataValue(src_ndv)
            #Need this to actually make src_ndv stick for dst_ds
            b.Fill(src_ndv)
    
    #Note: default maxerror=0.0
    #Shouldn't neet to specify srs?
    #result = gdal.ReprojectImage(src_ds, dst_ds, gra)
    result = gdal.ReprojectImage(src_ds, dst_ds, src_srs.ExportToWkt(), t_srs.ExportToWkt(), \
    gra, 0.0, 0.0, prog_func)

    #Return GDAL dataset object in memory
    if driver != mem_drv:
        dst_ds.FlushCache()
    return dst_ds

def memwarp_multi(src_ds_list, res='first', extent='intersection', t_srs='first', r='cubic'):
    #Type cast arguments as str for evaluation
    #Avoid path errors
    #res = str(res)
    #extent = str(extent)
    #t_srs = str(t_srs)

    if t_srs == 'first':
        t_srs = geolib.get_ds_srs(src_ds_list[0])
    elif t_srs == 'last':
        t_srs = geolib.get_ds_srs(src_ds_list[-1])
    elif isinstance(t_srs, gdal.Dataset):
        t_srs = geolib.get_ds_srs(t_srs)
    elif isinstance(t_srs, str) and os.path.exists(t_srs): 
        t_srs = geolib.get_ds_srs(gdal.Open(t_srs))
    else:
        temp = osr.SpatialReference()
        temp.ImportFromProj4(t_srs)
        t_srs = temp

    #Compute output resolution in t_srs
    #Returns min, max, mean, med
    res_stats = geolib.get_res_stats(src_ds_list, t_srs=t_srs)
    if res == 'first':
        res = geolib.get_res(src_ds_list[0], t_srs=t_srs, square=True)[0]
    elif res == 'last':
        res = geolib.get_res(src_ds_list[-1], t_srs=t_srs, square=True)[0]
    elif res == 'min':
        res = res_stats[0]
    elif res == 'max':
        res = res_stats[1]
    elif res == 'mean':
        res = res_stats[2]
    elif res == 'med':
        res = res_stats[3]
    elif res == 'source':
        res = None
    elif isinstance(res, gdal.Dataset):
        res = geolib.get_res(res, t_srs=t_srs, square=True)[0]
    elif isinstance(res, str) and os.path.exists(res): 
        res = geolib.get_res(gdal.Open(res), t_srs=t_srs, square=True)[0]
    else:
        res = float(res)

    #Extent returned is xmin, ymin, xmax, ymax
    if len(src_ds_list) == 1 and (extent == 'intersection' or extent == 'union'):
        extent = None
    elif extent == 'first':
        extent = geolib.get_extent_ds(src_ds_list[0], t_srs=t_srs)
    elif extent == 'last':
        extent = geolib.get_extent_ds(src_ds_list[-1], t_srs=t_srs)
    elif extent == 'intersection':
        #By default, compute_intersection takes ref_srs from ref_ds
        extent = geolib.compute_intersection_geom(*src_ds_list, t_srs=t_srs)
        if len(src_ds_list) > 1 and extent is None:
            #print "Input images do not intersect"
            sys.exit("Input images do not intersect")
    elif extent == 'union':
        #Need to clean up union t_srs handling
        extent = geolib.compute_union_geom(*src_ds_list, t_srs=t_srs)
        #extent = geolib.compute_union(*src_ds_list)
    elif extent == 'source':
        extent = None
    elif isinstance(extent, gdal.Dataset):
        extent = geolib.get_extent_ds(extent, t_srs=t_srs)
    elif isinstance(extent, str) and os.path.exists(extent): 
        extent = geolib.get_extent_ds(gdal.Open(extent), t_srs=t_srs)
    elif isinstance(extent, (list, tuple)):
        extent = list(extent)
    else:
        extent = [float(i) for i in extent.split(' ')]

    print
    print "Warping all inputs to the following:"
    print "Resolution: %s" % res
    print "Extent: %s" % str(extent)
    print "Projection: '%s'" % t_srs.ExportToProj4()
    print "Resampling alg: %s" % r  
    print

    out_ds_list = []
    for i, ds in enumerate(src_ds_list):
        fn_list = ds.GetFileList()
        fn = '[memory]'
        if fn_list is not None:
            fn = fn_list[0]
        print "%i of %i: %s" % (i+1, len(src_ds_list), fn)
        #Check each ds to see if warp is necessary
        ds_res = geolib.get_res(ds, square=True)[0]
        #Note: this rounds to %0.5f, extent_compare uses %0.6f
        #ds_extent = geolib.get_extent_ds(ds)
        ds_extent = geolib.get_extent(ds)
        ds_t_srs = geolib.get_ds_srs(ds)

        rescheck=False
        if res is None or geolib.res_compare(res, ds_res):
            rescheck=True
        extentcheck=False
        #Note: extent_compare is necessary to deal with rounding errors
        if extent is None or geolib.extent_compare(extent, ds_extent):
            extentcheck=True
        srscheck=False
        if (t_srs.IsSame(ds_t_srs)):
            srscheck=True
            
        #if geolib.res_compare(res, ds_res) and geolib.extent_compare(extent, ds_extent) and (t_srs.IsSame(ds_t_srs)):
        if rescheck and extentcheck and srscheck:
            out_ds_list.append(ds)
        else:
            dst_ds = memwarp(ds, res, extent, t_srs, r)
            out_ds_list.append(dst_ds)
            dst_ds = None

    return out_ds_list

def memwarp_multi_fn(src_fn_list, res='first', extent='intersection', t_srs='first', r='cubic'):
    #Should implement proper error handling here
    import genlib
    if not genlib.fn_list_check(src_fn_list):
        sys.exit('Missing input file(s)')
    src_ds_list = [gdal.Open(fn, gdal.GA_ReadOnly) for fn in src_fn_list]
    return memwarp_multi(src_ds_list, res, extent, t_srs, r)

def writeout(ds, outfn):
    print "Writing out %s" % outfn 
    #Use outfn extension to get driver
    #This may have issues if outfn already exists and the mem ds has different dimensions/res
    out_ds = gtif_drv.CreateCopy(outfn, ds, 0, options=gdal_opt)
    out_ds = None

def main():
    #Can't specify arbitrary fn, res when limiting choices
    tr_choices = ['first', 'last', 'min', 'max', 'mean', 'med', 'source', '"fn"', '"res"']
    te_choices = ['first', 'last', 'intersection', 'union', 'source', '"fn"', '"extent"']
    t_srs_choices = ['first', 'last', '"fn"', '"proj4str"']
    r_choices = ['near', 'bilinear', 'cubic', 'cubicspline', 'lanczos', 'average', 'mode']
    
    parser = argparse.ArgumentParser(description='Utility to warp stacks of rasters to the same res/extent/proj')
    #parser.add_argument('-tr', default='first', choices=tr_choices, help='Output resolution')
    parser.add_argument('-tr', default='first', help='Output resolution (default: %(default)s)')
    parser.add_argument('-te', default='intersection', help='Output extent (default: %(default)s)')
    parser.add_argument('-t_srs', default='first', help='Output projection (default: %(default)s)')
    parser.add_argument('-r', type=str, default='cubic', help='Resampling algorithm (default: %(default)s)', choices=r_choices)
    parser.add_argument('src_fn_list', nargs='+', help='Input filenames (img1.tif img2.tif ...)')
    args = parser.parse_args()

    #Can also provide filename for any of the tr, te, t_srs options

    #Automatically fix UL, LR coords?
    #Format is minx, miny, maxx, maxy
    #extent = [-197712.13, -2288712.83, -169650.72, -2253490.42] 
    #extent=[-195712.13, -2286712.83, -170650.72, -2256490.42]
    #res = 2.0

    print
    print "Input parameters"
    print "Resolution: %s" % str(args.tr)
    print "Extent: %s" % str(args.te)
    print "Projection: %s" % str(args.t_srs)
    print "Resampling alg: %s" % str(args.r)
    print
    
    ds_list = memwarp_multi_fn(args.src_fn_list, res=args.tr, extent=args.te, t_srs=args.t_srs, r=args.r)
   
    for i, ds in enumerate(ds_list):
        #Note: the following doesn't work for mem filenames
        #outfn = os.path.splitext(ds.GetFileList()[0])[0]+'_warp.tif'
        outfn = os.path.splitext(args.src_fn_list[i])[0]+'_warp.tif'
        writeout(ds, outfn)     

if __name__ == "__main__":
    main()
