#! /usr/bin/env python

#David Shean
#dshean@gmail.com
#3/20/12

#Library containing various functions for masked arrays

import sys
import os

import numpy as np
import gdal
import osr

gdal_opt = ['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=IF_SAFER']
#gdal_opt += ['BLOCKXSIZE=1024', 'BLOCKYSIZE=1024']

#Notes on geoma
#Note: Need better init overloading
#http://stackoverflow.com/questions/141545/overloading-init-in-python

#Might make more sense to create ma subclass, and add gdal ds as new object
#http://stackoverflow.com/questions/12597827/how-to-subclass-numpy-ma-core-masked-array
#http://docs.scipy.org/doc/numpy/user/basics.subclassing.html

#Want to implement basic array indexing with map coordinates, display in plt

#See pyresample
#https://github.com/talltom/PyRaster/blob/master/rasterIO.py
#http://www2-pcmdi.llnl.gov/cdat/tutorials/training/cdat_2004/06-arrays-variables-etc.pdf
#http://docs.scipy.org/doc/numpy/user/basics.subclassing.html

#Class bundling a gdal Dataset and numpy masked array
#Needed to preserve original filename and other dataset properties (e.g. projection)
#Need to worry about functions that modify either the ds or ma in place
#Destruction - write out ds with latest ma?
#class GeoMA(np.ma.MaskedArray, gdal.Dataset):
class GeoMA():
    def __init__(self, a):
        #Check to see if input is gdal datset
        if isinstance(a, gdal.Dataset):
            self.ds = a 
            #Note: this returns a list, let's force a single filename for now
            self.fn = a.GetFileList()[0]
        #Check to see if input is filename
        elif isinstance(a, basestring):
            self.fn = a
            self.ds = gdal.Open(a)
        #else:
        #   fail gracefully
        
        #Want to check to make sure ds is not none
        #Now open a masked array using getma method - basically a wrapper for masked_equal with ndv
        self.ma = ds_getma(self.ds)

#Want to add error attributes
#Want to make consistent with stack_count vs count keywords/attributes
class DEMStack:
    def __init__(self, fn_list, res=None, extent=None, trend=False, med=False, stats=True):
        self.fn_list = fn_list
        self.res = res
        self.extent = extent
        self.trend = trend
        self.med = med
        self.stats = stats
        if len(fn_list) == 0:
            raise ValueError('Input filename list empty')
        #Determine appropriate stack filename
        self.stack_fn = os.path.splitext(fn_list[0])[0] + '_' \
                        + os.path.splitext(os.path.split(fn_list[-1])[1])[0] \
                        + '_stack_%i' % len(self.fn_list) + '.npz'
        self.get_date_list() 
        if os.path.exists(self.stack_fn):
            self.loadstack()
            #This was an attempt to ensure that input fn/res/extent is consistent with saved npz
            #Should really check to see if new extent falls within existing extent
            #Only regenerate if new extent is larger
            #Res check is not working correctly
            #print extent, self.extent
            #print res, self.res
            #if (self.extent != extent) or (self.res != res):
                #self.res = res
                #self.extent = extent
                #self.makestack()
                #if self.stats: 
                    #self.compute_stats()
                    #self.write_stats()
                #self.savestack()
        else:
            self.makestack()
            if self.stats: 
                self.compute_stats()
                self.write_stats()
            #self.mean_hillshade()
            self.savestack()
        self.name = None

    #p = os.path.join(topdir, d, '*_warp.tif'))
    def get_fn_list(p):
        fn_list = glob.glob(p)
        return fn_list

    def get_res(self):
        self.res = np.mean([self.gt[1], np.abs(self.gt[5])])

    def get_extent(self):
        import geolib
        self.extent = geolib.get_extent_gt(self.gt, self.ma_stack.shape[1], self.ma_stack.shape[0])

    def makestack(self):
        import warplib
        print "Creating stack of %i files" % len(self.fn_list)
        #Jako front 
        #res = 16
        res = 'min'
        if self.res is not None:
            res=self.res
        #extent = '-195705.297256 -2286746.61662 -170642.601955 -2256442.61662'
        #extent = 'intersection'
        extent = 'union'
        if self.extent is not None:
            extent = self.extent
        ds_list = warplib.memwarp_multi_fn(self.fn_list, res=res, extent=extent, t_srs='first')
        #Note: might not need ma here in the 0 axis - shouldn't be any missing data
        self.ma_stack = np.ma.array([ds_getma(ds) for ds in ds_list])
        #Might want to convert to proj4
        self.proj = ds.GetProjectionRef()
        self.gt = ds.GetGeoTransform()
        self.res = self.get_res()

    #Note: want ot change the variable names min/max here
    #Might be better to save out as multiband GTiff here
    def savestack(self):
        print "Saving stack to: %s" % self.stack_fn
        if self.stats:
            out_args = dict(ma_stack_full=self.ma_stack.filled(np.nan), \
            count=self.stack_count.filled(0), mean=self.stack_mean.filled(np.nan), \
            min=self.stack_min.filled(np.nan), max=self.stack_max.filled(np.nan), \
            std=self.stack_std.filled(np.nan), proj=str(self.proj), gt=self.gt) 
        else:
            out_args = dict(ma_stack_full=self.ma_stack.filled(np.nan), \
            proj=str(self.proj), gt=self.gt)
        if self.trend:
            out_args['trend'] = self.stack_trend.filled(np.nan)
        if self.med:
            out_args['med'] = self.stack_med.filled(np.nan)
        #out_args['res'] = str(self.res)
        #out_args['extent'] = str(self.extent)
        np.savez_compressed(self.stack_fn, **out_args)

    def loadstack(self):
        print "Loading stack from: %s" % self.stack_fn
        data = np.load(self.stack_fn)
        self.ma_stack = np.ma.fix_invalid(data['ma_stack_full'])
        #Note: the str is an intermediate fix - all new stacks should have str written
        self.proj = str(data['proj'])
        self.gt = data['gt']
        #Load originally specified res
        #if 'res' in data:
        #    self.res = data['res']
        #else:
        self.res = self.get_res()
        #Load originally specified extent
        #if 'extent' in data:
        #    self.extent = data['extent']
        #else:
        #    self.extent = self.get_extent() 
        if self.stats:
            #Could do this individually to save time
            statlist = ['count', 'mean', 'std', 'min', 'max']
            if self.trend:
                statlist.append('trend')
            if self.med:
                statlist.append('med')
            if all([s in data for s in statlist]):
                self.stack_count = np.ma.masked_equal(data['count'], 0)
                self.stack_mean = np.ma.fix_invalid(data['mean'], fill_value=-9999)
                self.stack_std = np.ma.fix_invalid(data['std'], fill_value=-9999)
                self.stack_min = np.ma.fix_invalid(data['min'], fill_value=-9999)
                self.stack_max = np.ma.fix_invalid(data['max'], fill_value=-9999)
                if self.trend:
                    self.stack_trend = np.ma.fix_invalid(data['trend'], fill_value=-9999)
                if self.med:
                    self.stack_med = np.ma.fix_invalid(data['med'], fill_value=-9999)
            else:
                self.compute_stats()
                self.write_stats()
                self.savestack()
        data.close()

    #This needs some work - will break with nonstandard filenames
    def get_date_list(self):
        import dateutil
        import timelib
        import matplotlib.dates
        from datetime import datetime
        #self.date_list = np.ma.array([dateutil.parser.parse(os.path.split(fn)[1][0:13], fuzzy=True) for fn in self.fn_list])
        #This will return None if no date in fn
        date_list = [timelib.fn_getdatetime(os.path.split(fn)[-1]) for fn in self.fn_list]
        self.date_list = np.ma.masked_equal([datetime(1,1,1) if d is None else d for d in date_list], datetime(1,1,1))
        
        #self.date_list = np.ma.array([dateutil.parser.parse(os.path.split(fn)[1][3:12], fuzzy=True) for fn in self.fn_list])
        self.date_list_o = np.ma.array([matplotlib.dates.date2num(d) for d in self.date_list.filled()], mask=self.date_list.mask)

    def compute_stats(self):
        print "Compute stack count"
        self.stack_count = np.ma.masked_equal(self.ma_stack.count(axis=0), 0).astype(np.uint16)
        self.stack_count.set_fill_value(0)
        print "Compute stack mean"
        self.stack_mean = self.ma_stack.mean(axis=0).astype(np.float32)
        self.stack_mean.set_fill_value(-9999)
        print "Compute stack std"
        self.stack_std = self.ma_stack.std(axis=0).astype(np.float32)
        #Only want to preserve values where count > 1
        self.stack_std.mask = (self.stack_count <= 1) 
        self.stack_std.set_fill_value(-9999)
        print "Compute stack min"
        self.stack_min = self.ma_stack.min(axis=0).astype(np.float32)
        self.stack_min.set_fill_value(-9999)
        print "Compute stack max"
        self.stack_max = self.ma_stack.max(axis=0).astype(np.float32)
        self.stack_max.set_fill_value(-9999)
        if self.trend:
            print "Compute stack linear trend"
            self.linreg()
            self.stack_trend.set_fill_value(-9999)
        if self.med:
            print "Compute stack med"
            self.stack_med = np.ma.median(self.ma_stack, axis=0).astype(np.float32)
            self.stack_med.set_fill_value(-9999)

    def write_stats(self):
        #if not hasattr(self, 'stack_count'):
        stat_list = ['stack_count', 'stack_mean', 'stack_std', 'stack_min', 'stack_max']
        if self.trend:
            stat_list.append('stack_trend')
        if self.med:
            stat_list.append('stack_med')
        if any([not hasattr(self, i) for i in stat_list]):
            self.compute_stats()
        print "Writing out stats"
        #Create dummy ds - might want to use vrt here instead
        driver = gdal.GetDriverByName("MEM")
        ds = driver.Create('', self.ma_stack.shape[2], self.ma_stack.shape[1], 1, gdal.GDT_Float32)
        ds.SetGeoTransform(self.gt)
        ds.SetProjection(self.proj)
        #Write out with malib, should preserve ma type
        out_prefix = os.path.splitext(self.stack_fn)[0]
        writeGTiff(self.stack_count, out_prefix+'_count.tif', ds)
        writeGTiff(self.stack_mean, out_prefix+'_mean.tif', ds)
        writeGTiff(self.stack_std, out_prefix+'_std.tif', ds)
        writeGTiff(self.stack_min, out_prefix+'_min.tif', ds)
        writeGTiff(self.stack_max, out_prefix+'_max.tif', ds)
        if self.trend:
            writeGTiff(self.stack_trend, out_prefix+'_trend.tif', ds)
        if self.med:
            writeGTiff(self.stack_med, out_prefix+'_med.tif', ds)
    
    """
    def linreg_mstats(self):
        #import scipy.mstats
        from scipy.mstats import linregress
        x = self.date_list_o
        out = np.zeros_like(self.ma_stack[0])
        for row in np.arange(self.ma_stack.shape[1]):
            print '%i of %i rows' % (row, self.ma_stack.shape[1])
            for col in np.arange(self.ma_stack.shape[2])
                out[row, col] = linregress(x, self.ma_stack[:, row, col])[0]
        #out is in units of m/day since x is ordinal date
        out *= 365.25
        self.stack_trend = out
    """

    #Compute linear regression for every pixel in stack
    def linreg(self):
        from numpy.linalg import solve
        n_min = 3 
        #Only compute where we have n_min unmasked values in time
        valid_idx = self.stack_count.data >= n_min
        y = self.ma_stack[:, valid_idx]
        #Reshape to 2D
        #origshape = self.ma_stack.shape
        #newshape = (origshape[0], origshape[1] * origshape[2])
        #y = self.ma_stack.reshape(newshape)
        #Extract mask for axis 0 - invert, True where data is available
        mask = ~y.mask
        #Remove masks, fills with fill_value
        y = y.data
        #Independent variable is time ordinal
        x = self.date_list_o
        x = x.data
        #Prepare matrices
        X = np.c_[x, np.ones_like(x)]
        a = np.swapaxes(np.dot(X.T, (X[None, :, :] * mask.T[:, :, None])), 0, 1)
        b = np.dot(X.T, (mask*y))
        #Solve for slope/intercept
        r = np.linalg.solve(a, b.T)
        #Reshape to original dimensions
        #r = r.reshape(origshape[1], origshape[2], 2)
        #Create output grid with original dimensions
        out = np.ma.masked_all_like(self.ma_stack[0])
        #Fill in the valid indices
        out[valid_idx] = r[:,0]
        #out is in units of m/day since x is ordinal date
        out *= 365.25
        #Filter out clearly bogus values here?
        self.stack_trend = np.ma.array(out)
    
    def mean_hillshade(self):
        if hasattr(self, 'stack_med'):
            in_fn = os.path.splitext(self.stack_fn)[0]+'_med.tif'
        elif hasattr(self, 'stack_mean'):
            in_fn = os.path.splitext(self.stack_fn)[0]+'_mean.tif'
        else:
            self.compute_stats()
            in_fn = os.path.splitext(self.stack_fn)[0]+'_mean.tif'
        if not os.path.exists(in_fn):
            self.write_stats()
        import geolib
        print "Generate shaded relief from mean"
        self.stack_mean_hs = geolib.gdaldem_wrapper(in_fn, 'hs')

def get_edges(a, convex=False):
    a = checkma(a)
    #Need to deal with RGB images here
    #Need to be careful, probably want to take minimum value from masks
    if a.ndim == 3:
        #Assume that the same mask applies to each band
        #Only create edges for one band
        b = a[:,:,0]
        #Could convert to HSL and do edges for L channel
        #b = a[:,:,0].mask + a[:,:,1].mask + a[:,:,2].mask
    else:
        b = a

    #Compute edges along both axes, need both to handle undercuts
    #These are inclusive, indices indicate position of unmasked data
    edges0 = np.ma.notmasked_edges(b, axis=0)
    edges1 = np.ma.notmasked_edges(b, axis=1)
    edges = np.array([np.concatenate([edges0[0][0], edges0[1][0], edges1[1][0], edges1[0][0]]), np.concatenate([edges0[0][1], edges0[1][1], edges1[1][1], edges1[0][1]])])

    #This is a rough outline - needs testing
    if convex:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(edges.T)
        #edges = edges.T[hull.simplices]
        #This is in scipy v0.14
        #edges0 = edges1 = hull.vertices

    return edges0, edges1, edges

def get_edgemask(a, edge_env=False, convex=False, dilate=False):
    a = checkma(a)
    #Need to deal with RGB images here
    #Need to be careful, probably want to take minimum value from masks
    if a.ndim == 3:
        #Assume that the same mask applies to each band
        #Only create edges for one band
        b = a[:,:,0]
        #Could convert to HSL and do edges for L channel
        #b = a[:,:,0].mask + a[:,:,1].mask + a[:,:,2].mask
    else:
        b = a
    
    #Get pixel locations of edges
    edges0, edges1, edges = get_edges(b, convex)

    #Compute min/max indices
    #minrow, maxrow, mincol, maxcol
    #edge_bbox = [edges0[0][0].min(), edges0[1][0].max() + 1, edges0[0][1].min(), edges0[1][1].max() + 1]
    edge_bbox = [edges[0].min(), edges[0].max() + 1, edges[1].min(), edges[1].max() + 1]

    #Initialize new mask arrays
    #Need both to deal with undercuts
    colmask = np.empty_like(b.mask)
    colmask[:] = True
    rowmask = np.empty_like(b.mask)
    rowmask[:] = True

    #Loop through each item in the edge list
    #i is index, col is the column number listed at index i
    #unmask pixels between specified row numbers at index i
    for i,col in enumerate(edges0[0][1]):
        colmask[edges0[0][0][i]:edges0[1][0][i], col] = False

    for j,row in enumerate(edges1[0][0]):
        rowmask[row, edges1[0][1][j]:edges1[1][1][j]] = False

    #Combine the two masks with or operator
    newmask = np.logical_or(colmask, rowmask)

    if dilate:
        print "Dilating edgemask"
        import scipy.ndimage as ndimage 
        n = 3
        #Note: this results in unmasked elements near image corners
        #This will erode True values, which correspond to masked elements
        #Effect is to expand the unmasked area
        #newmask = ndimage.morphology.binary_erosion(newmask, iterations=n)
        #Now dilate to return to original size
        #newmask = ndimage.morphology.binary_dilation(newmask, iterations=n)
        #This is a more logical approach, dilating unmasked areas
        newmask = ~ndimage.morphology.binary_dilation(~newmask, iterations=n)
        newmask = ~ndimage.morphology.binary_erosion(~newmask, iterations=n)

    if edge_env:
        return newmask, edge_bbox
    else:
        return newmask

#This will update the mask to remove interior holes from unmasked region
def apply_edgemask(a, trim=False):
    newmask, edge_bbox = get_edgemask(a, edge_env=True)

    if a.ndim > 2:
        for i in range(a.ndim):
            a.mask[:,:,i] = newmask
    else:
        a.mask[:] = newmask
    
    #Option to return a trimmed array
    if trim:
        return a[edge_bbox[0]:edge_bbox[1], edge_bbox[2]:edge_bbox[3]]
    else:
        return a

#This will trim an array to the unmasked region
#Want to update gt here as well - see ndvtrim.py
def masktrim(a):
    a = checkma(a)
    idx = np.ma.notmasked_edges(a, axis=0)
    minrow = idx[0][0].min()
    maxrow = idx[1][0].max()
    mincol = idx[0][1].min()
    maxcol = idx[1][1].max()
    #print minrow, maxrow, mincol, maxcol
    #Want to make sure these are inclusive
    return a[minrow:maxrow, mincol:maxcol]

#Return a common mask for a set of input ma
def common_mask(ma_list, apply=False):
    if type(ma_list) is not list:
        print "Input must be list of masked arrays"
        return None
    #Note: a.mask will return single False if all elements are False
    #np.ma.getmaskarray(a) will return full array of False
    #ma_list = [np.ma.array(a, mask=np.ma.getmaskarray(a), shrink=False) for a in ma_list]
    a = np.ma.array(ma_list, shrink=False)
    #Check array dimensions
    #Check dtype = bool
    #Masked values are listed as true, so want to return any()
    #a+b+c - OR (any)
    mask = a.mask.any(axis=0)
    #a*b*c - AND (all)
    #return a.all(axis=0)
    if apply:
        return [np.ma.array(b, mask=mask) for b in ma_list] 
    else:
        return mask

#This will attempt to remove islands
def mask_islands(a, iterations=3):
    import scipy.ndimage as ndimage 
    a = checkma(a)
    #newmask = a.mask
    newmask = np.ma.getmaskarray(a)
    newmask = ndimage.morphology.binary_dilation(newmask, iterations=iterations)
    newmask = ndimage.morphology.binary_erosion(newmask, iterations=iterations)
    return np.ma.array(a, mask=newmask)

#This will fill internal holes in the mask
#This should be used to mask outer edges after inpainting or gdal_fillnodata
def maskfill(a, iterations=1, erode=False):
    import scipy.ndimage as ndimage 
    a = checkma(a)
    if erode:
        a = malib.mask_islands(a, iterations)
    #bmask = (~a.mask)
    bmask = (~np.ma.getmaskarray(a))
    bmask_filled = ndimage.morphology.binary_fill_holes(bmask)
    #This will create valid values with a.filled in the original ma
    #a_erode.mask[:] = ~bmask_filled
    #return a_erode
    return ~bmask_filled

#This is an alternative to the ma.notmasked_edges
#Note: probably faster/simpler to contour the mask
def contour_edges(a):
    import matplotlib.pyplot as plt
    a = checkma(a)
    #Contour nodata value
    levels = [a.fill_value]
    #kw = {'antialiased':True, 'colors':'r', 'linestyles':'-'}
    kw = {'antialiased':True}
    #Generate contours around nodata
    cs = plt.contour(a.filled(), levels, **kw)
    #This returns a list of numpy arrays
    #allpts = np.vstack(cs.allsegs[0])
    #Extract paths
    p = cs.collections[0].get_paths()
    #Sort by number of vertices in each path
    p_len = [i.vertices.shape[0] for i in p]
    p_sort = [x for (y,x) in sorted(zip(p_len,p), reverse=True)]    
    #cp = p[0].make_compound_path(*p)
    return p_sort

#Fill masked areas with random noise
#This is needed for any fft-based operations
def randomfill(a):
    a = checkma(a)
    #For data that have already been normalized,
    #This provides a proper normal distribution with mean=0 and std=1
    #a = (a - a.mean()) / a.std()
    #noise = a.mask * (np.random.randn(*a.shape))
    noise = a.mask * np.random.normal(a.mean(), a.std(), a.shape)
    #Add the noise
    b = a.filled(0) + noise
    return b

#Wrapper for functions that can't handle ma (e.g. scipy.ndimage)
#Substitutes np.nan
#This will force filters to ignore nan, but causes adjacent pixels to be set to nan as well
#http://projects.scipy.org/scipy/ticket/1155 
def nanfill(a, f_a, *args, **kwargs):
    a = checkma(a)
    ndv = a.fill_value  
    #Note: The following fails for arrays that are not float (np.nan is float)
    b = f_a(a.filled(np.nan), *args, **kwargs)
    #the fix_invalid fill_value parameter doesn't seem to work
    out = np.ma.fix_invalid(b, copy=False)
    out.set_fill_value(ndv)
    return out

#Compute median absolute difference
def mad(a, c=0.6745):
    a = checkma(a)
    return np.ma.median(np.fabs(a - np.ma.median(a))) / c

#Percentile values
def calcperc(b, perc=(0.1,99.9)):
    from scipy.stats import scoreatpercentile
    b = checkma(b)
    if b.count() > 0:
        low = scoreatpercentile(b.compressed(), perc[0])
        high = scoreatpercentile(b.compressed(), perc[1])
    else:
        low = 0
        high = 0

    #low = scipy.stats.mstats.scoreatpercentile(b, perc[0])
    #high = scipy.stats.mstats.scoreatpercentile(b, perc[1])

    #This approach can be used for unmasked array, but values less than 0 are problematic
    #bma_low = b.min()
    #bma_high = b.max()
    #low = scipy.stats.scoreatpercentile(b.data.flatten(), perc[0], (bma_low, bma_high))
    #high = scipy.stats.scoreatpercentile(b.data.flatten(), perc[1], (bma_low, bma_high))
    return low, high

def iqr(b, perc=(25, 75)):
    b = checkma(b)
    low, high = calcperc(b, perc)
    return low, high, high - low

def robust_spread(b, perc=(16,84)):
    p16, p84 = calcperc(b, perc)
    spread = np.abs((p84 - p16)/2)
    return p16, p84, spread

def robust_spread_idx(b, sigma=3):
    b = checkma(b)
    med = np.ma.median(b)
    p16, p84, spread = robust_spread(b)
    min = med - sigma*spread
    max = med + sigma*spread
    b_idx = (b > min) & (b < max)
    return b_idx

def robust_spread_fltr(b, sigma=3):
    #Should compare runtime w/ the ma solution
    b_idx = robust_spread_idx(b, sigma)
    newmask = np.ones_like(b, dtype=bool)
    newmask[b_idx] = False 
    return np.ma.array(b, mask=newmask)

#Note: need to specify dtype for mean/std calculations 
#full and fast give very different results
#Want to add 25/75 IQR - scipy.stats.mstats.scoreatpercentile
#Convert output to dictionary
#Need to convert stats to float before json.dumps
#The a.mean(dtype='float64') is needed for accuracte calculation
def print_stats(a, full=False):
    from scipy.stats.mstats import mode 
    a = checkma(a)
    thresh = 4E6
    if full or a.count() < thresh:
        q = (iqr(a))
        p16, p84, spread = robust_spread(a)
        #There has to be a better way to compute the mode for a ma
        #mstats.mode returns tuple of (array[mode], array[count])
        a_mode = float(mode(a, axis=None)[0])
        stats = (a.count(), a.min(), a.max(), a.mean(dtype='float64'), a.std(dtype='float64'), np.ma.median(a), mad(a), q[0], q[1], q[2], a_mode, p16, p84, spread) 
    else:
        ac = a.compressed()
        stride = np.around(ac.size / thresh)
        ac = np.ma.array(ac[::stride])
        #idx = np.random.permutation(ac.size)
        #Note: need the ma cast here b/c of a.count() below
        #ac = np.ma.array(ac[idx[::stride]])
        q = (iqr(ac))
        p16, p84, spread = robust_spread(ac)
        ac_mode = float(mode(ac, axis=None)[0])
        stats = (a.count(), a.min(), a.max(), a.mean(dtype='float64'), a.std(dtype='float64'), np.ma.median(ac), mad(ac), q[0], q[1], q[2], ac_mode, p16, p84, spread) 
    print "count: %i min: %0.2f max: %0.2f mean: %0.2f std: %0.2f med: %0.2f mad: %0.2f q1: %0.2f q2: %0.2f iqr: %0.2f mode: %0.2f p16: %0.2f p84: %0.2f spread: %0.2f" % stats
    return stats

def rmse(a):
    ac = checkma(a).compressed()
    rmse = np.sqrt(np.sum(a**2)/a.size)
    return rmse

#Replace nodata value in GDAL band
def replace_ndv(b, new_ndv):
    b_ndv = get_ndv_b(b)    
    bma = np.ma.masked_equal(b.ReadAsArray(), b_ndv)
    bma.set_fill_value(new_ndv)
    b.WriteArray(bma.filled())
    b.SetNoDataValue(new_ndv)
    return b

def set_ndv(dst_fn, ndv):
    dst_ds = gdal.Open(dst_fn, gdal.GA_Update)
    for n in range(1, dst_ds.RasterCount+1):
        b = dst_ds.GetRasterBand(1)
        b.SetNoDataValue(ndv)
    dst_ds = None

#Should overload these functions to handle fn, ds, or b
#Perhaps abstract, as many functions will need this functionality
def get_ndv_fn(fn):
    ds = gdal.Open(fn, gdal.GA_ReadOnly)
    return get_ndv_ds(ds)

#Want to modify to handle multi-band images and return list of ndv
def get_ndv_ds(ds, bnum=1):
    b = ds.GetRasterBand(bnum)
    return get_ndv_b(b)

#Return nodata value for GDAL band
def get_ndv_b(b):
    b_ndv = b.GetNoDataValue()
    if b_ndv is None:
        #Check ul pixel for ndv
        ns = b.XSize
        nl = b.YSize
        ul = float(b.ReadAsArray(0, 0, 1, 1))
        ur = float(b.ReadAsArray(ns-1, 0, 1, 1))
        lr = float(b.ReadAsArray(ns-1, nl-1, 1, 1))
        ll = float(b.ReadAsArray(0, nl-1, 1, 1))
        #Probably better to use 3/4 corner criterion
        #if ul == ur == lr == ll:
        if np.isnan(ul) or ul == lr:
            b_ndv = ul
        else:
            #Assume ndv is 0
            b_ndv = 0
    return b_ndv

#Want to define methods to load ma from OpenCV, PIL, etc.
#These formats should be directly readable as np arrays

#Note: want to modify this to import all bands as separate arrays in ndarray
#Unless the user requests a single band, or range of bands

#Given input filename, return a masked array for specified band
def fn_getma(fn, bnum=1):
    #Add check for filename existence
    ds = gdal.Open(fn, gdal.GA_ReadOnly)    
    return ds_getma(ds, bnum=bnum)

#Given input dataset, return a masked array for the input band
def ds_getma(ds, bnum=1):
    b = ds.GetRasterBand(bnum)
    return b_getma(b)

#Given input band, return a masked array
def b_getma(b):
    b_ndv = get_ndv_b(b)
    bma = np.ma.masked_equal(b.ReadAsArray(), b_ndv)
    bma = checkma(bma)
    return bma

def get_sub_dim(src_ds, scale=None, maxdim=1024.):
    ns = src_ds.RasterXSize
    nl = src_ds.RasterYSize
    if scale is None:
        scale_ns = ns/maxdim
        scale_nl = nl/maxdim
        scale = max(scale_ns, scale_nl)
    #Need to check to make sure scale is positive real 
    if scale > 1:
        ns = int(round(ns/scale))
        nl = int(round(nl/scale))
    return ns, nl

#Load a subsampled array
#Can specify scale factor or max dimension
#No need to load the entire dataset for stats computation
def gdal_getma_sub(src_ds, bnum=1, scale=None, maxdim=1024.):    
    #print src_ds.GetFileList()[0]
    b = src_ds.GetRasterBand(bnum)
    b_ndv = get_ndv_b(b)
    ns, nl = get_sub_dim(src_ds, scale, maxdim)
    #The buf_size parameters determine the final array dimensions
    b_array = b.ReadAsArray(buf_xsize=ns, buf_ysize=nl)
    bma = np.ma.masked_equal(b_array, b_ndv)
    bma = checkma(bma)
    return bma

#Note: need to consolidate with warplib.writeout (takes ds, not ma)
#Add option to build overviews when writing GTiff
#Input proj must be WKT
def writeGTiff(a, dst_fn, src_ds=None, bnum=1, ndv=None, gt=None, proj=None, create=False, sparse=False):
    #If input is not np.ma, this creates a new ma, which has default filL_value of 1E20
    #Must manually override with ndv
    a = checkma(a, fix=False)
    #Want to preserve fill_value if already specified
    if ndv is not None:
        a.set_fill_value(ndv)
    format = "GTiff"
    driver = gdal.GetDriverByName(format)
    #Currently only support writing singleband rasters
    #if a.ndim > 2:
    #   np_nbands = a.shape[2]
    #   if src_ds.RasterCount np_nbands: 
    #      for bnum in np_nbands:
    nbands = 1
    np_dt = a.dtype.name
    if src_ds is not None:
        src_dt = gdal.GetDataTypeName(src_ds.GetRasterBand(bnum).DataType)
        src_gt = src_ds.GetGeoTransform()
        #This is WKT
        src_proj = src_ds.GetProjection()
        src_srs = osr.SpatialReference()  
        src_srs.ImportFromWkt(src_ds.GetProjectionRef())

    #Probably a cleaner way to handle this
    if gt is None:
        gt = src_gt
    if proj is None:
        proj = src_proj

    #Need to create a new copy of the default options
    opt = list(gdal_opt)
    
    #Note: packbits is better for sparse data
    if sparse:
        opt.remove('COMPRESS=LZW')
        opt.append('COMPRESS=PACKBITS')
        #Not sure if VW can handle sparse tif
        #opt.append('SPARSE_OK=TRUE')

    #Use predictor=3 for floating point data
    if 'float' in np_dt.lower() and 'COMPRESS=LZW' in opt: 
        opt.append('PREDICTOR=3')

    #If input ma is same as src_ds, write out array using CreateCopy from existing dataset
    #if not create and (src_ds is not None) and ((a.shape[0] == src_ds.RasterYSize) and (a.shape[1] == src_ds.RasterXSize) and (np_dt.lower() == src_dt.lower())): 
    #Should compare srs.IsSame(src_srs)
    if not create and (src_ds is not None) and ((a.shape[0] == src_ds.RasterYSize) and (a.shape[1] == src_ds.RasterXSize) and (np_dt.lower() == src_dt.lower())) and (src_gt == gt) and (src_proj == proj):
        #Note: third option is strict flag, set to false
        dst_ds = driver.CreateCopy(dst_fn, src_ds, 0, options=opt)
    #Otherwise, use Create
    else:
        from osgeo import gdal_array
        dt_dict = gdal_array.codes        
        a_dt = a.dtype
        if a_dt.name == 'int8':
            gdal_dt = 1
        elif a_dt.name == 'bool':
            #Write out as Byte
            gdal_dt = 1 
            #Set ndv to 0
            a.fill_value = False
            opt.remove('COMPRESS=LZW')
            opt.append('COMPRESS=DEFLATE')
            #opt.append('NBITS=1')
        else:
            gdal_dt = dt_dict.keys()[dt_dict.values().index(a.dtype)]
        #Create(fn, nx, ny, nbands, dtype, opt)
        dst_ds = driver.Create(dst_fn, a.shape[1], a.shape[0], nbands, gdal_dt, options=opt)
        #Note: Need GeoMA here to make this work, or accept gt as argument
        #Could also do ds creation in calling script
        if gt is not None:
            dst_ds.SetGeoTransform(gt)
        if proj is not None:
            dst_ds.SetProjection(proj)
    
    dst_ds.GetRasterBand(bnum).WriteArray(a.filled())
    dst_ds.GetRasterBand(bnum).SetNoDataValue(float(a.fill_value))
    dst_ds = None

#Check that input is a masked array
def checkma(a, fix=True):
    #isinstance(a, np.ma.MaskedArray)
    if np.ma.is_masked(a):
        out=a
    else:
        out=np.ma.array(a)
    #Fix invalid values
    #This is not necessarily desirable for writing
    if fix:
        #Note: this fails for datetime arrays! Treated as objects.
        #Note: datetime ma returns '?' for fill value
        from datetime import datetime
        if isinstance(a[0], datetime):
            print "Input array appears to be datetime.  Skipping fix"
        else:
            out=np.ma.fix_invalid(out, copy=False)
    return out

#Image view of ma
#Should probably move this to imview.py
def iv(b, **kwargs):
    import matplotlib.pyplot as plt
    import imview 
    b = checkma(b)
    #if hasattr(kwargs,'imshow_kwargs'):
    #    kwargs['imshow_kwargs']['interpolation'] = 'bicubic'
    #else:
    #    kwargs['imshow_kwargs'] = {'interpolation': 'bicubic'}
    #bma_fig(fig, bma, cmap='gist_rainbow_r', clim=None, bg=None, n_subplt=1, subplt=1, label=None, **imshow_kwargs)
    fig = plt.figure()
    imview.bma_fig(fig, b)
    plt.show()
    return fig

#Simple gridding
#x, y, z = malib.get_xyz(a)
#xi, yi = np.indices(a_fp.shape[::-1])
#zi = sp.interpolate.griddata((x,y), z, (xi,yi), method='cubic').T
#zi = np.ma.fix_invalid(zi, fill_value=a.fill_value)

def get_xyz(a):
    a = checkma(a)
    y, x = np.indices(a.shape)
    y = np.ma.array(y, mask=np.ma.getmaskarray(a))
    x = np.ma.array(x, mask=np.ma.getmaskarray(a))
    return np.array([x.compressed(), y.compressed(), a.compressed()])

#Efficient solution for 1D
#Note, doesn't include edges
#http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html
def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

#This provides non-overlapping blocks of given size
#Could probably be modified for a stride of (1,1) to overlap
#From http://stackoverflow.com/questions/5073767/how-can-i-efficiently-process-a-numpy-array-in-blocks-similar-to-matlabs-blkpro
def block_view(A, block=(3, 3)):
    from numpy.lib.stride_tricks import as_strided as ast
    """Provide a 2D block view to 2D array. No error checking made.
    Therefore meaningful (as implemented) only for blocks strictly
    compatible with the shape of A."""
    # simple shape and strides computations may seem at first strange
    # unless one is able to recognize the 'tuple additions' involved ;-)
    shape= (A.shape[0]/ block[0], A.shape[1]/ block[1])+ block
    strides= (block[0]* A.strides[0], block[1]* A.strides[1])+ A.strides
    return ast(A, shape= shape, strides= strides)

def sliding_window_padded(a, ws, ss=(1,1), flatten=True):
    colpad = ws[0]/2
    col_a = np.empty((a.shape[0],colpad))
    col_a[:] = np.nan
    a = np.column_stack([col_a, a, col_a])
    rowpad = ws[1]/2
    row_a = np.empty((rowpad, a.shape[1]))
    row_a[:] = np.nan
    a = np.row_stack([row_a, a, row_a])
    return sliding_window(a, ws, ss, flatten) 

#From http://www.johnvinyard.com/blog/?p=268 
def sliding_window(a, ws, ss=None, flatten=True):
    from numpy.lib.stride_tricks import as_strided as ast
    '''
    Return a sliding window over a in any number of dimensions
     
    Parameters:
        a  - an n-dimensional numpy array
        ws - an int (a is 1D) or tuple (a is 2D or greater) representing the size
             of each dimension of the window
        ss - an int (a is 1D) or tuple (a is 2D or greater) representing the
             amount to slide the window in each dimension. If not specified, it
             defaults to ws.
        flatten - if True, all slices are flattened, otherwise, there is an
                  extra dimension for each dimension of the input.
     
    Returns
        an array containing each n-dimensional window from a
    '''
     
    if None is ss:
        # ss was not provided. the windows will not overlap in any direction.
        ss = ws
    ws = norm_shape(ws)
    ss = norm_shape(ss)
     
    # convert ws, ss, and a.shape to numpy arrays so that we can do math in every
    # dimension at once.
    ws = np.array(ws)
    ss = np.array(ss)
    shape = np.array(a.shape)
     
    # ensure that ws, ss, and a.shape all have the same number of dimensions
    ls = [len(shape),len(ws),len(ss)]
    if 1 != len(set(ls)):
        raise ValueError(\
        'a.shape, ws and ss must all have the same length. They were %s' % str(ls))
     
    # ensure that ws is smaller than a in every dimension
    if np.any(ws > shape):
        raise ValueError(\
        'ws cannot be larger than a in any dimension.\
        a.shape was %s and ws was %s' % (str(a.shape),str(ws)))
     
    # how many slices will there be in each dimension?
    newshape = norm_shape(((shape - ws) // ss) + 1)
    # the shape of the strided array will be the number of slices in each dimension
    # plus the shape of the window (tuple addition)
    newshape += norm_shape(ws)
    # the strides tuple will be the array's strides multiplied by step size, plus
    # the array's strides (tuple addition)
    newstrides = norm_shape(np.array(a.strides) * ss) + a.strides
    strided = ast(a,shape = newshape,strides = newstrides)
    if not flatten:
        return strided
     
    # Collapse strided so that it has one more dimension than the window.  I.e.,
    # the new array is a flat list of slices.
    meat = len(ws) if ws.shape else 0
    firstdim = (np.product(newshape[:-meat]),) if ws.shape else ()
    dim = firstdim + (newshape[-meat:])
    # remove any dimensions with size 1
    dim = filter(lambda i : i != 1,dim)
    return strided.reshape(dim)

def norm_shape(shape):
    '''
    Normalize numpy array shapes so they're always expressed as a tuple,
    even for one-dimensional shapes.
     
    Parameters
        shape - an int, or a tuple of ints
     
    Returns
        a shape tuple
    '''
    try:
        i = int(shape)
        return (i,)
    except TypeError:
        # shape was not a number
        pass
 
    try:
        t = tuple(shape)
        return t
    except TypeError:
        # shape was not iterable
        pass
     
    raise TypeError('shape must be an int, or a tuple of ints')
