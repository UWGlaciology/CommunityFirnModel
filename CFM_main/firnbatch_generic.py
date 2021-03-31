#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
This is what I use to batch CFM runs (in combination with GNU parallel).

But, it can also be used a standalone access to the CFM. Previously (and still)

the CFM is called using main.py and a .json file. 


This file takes a different direction - its inputs from the command line are the
lat/lon pair you want to model. Then it calls spin_up_generator_CFM to make 
forcing files and runs the CFM.

Model-specific parameters (i.e. things in the .json file) are changed here, so
you can take a (semi) generic .json file and set everything in here.

to run:
>>> python firnbatch_generic.py '72.5 -38.5'

or to batch some runs, e.g. make a text file called CFMsites.txt that has 3 lines:

72.5 -38.5
75.3 -34.8
77.8 -43.2

Then, you can run:
>>> parallel -a CFMsites.txt python firnbatch_generic.py

(you need to have GNU parallel installed)

@author: maxstev
'''

import os
import numpy as np
import sys
# from shutil import copyfile
import shutil
import sys
import subprocess
# from string import join
from firn_density_spin import FirnDensitySpin
from firn_density_nospin import FirnDensityNoSpin
import time
import json
# from multiprocessing import Pool
from functools import partial
from siteClimate_from_RCM import getClimate
from RCMpkl_to_spin import makeSpinFiles
import xarray as xr
import pandas as pd

def run_CFM(LLpair, json_base, timeres = '1D', Tinterp = 'mean', MELT= True, runtype = 'local', datatype = 'MERRA', movefiles=False,RCdrive = None):

    '''
    Run the CFM for a lat/lon pair. This is a bit of a hack because there are
    many things to go through and customize for a particular model run.

    It can optionally transfer your results using rclone.

    Parameters
    ----------
    LLpair: string with lat lon pair
        Latitude, Longitude (Both in Degrees)
        LLpair is a string like '72.5 -38.5', which is then parsed into the 
        constitutive parts.
    timeres: pandas Timedelta (string)
        Resampling frequency, e.g. '1D' is 1 day; '1M' for 1 month.
    Tinterp: 'mean', 'effective', or 'weighted'
        how to resample the temperature; mean is regular mean, 'effective' is 
        Arrhenius mean; 'weighted' is accumulation-weighted mean
    RCdrive: string
        name of the rclone drive to move files to
    runtype: 'local' or 'remote'
        whether you are running locally or remotely (you can customize these)
    datatype: 'MAR' or 'MERRA'
        which RCM product you want to use
    movefiles: boolean
        whether or not you want to move the CFM outputs somewhere using rclone
    json_base: string
        path and name of the 'base' .json that will be used for the CFM run;
        it will be edited a bit by this script.
    '''

    tnow = time.time()
    npa = np.fromstring(LLpair,dtype =float, sep=' ')   
    
    lat_int = npa[0]
    lon_int = npa[1]

    print('lat_int: ',lat_int)
    print('lon_int: ',lon_int)

    if datatype == 'MAR':
        dsource = 'ERA6k' # [ERA10k, ERA6k, 'NCEP20k']
        dwriter = datatype+'_'+dsource
    else:
        dsource = None
        dwriter = datatype

    df_CLIM = getClimate(lat_int,lon_int, writer = True, runtype = runtype, datatype=datatype, melt=MELT, dsource = dsource)

    Cd, StpsPerYr, depth_S1, depth_S2, grid_bottom = makeSpinFiles(df_CLIM,timeres='1D',Tinterp='mean',spin_date_st = 1980.0, spin_date_end = 1995.0)

    # Cd, StpsPerYr, depth_S1, depth_S2, grid_bottom = makeSpinFiles(lat_int,lon_int,timeres=timeres,writer = True,Tinterp = Tinterp, runtype = runtype,datatype=datatype,dsource = dsource)


    writeCFMinputs = True
    if writeCFMinputs:
        # This is a bit redundant from spin up generator, but that pickle is the raw data from the RCM;
        # This is the resampled (i.e. what gets fed to CFM)
        CFM_df = pd.DataFrame(Cd).set_index('time')

        if os.path.exists('CFMpickle/'):
            pass
        else:
            os.makedirs('CFMpickle/')
        CFM_df.to_pickle('CFMpickle/CFM_df_{}_{}_{}_{}_{}.pkl'.format(lat_int,lon_int,timeres,Tinterp,dwriter))


    print('grid_bottom: ', grid_bottom)
    telap = (time.time()-tnow)/60
    print('part 1 done, {} minutes'.format(telap))

    gridtype = 'Double' # If you want to use the regrid feature to save computing time

    if os.path.exists('json/'):
        pass
    else:
        os.makedirs('json/')

    connm='json/example_{}_{}_{}_{}_{}.json'.format(lat_int,lon_int,timeres,Tinterp,dwriter)
    shutil.copyfile(json_base, connm)

    jsonFile = open(connm, "r")
    data = json.load(jsonFile)
    jsonFile.close()

    ### Configure these as much as you would like.
    re = 'CFMresults/example_{}_{}_{}_{}_{}'.format(lat_int,lon_int,timeres,Tinterp,dwriter)
    print('results folder:', re)
    data["resultsFolder"] = re
    data["stpsPerYear"] = float('%.2f' % (StpsPerYr))
    data["stpsPerYearSpin"] = float('%.2f' % (StpsPerYr))
    data["grid1bottom"] = float('%.1f' %(depth_S1))
    data["grid2bottom"] = float('%.1f' %(depth_S2))
    data["HbaseSpin"] = float('%.1f' %(3000 - grid_bottom))
    if gridtype == "Single":
        data["doublegrid"] = False
    elif gridtype == "Double":
        data["doublegrid"] = True
    if np.mean(Cd['BDOT'])<0.05:
        data["nodestocombine"] = 60
        data["multnodestocombine"] = 24
    else:        
        data["nodestocombine"] = 30
        data["multnodestocombine"] = 12
    data["conductivity"] = 'Calonne2019'
    data["DIPhorizon"] = np.floor(0.8*grid_bottom)

    if 'SMELT' in Cd.keys(): #snowmelt; use SMELT rather than MELT b/c pandas has a function called melt
        data['MELT'] = True
        data['liquid'] = 'bucketVV'
    if 'RAIN' in Cd.keys():
        data['RAIN'] = True
    # if 'subl' in Cd.keys():
    #     data[]

    jsonFile = open(connm, "w+")
    jsonFile.write(json.dumps(data,sort_keys=True, indent=4, separators=(',', ': ')))
    jsonFile.close()
    configName = connm

    firn = FirnDensityNoSpin(configName,climateTS = Cd, NewSpin = True)
    firn.time_evolve()
    telap = (time.time()-tnow)/60
    print('main done, {} minutes'.format(telap))

    shutil.move(configName,re)

    if movefiles:
        subprocess.call(['rclone','move', re, RCdrive +'/'+re])
        try:
            os.rmdir(re)
        except Exception:
            pass
    
        
if __name__ == '__main__':
    gridlist = sys.argv[1]
    #############################
    ### THINGS TO CHANGE HERE ###
    ### See docstring ###########
    timeres = '5D' 
    Tinterp = 'mean' # [mean, effective, weighted]
    runtype = 'local' 
    datatype = 'MAR' 
    movefiles = False 
    RCdrive = 'drive_name:' 
    json_base = 'example.json' 
    MELT = True
    #############################
    ### THINGS TO CHANGE HERE ###
    
    run_CFM(gridlist, json_base, timeres = timeres, Tinterp = Tinterp, MELT= MELT, runtype = runtype, datatype = datatype, movefiles=movefiles,RCdrive=RCdrive)

