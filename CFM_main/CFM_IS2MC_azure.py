#!/usr/bin/env python

import h5py
import numpy as np 
import scipy
from scipy.spatial import cKDTree
import pandas as pd
import xarray as xr
# import netCDF4
import glob
# import boto3
# import botocore
import io
from pathlib import Path
# import s3fs
# import fsspec

"""
CFM_hpc.py
=======
This file configures and runs the CFM.
- On SMCE, it goes in CFM_main 
- input climate data come from a zarr with daily MERRA-2
- (the merra-2 data are  preprocessed for this purpose)
- this script creates a .json configuration file
- CFM runs with the climate input and config file, results get put
  in the specified directory.

to test on demand: 
srun -N1 -n1 -c1 --exclusive --partition=hpc-demand-36 singularity run -B /efs/maxstev/CFM/CommunityFirnModel/CFM_main,/efs/maxstev/CFMresultsSNOWPACK_varrho,/efs/maxstev/ ~/containers/ilab-cfm-1.1.0.sif python /efs/maxstev/CFM/CommunityFirnModel/CFM_main/CFM_hpc_SEB_zarr.py 589 -9999

run this script using:
>>>python CFM_hpc.py "LAT,LON" 
where LAT,LON are your latitude and longitude, e.g. "72.5,-36.5"

"""

import sys
import os

### make sure to edit this to ensure CFM is in pypath
sys.path.insert(0, 'CommunityFirnModel/CFM_main')

# from firn_density_spin import FirnDensitySpin
from firn_density_nospin import FirnDensityNoSpin
import time
import json
import shutil
import RCMpkl_to_spin as RCM

def MERRA2_zarr_to_dataframe(y_int,x_int,zarr_source='blob'):
    '''
    Create a pandas dataframe for a site in Greenland
    returns:
    df_daily: a dataframe holding all of the forcing fields needed for CFM run
    
    input mass fluxes from zarr are in kg/m2/s.
    mass fluxes in df_daily are in kg/m2/timestep
    energy fluxes are W/m2
    
    in MERRA, ?positive? sublimation flux = sublimation, ?negative? = deposition
    '''
    def make_dataframe(dsZ,y_int,x_int):
            y_ll = dsZ.y.data
            x_ll = dsZ.x.data
            ii, y_val = min(enumerate(y_ll), key=lambda x: abs(x[1]-y_int))
            jj, x_val = min(enumerate(x_ll), key=lambda x: abs(x[1]-x_int))

            df_sub = pd.DataFrame(index=dsZ.time.values)

            varlist = [jj for jj in dsZ.variables if jj not in ['time','y','x']]
            
            if 'EMIS_eff' in varlist:
                df_sub['EMIS_eff'] = dsZ.isel(y=[ii],x=[jj])['EMIS_eff'].values.flatten()
            else:
                for vv in varlist:
                        df_sub[vv] = dsZ.isel(y=[ii],x=[jj])[vv].values.flatten()

                df_seconds = df_sub.index.diff().mean().total_seconds()

                df_sub['RAIN'] = (df_sub['PRECLS'] + df_sub['PRECCU']) * df_seconds
                df_sub['EVAP'] = -1 * df_sub['EVAP'] * df_seconds # multiply by -1 because of MERRA2 sign convention
                df_sub['PRECSN'] = df_sub['PRECSN'] * df_seconds
                df_sub['SMELT'] = df_sub['SMELT'] * df_seconds

                drn = {'T2M_i':'T2m','TS_i':'TSKIN','EVAP':'SUBLIM','HFLUX':'QH','EFLUX':'QL','SWGDN':'SW_d','LWGAB':'LW_d_M2','RAIN':'RAIN','PRECSN':'BDOT','ALBEDO':'ALBEDO_i','SMELT':'SMELT','SWGNT':'SW_n'}

                df_sub = df_sub[drn.keys()]
                df_sub.rename(mapper=drn,axis=1,inplace=True)

            return ii,jj,y_val,x_val,df_sub
        
    decades = [1980,1990,2000,2010,2020]
    df_dict = {}
    for decade in decades:
        if zarr_source=='discover':
            filename = f"/discover/nobackup/cdsteve2/climate/MERRA2/GrIS_IS2mc/zip/M2_GrIS_4h_IS2mc_{decade}.zarr.zip"
            fn_EE = '/discover/nobackup/cdsteve2/climate/MERRA2/GrIS_IS2mc/MERRA2_GrIS_4h_eff_emis_ALL.nc'
        elif zarr_source=='azure':
            filename = f'/mnt/firnadls/M2_GrIS_4h_IS2mc_{decade}.zarr.zip'
            fn_EE = '/mnt/firnadls/MERRA2_GrIS_4h_eff_emis_ALL.nc'

        with xr.open_dataset(filename,engine='zarr') as dsZ:
            ii,jj,y_val,x_val,df_sub= make_dataframe(dsZ,y_int,x_int)
            df_dict[decade] = df_sub

    with xr.open_dataset(fn_EE) as fE:
        _ii,_jj,_y_val,_x_val,df_EE = make_dataframe(fE,y_int,x_int)        
            
    df_out = pd.concat(df_dict.values())

    df_out['QL'] = -1 * df_out['QL']
    df_out['QH'] = -1 * df_out['QH']

    sigma = 5.670374419e-8
    df_out['LW_d_EE'] = df_EE['EMIS_eff'] * sigma * df_out['T2m']**4
    df_out['ALBEDO_post'] = 1 - (df_out['SW_n']/df_out['SW_d'])

    return ii,jj,y_val,x_val,df_out

#################################
#################################
#################################

if __name__ == '__main__':

    tic=time.time()

    ### sys.argv[1] is the run number, which corresponds to a x/y pair in on the IS2 grid
    ### DYE2:
    ### x_int = -54042
    ### y_int = -2579982
    ### point in IS2_icepixels.csv: 6304
    
    runloc = 'discover'
    seb = True
    LWdown_source = 'EMIS_eff'
    ALBEDO_source = 'post'
     
    if runloc =='discover':
        CFM_path = Path('/discover/nobackup/cdsteve2/ATL_masschange/CommunityFirnModel/CFM_main')
    elif runloc == 'azure':
        CFM_path = Path('/shared/home/cdsteve2/CommunityFirnModel/CFM_main/')

    config_in = Path(CFM_path,'HPC_config_default.json')

    with open(config_in, "r") as f:
        jsonString      = f.read()
        c          = json.loads(jsonString) 

    c['runloc'] = runloc

    if c['runloc'] == 'azure':
        zarr_source = 'azure'
        ll_list = np.genfromtxt(Path(CFM_path,'IS2_icepixels.csv'),delimiter=',',skip_header=1)
    
    elif c['runloc'] == 'discover':
        zarr_source = 'discover'
        ll_list = np.genfromtxt(Path(CFM_path,'IS2_icepixels.csv'),delimiter=',',skip_header=1)
    
    if c['runloc']=='local':
        x_int = c['x_val']
        y_int = c['y_val']
    else:
        print(f'sys.argv[1]: {sys.argv[1]}')
        dkey = int(sys.argv[1]) # this is the array value
        x_int = float(ll_list[dkey][0])
        y_int = float(ll_list[dkey][1])

    if np.isnan(y_int):
        print('y_int is nan')
        sys.exit()

    try:
        runid = float(sys.argv[1])
    except:
        runid=-9999
        
    ### Get climate data from zarr
    ii,jj,y_val,x_val,df_daily = MERRA2_zarr_to_dataframe(y_int,x_int,zarr_source=zarr_source)
    sds = 1980.0 #spin date start
    sde = 1995.0 #spin date end
    y_w = y_val
    x_w = x_val

    if LWdown_source == 'MERRA2':
        df_daily = df_daily.rename({'LW_d_M2':'LW_d'},axis=1).drop(['LW_d_EE'],axis=1)
    # elif LWdown_source == 'EMIS_eff':
    else:
        df_daily = df_daily.rename({'LW_d_EE':'LW_d'},axis=1).drop(['LW_d_M2'],axis=1)

    if ALBEDO_source == 'M2_interp':
        df_daily = df_daily.rename({'ALBEDO_i':'ALBEDO'},axis=1).drop(['ALBEDO_post','SW_n'],axis=1)
    else:
        df_daily = df_daily.rename({'ALBEDO_post':'ALBEDO'},axis=1).drop(['ALBEDO_i','SW_n'],axis=1)

    print(ii, jj, y_val, x_val)
    print(df_daily.head())
    #######

    c['SEB'] = True
    calc_melt = False
    c['MELT'] = True

    rf_po = f'/CFMresults_{dkey}' #results path

    if runloc == 'azure':
        c['resultspath'] = '/mnt/firnadls/CFM_outputs'
    elif runloc == 'discover':
        c['resultspath'] = '/discover/nobackup/cdsteve2/ATL_masschange/CFMoutputs'

    # c['resultsFolder'] = c['resultspath'] + c['results_ext'] + rf_po
    c['resultsFolder'] = str(Path(c['resultspath'], rf_po))
    
    c['y_int'] = float(y_int)
    c['x_int'] = float(x_int)
    c['y_val'] = float(y_val)
    c['x_val'] = float(x_val)
    c['runid'] = runid
    ##########

    climateTS, StpsPerYr, depth_S1, depth_S2, grid_bottom, SEBfluxes = (
        RCM.makeSpinFiles(df_daily,timeres=c['DFresample'],Tinterp='mean',spin_date_st = sds, 
        spin_date_end = sde,melt=c['MELT'],desired_depth = None,SEB=c['SEB'],rho_bottom=916,calc_melt=calc_melt))

    c["stpsPerYear"] = float('%.2f' % (StpsPerYr))
    c["stpsPerYearSpin"] = float('%.2f' % (StpsPerYr))
    c["grid1bottom"] = float('%.1f' %(depth_S1))
    c["grid2bottom"] = float('%.1f' %(depth_S2))
    c["HbaseSpin"] = float('%.1f' %(3000 - grid_bottom))
    
    c["NewSpin"] = True

    configName = f'CFMconfig_{y_w}_{x_w}.json'
    shutil.copyfile(config_in, configName)
    
    if os.path.isfile(os.path.join(c['resultsFolder'],configName)):
        CFMconfig = os.path.join(c['resultsFolder'],configName)
        if os.path.isfile(os.path.join(os.getcwd(), configName)):
            os.remove(os.path.join(os.getcwd(), configName))
        shutil.move(CFMconfig, os.getcwd())
    else:
        CFMconfig = configName     
    
    with open(CFMconfig,'w+') as fp:
        fp.write(json.dumps(c,sort_keys=True, indent=4, separators=(',', ': ')))

    if 'NewSpin' in c:
        NewSpin = c['NewSpin']
    else:
        NewSpin = True

    ### Create CFM instance by passing config file and forcing data, then run the model
    firn = FirnDensityNoSpin(CFMconfig, climateTS = climateTS, NewSpin = NewSpin, SEBfluxes = SEBfluxes)
    firn.time_evolve()
    ###

    shutil.move(configName,os.path.join(c['resultsFolder'],configName))

    print ("output folder = ", c['resultsFolder'])
    print('run time =' , time.time()-tic , 'seconds')