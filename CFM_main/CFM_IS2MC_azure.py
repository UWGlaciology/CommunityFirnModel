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
import socket

"""
CFM_IS2MC_azure.py
=======
This file configures and runs the CFM.
- It goes in CFM_main (changing this...) 
- input climate data come from a zarr with 4h MERRA-2
- (the merra-2 data are  preprocessed for this purpose)
- this script creates a .json configuration file
- CFM runs with the climate input and config file, results get put
  in the specified directory.

to test on demand: NEED TO UPDATE THIS
srun -N1 -n1 -c1 --exclusive --partition=hpc-demand-36 singularity run -B /efs/maxstev/CFM/CommunityFirnModel/CFM_main,/efs/maxstev/CFMresultsSNOWPACK_varrho,/efs/maxstev/ ~/containers/ilab-cfm-1.1.0.sif python /efs/maxstev/CFM/CommunityFirnModel/CFM_main/CFM_hpc_SEB_zarr.py 589 -9999

run this script using:
>>>python CFM_IS2MC_azure.py X (where X is an integer) 
"""

import sys
import os

hh = socket.gethostname()
global runloc
if 'disc' in hh:
    runloc = 'discover'
elif 'gs615-meltwater' in hh:
    runloc = 'local'
else:
    runloc = 'azure'
print(runloc)

### make sure to edit this to ensure CFM is in pypath
if runloc=='discover':
    sys.path.insert(0, '/discover/nobackup/cdsteve2/ATL_masschange/CommunityFirnModel/CFM_main')
elif runloc=='azure':
    sys.path.insert(0, '/shared/home/cdsteve2/CommunityFirnModel/CFM_main')
elif runloc=='local':
    sys.path.insert(0, '/Users/cdsteve2/research/firn/CommunityFirnModel/CFM_main')

# from firn_density_spin import FirnDensitySpin
from firn_density_nospin import FirnDensityNoSpin
import time
import json
import shutil
import RCMpkl_to_spin as RCM

def MERRA2_zarr_to_dataframe(y_int,x_int,icesheet,zarr_source=runloc):
    '''
    Create a pandas dataframe for a site in Greenland or Antarctica
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
            
            # if 'EMIS_eff' in varlist:
            #     df_sub['EMIS_eff'] = dsZ.isel(y=[ii],x=[jj])['EMIS_eff'].values.flatten()
            # else:
            for vv in varlist: # now includes EMIS_eff in main
                    df_sub[vv] = dsZ.isel(y=[ii],x=[jj])[vv].values.flatten()

            df_seconds = df_sub.index.diff().mean().total_seconds()

            df_sub['RAIN'] = (df_sub['PRECLS'] + df_sub['PRECCU']) * df_seconds
            df_sub['EVAP'] = -1 * df_sub['EVAP'] * df_seconds # multiply by -1 because of MERRA2 sign - this makes negative sign mean sublimation (mass loss)
            df_sub['PRECSN'] = df_sub['PRECSN'] * df_seconds
            df_sub['SMELT'] = df_sub['SMELT'] * df_seconds

            ### Below: use the T2M_i line to use the temps calculated from lapse rate; use the other for hte bilinear interp.
            # drn = {'T2M_i':'T2m','TS_i':'TSKIN','EVAP':'SUBLIM','HFLUX':'QH','EFLUX':'QL','SWGDN':'SW_d','LWGAB':'LW_d_M2','RAIN':'RAIN','PRECSN':'BDOT','ALBEDO':'ALBEDO_i','SMELT':'SMELT','SWGNT':'SW_n','EMIS_eff':'EMIS_eff'}
            drn = {'T2M':'T2m','TS':'TSKIN','EVAP':'SUBLIM','HFLUX':'QH','EFLUX':'QL','SWGDN':'SW_d','LWGAB':'LW_d_M2','RAIN':'RAIN','PRECSN':'BDOT','ALBEDO':'ALBEDO_i','SMELT':'SMELT','SWGNT':'SW_n','EMIS_eff':'EMIS_eff'}

            df_sub = df_sub[drn.keys()]
            df_sub.rename(mapper=drn,axis=1,inplace=True)

            return ii,jj,y_val,x_val,df_sub
        
    decades = [1980,1990,2000,2010,2020]
    df_dict = {}
    for decade in decades:
        if zarr_source=='discover':
            if icesheet=='GrIS':
                zarr_path = Path("/discover/nobackup/cdsteve2/climate/MERRA2/GrIS_emis/zarr/")
            elif icesheet=='AIS':
                zarr_path = Path("/discover/nobackup/projects/icesat2/firn/ATL_masschange/CFM_forcing/AIS/zarr")    
        elif zarr_source=='azure':
            if icesheet=='GrIS':
                zarr_path = Path("/shared/firndata/")
            elif icesheet=='AIS':
                zarr_path = Path("/shared/home/cdsteve2/firnadls/CFM_inputs/AIS/")
            
        filename = Path(zarr_path,f"M2_{icesheet}_4h_IS2mc_{decade}.zarr.zip")

        # with xr.open_dataset(filename,engine='zarr') as dsZ:
        with xr.open_zarr(filename) as dsZ:
            ii,jj,y_val,x_val,df_sub = make_dataframe(dsZ,y_int,x_int)
            df_dict[decade] = df_sub

    # with xr.open_dataset(fn_EE,engine='zarr') as fE:
    #     _ii,_jj,_y_val,_x_val,df_EE = make_dataframe(fE,y_int,x_int)        
            
    df_out = pd.concat(df_dict.values())

    df_out['QL'] = -1 * df_out['QL']
    df_out['QH'] = -1 * df_out['QH']

    sigma = 5.670374419e-8
    # df_out['LW_d_EE'] = df_EE['EMIS_eff'] * sigma * df_out['T2m']**4
    df_out['LW_d_EE'] = df_out['EMIS_eff'] * sigma * df_out['T2m']**4
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

    # hh = socket.gethostname()
    # if 'disc' in hh:
    #     runloc = 'discover'
    # else:
    #     runloc = 'azure'

    icesheet = 'AIS'
        
    seb = True
    LWdown_source = 'EMIS_eff' #EMIS_eff, MERRA2
    ALBEDO_source = 'M2_interp' #post, M2_interp
     
    if runloc =='discover':
        CFM_path = Path('/discover/nobackup/cdsteve2/ATL_masschange/CommunityFirnModel/CFM_main')
    elif runloc == 'azure':
        CFM_path = Path('/shared/home/cdsteve2/CommunityFirnModel/CFM_main/')
    elif runloc == 'local':
        CFM_path = Path('/Users/cdsteve2/research/firn/CommunityFirnModel/CFM_main')

    if runloc!='local':
        config_in = Path(CFM_path,'HPC_config_default.json')
    else:
        config_in = Path('/Users/cdsteve2/research/ATL_masschange','HPC_config_default.json')

    with open(config_in, "r") as f:
        jsonString      = f.read()
        c          = json.loads(jsonString) 

    c['runloc'] = runloc
    quad = 'A1'
    c['quad'] = quad 

    if c['runloc'] == 'azure':
        zarr_source = 'azure'
        # ll_list = np.genfromtxt(Path(CFM_path,f'IS2_icepixels_{icesheet}.csv'),delimiter=',',skip_header=1)
        ll_list = np.genfromtxt(Path(CFM_path,f'IS2_pixelstorun_{icesheet}_{quad}_full.csv'),delimiter=',',skip_header=1)
    
    elif c['runloc'] == 'discover':
        zarr_source = 'discover'
        pixel_path = Path('/discover/nobackup/cdsteve2/ATL_masschange/pixels_to_run')
        # ll_list = np.genfromtxt(Path(CFM_path,f'IS2_icepixels_{icesheet}.csv'),delimiter=',',skip_header=1)
        ll_list = np.genfromtxt(Path(pixel_path,f'IS2_pixelstorun_{icesheet}_{quad}_full.csv'),delimiter=',',skip_header=1)

    elif c['runloc'] == 'local':
        zarr_source = runloc
        pixel_path = Path('/Users/cdsteve2/research/ATL_masschange/pixels_to_run')
        # ll_list = np.genfromtxt(Path(CFM_path,f'IS2_icepixels_{icesheet}.csv'),delimiter=',',skip_header=1)
        ll_list = np.genfromtxt(Path(pixel_path,f'IS2_pixelstorun_{icesheet}_{quad}_add.csv'),delimiter=',',skip_header=1)
    
    # if c['runloc']=='local':
    #     x_int = c['x_val']
    #     y_int = c['y_val']
    # else:
    print(f'pixel number: {sys.argv[1]}')
    dkey = int(sys.argv[1]) # this is the pixel number
    x_int = float(ll_list[dkey][0])
    y_int = float(ll_list[dkey][1])

    if np.isnan(y_int):
        print('y_int is nan')
        sys.exit()

    try:
        runid = float(sys.argv[1])
    except:
        runid=-9999

    c['SEB'] = True
    calc_melt = False
    c['MELT'] = True
    
    c['physRho'] = "GSFC2020"
    c['spinUpdate'] = True

    rf_po = f'CFMresults_{quad}_{dkey}_{c["physRho"]}_LW-{LWdown_source}_ALB-{ALBEDO_source}' #results directory name

    if runloc == 'azure':
        # c['resultspath'] = '/shared/firndata/CFM_outputs' # previous GrIS outputs
        c['resultspath'] = f'/shared/home/cdsteve2/firnadls/CFM_outputs/{icesheet}_{quad}' # cheaper to put on alds
    elif runloc == 'discover':
        c['resultspath'] = f'/discover/nobackup/cdsteve2/ATL_masschange/CFMoutputs/{icesheet}_{quad}'
    elif runloc == 'local':
        c['resultspath'] = f'/Users/cdsteve2/research/ATL_masschange/CFMoutputs/{icesheet}_{quad}'

    Path(c['resultspath']).mkdir(parents=True, exist_ok=True)

    # c['resultsFolder'] = c['resultspath'] + c['results_ext'] + rf_po
    c['resultsFolder'] = str(Path(c['resultspath'], rf_po))

    if os.path.exists(Path(c['resultsFolder'],'CFMresults.hdf5')):
        rp_str = str(Path(c['resultsFolder'],'CFMresults.hdf5'))
        print(f'run has already completed at:')
        print(f'{rp_str}')
        print('exiting')
        sys.exit()

    ### Get climate data from zarr
    if runloc != 'local':
        ii,jj,y_val,x_val,df_daily = MERRA2_zarr_to_dataframe(y_int,x_int,icesheet,zarr_source=zarr_source)
        write_df = False
    else: 
        ### local requires first building the csv from a run on discover or azure (zarr too large for local)
        ii=-9999
        jj=-9999
        x_val=x_int
        y_val=y_int
        df_daily = pd.read_csv(Path('/Users/cdsteve2/research/ATL_masschange/CFMforcing',f'CFMforcing_df_{dkey}.csv'),index_col=0,parse_dates=True)
        write_df = False # This stays false

    if write_df:
        df_daily.to_csv(f'CFMforcing_df_{int(runid)}.csv')

    ### spin date end is inclusive, so e.g. if sde is 2019, it goes to 12/31/19
    if icesheet=='GrIS':
        sds = 1980.0 #spin date start
        sde = 1995.0 #spin date end
    elif icesheet=='AIS': 
        sds = 1980.0 #spin date start
        sde = 2019.0 #spin date end        
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
        
    df_daily = df_daily.drop('EMIS_eff',axis=1)

    print(ii, jj, y_val, x_val)
    print(df_daily.head())
    df_spy = 365.25*24*3600 / (df_daily.index.to_series().diff()).dt.total_seconds().mean()
    print(f'stepsperyear (az): {df_spy}')
    
    c['bdm_sublim'] = True
    
    if c['bdm_sublim']:
        bdot_mean = ((df_daily['BDOT']+(df_daily['SUBLIM']))*df_spy/917).mean()
    else:
        bdot_mean = (df_daily['BDOT']*df_spy/917).mean()

    print(f'bdot mean: {bdot_mean}')

    #######

    c['y_int'] = float(y_int)
    c['x_int'] = float(x_int)
    c['y_val'] = float(y_val)
    c['x_val'] = float(x_val)
    c['runid'] = runid
    ##########

    if bdot_mean>0.2:
        rho_bottom=916
    elif bdot_mean<0.1:
        rho_bottom=900
    else:
        rho_bottom=910
        
    climateTS, StpsPerYr, depth_S1, depth_S2, grid_bottom, SEBfluxes = (
        RCM.makeSpinFiles(df_daily,timeres=c['DFresample'],Tinterp='mean',spin_date_st = sds, 
        spin_date_end = sde,melt=c['MELT'],desired_depth = None,SEB=c['SEB'],rho_bottom=rho_bottom,calc_melt=calc_melt,bdm_sublim=c['bdm_sublim']))
    print('spin file made')

    # print(f'climateTS_keys:{climateTS.keys()}')
    # print(len(climateTS['time']))
    # print(f'sebf_keys:{SEBfluxes.keys()}')
    # print(len(SEBfluxes['time']))
    
    # if write_df:
    #     i_dec = np.where(SEBfluxes['time']>=sds)[0]
    #     df_daily['dectime'] = SEBfluxes['time'][i_dec]
    #     df_daily.to_csv(f'CFMforcing_df_{int(runid)}.csv')
    
    i1 = np.where(climateTS['time']==sds)[0][0]
    i2 = np.where(climateTS['time']==sde+1)[0][0]
    tsu = climateTS['time'][:i1]
    trep = climateTS['time'][i1:i2]
    num_reps = len(tsu)/(len(trep))
    climateTS['sds'] = sds
    climateTS['sde'] = sde
    climateTS['num_reps'] = num_reps

    c["stpsPerYear"] = float('%.2f' % (StpsPerYr))
    c["stpsPerYearSpin"] = float('%.2f' % (StpsPerYr))
    c["grid1bottom"] = float('%.1f' %(depth_S1))
    c["grid2bottom"] = float('%.1f' %(depth_S2))
    c["HbaseSpin"] = float('%.1f' %(3000 - grid_bottom))
    
    ####
    if bdot_mean>=0.15:
        print('one')
        pass
    elif ((bdot_mean>0.06) & (bdot_mean<0.15)):
        print('two')
        c["grid1bottom"] = min(3,depth_S1)
        c["grid2bottom"] = min(10,depth_S2)
        c['nodestocombine'] = 30 
        c['multnodestocombine'] = 12
    elif bdot_mean<0.03:
        print('three')
        c["grid1bottom"] = 1
        c["grid2bottom"] = 5
        c['nodestocombine'] = 180 
        c['multnodestocombine'] = 12
    else: # between 0.02 and 0.05
        print('four')
        c["grid1bottom"] = 2
        c["grid2bottom"] = 10
        c['nodestocombine'] = 90 
        c['multnodestocombine'] = 12
    ####
    print(f'depth 1: {c["grid1bottom"]}')
    print(f'depth 2: {c["grid2bottom"]}')
    
    c["NewSpin"] = False

    # configName = f'CFMconfig_{y_w}_{x_w}.json'
    configName = f'CFMconfig_{icesheet}_{dkey}_{c["physRho"]}_LW-{LWdown_source}_ALB-{ALBEDO_source}.json'
    configPath_in = Path(CFM_path,'json',configName)
    shutil.copyfile(config_in, configPath_in)
    
    if os.path.isfile(os.path.join(c['resultsFolder'],configName)):
        CFMconfig = os.path.join(c['resultsFolder'],configName)
        if os.path.isfile(os.path.join(os.getcwd(), configPath_in)):
            os.remove(os.path.join(os.getcwd(), configPath_in))
        shutil.move(CFMconfig, Path(CFM_path,'json'))
    else:
        CFMconfig = configPath_in     
    
    with open(configPath_in,'w+') as fp:
        fp.write(json.dumps(c,sort_keys=True, indent=4, separators=(',', ': ')))

    if runloc != 'local':
        if 'NewSpin' in c:
            NewSpin = c['NewSpin']
        else:
            NewSpin = False
    else: ### probably want new spin if debugging on local, or at least easy control here.
        NewSpin=True


    ### Create CFM instance by passing config file and forcing data, then run the model
    print('Configuring complete. Starting run.')
    firn = FirnDensityNoSpin(CFMconfig, climateTS = climateTS, NewSpin = NewSpin, SEBfluxes = SEBfluxes)
    firn.time_evolve()
    ###

    shutil.move(configPath_in,os.path.join(c['resultsFolder'],configName))

    print ("output folder = ", c['resultsFolder'])
    print('run time =' , time.time()-tic , 'seconds')
