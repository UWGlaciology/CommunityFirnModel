#!/usr/bin/env python

import numpy as np 
import pandas as pd
import xarray as xr
import glob
import io
import multiprocessing
import sys
import os
import time
import json
import shutil
from pathlib import Path

"""
run_CFM_example.py
=======
This script is very similar to run_CFM_example_notebook.ipynb. However, I consider
the notebook to be more of a user-friendly introduction to running the CFM for a 
specific site, whereas this script is more of an example of how I run the CFM in batch mode
(e.g., calling this script a bunch of times in a slurm job). As such, there are a number of 
lines of code in here that probably won't work and/or do anything for most users, but I've 
kept them in here for illustrative purposes. This script was more or less designed to 
run CFM when forced with MERRA-2 data. 

Feel free to email me (maxstev@umd.edu) if you have questions. 

=======

This file configures a CFM run and then runs it.
- it goes in CFM_main 
- it will work with the 2 provided pkl or csv files
    - CFM_example_72.5_-38.75.pkl/CFM_example_72.5_-38.75.csv (Summit)
    - CFM_example_66.5_-46.25.pkl/CFM_example_66.5_-46.25.csv (DYE-2)

- run this script using:
>>>python run_CFM_example.py "LAT,LON" 
where LAT,LON are your latitude and longitude, e.g. "66.5,-46.25" or "72.5,-38.75"
(include the quote marks when invoking the above call)

- Forcing data:
    - this example file uses a pickled dataframe with climate data (MERRA-2) for DYE-2 to force the CFM
    - these climate data come from a zarr store with MERRA-2 climate data
    - That zarr contains MERRA-2 climate fields that I have subsetted, getting just the variables
      from MERRA-2 that CFM needs for either Antarctica or Greenland.
    - the zarr store is here: https://zenodo.org/records/17317018. Note that the zarr is divided into 5 zip
      stores because of the file size. If you download them, put the 5 zips into a directory and update the 
      zarr_path below.

- this script creates a .json configuration file
- CFM runs with the climate input and config file, results get put
  in the specified directory.

But, it should be a bit more of an example of how CFM run(s) could be integrated as part
of a workflow by creating a CFM class instance.

A key bit to using the CFM with this or a similar script is that all of the forcing data goes into a 
pandas data frame with a datetime index.

In that dataframe, mass fluxes (ie., precip, sublimation) need to be in units kg/m2/dt, 
where dt is the the time delta of the dataframe. I.e., if your dataframe 
is daily resolution, the units are kg/m2/day.

Hopefully, using this example, if you want to you can 
make a working dataframe with your own climate data. 
Email me if you have questions or issues! maxstev@umd.edu

"""

### CHANGE THESE PATHS TO MATCH YOUR FILESYSTEM
### set paths:
# cfm_path = Path('/Path/To/CommunityFirnModel/CFM_main') #You can use this if you want this notebook in a directory other than CFM_main
# zarr_path = Path('/Path/To/zarr')

cfm_path = Path('/Users/cdsteve2/research/firn/CommunityFirnModel/CFM_main')
zarr_path = Path('/Users/cdsteve2/nobackup/RCMdata/MERRA2/GrIS/zarr')
###

sys.path.append(str(cfm_path))

from firn_density_nospin import FirnDensityNoSpin
import RCMpkl_to_spin as RCM

##########
### below function gets climate data from the zarr store
### will not work unless you get the zarr from me
def MERRA2_zarr_to_dataframe(lat_int,lon_int, zarr_path=None):
    '''
    Create a pandas dataframe for a site in Greenland
    returns:
    df_daily: a dataframe holding all of the forcing fields needed for CFM run
    
    input mass fluxes from zarr are in kg/m2/s.
    mass fluxes in df_daily are in kg/m2/timestep
    energy fluxes are W/m2
    
    in MERRA, ?positive? sublimation flux = sublimation, ?negative? = deposition
    '''
    def make_dataframe(dsZ,lat_int,lon_int):
            lat_ll = dsZ.lat.data
            lon_ll = dsZ.lon.data
            ii, lat_val = min(enumerate(lat_ll), key=lambda x: abs(x[1]-lat_int))
            jj, lon_val = min(enumerate(lon_ll), key=lambda x: abs(x[1]-lon_int))

            df_sub = pd.DataFrame(index=dsZ.time.values)

            varlist = [jj for jj in dsZ.variables if jj not in ['time','lat','lon']]

            for vv in varlist:
                df_sub[vv] = dsZ.isel(lat=[ii],lon=[jj])[vv].values.flatten()

            df_sub['RAIN'] = (df_sub['PRECLS']+df_sub['PRECCU'])*4*3600
            df_sub['EVAP'] = df_sub['EVAP']*4*3600
            df_sub['PRECSN'] = df_sub['PRECSN']*4*3600
            df_sub['SMELT'] = df_sub['SMELT']*4*3600

            drn = {'T2M':'T2m','TS':'TSKIN','EVAP':'SUBLIM','HFLUX':'QH','EFLUX':'QL','SWGDN':'SW_d','LWGAB':'LW_d','RAIN':'RAIN','PRECSN':'BDOT','ALBEDO':'ALBEDO','SMELT':'SMELT'}

            df_sub = df_sub[drn.keys()]
            df_sub.rename(mapper=drn,axis=1,inplace=True)

            return ii,jj,lat_val,lon_val,df_sub
        
    decades = [1980,1990,2000,2010,2020]
    df_dict = {}

    if zarr_path is None: # this is old behavior. Now zarr_path is set at beginning of script.
        zarr_path = Path('/Users/cdsteve2/nobackup/RCMdata/MERRA2/GrIS/zarr/')

    for decade in decades:
        file = Path(zarr_path,f'M2_GrIS_daily_IS2mc_{decade}.zarr.zip')
        with xr.open_dataset(file,engine='zarr') as dsZ:
            ii,jj,lat_val,lon_val,df_sub = make_dataframe(dsZ,lat_int,lon_int)
            df_dict[decade] = df_sub
            
    df_out = pd.concat(df_dict.values())

    df_out['QL'] = -1 * df_out['QL']
    df_out['QH'] = -1 * df_out['QH']

    return ii,jj,lat_val,lon_val,df_out

#################################
#################################
#################################

class M2_CFM():
    def __init__(self):
        pass

    def run_CFM(self, input_coords):

        lon_int = float(input_coords.split(",")[1])
        lat_int = float(input_coords.split(",")[0])

        runloc='local' #'local','SMCE','MAAP','loki' # if I am working on a HPC, I will code in a switch (if/else statements setting paths) to run locally for testing and on the remote
        
        seb = True
        RCMtype = 'MERRA2'
        variable_srho = False

        try:
            runid = float(sys.argv[2]) # option to assign an ID to the run to keep track of things. 
        except:
            runid=-9999

        climate_source='dataframe' # 'dataframe' (option 2 below) or 'zarr' (option 1 below)

        if climate_source=='zarr':
            #############################
            ### OPTION 1: GET CLIMATE FORCING FROM THE ZARR
            ### Below is to generate the forcing data from from zarr
            ### get zarr from: https://zenodo.org/records/17317018

            ### Call the function to create a pandas dataframe with climte forcing
            ### (need the zarr for this)
            ### optionally, save the dataframe as a pkl
            ### returns "lat_val" and "lon_val", which are "latitude value" 
            ### and "longitude value", i.e. the MERRA-2 grid point
            ### closest to lat_int and lon_int.

            ii,jj,lat_val,lon_val,df_daily = MERRA2_zarr_to_dataframe(lat_int,lon_int,zarr_path=zarr_path) # make sure zarr_path is set correctly (first cell)
            print(ii, jj, lat_val, lon_val)

            save_df=True
            save_format = '.pkl' # .csv or pkl

            if save_df:
                if save_format=='.pkl':
                    pd.to_pickle(df_daily,f'CFM_example_{lat_val}_{lon_val}.pkl')
                else: # csv
                    df_daily.to_csv(f'CFMinput_example/CFM_example_{lat_val}_{lon_val}.csv')

            ### END OPTION 1
            #############################

        elif climate_source=='dataframe':
            #############################
            ### OPTION 2: LOAD CLIMATE DATA FROM AN EXTANT DATAFRAME
            ### (which is saved in .pkl or .csv file)
            '''
            If you already have the forcing data for the lat/lon pair pickled,
            you can instead load it here.

            This example notebooke only works with the provided example files.

            But, you can also create your own dataframe with climate data
            and load it here
            If you are doing this you can mimic the format of the example files.
            '''

            lat_val = lat_int
            lon_val = lon_int

            ### e.g.:
            df_daily = pd.read_pickle(f'CFMinput_example/CFM_example_{lat_val}_{lon_val}.pkl') # if your dataframe is in a pkl
            # df_daily = pd.read_csv(f'CFMinput_example/CFM_example_{lat_val}_{lon_val}.csv') # if your dataframe is in a csv
            
            ### END OPTION 2
            #############################

        print(df_daily.head())
        df_spy = 365.25
        bdot_mean = (df_daily['BDOT']*df_spy/917).mean()
        print(f'bdot mean: {bdot_mean}')

        if RCMtype=='MERRA2':
            sds = 1980.0 #spin date start
            sde = 1985.0 #spin date end
            lat_w = lat_val
            lon_w = lon_val

        #######
        ### Prepare config .json (which is a dictionary called c within this python script) ###
        ### start by importing the example.json (ie., load defaults), and then edit just the ones you want to
        ### the edited json will be saved and used for the run.

        config_in = 'example_df.json' # This is an example .json with defaults loaded in it.

        with open(config_in, "r") as f:
            jsonString      = f.read()
            c          = json.loads(jsonString) 

        c['physRho'] = 'GSFC2020'
        c['runID'] = runid
        c['DFresample'] = '1d' # resolution of the model run, e.g. '1d' is 1 day.
        
        c['SEB'] = seb
        c['lat_int'] = float(lat_int)
        c['lon_int'] = float(lon_int)
        c['lat_val'] = float(lat_val)
        c['lon_val'] = float(lon_val)
        
        '''
        CFM regrids (merges) deeper nodes to save computation. There are 2 mergings
        nodestocombine and multnodestocombine should be adjusted based on the time resolution of the run
        e.g. if DFresample is '1d', nodestocombine = 30 will combine 30 layers at an intermediate depth, 
        and multnodestocombine = 12 will combine 12 of those layers at a greater depth (which in this case 
        will give 3 sections of firn - near the surface very thin layers, representing a day's accumulation,
        middle, which is a month's accumulation, and deep, that should be a year's accumulation. 
        e.g. if I am doing DFresample = '5d', I would set nodestocombine to 6 to still get layers that are a
        month's worth of accumulation.
        '''

        c['MELT'] = True #you can set to false to make things run a lot faster if don't care about results

        ### surface density
        if variable_srho:
            c['variable_srho'] = True
            c['srho_type'] = "noise"
        else:
            c['rhos0'] = 150.0 #e.g here you could change the surface density
            rhotype=f"rho{c['rhos0']}"
        #######

        c["physGrain"] =  False
        c["calcGrainSize"] = False

        ### Specify where results should go ###
    
        rf_pre = 'CFMoutput_example'
        rf_po = f'/CFMresults_{lat_w}_{lon_w}_{c["physRho"]}'

        c['resultsFolder'] = rf_pre + rf_po

        ##########

        climateTS, StpsPerYr, depth_S1, depth_S2, grid_bottom, SEBfluxes = (
            RCM.makeSpinFiles(df_daily,timeres=c['DFresample'],Tinterp='mean',spin_date_st = sds, 
            spin_date_end = sde,melt=c['MELT'],desired_depth = None,SEB=seb,rho_bottom=900))

        climateTS['SUBLIM'] = -1 * climateTS['SUBLIM'] #ADDED THIS FOR MERRA2 TO GET THE SIGN CORRECT.

        climateTS['forcing_data_start'] = sds

        c["stpsPerYear"] = float('%.2f' % (StpsPerYr))
        c["stpsPerYearSpin"] = float('%.2f' % (StpsPerYr))
        c["grid1bottom"] = float('%.1f' %(depth_S1))
        c["grid2bottom"] = float('%.1f' %(depth_S2))

        c['nodestocombine'] = 30 
        c['multnodestocombine'] = 12
        
        c["HbaseSpin"] = float('%.1f' %(3000 - grid_bottom))

        c['keep_firnthickness'] = True
        c['grid_outputs'] = True
        c['grid_output_res'] = 0.05

        print(f'grid_bottom={grid_bottom}')
        print(f'grid1bottom={c["grid1bottom"]}')
        print(f'grid2bottom={c["grid2bottom"]}')

        configName = f'CFMconfig_{lat_w}_{lon_w}_{c["physRho"]}_{RCMtype}_example.json'
        if os.path.isfile(os.path.join(c['resultsFolder'],configName)):
            CFMconfig = os.path.join(c['resultsFolder'],configName)
            if os.path.isfile(os.path.join(os.getcwd(), configName)):
                os.remove(os.path.join(os.getcwd(), configName))
            shutil.move(CFMconfig, os.getcwd())
        else:
            CFMconfig = configName

        with open(CFMconfig,'w') as fp:
            fp.write(json.dumps(c,sort_keys=True, indent=4, separators=(',', ': ')))

        if 'NewSpin' in c:
            NewSpin = c['NewSpin']
        else:
            NewSpin = False

        ### Create CFM instance by passing config file and forcing data, then run the model
        firn = FirnDensityNoSpin(CFMconfig, climateTS = climateTS, NewSpin = NewSpin, SEBfluxes = SEBfluxes)
        firn.time_evolve()
        ###

        shutil.move(configName,os.path.join(c['resultsFolder'],configName))

        print ("output folder = ", c['resultsFolder'])
        # print('run time =' , time.time()-tic , 'seconds')

def unwrap_fun(ll):
    CFM = M2_CFM()
    return CFM.run_CFM(ll)
    # print(ll)

if __name__ == '__main__':

        ### a bit of code below that can but used with multiprocessing

        ### all of the input lat/lon pairs need to be strings in a list.
    # with open('M2_tr.txt') as file:
    #     ll_list = [line.rstrip() for line in file]

    # NUM_PROCESSORS = 1 #This is the number of processors to use. set it here manually

    # ### or set it automatically to be 2 fewer than the number on the server:
    # # NUM_PROCESSORS = psutil.cpu_count(logical=False) - 2

    # with multiprocessing.Pool(processes=NUM_PROCESSORS) as p:
    #     p.map(unwrap_fun,ll_list)
    #     p.close()
    #     p.join()

    tic=time.time()

    CFM = M2_CFM()
    CFM.run_CFM(sys.argv[1])
    print('run time =' , time.time()-tic , 'seconds')