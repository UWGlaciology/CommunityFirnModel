#!/usr/bin/env python
'''
writer.py
=========

Functions for writing model outputs.
'''

import csv
import os
import numpy as np
import h5py
from constants import *

def write_nospin_hdf5(self,Mout_dict):
    '''
    Write the results fromt the main model run to hdf file.

    Parameters
    ----------
    Mout_dict: dict
        contains all of the model outputs; each key is the name of the output 
    '''

    f4 = h5py.File(os.path.join(self.c['resultsFolder'], self.c['resultsFileName']),'w')

    for VW in Mout_dict.keys():

        if VW == 'rho': 
            wn = 'density'
        elif VW == 'Tz':
            wn = 'temperature'
        elif VW == 'z':
            wn = 'depth'
        elif VW == 'age':
            Mout_dict[VW] = Mout_dict[VW]/S_PER_YEAR
            wn = 'age'
        elif VW == 'climate':
            wn = 'Modelclimate'
        elif VW == 'Hx':
            wn = 'temp_Hx'
        else:
            wn = VW

        f4.create_dataset(wn, data = Mout_dict[VW])

    f4.close()


def write_spin_hdf5(self):
    '''
    Write the model outputs to hdf file at the end of spin up.
    '''

    f5 = h5py.File(os.path.join(self.c['resultsFolder'], self.c['spinFileName']), 'w')

    f5.create_dataset('densitySpin', data = self.rho_time)
    f5.create_dataset('tempSpin', data = self.Tz_time)
    f5.create_dataset('ageSpin', data = self.age_time)
    f5.create_dataset('depthSpin', data = self.z_time)
    if self.c['physGrain']:
        f5.create_dataset('r2Spin', data = self.r2_time)
    if self.THist:
        f5.create_dataset('HxSpin', data = self.Hx_time)
    if self.c['isoDiff']:
        for isotope in self.c['iso']:
            f5.create_dataset('IsoSpin_{}'.format(isotope), data = self.iso_out[isotope])
            f5.create_dataset('iso_sig2_{}'.format(isotope), data = self.iso_sig2_out[isotope])
    if self.doublegrid:
        f5.create_dataset('gridSpin', data = self.grid_time)
    if self.c['MELT']: #VV
        f5.create_dataset('LWCSpin', data = self.LWC_time)
    # if self.write_bdot:
        # f5.create_dataset('bdot_meanSpin', data = self.bdot_mean_time)
    f5.close()

def SpinUpdate(self,mtime):
    '''
    Overwrite the variables in the spin file to whatever they are 
    at time = mtime

    Parameters
    ----------
    mtime: float
        Time (model time) at which the 
    '''

    spin_results = h5py.File(os.path.join(self.c['resultsFolder'], self.c['spinFileName']),'r+')
    
    spin_results['densitySpin'][:] = np.append(mtime,self.rho)
    spin_results['tempSpin'][:]    = np.append(mtime,self.Tz)
    spin_results['ageSpin'][:]     = np.append(mtime,self.age)
    spin_results['depthSpin'][:]   = np.append(mtime,self.z)

    try:
        spin_results.create_dataset('bdot_meanSpin',data = np.append(mtime,self.bdot_mean))
    except:
        spin_results['bdot_meanSpin'][:] = np.append(mtime,self.bdot_mean)

    if self.c['MELT']:
        try:
            spin_results.create_dataset('LWCSpin',data = np.append(mtime,self.LWC))
        except:
            spin_results['LWCSpin'][:] = np.append(mtime,self.LWC)

    if self.c['physGrain']:
        try:
            spin_results['r2Spin'][:] = np.append(mtime,self.r2)
        except:
            pass

    if self.c['physRho']=='Morris2014':
        spin_results['HxSpin'][:] = np.append(mtime,self.Hx)

    if self.c['isoDiff']:
        for isotope in self.c['iso']:
            spin_results['IsoSpin_{}'.format(isotope)][:]  = np.append(mtime, self.Isoz[isotope])
            spin_results['iso_sig2_{}'.format(isotope)][:] = np.append(mtime, self.Iso_sig2_z[isotope])

    if self.doublegrid:
        spin_results['gridSpin'][:] = np.append(mtime,self.gridtrack)

    spin_results.close()