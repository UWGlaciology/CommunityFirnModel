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
    Write the results fromt the main model run
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
    
    # if 'rho' in self.output_list:
    #     f4.create_dataset('density', data = Mout_dict['rho'])
    # if 'Tz' in self.output_list:
    #     f4.create_dataset('temperature', data = Mout_dict['Tz'])
    # if 'age' in self.output_list:
    #     f4.create_dataset('age', data = Mout_dict['age']) # use this one if you want a matrix of ages
    #     # f4.create_dataset('age', data = Mout_dict['age'][-1,:]) # use this one if you want just the last row 
    # if 'z' in self.output_list:    
    #     f4.create_dataset('depth', data = Mout_dict['z'])
    # if 'Dcon' in self.output_list:    
    #     f4.create_dataset('Dcon', data = Mout_dict['Dcon'])
    # if 'bdot_mean' in self.output_list:    
    #     f4.create_dataset('bdot_mean', data = Mout_dict['bdot_mean'])
    # if 'climate' in self.output_list:    
    #     f4.create_dataset('Modelclimate', data = Mout_dict['climate'])
    # if 'compaction' in self.output_list:    
    #     f4.create_dataset('compaction', data = Mout_dict['compaction'])

    # if self.c['FirnAir']:
    #     for gas in self.cg['gaschoice']:       
    #         f4.create_dataset(gas, data = Mout_dict[gas])
    #     if "diffusivity" in self.cg['outputs']:
    #         f4.create_dataset('diffusivity', data = Mout_dict['diffusivity'])
    #     if "gas_age" in self.cg['outputs']:
    #         f4.create_dataset('gas_age', data = Mout_dict['gas_age'])
    #     if "advection_rate" in self.cg['outputs']:
    #         f4.create_dataset('w_air', data = self.w_air_out)
    #         f4.create_dataset('w_firn', data = self.w_firn_out)
    # if 'grainsize' in self.output_list:
    #     f4.create_dataset('r2', data = Mout_dict['r2'])
    #     f4.create_dataset('dr2_dt', data = Mout_dict['dr2_dt'])
    # if 'temp_Hx' in self.output_list:    
    #     f4.create_dataset('temp_Hx',data = Mout_dict['Hx'])
    # if self.c['isoDiff']:
    #     for isotope in self.c['iso']:
    #         f4.create_dataset('isotopes_{}'.format(isotope), data = self.iso_out[isotope])
    #         f4.create_dataset('iso_sig2_{}'.format(isotope), data = self.iso_sig2_out[isotope])
    # if 'LWC' in self.output_list:
    #     f4.create_dataset('LWC',data = Mout_dict['LWC'])
    # if 'PLWC_mem' in self.output_list:
    #     f4.create_dataset('PLWC_mem',data = Mout_dict['PLWC_mem'])
    # if 'DIP' in self.output_list:
    #     f4.create_dataset('DIP',data = Mout_dict['DIP'])
    # if 'DIPc' in self.output_list:
    #     f4.create_dataset('DIPc',data = Mout_dict['DIPc'])
    # if 'BCO' in self.output_list:
    #     f4.create_dataset('BCO',data = Mout_dict['BCO']) 
    # if 'viscosity' in self.output_list:
    #     f4.create_dataset('viscosity',data = Mout_dict['viscosity'])
    # if 'meltoutputs' in self.output_list:
    #     f4.create_dataset('runoff',data = Mout_dict['runoff'])
    #     f4.create_dataset('refrozen',data = Mout_dict['refrozen'])
    #     # f4.create_dataset('icecon',data = self.icecon_out)
    #     # f4.create_dataset('trfrz',data = self.trfrz_out)
    #     # f4.create_dataset('tfac',data = self.tfac_out)
    #     # f4.create_dataset('tlwc',data = self.tlwc_out) 
    #     # f4.create_dataset('totcumrunoff',data = self.totcumrunoff_out)
    #     # f4.create_dataset('cumrefrozen',data = self.cumrefrozen_out)
    f4.close()


def write_spin_hdf5(self):
    '''
    Write the model outputs at the end of spin up
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
    # tind = np.where(self.rho_out[:,0]>=self.c['spinUpdateDate'])[0][0]

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
