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
from reader import expand_rci, RCI_KEYS, FORCING_FILENAME

def write_nospin_hdf5(self,Mout_dict,forcing_dict=None):
    '''
    Write the results fromt the main model run to hdf file.

    Parameters
    ----------
    Mout_dict: dict
        contains all of the model outputs; each key is the name of the output 
    '''

    f4 = h5py.File(os.path.join(self.c['resultsFolder'], self.c['resultsFileName']),'w')

    if 'truncate_outputs' in self.c:
        truncate_outputs = self.c['truncate_outputs']
    else:
        truncate_outputs = False
        print('truncate_outputs not in .json. Setting to false')

    for VW in Mout_dict.keys():

        if VW == 'rho': 
            wn = 'density'
        elif VW == 'Tz':
            wn = 'temperature'
        elif VW == 'z':
            wn = 'depth'
        elif VW == 'age':
            Mout_dict[VW][:,1:] = Mout_dict[VW][:,1:] / S_PER_YEAR
            wn = 'age'
        elif VW == 'climate':
            wn = 'Modelclimate'
        elif VW == 'Hx':
            wn = 'temp_Hx'
        else:
            wn = VW

        subvars = ['rho','Tz','LWC','age']

        if ((VW in subvars) and (truncate_outputs)):
            # data_out = np.column_stack((Mout_dict[VW][:,0],Mout_dict[VW][:,1::5]))
            data_out = np.vstack((Mout_dict[VW][0,:],Mout_dict[VW][1::5,:]))

        else:
            data_out = Mout_dict[VW]

        f4.create_dataset(wn, data = data_out)

    # if forcing_dict:
    #     ks = list(forcing_dict)
    #     ll = len(forcing_dict[ks[0]])
    #     forcing_out = np.zeros([ll,6])
    #     forcing_out[:,0] = forcing_dict['dectime']
    #     forcing_out[:,1] = forcing_dict['TSKIN']
    #     forcing_out[:,2] = forcing_dict['BDOT']
    #     try:
    #         forcing_out[:,3] = forcing_dict['SMELT']
    #     except:
    #         forcing_out[:,3] = -9999* np.ones_like(forcing_dict['dectime'])
    #     try:
    #         forcing_out[:,4] = forcing_dict['RAIN']
    #     except:
    #         forcing_out[:,4] = -9999* np.ones_like(forcing_dict['dectime'])
    #     try:
    #         forcing_out[:,5] = forcing_dict['SUBLIM']
    #     except:
    #         forcing_out[:,5] = -9999* np.ones_like(forcing_dict['dectime'])
    #     f4.create_dataset('forcing',data=forcing_out,dtype='float64')

    f4.close()

def write_spin_hdf5(self):
    '''
    Write the model outputs to hdf file at the end of spin up.
    '''

    f5 = h5py.File(os.path.join(self.c['resultsFolder'], self.c['spinFileName']), 'w')

    vec_len = len(self.rho_time)

    f5.create_dataset('densitySpin', data = self.rho_time[...,None],maxshape=(vec_len,None,))
    f5.create_dataset('tempSpin', data = self.Tz_time[...,None],maxshape=(vec_len,None,))
    f5.create_dataset('ageSpin', data = self.age_time[...,None],maxshape=(vec_len,None,))
    f5.create_dataset('depthSpin', data = self.z_time[...,None],maxshape=(vec_len,None,))
    if self.c['physGrain']:
        f5.create_dataset('r2Spin', data = self.r2_time[...,None],maxshape=(vec_len,None,))
    if self.THist:
        f5.create_dataset('HxSpin', data = self.Hx_time[...,None],maxshape=(vec_len,None,))
    if self.c['isoDiff']:
        for isotope in self.c['iso']:
            f5.create_dataset(f'IsoSpin_{isotope}', data = self.iso_out[isotope][...,None],maxshape=(vec_len,None,))
            f5.create_dataset(f'iso_sig2_{isotope}', data = self.iso_sig2_out[isotope][...,None],maxshape=(vec_len,None,))
    if self.doublegrid:
        f5.create_dataset('gridSpin', data = self.grid_time[...,None],maxshape=(vec_len,None,))
    if self.c['MELT']: #VV
        f5.create_dataset('LWCSpin', data = self.LWC_time[...,None],maxshape=(vec_len,None,))
    # if self.write_bdot:
        # f5.create_dataset('bdot_meanSpin', data = self.bdot_mean_time)
    f5.close()

# def write_spin_hdf5(self):
#     '''
#     older version without maxshape
#     Write the model outputs to hdf file at the end of spin up.
#     '''

#     f5 = h5py.File(os.path.join(self.c['resultsFolder'], self.c['spinFileName']), 'w')

#     f5.create_dataset('densitySpin', data = self.rho_time)
#     f5.create_dataset('tempSpin', data = self.Tz_time)
#     f5.create_dataset('ageSpin', data = self.age_time)
#     f5.create_dataset('depthSpin', data = self.z_time)
#     if self.c['physGrain']:
#         f5.create_dataset('r2Spin', data = self.r2_time)
#     if self.THist:
#         f5.create_dataset('HxSpin', data = self.Hx_time)
#     if self.c['isoDiff']:
#         for isotope in self.c['iso']:
#             f5.create_dataset(f'IsoSpin_{isotope}', data = self.iso_out[isotope])
#             f5.create_dataset(f'iso_sig2_{isotope}', data = self.iso_sig2_out[isotope])
#     if self.doublegrid:
#         f5.create_dataset('gridSpin', data = self.grid_time)
#     if self.c['MELT']: #VV
#         f5.create_dataset('LWCSpin', data = self.LWC_time)
#     # if self.write_bdot:
#         # f5.create_dataset('bdot_meanSpin', data = self.bdot_mean_time)
#     f5.close()

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
            spin_results[f'IsoSpin_{isotope}'][:]  = np.append(mtime, self.Isoz[isotope])
            spin_results[f'iso_sig2_{isotope}'][:] = np.append(mtime, self.Iso_sig2_z[isotope])

    if self.doublegrid:
        spin_results['gridSpin'][:] = np.append(mtime,self.gridtrack)

    spin_results.close()

def SpinUpdate_res(self,mtime):
    '''
    keep model snapshots in the spin file (restarts)
    at time = mtime

    Parameters
    ----------
    mtime: float
        Time (model time) at which the 
    '''

    spin_results = h5py.File(os.path.join(self.c['resultsFolder'], self.c['spinFileName']),'a')

    vec_len = spin_results['densitySpin'].shape[0]
    res_dim = spin_results['densitySpin'].shape[1]

    spin_results['densitySpin'].resize((vec_len,res_dim+1))
    spin_results['tempSpin'].resize((vec_len,res_dim+1))
    spin_results['ageSpin'].resize((vec_len,res_dim+1))
    spin_results['depthSpin'].resize((vec_len,res_dim+1))

    spin_results['densitySpin'][:,-1] = np.append(mtime,self.rho)
    spin_results['tempSpin'][:,-1]    = np.append(mtime,self.Tz)
    spin_results['ageSpin'][:,-1]     = np.append(mtime,self.age)
    spin_results['depthSpin'][:,-1]   = np.append(mtime,self.z)

    if 'bdot_meanSpin' not in spin_results.keys():
        b_ms = np.append(mtime,self.bdot_mean)
        spin_results.create_dataset('bdot_meanSpin', data = b_ms[...,None], maxshape=(vec_len,None,))
    else:
        spin_results['bdot_meanSpin'].resize((vec_len,res_dim+1))
        spin_results['bdot_meanSpin'][:,-1] = np.append(mtime,self.bdot_mean)

    if self.c['MELT']:
        # try:
        #     spin_results.create_dataset('LWCSpin',data = np.append(mtime,self.LWC))
        # except:
        spin_results['LWCSpin'].resize((vec_len,res_dim+1))
        spin_results['LWCSpin'][:,-1] = np.append(mtime,self.LWC)

    if self.c['physGrain']:
        try:
            spin_results['r2Spin'].resize((vec_len,res_dim+1))
            spin_results['r2Spin'][:,-1] = np.append(mtime,self.r2)
        except:
            pass

    if self.c['physRho']=='Morris2014':
        spin_results['HxSpin'].resize((vec_len,res_dim+1))
        spin_results['HxSpin'][:,-1] = np.append(mtime,self.Hx)

    if self.c['isoDiff']:
        for isotope in self.c['iso']:
            spin_results[f'IsoSpin_{isotope}'].resize((vec_len,res_dim+1))
            spin_results[f'iso_sig2_{isotope}'].resize((vec_len,res_dim+1))
            spin_results[f'IsoSpin_{isotope}'][:,-1]  = np.append(mtime, self.Isoz[isotope])
            spin_results[f'iso_sig2_{isotope}'][:,-1] = np.append(mtime, self.Iso_sig2_z[isotope])

    if self.doublegrid:
        spin_results['gridSpin'].resize((vec_len,res_dim+1))
        spin_results['gridSpin'][:,-1] = np.append(mtime,self.gridtrack)

    spin_results.close()

MAX_PATCH = 10 # the most elements per variable that the RCI repeats may fail to reproduce (see _trim_rci)

def _trim_rci(d, sds, sde, num_reps):
    '''
    Remove the repeats of the reference climate interval from a forcing series.

    The record is what follows the repeats (it starts at sds). The record is only used if reader.expand_rci
    gives back the series d from it. A few elements may differ: for the mass fluxes of the series without
    the surface energy balance, makeSpinFiles gives the first time step the mean step length (flux_to_rate),
    and that is not what the repeat of the first year of the record has. Those elements are returned
    as a patch, so that the series is still exactly the same.

    :return: (record, patch), where patch is {key: (indices, values)}; or None if the series cannot be
        reproduced from its record (then the whole series has to be stored)
    '''
    t = d['time']
    if any(np.ndim(v) > 0 and len(v) != len(t) for v in d.values()):
        return None
    i0 = np.searchsorted(t, sds)
    record = {k: (v[i0:] if np.ndim(v) > 0 else v) for k, v in d.items()}
    rebuilt = expand_rci(record, sds, sde, num_reps)

    patch = {}
    for k, v in d.items():
        if np.ndim(v) == 0:
            continue
        r = rebuilt[k]
        if r.shape != v.shape or r.dtype != v.dtype:
            return None
        same = (r == v)
        if np.issubdtype(v.dtype, np.floating):
            same = same | (np.isnan(r) & np.isnan(v))
        bad = np.flatnonzero(~same)
        if len(bad) > MAX_PATCH:
            return None
        if len(bad) > 0:
            patch[k] = (bad, v[bad])
    return record, patch


def _write_series(f, name, d, rci):
    '''
    Write one forcing series (a dictionary of arrays and scalars) as the group name of the file f. If rci
    is (sds, sde, num_reps), the repeats of the reference climate interval are left out when that is possible.
    '''
    got = _trim_rci(d, *rci) if rci is not None else None
    if rci is not None and got is None:
        print(f'forcing_writer: the repeats of the reference climate interval cannot be left out of "{name}" '
              f'(it does not start at sds, or it was changed). The whole series is stored.')
    record, patch = got if got is not None else (d, {})

    g = f.create_group(name)
    for k, v in record.items():
        g.create_dataset(k, data=v)
    g.attrs['rci_trimmed'] = got is not None
    if got is not None:
        for k, v in zip(RCI_KEYS, rci):
            g.create_dataset(k, data=v)
        for k, (idx, val) in patch.items():
            g.create_dataset(f'patch/{k}/index', data=idx)
            g.create_dataset(f'patch/{k}/value', data=val)


def forcing_writer(self, climateTS, SEBfluxes = None, strain_rates = None):

    '''
    write forcing data to its own hdf5 file
    units on mass fluxes (e.g., bdot, rain, etc.) are m ice eq. per year.

    Together with the configuration file, this file is what is needed to repeat the run (reader.read_forcing
    reads it). The forcing of a run with a reference climate interval (RCI) spin-up is the RCI repeated
    num_reps times and then the record. The repeats are not stored: the file has the record, and sds, sde
    and num_reps (which makeSpinFiles puts in climateTS) to make the repeats again. If that is not
    possible without changing the series, the whole series is written.

    Groups: main (climateTS), SEB (SEBfluxes, if given), strain_rates (the strain_rates argument, if given).
    '''

    try:
        forcing_filename = self.c['forcingFileName']
    except:
        print('forcing_filename not in json. Defaulting to CFMforcing.hdf5')
        forcing_filename = FORCING_FILENAME

    rci = None
    if all(k in climateTS for k in RCI_KEYS):
        rci = (float(climateTS['sds']), float(climateTS['sde']), int(climateTS['num_reps']))
    series = {k: v for k, v in climateTS.items() if k not in RCI_KEYS}

    with h5py.File(os.path.join(self.c['resultsFolder'], forcing_filename), 'w') as f6:
        _write_series(f6, 'main', series, rci)
        if SEBfluxes is not None:
            _write_series(f6, 'SEB', SEBfluxes, rci)
        if strain_rates is not None:
            if isinstance(strain_rates, dict):
                for k, v in strain_rates.items():
                    f6.create_dataset(f'strain_rates/{k}', data=np.asarray(v))
            else:
                f6.create_dataset('strain_rates', data=np.asarray(strain_rates, dtype=float))
